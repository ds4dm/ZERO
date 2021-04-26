/* #############################################
 *             This file is part of
 *                    ZERO
 *
 *             Copyright (c) 2020
 *     Released under the Creative Commons
 *        Zero v1.0 Universal License
 *
 *              Find out more at
 *        https://github.com/ds4dm/ZERO
 * #############################################*/


#include <armadillo>
#include <games/ipg.h>
#include <interfaces/ipg_models.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <zero.h>


int main(int argc, char **argv) {

  loguru::g_stderr_verbosity                   = -5;
  std::vector<unsigned int>          ItemSizes = {10, 20, 40, 80, 100};
  std::vector<unsigned int>          Players   = {2, 3};
  unsigned int                       insUB     = 10;
  bool                               print     = false;
  GRBEnv                             test;
  std::string                        BasePath      = "dat/KPMarg/dest/";
  std::string                        MargInstances = "dat/KPMarg/source/";
  std::string                        Logfile       = "dat/KPMarg/resultsMIPCuts.csv";
  std::random_device                 random;
  std::mt19937                       Oracle(random());
  std::uniform_int_distribution<int> hRandom(-100, 100);
  std::uniform_int_distribution<int> fRandom(-50, 50);

  for (auto &&m : Players) {
	 for (unsigned int n : ItemSizes) {
		for (unsigned int ins = 0; ins < insUB; ++ins) {

		  std::string InstanceName =
				"Instance_" + std::to_string(m) + "-" + std::to_string(n) + "-" + std::to_string(ins);
		  {
			 try {

				// InstanceName      = "Instance_2-80-9";
				std::string filec = BasePath + InstanceName + ".json";
				std::cout << filec << std::endl;
				if (access(filec.c_str(), F_OK) == 0) {

				  std::vector<Data::LCP::Algorithms> algos = {Data::LCP::Algorithms::MIP,
																			 // Data::LCP::Algorithms::MINLP,
																			 Data::LCP::Algorithms::PATH};

				  std::vector<Data::IPG::Objectives> objectives = {Data::IPG::Objectives::Quadratic,
																					// Data::IPG::Objectives::Linear,
																					Data::IPG::Objectives::Feasibility};

				  std::vector<Data::IPG::CutsAggressiveness> cuts = {
						Data::IPG::CutsAggressiveness::KeepItCool,
						Data::IPG::CutsAggressiveness::Truculent};

				  Models::IPG::IPGInstance Instance;
				  Instance.load(BasePath + InstanceName);
				  for (auto algo : algos) {
					 for (auto obj : objectives) {
						for (auto cut : cuts) {
						  if (!((obj == Data::IPG::Objectives::Quadratic ||
									obj == Data::IPG::Objectives::Linear) &&
								  algo == Data::LCP::Algorithms::PATH)) {
							 std::cout << "\t Testing: Algorithm " << std::to_string(algo)
										  << " with objective " << std::to_string(obj) << ", CutsLevel"
										  << std::to_string(cut) << "\n";
							 Models::IPG::IPG Test(&test, Instance);
							 Test.setAlgorithm(Data::IPG::Algorithms::Oracle);
							 Test.setDeviationTolerance(3e-4);
							 Test.setNumThreads(8);
							 Test.setCutsAggressiveness(cut);
							 Test.setTimeLimit(300);
							 Test.setLCPAlgorithm(algo);
							 Test.setGameObjective(obj);
							 Test.finalize();
							 std::chrono::steady_clock::time_point begin =
								  std::chrono::steady_clock::now();
							 Test.findNashEq();
							 std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

							 // std::cout <<
							 // std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count()/1000000.0
							 // <<"s";

							 auto stat = Test.getStatistics();

							 // std::cout << "\n" <<stat.WallClockTime.get();
							 std::ifstream existCheck(Logfile);
							 std::ofstream results(Logfile, std::ios::app);

							 auto cuts      = stat.AlgorithmData.Cuts.get();
							 int  totalCuts = 0;
							 for (auto cut : cuts)
								totalCuts += cut.second;

							 if (!existCheck.good()) {
								results << "Instance;Algorithm;LCPAlgo;ObjType;Cuts;Players;isPureNE;"
											  "Status;"
											  "NumVar;ClockTime(s);Threads;numIterations;ValueCuts;VPolyCuts;"
											  "MIRCuts;GMICuts;KPCuts;TotalCuts;\n";
							 }
							 existCheck.close();


							 results << InstanceName << ";"
										<< std::to_string(stat.AlgorithmData.Algorithm.get()) << ";"
										<< std::to_string(algo) << ";" << std::to_string(obj) << ";"
										<< std::to_string(cut) << ";" << Test.getNumPlayers() << ";"
										<< std::to_string(stat.PureNashEquilibrium.get()) << ";"
										<< std::to_string(stat.Status.get()) << ";" << stat.NumVar.get()
										<< ";" << std::to_string(stat.WallClockTime.get()) << ";"
										<< std::to_string(stat.AlgorithmData.Threads.get()) << ";"
										<< stat.NumIterations.get() << ";" << cuts.at(0).second << ";"
										<< cuts.at(1).second << ";" << cuts.at(2).second << ";"
										<< cuts.at(3).second << ";" << cuts.at(4).second << ";" << totalCuts
										<< "\n";



							 results.close();
						  } // endif objectives and PATH
						}   // algo
					 }     // cuts
				  }
				}
			 } catch (ZEROException e) {
				std::cout << "ZEROException: " << std::to_string(e.which()) << " | " << e.more();
			 }
		  }
		}
	 }
  }
}
