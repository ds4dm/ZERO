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

  loguru::g_stderr_verbosity                   = 0;
  std::vector<unsigned int>          ItemSizes = {10, 20, 40, 80, 100};
  std::vector<unsigned int>          Players   = {2, 3};
  unsigned int                       insUB     = 10;
  bool                               print     = false;
  GRBEnv                             test;
  std::string                        BasePath      = "dat/KPMarg/dest/";
  std::string                        MargInstances = "dat/KPMarg/source/";
  std::string                        Logfile       = "dat/KPMarg/results.csv";
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


				std::string filec = BasePath + InstanceName + ".json";
				std::cout << "-------------------" << filec << "------------------\n\n";
				if (access(filec.c_str(), F_OK) == 0) {

				  std::vector<Data::LCP::Algorithms> algos = {
						//Data::LCP::Algorithms::MIP,
						//Data::LCP::Algorithms::MINLP,
						Data::LCP::Algorithms::PATH
				  };

				  std::vector<Data::IPG::Objectives> objectives = {
						//Data::IPG::Objectives::Quadratic,
						//Data::IPG::Objectives::Linear,
						Data::IPG::Objectives::Feasibility
				  };

				  std::vector<Data::IPG::CutsAggressiveness> cuts = {
						Data::IPG::CutsAggressiveness::NoThanks,
						Data::IPG::CutsAggressiveness::KeepItCool,
						Data::IPG::CutsAggressiveness::Truculent
				  };


				  //InstanceName      = "Instance_2-100-8";

				  Models::IPG::IPGInstance Instance;
				  Instance.load(BasePath + InstanceName);
				  for (auto algo : algos) {
					 for (auto obj : objectives) {
						for (auto cut : cuts) {
						  std::cout << "\tTesting:" << std::to_string(algo) <<"-"<<std::to_string(obj)<<"-"<<std::to_string(cut)<<"\n";
						  auto realobj = obj;
						  if (algo == Data::LCP::Algorithms::PATH &&
						      obj == Data::IPG::Objectives::Quadratic)
							 realobj = Data::IPG::Objectives::Feasibility;
						  Models::IPG::IPG Test(&test, Instance);
						  Test.setAlgorithm(Data::IPG::Algorithms::Oracle);
						  Test.setDeviationTolerance(3e-4);
						  Test.setNumThreads(2);
						  Test.setCutsAggressiveness(cut);
						  Test.setTimeLimit(1000);
						  Test.setLCPAlgorithm(algo);
						  Test.setGameObjective(realobj);
						  Test.finalize();
						  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
						  Test.findNashEq();
						  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

						  // std::cout <<
						  // std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count()/1000000.0
						  // <<"s";

						  auto stat = Test.getStatistics();
						  double time = stat.WallClockTime.get();
						  std::cout << "\t\t Time:" <<std::to_string(time) << " Status:" << std::to_string(stat.Status.get())<<std::endl;

						  std::ifstream existCheck(Logfile);
						  std::ofstream results(Logfile, std::ios::app);

						  auto theCuts      = stat.AlgorithmData.Cuts.get();
						  int  totalCuts = 0;
						  for (const auto& thecut : theCuts)
							 totalCuts += thecut.second;

						  if (!existCheck.good()) {
							 results << "Instance;m;n;ins;Algorithm;LCPAlgo;ObjType;Cuts;Players;isPureNE;"
							            "Status;"
							            "NumVar;Time;Threads;numIterations;ValueCuts;VPolyCuts;"
							            "MIRCuts;GMICuts;KPCuts;TotalCuts;\n";
						  }
						  existCheck.close();


						  results << InstanceName << ";" << std::to_string(m) << ";" << std::to_string(n)
						          << ";" << std::to_string(ins) << ";"
						          << std::to_string(stat.AlgorithmData.Algorithm.get()) << ";"
						          << std::to_string(algo) << ";" << std::to_string(obj) << ";"
						          << std::to_string(cut) << ";" << Test.getNumPlayers() << ";"
						          << std::to_string(stat.PureNashEquilibrium.get()) << ";"
						          << std::to_string(stat.Status.get()) << ";" << stat.NumVar.get() << ";"
						          << std::to_string(stat.WallClockTime.get()) << ";"
						          << std::to_string(stat.AlgorithmData.Threads.get()) << ";"
						          << stat.NumIterations.get() << ";" << theCuts.at(0).second << ";"
						          << theCuts.at(1).second << ";" << theCuts.at(2).second << ";"
						          << theCuts.at(3).second << ";" << theCuts.at(4).second << ";" << totalCuts
						          << "\n";



						  results.close();
						} // algo
					 }   // cuts
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
