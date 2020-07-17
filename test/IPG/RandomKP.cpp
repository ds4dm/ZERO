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
#include <interfaces/epec_models.h>
#include <interfaces/ipg_models.h>
#include <iostream>
#include <string>


int main(int argc, char **argv) {


  std::vector<unsigned int>          ItemSizes = {20, 40, 80, 100};
  std::vector<unsigned int>          Players   = {2, 3, 4};
  unsigned int                       insUB     = 10;
  bool                               print     = false;
  bool                               generate  = false;
  GRBEnv                             test;
  std::string                        BasePath = "dat/KP/";
  std::string                        Logfile  = "dat/KP/results.csv";
  std::random_device                 random;
  std::mt19937                       Oracle(random());
  std::uniform_int_distribution<int> hRandom(-100, 100);
  std::uniform_int_distribution<int> fRandom(-50, 50);

  for (unsigned int n : ItemSizes) {
	 for (unsigned int ins = 0; ins < insUB; ++ins) {
		for (auto &&m : Players) {

		  std::string InstanceName =
				"Instance_" + std::to_string(m) + "-" + std::to_string(n) + "-" + std::to_string(ins);
		  if (generate) {
			 Models::IPG::IPGInstance Instance;
			 for (unsigned int j = 0; j < m; ++j) {
				arma::vec _c(n);
				arma::mat _Q(n, n * (m - 1));
				arma::mat _A(1, n);

				// Bounds on integer variables to one
				arma::mat _Abounds(n, n, arma::fill::eye);
				arma::vec _ones(n, arma::fill::ones);
				arma::vec _integers(n);

				for (unsigned int i = 0; i < n; ++i)
				  _integers.at(i) = i;

				for (auto &element : _Q)
				  element = hRandom(Oracle);

				for (auto &element : _A)
				  element = hRandom(Oracle);

				for (auto &element : _c) {
				  element = fRandom(Oracle);
				}

				double b = arma::as_scalar(arma::accu(_A)) * ins / 11;

				// Explicitly join the upper bounds
				_A           = arma::join_cols(_A, _Abounds);
				arma::vec _b = arma::join_cols(arma::vec{b}, _ones);

				if (print) {
				  _c.print("c:");
				  _b.print("b:");
				  _Q.print("Q");
				  _A.print("A");
				}

				arma::sp_mat      _Q2 = arma::sp_mat{_Q};
				arma::sp_mat      _A2 = arma::sp_mat{_A};
				MathOpt::IP_Param ipParam(_Q2, _A2, _b, _c, _integers, &test);
				Instance.addIPParam(ipParam,
										  BasePath + "raw/" + InstanceName + "-p" + std::to_string(j));
			 }
			 Instance.save(BasePath + InstanceName);
			 std::cout << "\tSaving instance " << InstanceName << "\n";
			 // Save the instance
		  } else {
			 try {
				Models::IPG::IPG Test(&test, BasePath + InstanceName);
				Test.setAlgorithm(Data::IPG::Algorithms::Oracle);
				Test.setDeviationTolerance(3e-4);
				Test.setIndicators(true);
				Test.setNumThreads(8);
				Test.finalize();
				Test.findNashEq();

				auto          stat = Test.getStatistics();
				std::ifstream existCheck(Logfile);
				std::ofstream results(Logfile, std::ios::app);

				if (!existCheck.good()) {
				  results << "instance;Algorithm;Players;isPureNE;RequiredPureNE;"
								 "Status;"
								 "NumVar;ClockTime(s);Threads;Indicators;numInnerIterations\n";
				}
				existCheck.close();


				results << InstanceName << ";" << std::to_string(stat.AlgorithmData.Algorithm.get())
						  << ";" << Test.getNumPlayers() << ";"
						  << std::to_string(stat.PureNashEquilibrium.get()) << ";"
						  << "0"
						  << ";" << std::to_string(stat.Status.get()) << ";" << stat.NumVar.get()
						  << ";inf;" << std::to_string(stat.AlgorithmData.Threads.get()) << ";"
						  << stat.NumIterations.get() << "\n";

				results.close();
			 } catch (ZEROException e) {
				std::cout << "ZEROException: " << std::to_string(e.which()) << " | " << e.more();
			 }
		  }
		}
	 }
  }
}
