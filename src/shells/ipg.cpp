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


#include "interfaces/ipg_models.h"
#include "zero.h"
#include <armadillo>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/program_options.hpp>
#include <chrono>
#include <cstdlib>
#include <gurobi_c++.h>
#include <iostream>
#include <iterator>
#include <math.h>

using namespace std;
namespace logging = boost::log;
using namespace boost::program_options;
namespace po = boost::program_options;

int main(int argc, char **argv) {
  string resFile, instanceFile, logFile;
  int    nThreads = 0, verbosity = 0, algorithm = 0;
  double timeLimit = NAN, devtol = NAN;
  bool   pure = 0;

  po::options_description desc("ZERO-IPG: Allowed options");
  desc.add_options()("help,h", "Shows this help message")("version,v", "Shows ZERO version")(
		"input,i",
		po::value<string>(&instanceFile),
		"Sets the input path/filename of the instance file (.json appended "
		"automatically)")("pure,p",
								po::value<bool>(&pure)->default_value(false),
								"Controls whether the Algorithm should seek for a pure NE or not")(
		"Algorithm,a", po::value<int>(&algorithm), "Sets the Algorithm. 0:Oracle")(
		"solution,s",
		po::value<string>(&resFile)->default_value("dat/Solution"),
		"Sets the output path/filename of the solution file (.json appended "
		"automatically)")("log,l",
								po::value<string>(&logFile)->default_value("dat/Results.csv"),
								"Sets the output path/filename of the csv log file")(
		"timelimit,tl", po::value<double>(&timeLimit)->default_value(-1.0), "Sets the timelimit")(
		"message,m",
		po::value<int>(&verbosity)->default_value(0),
		"Sets the verbosity level for info and warning messages. 0: "
		"warning and critical. 1: info. 2: debug. 3: trace")(
		"Threads,t",
		po::value<int>(&nThreads)->default_value(1),
		"Sets the number of Threads for Gurobi. (int): number of Threads. 0: "
		"auto (number of processors)")("devtol,dt",
												 po::value<double>(&devtol)->default_value(-1.0),
												 "Sets the deviation tolerance.");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
	 cout << desc;
	 return EXIT_SUCCESS;
  }
  if (vm.count("version") || verbosity >= 2) {
	 arma::arma_version ver;
	 int                major = 0, minor = 0, technical = 0;
	 string             M, m, p;
	 ZEROVersion(M, m, p);
	 BOOST_LOG_TRIVIAL(info) << "ZERO Version: " << M << "." << m << "." << p;
	 if (vm.count("version"))
		return EXIT_SUCCESS;
  }

  if (instanceFile == "") {
	 cout << "-i [--input] option missing.\n Use with --help for help on list "
				"of arguments\n";
	 return EXIT_SUCCESS;
  }
  switch (verbosity) {
  case 0:
	 logging::core::get()->set_filter(logging::trivial::severity > logging::trivial::info);
	 break;
  case 1:
	 logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::info);
	 break;
  case 2:
	 logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::debug);
	 break;
  case 3:
	 logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::trace);
	 break;
  default:
	 BOOST_LOG_TRIVIAL(warning) << "Invalid option for --message (-m). Setting default value: 0";
	 verbosity = 0;
	 logging::core::get()->set_filter(logging::trivial::severity > logging::trivial::info);
	 break;
  }

  // --------------------------------
  // START
  // --------------------------------
  auto timeStart = std::chrono::high_resolution_clock::now();
  try {
	 GRBEnv env = GRBEnv();

	 // OPTIONS
	 //------------
	 Models::IPG::IPG ipg(&env, instanceFile);
	 // Num Threads
	 if (nThreads != 0)
		ipg.setNumThreads(nThreads);
	 // Pure NE
	 if (pure)
		ipg.setPureNashEquilibrium(true);
	 // TimeLimit
	 ipg.setTimeLimit(timeLimit);
	 if (devtol > 0)
		ipg.setDeviationTolerance(devtol);

	 // Algorithm

	 switch (algorithm) {
	 default:
		ipg.setAlgorithm(Data::IPG::Algorithms::Oracle);
	 }

	 ipg.findNashEq();

	 auto                          timeStop      = std::chrono::high_resolution_clock::now();
	 std::chrono::duration<double> timeDiff      = timeStop - timeStart;
	 double                        wallClockTime = timeDiff.count();
	 int realThreads = nThreads > 0 ? env.get(GRB_IntParam_Threads) : nThreads;

	 // --------------------------------
	 // WRITING STATISTICS AND SOLUTION
	 // --------------------------------
	 auto          stat = ipg.getStatistics();
	 ifstream      existCheck(logFile);
	 std::ofstream results(logFile, ios::app);

	 if (!existCheck.good()) {
		results << "instance;Algorithm;Players;isPureNE;RequiredPureNE;"
					  "Status;"
					  "NumVar;ClockTime(s);Threads;Indicators;numInnerIterations\n";
	 }
	 existCheck.close();


	 results << instanceFile << ";" << to_string(stat.AlgorithmData.Algorithm.get()) << ";"
				<< ipg.getNumPlayers() << ";"
				<< to_string(ipg.getStatistics().PureNashEquilibrium.get()) << ";" << to_string(pure)
				<< ";" << to_string(stat.Status.get()) << ";" << stat.NumVar.get() << ";"
				<< wallClockTime << ";" << realThreads << ";" << stat.NumIterations.get() << "\n";

	 results.close();
  } catch (ZEROException &e) {
	 std::cerr << "" << e.what() << "--" << e.more();
  }

  return EXIT_SUCCESS;
}
