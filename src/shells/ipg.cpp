/* #############################################
 *             This file is part of
 *                    ZERO
 *
 *             Copyright (c) 2020
 *     Released under the Creative Commons
 *         CC BY-NC-SA 4.0 License
 *
 *              Find out more at
 *        https://github.com/ds4dm/ZERO
 * #############################################*/


#include "zero.h"
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv) {
  string resFile, instanceFile, logFile;
  int    nThreads = 0, verbosity = 0, algorithm = 0, aggressiveness = 0, objective = 0, LCPalgo = 0;
  double timeLimit = -1, devtol = -1;

  po::options_description desc("ZERO-IPG: Allowed options");
  desc
		.add_options()(
			 "help,h", "Shows this help message")("version,v",
															  "Shows ZERO version")("input,i",
																							po::value<string>(
																								 &instanceFile),
																							"Sets the input "
																							"path/filename of the "
																							"instance file (.json "
																							"appended "
																							"automatically)")("algorithm"
																													"m,a",
																													po::value<
																														 int>(
																														 &algorithm),
																													"Sets "
																													"the "
																													"Algorithm"
																													"m. "
																													"0:"
																													"CutAndPlay")("s"
																																 "o"
																																 "l"
																																 "u"
																																 "t"
																																 "i"
																																 "o"
																																 "n"
																																 ","
																																 "s",
																																 po::value<
																																	  string>(
																																	  &resFile)
																																	  ->default_value(
																																			"dat/Solution"),
																																 "S"
																																 "e"
																																 "t"
																																 "s"
																																 " "
																																 "t"
																																 "h"
																																 "e"
																																 " "
																																 "o"
																																 "u"
																																 "t"
																																 "p"
																																 "u"
																																 "t"
																																 " "
																																 "p"
																																 "a"
																																 "t"
																																 "h"
																																 "/"
																																 "f"
																																 "i"
																																 "l"
																																 "e"
																																 "n"
																																 "a"
																																 "m"
																																 "e"
																																 " "
																																 "o"
																																 "f"
																																 " "
																																 "t"
																																 "h"
																																 "e"
																																 " "
																																 "s"
																																 "o"
																																 "l"
																																 "u"
																																 "t"
																																 "i"
																																 "o"
																																 "n"
																																 " "
																																 "f"
																																 "i"
																																 "l"
																																 "e"
																																 " "
																																 "("
																																 "."
																																 "j"
																																 "s"
																																 "o"
																																 "n"
																																 " "
																																 "a"
																																 "p"
																																 "p"
																																 "e"
																																 "n"
																																 "d"
																																 "e"
																																 "d"
																																 " "
																																 "a"
																																 "u"
																																 "t"
																																 "o"
																																 "m"
																																 "a"
																																 "t"
																																 "i"
																																 "c"
																																 "a"
																																 "l"
																																 "l"
																																 "y"
																																 ")")("log,l",
																																		po::value<
																																			 string>(
																																			 &logFile)
																																			 ->default_value(
																																				  "dat/"
																																				  "Results."
																																				  "csv"),
																																		"Sets "
																																		"the "
																																		"output "
																																		"path/"
																																		"filenam"
																																		"e of "
																																		"the "
																																		"csv "
																																		"log "
																																		"file")("timelimit,"
																																				  "tl",
																																				  po::value<
																																						double>(
																																						&timeLimit)
																																						->default_value(
																																							 -1.0),
																																				  "Sets the "
																																				  "timelimit")("verbosity,ve",
																																									po::value<
																																										 int>(
																																										 &verbosity)
																																										 ->default_value(
																																											  0),
																																									"Sets the verbosity level for info and warning messages. 0: "
																																									"warning and critical. 1: info. 2: debug. 3: trace")("threads,t",
																																																										  po::value<
																																																												int>(
																																																												&nThreads)
																																																												->default_value(
																																																													 1),
																																																										  "Sets the number of Threads for Gurobi. (int): number of Threads. 0: "
																																																										  "auto (number of processors)")("devtol,dt",
																																																																					po::value<
																																																																						 double>(
																																																																						 &devtol)
																																																																						 ->default_value(
																																																																							  3e-4),
																																																																					"Sets the deviation tolerance.")("aggressiveness,aggr",
																																																																																po::value<
																																																																																	 int>(
																																																																																	 &aggressiveness)
																																																																																	 ->default_value(
																																																																																		  0),
																																																																																"Cutting planes aggressiveness (int): -1: NotEvenTry; 0: NoThanks; 1: KeepItCool; 2: Truculent")("objective,o",
																																																																																																																 po::value<
																																																																																																																	  int>(&objective)
																																																																																																																	  ->default_value(
																																																																																																																			0),
																																																																																																																 "LCP Objective type (int): 0: Feasibility; 1: Linear; 2: Quadratic")("lcpalgo,l",
																																																																																																																																							 po::value<
																																																																																																																																								  int>(
																																																																																																																																								  &LCPalgo)
																																																																																																																																								  ->default_value(
																																																																																																																																										0),
																																																																																																																																							 "LCP Algorithm type (int): 0: MIP; 1: MINLP; 2: PATH");

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
	 LOG_S(INFO) << "ZERO Version: " << M << "." << m << "." << p;
	 if (vm.count("version"))
		return EXIT_SUCCESS;
  }

  if (instanceFile.empty()) {
	 cout << "-i [--input] option missing.\n Use with --help for help on list "
				"of arguments\n";
	 return EXIT_SUCCESS;
  }

  loguru::g_stderr_verbosity = verbosity;


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
	 // TimeLimit
	 ipg.setTimeLimit(timeLimit);
	 if (devtol > 0)
		ipg.setDeviationTolerance(devtol);

	 // Algorithm
	 // So far, just one...
	 ipg.setAlgorithm(Data::IPG::Algorithms::CutAndPlay);
	 switch (aggressiveness) {
	 case -1:
		ipg.setCutsAggressiveness(Data::IPG::CutsAggressiveness::NotEvenTry);
		break;
	 case 0:
		ipg.setCutsAggressiveness(Data::IPG::CutsAggressiveness::NoThanks);
		break;
	 case 1:
		ipg.setCutsAggressiveness(Data::IPG::CutsAggressiveness::KeepItCool);
		break;
	 default:
		ipg.setCutsAggressiveness(Data::IPG::CutsAggressiveness::Truculent);
	 }

	 switch (LCPalgo) {
	 case 0:
		ipg.setLCPAlgorithm(Data::LCP::Algorithms::MIP);
		break;
	 case 1:
		ipg.setLCPAlgorithm(Data::LCP::Algorithms::MINLP);
		break;
	 default:
		ipg.setLCPAlgorithm(Data::LCP::Algorithms::PATH);
	 }

	 if (LCPalgo != 2) {
		switch (objective) {
		case 0:
		  ipg.setGameObjective(Data::IPG::Objectives::Feasibility);
		  break;
		case 1:
		  ipg.setGameObjective(Data::IPG::Objectives::Linear);
		  break;
		default:
		  ipg.setGameObjective(Data::IPG::Objectives::Quadratic);
		}
	 } else
		ipg.setGameObjective(Data::IPG::Objectives::Feasibility);


	 ipg.findNashEq();

	 auto                          timeStop      = std::chrono::high_resolution_clock::now();
	 std::chrono::duration<double> timeDiff      = timeStop - timeStart;
	 double                        wallClockTime = timeDiff.count();
	 int realThreads = nThreads > 0 ? env.get(GRB_IntParam_Threads) : nThreads;

	 // --------------------------------
	 // WRITING STATISTICS AND SOLUTION
	 // --------------------------------
	 auto stat = ipg.getStatistics();
	 if (!logFile.empty()) {
		ifstream      existCheck(logFile);
		std::ofstream results(logFile, ios::app);

		if (!existCheck.good()) {
		  results << "Instance;Algorithm;LCPAlgo;ObjType;Cuts;Players;isPureNE;"
						 "Status;"
						 "NumVar;Time;SocialWelfare;Threads;numIterations;ValueCuts;VPolyCuts;"
						 "MIPCuts;TotalCuts\n";
		}
		existCheck.close();

		auto cuts      = stat.AlgorithmData.Cuts.get();
		int  totalCuts = 0;
		for (const auto &cut : cuts)
		  totalCuts += cut.second;
		results << instanceFile << ";" << std::to_string(stat.AlgorithmData.Algorithm.get()) << ";"
				  << std::to_string(stat.AlgorithmData.LCPSolver.get()) << ";"
				  << std::to_string(stat.AlgorithmData.Objective.get()) << ";"
				  << std::to_string(stat.AlgorithmData.CutAggressiveness.get()) << ";"
				  << ipg.getNumPlayers() << ";" << std::to_string(stat.PureNashEquilibrium.get()) << ";"
				  << std::to_string(stat.Status.get()) << ";" << stat.NumVar.get() << ";"
				  << std::to_string(stat.WallClockTime.get()) << ";"
				  << std::to_string(ipg.getSocialWelfare()) << ";"
				  << std::to_string(stat.AlgorithmData.Threads.get()) << ";" << stat.NumIterations.get()
				  << ";" << cuts.at(0).second << ";" << cuts.at(1).second << ";" << cuts.at(2).second
				  << ";" << totalCuts << "\n";


		results.close();
	 }
  } catch (ZEROException &e) {
	 std::cerr << "" << e.what() << "--" << e.more();
  }

  return EXIT_SUCCESS;
}
