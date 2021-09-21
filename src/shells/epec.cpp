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


#include "boost/program_options.hpp"
#include "zero.h"


using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv) {
  string resFile, instanceFile, logFile;
  int    writeLevel = 0, nThreads = 0, verbosity = 0, algorithm = 0, aggressiveness = 0, add{0},
		recover = 0, branching = 0;
  double timeLimit = NAN, boundBigM = NAN, devtol = NAN;
  bool   bound = false, pure = false;

  po::options_description desc("ZERO-EPEC: Allowed options");
  desc.add_options()(
		"help,h", "Shows this help message")("version,v",
														 "Shows ZERO version")("input,i",
																					  po::value<string>(&instanceFile),
																					  "Sets the input path/filename of "
																					  "the instance file (.json "
																					  "appended "
																					  "automatically)")("pure,p",
																											  po::value<bool>(
																													&pure)
																													->default_value(
																														 false),
																											  "Controls "
																											  "whether the "
																											  "Algorithm "
																											  "should "
																											  "seek for a "
																											  "pure NE or "
																											  "not. If "
																											  "Algorithm is "
																											  "CombinatorialPN"
																											  "E, this is "
																											  "automatically "
																											  "true.")("recove"
																														  "r,r",
																														  po::value<
																																int>(
																																&recover)
																																->default_value(
																																	 0),
																														  "If "
																														  "InnerA"
																														  "pproxi"
																														  "mati"
																														  "on is "
																														  "used "
																														  "along "
																														  "with "
																														  "PureNa"
																														  "shEqui"
																														  "libr"
																														  "ium, "
																														  "which "
																														  "strate"
																														  "gy "
																														  "should"
																														  " "
																														  "be "
																														  "used "
																														  "to "
																														  "retrie"
																														  "ve a "
																														  "pure "
																														  "NE. "
																														  "0: "
																														  "Increm"
																														  "entalE"
																														  "nume"
																														  "ration"
																														  ", "
																														  "1:"
																														  "Combin"
																														  "atoria"
																														  "lPN"
																														  "E")("A"
																																 "l"
																																 "g"
																																 "o"
																																 "r"
																																 "i"
																																 "t"
																																 "h"
																																 "m"
																																 ","
																																 "a",
																																 po::value<
																																	  int>(
																																	  &algorithm),
																																 "S"
																																 "e"
																																 "t"
																																 "s"
																																 " "
																																 "t"
																																 "h"
																																 "e"
																																 " "
																																 "A"
																																 "l"
																																 "g"
																																 "o"
																																 "r"
																																 "i"
																																 "t"
																																 "h"
																																 "m"
																																 "."
																																 " "
																																 "0"
																																 ":"
																																 "F"
																																 "u"
																																 "l"
																																 "l"
																																 "E"
																																 "n"
																																 "u"
																																 "m"
																																 "e"
																																 "r"
																																 "a"
																																 "t"
																																 "i"
																																 "o"
																																 "n"
																																 ","
																																 " "
																																 "1"
																																 ":"
																																 "I"
																																 "n"
																																 "n"
																																 "e"
																																 "r"
																																 "A"
																																 "p"
																																 "p"
																																 "r"
																																 "o"
																																 "x"
																																 "i"
																																 "m"
																																 "a"
																																 "t"
																																 "i"
																																 "o"
																																 "n"
																																 ","
																																 " "
																																 "2"
																																 ":"
																																 "C"
																																 "o"
																																 "m"
																																 "b"
																																 "i"
																																 "n"
																																 "a"
																																 "t"
																																 "o"
																																 "r"
																																 "i"
																																 "a"
																																 "l"
																																 "P"
																																 "N"
																																 "E"
																																 ","
																																 " "
																																 "3"
																																 ":"
																																 "O"
																																 "u"
																																 "t"
																																 "e"
																																 "r"
																																 "A"
																																 "p"
																																 "p"
																																 "r"
																																 "o"
																																 "x"
																																 "i"
																																 "m"
																																 "a"
																																 "t"
																																 "i"
																																 "o"
																																 "n")("solution,s",
																																		po::value<
																																			 string>(
																																			 &resFile)
																																			 ->default_value(
																																				  "dat/Solution"),
																																		"Sets the output path/filename of the solution file (.json appended "
																																		"automatically)")("log,l",
																																								po::value<
																																									 string>(
																																									 &logFile)
																																									 ->default_value(
																																										  "dat/Results.csv"),
																																								"Sets the output path/filename of the csv log file")("timelimit,tl",
																																																									  po::value<
																																																											double>(
																																																											&timeLimit)
																																																											->default_value(
																																																												 -1.0),
																																																									  "Sets the timelimit for solving the Nash Equilibrium model")("writelevel,w",
																																																																														po::value<
																																																																															 int>(
																																																																															 &writeLevel)
																																																																															 ->default_value(
																																																																																  0),
																																																																														"Sets the writeLevel param. 0: only Json. 1: only human-readable. 2: "
																																																																														"both")("message,m",
																																																																																  po::value<
																																																																																		int>(
																																																																																		&verbosity)
																																																																																		->default_value(
																																																																																			 0),
																																																																																  "Sets the verbosity level for info and warning messages. 0: "
																																																																																  "warning and critical. 1: info. 2: debug. 3: trace")("Threads,t",
																																																																																																		 po::value<
																																																																																																			  int>(
																																																																																																			  &nThreads)
																																																																																																			  ->default_value(
																																																																																																					1),
																																																																																																		 "Sets the number of Threads for Gurobi. (int): number of Threads. 0: "
																																																																																																		 "auto (number of processors)")("aggr,ag",
																																																																																																												  po::value<
																																																																																																														int>(
																																																																																																														&aggressiveness)
																																																																																																														->default_value(
																																																																																																															 1),
																																																																																																												  "Sets the Aggressiveness for the InnerApproximation, namely the number "
																																																																																																												  "of Random polyhedra added if no deviation is found. (int)")("bound,bo",
																																																																																																																																	po::value<
																																																																																																																																		 bool>(
																																																																																																																																		 &bound)
																																																																																																																																		 ->default_value(
																																																																																																																																			  false),
																																																																																																																																	"Decides whether primal variables should be bounded or not.")("devtol,dt",
																																																																																																																																																					  po::value<
																																																																																																																																																							double>(
																																																																																																																																																							&devtol)
																																																																																																																																																							->default_value(
																																																																																																																																																								 -1.0),
																																																																																																																																																					  "Sets the deviation tolerance.")("BoundBigM,bbm",
																																																																																																																																																																  po::value<
																																																																																																																																																																		double>(
																																																																																																																																																																		&boundBigM)
																																																																																																																																																																		->default_value(
																																																																																																																																																																			 1e5),
																																																																																																																																																																  "Set the bounding BigM related to the parameter --bound")("add,ad",
																																																																																																																																																																																				po::
																																																																																																																																																																																					 value<
																																																																																																																																																																																						  int>(
																																																																																																																																																																																						  &add)
																																																																																																																																																																																						  ->default_value(0),
																																																																																																																																																																																				"Sets the Game::EPECAddPolyMethod for the InnerApproximation. 0: "
																																																																																																																																																																																				"Sequential. "
																																																																																																																																																																																				"1: ReverseSequential. 2:Random.")("branching,br",
																																																																																																																																																																																															  po::value<
																																																																																																																																																																																																	int>(
																																																																																																																																																																																																	&branching)
																																																																																																																																																																																																	->default_value(
																																																																																																																																																																																																		 0),
																																																																																																																																																																																															  "Sets the Branching Strategy for the CutAndPlay. 0: "
																																																																																																																																																																																															  "Hybrid. "
																																																																																																																																																																																															  "1: Deviation.");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
	 cout << desc;
	 return EXIT_SUCCESS;
  }
  if (vm.count("version") || verbosity >= 2) {
	 string M, m, p;
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
  // LOADING INSTANCE
  // --------------------------------
  Models::EPEC::EPECInstance instance(instanceFile);
  if (instance.Countries.empty()) {
	 cerr << "Error: instance is empty\n";
	 return 1;
  }

  // --------------------------------
  // TEST STARTS
  // --------------------------------
  auto timeStart = std::chrono::high_resolution_clock::now();
  try {
	 GRBEnv env = GRBEnv();
	 //As to have consistency among results
	 env.set(GRB_IntParam_Seed,420);

	 // OPTIONS
	 //------------
	 Models::EPEC::EPEC epec(&env);
	 // Num Threads
	 if (nThreads != 0)
		epec.setNumThreads(nThreads);
	 // Pure NE
	 if (pure)
		epec.setPureNashEquilibrium(true);
	 // TimeLimit
	 epec.setTimeLimit(timeLimit);
	 // bound QPs
	 if (devtol > 0)
		epec.setDeviationTolerance(devtol);

	 // Algorithm

	 switch (algorithm) {
	 case 1: {
		epec.setAlgorithm(Data::EPEC::Algorithms::InnerApproximation);
		if (aggressiveness != 1)
		  epec.setAggressiveness(aggressiveness);
		switch (add) {
		case 1:
		  epec.setAddPolyMethod(Data::LCP::PolyhedraStrategy::ReverseSequential);
		  break;
		case 2:
		  epec.setAddPolyMethod(Data::LCP::PolyhedraStrategy::Random);
		  break;
		default:
		  epec.setAddPolyMethod(Data::LCP::PolyhedraStrategy::Sequential);
		}
		if (recover != 0)
		  epec.setRecoverStrategy(Data::EPEC::RecoverStrategy::Combinatorial);
		break;
	 }
	 case 2: {
		epec.setAlgorithm(Data::EPEC::Algorithms::CombinatorialPne);
		break;
	 }
	 case 3: {
		epec.setAlgorithm(Data::EPEC::Algorithms::OuterApproximation);
		switch (branching) {
		case 0:
		  epec.setBranchingStrategy(Data::EPEC::BranchingStrategy::HybridBranching);
		  break;
		default:
		  epec.setBranchingStrategy(Data::EPEC::BranchingStrategy::DeviationBranching);
		}

		break;
	 }
	 default:
		epec.setAlgorithm(Data::EPEC::Algorithms::FullEnumeration);
	 }

	 //------------

	 for (auto &Countrie : instance.Countries)
		epec.addCountry(Countrie);
	 epec.addTranspCosts(instance.TransportationCosts);
	 epec.finalize();
	 epec.findNashEq();

	 auto                          timeStop      = std::chrono::high_resolution_clock::now();
	 std::chrono::duration<double> timeDiff      = timeStop - timeStart;
	 double                        wallClockTime = timeDiff.count();
	 int realThreads = nThreads > 0 ? env.get(GRB_IntParam_Threads) : nThreads;

	 // --------------------------------
	 // WRITING STATISTICS AND SOLUTION
	 // --------------------------------
	 auto stat = epec.getStatistics();
	 if (stat.Status.get() == ZEROStatus::NashEqFound)
		epec.writeSolution(writeLevel, resFile);
	 ifstream      existCheck(logFile);
	 std::ofstream results(logFile, ios::app);

	 if (!existCheck.good()) {
		results << "instance;Algorithm;Countries;Followers;isPureNE;RequiredPureNE;"
					  "Status;"
					  "numFeasiblePolyhedra/Complementarities;"
					  "NumVar;NumConstraints;NumNonZero;ClockTime"
					  "(s);Threads;numInnerIterations;LostIntermediateEq;"
					  "Aggressiveness;"
					  "AddPolyMethod;NumericalIssues;"
					  "recoveryStrategy;branchingStrategy\n";
	 }
	 existCheck.close();

	 std::stringstream         polyT;
	 std::vector<unsigned int> target;
	 if (stat.AlgorithmData.Algorithm.get() != Data::EPEC::Algorithms::OuterApproximation)
		target = stat.AlgorithmData.FeasiblePolyhedra.get();
	 else
		target = stat.AlgorithmData.OuterComplementarities.get();

	 copy(target.begin(), target.end(), ostream_iterator<int>(polyT, " "));
	 auto test2 = polyT.str();


	 results << instanceFile << ";" << to_string(stat.AlgorithmData.Algorithm.get()) << ";"
				<< instance.Countries.size() << ";[";
	 for (auto &countrie : instance.Countries)
		results << " " << countrie.n_followers;

	 results << " ];" << to_string(epec.getStatistics().PureNashEquilibrium.get()) << ";"
				<< to_string(pure) << ";" << to_string(stat.Status.get()) << ";[ " << polyT.str()
				<< "];" << stat.NumVar.get() << ";" << stat.NumConstraints.get() << ";"
				<< stat.NumNonZero.get() << ";" << wallClockTime << ";" << realThreads;
	 if (stat.AlgorithmData.Algorithm.get() == Data::EPEC::Algorithms::InnerApproximation) {
		results << ";" << stat.NumIterations.get() << ";"
				  << epec.getStatistics().AlgorithmData.LostIntermediateEq.get() << ";"
				  << stat.AlgorithmData.Aggressiveness.get() << ";"
				  << to_string(stat.AlgorithmData.PolyhedraStrategy.get()) << ";"
				  << stat.NumericalIssues.get() << ";"
				  << to_string(stat.AlgorithmData.RecoverStrategy.get()) << ";-";
	 } else {
		results << ";-;-;-;-;-;-;" << std::to_string(stat.AlgorithmData.BranchingStrategy.get());
	 }
	 results << "\n";
	 results.close();
  } catch (ZEROException &e) {
	 std::cerr << "" << e.what() << "--" << e.more();
	 throw ZEROException(e);
  } catch (GRBException &e) {
	 std::cerr << "" << e.getErrorCode() << "--" << e.getMessage();
	 throw;
  }

  return EXIT_SUCCESS;
}
