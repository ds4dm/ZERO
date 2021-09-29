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

#include "include/settings.h"
#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(IPG_Tests)

/**
 * @brief A function creating an IP_Param with 2 players, as in examples/IPG.cpp
 * @param one Whether it is the first or second player
 * @param e The gurobi environment
 * @return A new IP_Param
 */
MathOpt::IP_Param Player(bool one, GRBEnv *e) {
  arma::vec      c(2);                         // Profits c in the objective
  arma::sp_mat   C(2, 2);                      // C terms in the objective
  arma::sp_mat   a(1, 2);                      // LHS for Knapsack constraint
  arma::vec      b(1);                         // RHS for constraints
  arma::vec      IntegerIndexes(2);            // The index of the integer variables
  VariableBounds VarBounds = {{0, 1}, {0, 1}}; // Implicit bounds on variables

  // Fill the values in the paramterized integer problem
  for (unsigned int i = 0; i < 2; ++i)
	 IntegerIndexes.at(i) = i;
  b(0) = 5; // Knapsack Capacity
  if (one) {
	 C(0, 0) = 2; // C terms in the objective
	 C(1, 1) = 3;
	 a(0, 0) = 3; // Knapsack Constraints
	 a(0, 1) = 4;
	 c(0)    = -1; // The standard is minimization, hence minus
	 c(1)    = -2;
  } else {
	 C(0, 0) = 5;
	 C(1, 1) = 4;
	 a(0, 0) = 2;
	 a(0, 1) = 5;
	 c(0)    = -3;
	 c(1)    = -5;
  }


  return MathOpt::IP_Param(C, a, b, c, IntegerIndexes, VarBounds, e);
}

BOOST_AUTO_TEST_CASE(IPGModelInstance_test) {
  BOOST_TEST_MESSAGE("\n\n");
  BOOST_TEST_MESSAGE("Testing Models::IPG::IPG");

  GRBEnv env       = GRBEnv();
  auto   PlayerOne = Player(true, &env);
  auto   PlayerTwo = Player(false, &env);
  // Constructor
  BOOST_TEST_MESSAGE("Constructor tests");
  Models::IPG::IPGInstance IPG_Instance1;
  BOOST_CHECK_NO_THROW(IPG_Instance1.addIPParam(PlayerOne, "../test/run_data/PlayerOne"));
  BOOST_CHECK_NO_THROW(IPG_Instance1.addIPParam(PlayerTwo, "../test/run_data/PlayerTwo"));

  BOOST_TEST_MESSAGE("Write/save tests");
  BOOST_CHECK_NO_THROW(IPG_Instance1.save("../test/run_data/IPG1"));
  Models::IPG::IPGInstance IPG_Instance2;
  BOOST_CHECK_NO_THROW(IPG_Instance2.load("../test/run_data/IPG1"));

  BOOST_TEST_MESSAGE("Loading into a game");
  Models::IPG::IPG IPG1(&env, IPG_Instance1);
  Models::IPG::IPG IPG2(&env, IPG_Instance2);


  BOOST_CHECK_MESSAGE(*IPG1.getIPParam(0) == *IPG2.getIPParam(0),
							 "IPG1::IP_Param1 is equal to IPG2::IP_Param1");
  BOOST_CHECK_MESSAGE(*IPG1.getIPParam(1) == *IPG2.getIPParam(1),
							 "IPG1::IP_Param2 is equal to IPG2::IP_Param2");
}

BOOST_AUTO_TEST_CASE(IPGRandomKP_Test) {
  BOOST_TEST_MESSAGE("\n\n");
  BOOST_TEST_MESSAGE("Testing RandomKP Game");

  GRBEnv env       = GRBEnv();
  auto   PlayerOne = Player(true, &env);
  auto   PlayerTwo = Player(false, &env);
  // Constructor
  BOOST_TEST_MESSAGE("Constructor tests");
  Models::IPG::IPGInstance IPG_Instance1;
  BOOST_CHECK_NO_THROW(IPG_Instance1.addIPParam(PlayerOne, "../test/run_data/PlayerOne"));
  BOOST_CHECK_NO_THROW(IPG_Instance1.addIPParam(PlayerTwo, "../test/run_data/PlayerTwo"));

  std::vector<arma::vec> solutions;
  solutions.push_back(arma::vec{0, 1, 1, 0});
  solutions.push_back(arma::vec{1, 0, 0, 1});
  solutions.push_back(arma::vec{2.0 / 9, 7.0 / 9, 2.0 / 5, 3.0 / 5});


  std::vector<Data::LCP::Algorithms> algos = {
		Data::LCP::Algorithms::MIP, Data::LCP::Algorithms::MINLP, Data::LCP::Algorithms::PATH};

  std::vector<Data::IPG::Objectives> objectives = {Data::IPG::Objectives::Quadratic,
																	Data::IPG::Objectives::Linear,
																	Data::IPG::Objectives::Feasibility};

  std::vector<Data::IPG::CutsAggressiveness> cuts = {Data::IPG::CutsAggressiveness::NotEvenTry,
																	  Data::IPG::CutsAggressiveness::NoThanks,
																	  Data::IPG::CutsAggressiveness::KeepItCool,
																	  Data::IPG::CutsAggressiveness::Truculent};

  for (auto algo : algos) {
	 for (auto obj : objectives) {
		for (auto cut : cuts) {
		  if (!((obj == Data::IPG::Objectives::Quadratic || obj == Data::IPG::Objectives::Linear) &&
				  algo == Data::LCP::Algorithms::PATH)) {
			 BOOST_TEST_MESSAGE(" Testing: Algorithm " << std::to_string(algo) << " with objective "
																	 << std::to_string(obj) << " and Cuts"
																	 << std::to_string(cut));
			 Models::IPG::IPG Test(&env, IPG_Instance1);
			 BOOST_CHECK_NO_THROW(Test.setAlgorithm(Data::IPG::Algorithms::CutAndPlay));
			 BOOST_CHECK_NO_THROW(Test.setDeviationTolerance(3e-4));
			 BOOST_CHECK_NO_THROW(Test.setNumThreads(TEST_NUM_THREADS));
			 BOOST_CHECK_NO_THROW(Test.setCutsAggressiveness(cut));
			 BOOST_CHECK_NO_THROW(Test.setTimeLimit(10));
			 BOOST_CHECK_NO_THROW(Test.setLCPAlgorithm(algo));
			 BOOST_CHECK_NO_THROW(Test.setGameObjective(obj));
			 BOOST_CHECK_NO_THROW(Test.finalize());
			 BOOST_CHECK_NO_THROW(Test.findNashEq());
			 BOOST_CHECK_MESSAGE(Test.isSolved(),
										"Invoking isSolved method for " + std::to_string(algo));

			 bool found = false;
			 auto x1    = Test.getX().at(0);
			 auto x2    = Test.getX().at(1);
			 for (unsigned int s = 0; s < solutions.size(); ++s) {
				if (Utils::isEqual(x1.at(0), solutions.at(s).at(0)) &&
					 Utils::isEqual(x1.at(1), solutions.at(s).at(1)) &&
					 Utils::isEqual(x2.at(0), solutions.at(s).at(2)) &&
					 Utils::isEqual(x2.at(1), solutions.at(s).at(3)))
				  found = true;
			 }
			 BOOST_CHECK_MESSAGE(found, "The instance was solved correctly");


		  } // endif objectives and PATH
		}   // objectives
	 }     // algo
  }       // cuts
}

BOOST_AUTO_TEST_SUITE_END()
