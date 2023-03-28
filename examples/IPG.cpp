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


#include <zero.h>

int main(int argc, char **argv) {

  loguru::g_stderr_verbosity = 10;
  GRBEnv GurobiEnv;
  try {
	 Models::IPG::IPGInstance IPG_Instance; // The IPG Instance
	 int                      numItems = 2, numPlayers = 2;

	 arma::vec      c(numItems);                              // Profits c in the objective
	 arma::sp_mat   C(numItems * (numPlayers - 1), numItems); // C terms in the objective
	 arma::sp_mat   a(1, numItems);                           // LHS for Knapsack constraint
	 arma::vec      b(1), d(2,arma::fill::zeros);                                     // RHS for constraints
	 arma::vec      IntegerIndexes(numItems);                 // The index of the integer variables
	 VariableBounds VarBounds = {{0, 1}, {0, 1}};             // Implicit bounds on variables

	 // Fill the values in the paramterized integer problem
	 for (unsigned int i = 0; i < numItems; ++i)
		IntegerIndexes.at(i) = i;

	 C(0, 0) = 2; // C terms in the objective
	 C(1, 1) = 3;
	 a(0, 0) = 3; // Knapsack Constraints
	 a(0, 1) = 4;
	 b(0)    = 5;  // Knapsack Capacity
	 c(0)    = -1; // The standard is minimization, hence minus
	 c(1)    = -2;

	 // Create a parametrized Integer Program
	 MathOpt::IP_Param PlayerOne(C, a, b, c, d, IntegerIndexes, VarBounds, &GurobiEnv);

	 // Parametrized Integer Program for the second player.
	 C(0, 0) = 5;
	 C(1, 1) = 4;
	 a(0, 0) = 2;
	 a(0, 1) = 5;
	 c(0)    = -3;
	 c(1)    = -5;

	 MathOpt::IP_Param PlayerTwo(C, a, b, c, d, IntegerIndexes, VarBounds, &GurobiEnv);

	 // Add the players to the instance. We can also specify a file path to write the instance
	 IPG_Instance.addIPParam(PlayerOne, "A_Parametrized_KnapsackProblem1");
	 IPG_Instance.addIPParam(PlayerTwo, "A_Parametrized_KnapsackProblem2");
	 IPG_Instance.save("A_Knapsack_Game"); // Save the instance with the standardize format
	 Models::IPG::IPG KnapsackGame(&GurobiEnv, IPG_Instance); // Create a model from the instance
	 // A few optional settings
	 KnapsackGame.setNumThreads(4);            // How many threads, if supported by the solver?
	 KnapsackGame.setTimeLimit(5);             // Time limit in second
	 KnapsackGame.finalize();                  // Lock the model
	 KnapsackGame.setDeviationTolerance(3e-4); // Numerical tolerance
	 // Run and get the results


	 // Cut and Play
	 KnapsackGame.setAlgorithm(Data::IPG::Algorithms::CutAndPlay);
	 KnapsackGame.setLCPAlgorithm(Data::LCP::Algorithms::MIP); // How do we solve the LCPs?
	 KnapsackGame.findNashEq();
	 std::cout << "The Cut-and-Play solution" << std::endl;
	 KnapsackGame.getX().at(0).print("Player 1:"); // Print the solution
	 KnapsackGame.getX().at(1).print("\n Player 2:");

	 // Zero Regrets
	 KnapsackGame.setAlgorithm(Data::IPG::Algorithms::ZERORegrets);
	 KnapsackGame.setGameObjective(Data::IPG::Objectives::ZERORegrets_PlayerOne);
	 KnapsackGame.findNashEq();
	 std::cout << "The ZERO Regrets solution" << std::endl;
	 KnapsackGame.getX().at(0).print("Player 1:"); // Print the solution
	 KnapsackGame.getX().at(1).print("\n Player 2:");

  } catch (ZEROException &e) {
	 throw ZEROException(e);
  }
}