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

int main() {
  GRBEnv GurobiEnv;
  try {
	 Models::IPG::IPGInstance IPG_Instance; // The IPG Instance


	 std::vector<double> a_in       = {10, 45, 35};
	 std::vector<double> b_in       = {1000, 100, 2000};
	 std::vector<double> c_in       = {0.05, 0.1, 0.002};
	 std::vector<double> UB_in      = {600, 250, 500};
	 std::vector<double> LB_in      = {400, 200, 300};
	 int                 numPlayers = 3;

	 // 3 vars (in order) : --- x^i_c --- x^i_b  ---- z^i
	 int    numVars   = 3;
	 int    numConstr = 5;
	 double alpha     = 200;
	 double beta      = 0.2;


	 for (int i = 0; i < numPlayers; ++i) {
		arma::vec      c(numVars);                             // Linear profit in the objective
		arma::sp_mat   C(numVars, numVars * (numPlayers - 1)); // C terms in the objective
		arma::sp_mat   A(numConstr, numVars);                  // Constraints
		arma::vec      b(numConstr);                           // RHS for constraints
		arma::vec      IntegerIndexes = {1};                   // The index of the integer variables
		VariableBounds VarBounds      = {{LB_in.at(i), UB_in.at(i)},
													{0, 1},
													{0, UB_in.at(i) * UB_in.at(i)}}; // Bounds on variables

		// Linear coefficients of x^i_c
		c.at(0) = a_in.at(i) - alpha;
		// Linear coefficients of x^i_b
		c.at(1) = b_in.at(i);
		// Linear coefficients of z^i
		c.at(2) = 0.5 * c_in.at(i) + beta;

		for (unsigned int o = 0, index = 0; o < numPlayers; ++o) {
		  if (i != o) {
			 C.at(0, index * numVars) = beta;
			 index++;
		  }
		}

		// Lower bound
		A.at(0, 1) = LB_in.at(i);
		A.at(0, 0) = -1;
		b.at(0)    = 0;
		// Upper bound
		A.at(1, 1) = -UB_in.at(i);
		A.at(1, 0) = 1;
		b.at(1)    = 0;
		// McCormick 1
		A.at(2, 2) = 1;
		A.at(2, 0) = -(LB_in.at(i) + UB_in.at(i));
		b.at(2)    = -LB_in.at(i) * UB_in.at(i);
		// McCormick 2
		A.at(3, 2) = -1;
		A.at(3, 0) = 2 * LB_in.at(i);
		b.at(3)    = LB_in.at(i) * LB_in.at(i);
		// McCormick 3
		A.at(4, 2) = -1;
		A.at(4, 0) = 2 * UB_in.at(i);
		b.at(4)    = UB_in.at(i) * UB_in.at(i);


		MathOpt::IP_Param Player(C, A, b, c, IntegerIndexes, VarBounds, &GurobiEnv);
		IPG_Instance.addIPParam(Player, "UCP_Player" + std::to_string(i));

		arma::vec xTest(6, arma::fill::zeros);
		Player.getIPModel(xTest, false)->write("test.lp");
	 }



	 IPG_Instance.save("UCP_Game"); // Save the instance with the standardized format
	 Models::IPG::IPG UCPGame(&GurobiEnv, IPG_Instance); // Create a model from the instance
	 // Select the equilibrium to compute a Nash Equilibrium
	 UCPGame.setAlgorithm(Data::IPG::Algorithms::CutAndPlay);
	 // A few optional settings
	 UCPGame.setDeviationTolerance(3e-4); // Numerical tolerance
	 UCPGame.setNumThreads(4);            // How many threads, if supported by the solver?
	 UCPGame.setLCPAlgorithm(Data::LCP::Algorithms::MIP); // How do we solve the LCPs?
	 UCPGame.setTimeLimit(5);                             // Time limit in second
	 UCPGame.finalize();                                  // Lock the model
	 // Run and get the results
	 UCPGame.findNashEq();
	 UCPGame.getX().at(0).print("Player 1:"); // Print the solution
	 UCPGame.getX().at(1).print("\n Player 2:");

  } catch (ZEROException &e) {
	 throw ZEROException(e);
  }
}