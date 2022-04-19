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
  GRBEnv GurobiEnv;
  try {
	 Models::IPG::IPGInstance IPG_Instance; // The IPG Instance
	 int                      numVars = 20, numPlayers = 2;
	 int                      LB = -500, UB = 500, boundCorrector = 0;

	 int  distrBound = 5;
	 int  size       = numVars * numPlayers;
	 bool test       = true;


	 arma::mat JF;
	 while (test) {
		// Random matrix
		JF = arma::conv_to<arma::mat>::from(
			 arma::randi(size, size, arma::distr_param(-distrBound, +distrBound)));
		// Make it PSD
		JF   =  (0.5 * (JF + JF.t()) + size * arma::eye(size, size));
		test = !JF.is_sympd();
	 }

	 double M_max = std::max(JF.max(), abs(JF.min()));
	 double asym  = round(distrBound / 10);


	 for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
		  if (i < j) {
			 int v    = arma::randi<int>(arma::distr_param(-floor(M_max * asym), ceil(M_max * asym)));
			 JF[i, j] = JF[i, j] + v;
			 JF[j, i] = JF[j, i] - v;
		  }
		}
	 }
	 JF.print("JF");



	 for (unsigned int i = 0; i < numPlayers; ++i) {


		arma::vec c = arma::conv_to<arma::vec>::from(arma::randi(
						  numVars, arma::distr_param(0, +distrBound))); // Profits c in the objective
		arma::sp_mat a(0, numVars);                                 // LHS for Knapsack constraint
		arma::vec    b;                                             // RHS for constraints
		arma::vec    IntegerIndexes(numVars); // The index of the integer variables

		VariableBounds VarBounds = {{LB, UB}, {LB, UB}}; // Implicit bounds on variables


		// Fill the values in the paramterized integer problem
		for (unsigned int j = 0; j < numVars; ++j)
		  IntegerIndexes.at(j) = j;


		arma::mat C = JF.submat(numVars * i, 0, numVars * (i + 1) - 1, JF.n_cols - 1);
		C.shed_cols(numVars * i, numVars * (i + 1) - 1);
		arma::sp_mat Ct(C);
		assert(Ct.n_cols == numVars * (numPlayers - 1) && Ct.n_rows == numVars);

		MathOpt::IP_Param Player(-Ct, a, b, -c, IntegerIndexes, VarBounds, &GurobiEnv);
		IPG_Instance.addIPParam(Player, "Parametrized_" + std::to_string(i));
		Ct.print_dense("Ct");
		c.print("c");
	 }
	 std::cout << "Generated!"<<std::endl;

	 IPG_Instance.save("AnInstance"); // Save the instance with the standardize format
	 Models::IPG::IPG KnapsackGame(&GurobiEnv, IPG_Instance); // Create a model from the instance
	 // Select the equilibrium to compute a Nash Equilibrium
	 KnapsackGame.setAlgorithm(Data::IPG::Algorithms::CutAndPlay);
	 // A few optional settings
	 KnapsackGame.setDeviationTolerance(3e-4); // Numerical tolerance
	 KnapsackGame.setNumThreads(4);            // How many threads, if supported by the solver?
	 KnapsackGame.setLCPAlgorithm(Data::LCP::Algorithms::MIP); // How do we solve the LCPs?
	 KnapsackGame.setCutsAggressiveness(Data::IPG::CutsAggressiveness::NotEvenTry);
	 KnapsackGame.setTimeLimit(10000); // Time limit in second
	 KnapsackGame.setGameObjective(Data::IPG::Objectives::Feasibility);
	 KnapsackGame.finalize();      // Lock the model
	 // Run and get the results
	 KnapsackGame.findNashEq();
	 KnapsackGame.getX().at(0).print("Player 1:"); // Print the solution
	 KnapsackGame.getX().at(1).print("\n Player 2:");

  } catch (ZEROException &e) {
	 throw ZEROException(e);
  }
}