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

double bestResponse(GRBEnv *env, double a, double b, double c, double LB, double UB, double price) {

  GRBModel brModel(*env);

  GRBVar xContinuous = brModel.addVar(0.0, UB, 0, GRB_CONTINUOUS, "xCont");
  GRBVar xBinary     = brModel.addVar(0.0, 1, 0, GRB_BINARY, "xBin");

  brModel.addConstr(LB * xBinary <= xContinuous);
  brModel.addConstr(xContinuous <= UB * xBinary);

  brModel.setObjective(a * xContinuous + 0.5 * c * xContinuous * xContinuous + b * xBinary -
									price * xContinuous,
							  GRB_MINIMIZE);

  brModel.set(GRB_IntParam_OutputFlag, 0);
  brModel.set(GRB_IntParam_NonConvex, 2);
  brModel.optimize();

  return brModel.getObjective().getValue();
}


int main() {
  loguru::g_stderr_verbosity = 0;
  GRBEnv GurobiEnv;
  try {
	 Models::IPG::IPGInstance IPG_Instance; // The IPG Instance


	 std::vector<double> a_in         = {10, 45, 35, 20};
	 std::vector<double> c_in         = {0.05, 0.1, 0.002, 0.05};
	 std::vector<double> b_in         = {4000, 100, 2000, 1500};
	 std::vector<double> LB_in        = {400, 200, 300, 500};
	 std::vector<double> UB_in        = {600, 250, 500, 800};
	 std::vector<double> solutionSeed = {502.5, 0, 300};
	 double              alpha        = 500;
	 double              beta         = 0.2;

	 int  numPlayers = a_in.size();
	 bool solution   = false;

	 // 3 vars (in order) : --- x^i_c --- x^i_b  ---- z^i
	 int numVars   = 3;
	 int numConstr = 5 + (solution ? 2 : 0);



	 for (int i = 0; i < numPlayers; ++i) {
		arma::vec      c(numVars); // Linear profit in the objective
		arma::vec      d(numVars * (numPlayers - 1), arma::fill::zeros);
		arma::sp_mat   C(numVars, numVars * (numPlayers - 1)); // C terms in the objective
		arma::sp_mat   A(numConstr, numVars);                  // Constraints
		arma::vec      b(numConstr);                           // RHS for constraints
		arma::vec      IntegerIndexes = {1};                   // The index of the integer variables
		VariableBounds VarBounds      = {
          {0, UB_in.at(i)}, {0, 1}, {0, UB_in.at(i) * UB_in.at(i) + 10}}; // Bounds on variables

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
		A.at(0, 0) = -1;
		A.at(0, 1) = LB_in.at(i);
		b.at(0)    = 0;
		// Upper bound
		A.at(1, 0) = 1;
		A.at(1, 1) = -UB_in.at(i);
		b.at(1)    = 0;


		// McCormick 1
		A.at(2, 2) = 1;
		A.at(2, 0) = -(LB_in.at(i) + UB_in.at(i));
		A.at(2, 1) = LB_in.at(i) * UB_in.at(i);
		b.at(2)    = 0;
		// McCormick 2
		A.at(3, 2) = -1;
		A.at(3, 0) = 2 * LB_in.at(i);
		A.at(3, 1) = -LB_in.at(i) * LB_in.at(i);
		b.at(3)    = 0;
		//  McCormick 3
		A.at(4, 2) = -1;
		A.at(4, 0) = 2 * UB_in.at(i);
		A.at(4, 1) = -UB_in.at(i) * UB_in.at(i);
		b.at(4)    = 0;

		if (solution) {
		  A.at(5, 0) = 1;
		  b.at(5)    = solutionSeed.at(i);
		  A.at(6, 0) = -1;
		  b.at(6)    = -solutionSeed.at(i);
		}


		MathOpt::IP_Param Player(C, A, b, c, d, IntegerIndexes, VarBounds, &GurobiEnv);
		IPG_Instance.addIPParam(Player, "UCP_Player" + std::to_string(i));


		// arma::vec xTest(numVars * (numPlayers - 1), arma::fill::zeros);
		// Player.getIPModel(xTest, false)->write("test" + std::to_string(i) + ".lp");
	 }



	 IPG_Instance.save("UCP_Game"); // Save the instance with the standardized format
	 Models::IPG::IPG UCPGame(&GurobiEnv, IPG_Instance); // Create a model from the instance
	 // Select the equilibrium to compute a Nash Equilibrium
	 UCPGame.setAlgorithm(Data::IPG::Algorithms::CutAndPlay);
	 // A few optional settings
	 UCPGame.setDeviationTolerance(1e-6); // Numerical tolerance
	 UCPGame.setNumThreads(8);            // How many threads, if supported by the solver?
	 UCPGame.setLCPAlgorithm(Data::LCP::Algorithms::PATH); // How do we solve the LCPs?
	 UCPGame.setGameObjective(Data::IPG::Objectives::Quadratic);
	 UCPGame.setTimeLimit(1000000); // Time limit in second
	 UCPGame.finalize();            // Lock the model
	 UCPGame.setPresolve(false);
	 // Run and get the results
	 UCPGame.findNashEq();

	 double quantity = 0;
	 for (int i = 0; i < numPlayers; ++i) {
		UCPGame.getX().at(i).print("Player " + std::to_string(i) + ":");
		quantity += UCPGame.getX().at(i).at(0);
	 }

	 double price      = alpha - beta * quantity;
	 double maxEpsilon = 0;

	 std::cout << "\n\nThe social welfare is: \t\t" << UCPGame.getSocialWelfare() << std::endl;
	 std::cout << "The Price is: \t\t" << price << std::endl;
	 std::cout << "The Quantity is: \t\t" << quantity << std::endl;
	 for (int i = 0; i < numPlayers; ++i) {
		double x = UCPGame.getX().at(i).at(0);
		double poff =
			 (fabs(x) < 1e-5 ? 0 : a_in.at(i) * x + 0.5 * c_in.at(i) * x * x + b_in.at(i) - price * x);
		std::cout << "\tPlayer " << i << " produces " << x << " and has a (true) payoff of " << poff
					 << std::endl;

		double epsilon =
			 poff -
			 bestResponse(
				  &GurobiEnv, a_in.at(i), b_in.at(i), c_in.at(i), LB_in.at(i), UB_in.at(i), price);
		maxEpsilon = epsilon > maxEpsilon ? epsilon : maxEpsilon;
	 }
	 std::cout << "The MaxEpsilon is: \t\t" << maxEpsilon << std::endl;



	 GRBModel potentialModel(GurobiEnv);

	 GRBVar xContinuous[numPlayers], xBinary[numPlayers], p, payoff[numPlayers];
	 p = potentialModel.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "price");
	 GRBLinExpr sumQuantities = 0;

	 for (int i = 0; i < numPlayers; ++i) {
		xContinuous[i] =
			 potentialModel.addVar(0.0, UB_in.at(i), 0, GRB_CONTINUOUS, "xCont_" + std::to_string(i));
		sumQuantities += xContinuous[i];
		xBinary[i] = potentialModel.addVar(0.0, 1, 0, GRB_BINARY, "xBin_" + std::to_string(i));
		payoff[i]  = potentialModel.addVar(
          -GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "payoff_" + std::to_string(i));

		potentialModel.addConstr(LB_in.at(i) * xBinary[i] <= xContinuous[i]);
		potentialModel.addConstr(xContinuous[i] <= UB_in.at(i) * xBinary[i]);
	 }

	 potentialModel.addConstr(p == alpha - beta * sumQuantities);

	 GRBQuadExpr objective;
	 for (int i = 0; i < numPlayers; ++i) {
		potentialModel.addQConstr(a_in.at(i) * xContinuous[i] +
												0.5 * c_in.at(i) * xContinuous[i] * xContinuous[i] +
												b_in.at(i) * xBinary[i] - p * xContinuous[i] ==
										  payoff[i]);
		objective += payoff[i];
		GRBLinExpr sumOtherQuantities = 0;
		for (int o = 0; o < numPlayers; ++o) {
		  if (i != o)
			 sumOtherQuantities += xContinuous[o];
		}
		objective += 0.5 * xContinuous[i] * sumOtherQuantities;
	 }
	 potentialModel.set(GRB_IntParam_OutputFlag, 0);
	 potentialModel.set(GRB_IntParam_NonConvex, 2);
	 potentialModel.setObjective(objective, GRB_MAXIMIZE);
	 potentialModel.optimize();
	 potentialModel.write("dat/potential.lp");

	 double sw = 0;
	 quantity  = 0;
	 for (int i = 0; i < numPlayers; ++i) {
		sw += payoff[i].get(GRB_DoubleAttr_X);
		quantity += xContinuous[i].get(GRB_DoubleAttr_X);
	 }

	 price      = alpha - beta * quantity;
	 maxEpsilon = 0;

	 std::cout << "\n\n------Exact Potential Function\nThe social welfare is: \t\t" << sw
				  << std::endl;
	 std::cout << "The Price is: \t\t" << p.get(GRB_DoubleAttr_X) << std::endl;
	 std::cout << "The Quantity is: \t\t" << sumQuantities.getValue() << std::endl;
	 for (int i = 0; i < numPlayers; ++i) {
		double x = xContinuous[i].get(GRB_DoubleAttr_X);
		double poff =
			 (fabs(x) < 1e-5 ? 0 : a_in.at(i) * x + 0.5 * c_in.at(i) * x * x + b_in.at(i) - price * x);
		std::cout << "\tPlayer " << i << " produces " << x << " and has a (true) payoff of " << poff
					 << std::endl;
		double guPayoff = payoff[i].get(GRB_DoubleAttr_X);
		assert(fabs(poff - guPayoff) < 1);
		double epsilon =
			 poff -
			 bestResponse(
				  &GurobiEnv, a_in.at(i), b_in.at(i), c_in.at(i), LB_in.at(i), UB_in.at(i), price);
		maxEpsilon = epsilon > maxEpsilon ? epsilon : maxEpsilon;
	 }
	 std::cout << "The MaxEpsilon is: \t\t" << maxEpsilon << std::endl;
	 assert(maxEpsilon == 0);



  } catch (ZEROException &e) {
	 throw ZEROException(e);
  }
}