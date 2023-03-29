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
  GRBEnv        GurobiEnv;
  std::string   logFile = "CNG_Results.csv";
  std::ifstream existCheck(logFile);
  std::ofstream results(logFile, std::ios::app);

  if (!existCheck.good()) {
	 results << "Instance;Algorithm;isPureNE;Status;Time;SocialWelfare;numIterations\n";
  }
  existCheck.close();
  try {
	 std::vector<int> node_list = {10, 20, 25, 50};
	 for (int n = 0; n < node_list.size(); ++n) {
		int nodes = node_list.at(n);
		for (unsigned instance = 1; instance < 21; instance++) {

		  Models::IPG::IPGInstance IPG_Instance; // The IPG Instance

		  arma::vec      c_attacker(nodes), d_attacker(nodes); // Profits c in the objective
		  arma::sp_mat   C_attacker(nodes, nodes);             // C terms in the objective
		  arma::sp_mat   a_attacker(1, nodes);                 // LHS for Knapsack constraint
		  arma::vec      b_attacker(1);                        // RHS for constraints
		  arma::vec      c_defender(nodes), d_defender(nodes); // Profits c in the objective
		  arma::sp_mat   C_defender(nodes, nodes);             // C terms in the objective
		  arma::sp_mat   a_defender(1, nodes);                 // LHS for Knapsack constraint
		  arma::vec      b_defender(1);                        // RHS for constraints
		  arma::vec      IntegerIndexes(nodes);                // The index of the integer variables
		  VariableBounds VarBounds;                            // Implicit bounds on variables
		  auto           data = arma::mat();
		  std::string    instanceName =
				"dat/CNG/Instance-" + std::to_string(nodes) + "_" + std::to_string(instance);
		  data.load(instanceName);
		  std::cout << "###INSTANCE " << instanceName << std::endl;
		  // data.load("dat/CNG/Toy");


		  data.print("Instance Data");
		  // Row 0: traffic
		  // Row 1: Costs attacker
		  // Row 2: Costs defender
		  // Row 3: Weights attacker
		  // Row 4: Weights defender
		  // Row 5: [d, a, eta, epsilon, gamma, delta, 0 if payoff ==OBJ.Attacker else 1]


		  double eta     = data.at(5, 2);
		  double epsilon = data.at(5, 3);
		  double gamma   = data.at(5, 4);
		  double delta   = data.at(5, 5);

		  // Fill the values in the paramterized integer problem
		  double A = 0, D = 0;
		  for (unsigned int i = 0; i < nodes; ++i) {
			 IntegerIndexes.at(i) = i;
			 VarBounds.emplace_back(0, 1);

			 double weightA = data.at(3, i);
			 double costA   = data.at(1, i);

			 c_attacker.at(i)    = -costA * (1 + gamma);
			 d_attacker.at(i)    = -costA * gamma;
			 C_attacker.at(i, i) = -costA * (-gamma - eta);
			 a_attacker.at(i)    = weightA;
			 A += weightA;

			 double weightD      = data.at(4, i);
			 double costD        = data.at(2, i);
			 c_defender.at(i)    = -costD * (epsilon - 1);
			 d_defender.at(i)    = -costD * (delta - 1);
			 C_defender.at(i, i) = -costD * ((1 + eta - epsilon - delta));
			 a_defender.at(i)    = weightD;
			 D += weightD;
		  }

		  b_attacker.at(0) = floor(A * data.at(5, 1));
		  b_defender.at(0) = floor(D * data.at(5, 0));

		  // Create a parametrized Integer Program
		  MathOpt::IP_Param Attacker(C_attacker,
											  a_attacker,
											  b_attacker,
											  c_attacker,
											  d_attacker,
											  IntegerIndexes,
											  VarBounds,
											  &GurobiEnv);
		  MathOpt::IP_Param Defender(C_defender,
											  a_defender,
											  b_defender,
											  c_defender,
											  d_defender,
											  IntegerIndexes,
											  VarBounds,
											  &GurobiEnv);

		  // Add the players to the instance. We can also specify a file path to write the instance
		  IPG_Instance.addIPParam(Defender, "Defender_Problem");
		  IPG_Instance.addIPParam(Attacker, "Attacker_Problem");
		  IPG_Instance.save("CriticalNodeGame"); // Save the instance with the standardize format

		  Models::IPG::IPG CNG(&GurobiEnv, IPG_Instance); // Create a model from the instance
		  // A few optional settings
		  CNG.setNumThreads(8);            // How many threads, if supported by the solver?
		  CNG.setTimeLimit(50);            // Time limit in second
		  CNG.finalize();                  // Lock the model
		  CNG.setDeviationTolerance(3e-4); // Numerical tolerance
		  // Cut and Play
		  CNG.setAlgorithm(Data::IPG::Algorithms::CutAndPlay);
		  CNG.setGameObjective(Data::IPG::Objectives::Feasibility);
		  CNG.setLCPAlgorithm(Data::LCP::Algorithms::PATH); // How do we solve the LCPs?
		  CNG.findNashEq();
		  std::cout << "The Cut-and-Play solution" << std::endl;
		  CNG.getX().at(0).print("Player 1:"); // Print the solution
		  CNG.getX().at(1).print("\n Player 2:");
		  auto stat = CNG.getStatistics();

		  results << std::to_string(nodes) + "_" + std::to_string(instance) << ";CutAndPlay;"
					 << std::to_string(stat.PureNashEquilibrium.get()) << ";"
					 << std::to_string(stat.Status.get()) << ";"
					 << std::to_string(stat.WallClockTime.get()) << ";" << CNG.getSocialWelfare() << ";"
					 << stat.NumIterations.get() << "\n";

		  Models::IPG::IPG CNG2(&GurobiEnv, IPG_Instance); // Create a model from the instance
		  CNG2.setNumThreads(8);            // How many threads, if supported by the solver?
		  CNG2.setTimeLimit(50);            // Time limit in second
		  CNG2.finalize();                  // Lock the model
		  CNG2.setDeviationTolerance(3e-4); // Numerical tolerance
		  // Zero Regrets
		  CNG2.setAlgorithm(Data::IPG::Algorithms::ZERORegrets);
		  CNG2.setGameObjective(Data::IPG::Objectives::ZERORegrets_PlayerOne);
		  CNG2.findNashEq();
		  std::cout << "The ZERO Regrets solution" << std::endl;
		  CNG2.getX().at(0).print("Player 1:"); // Print the solution
		  CNG2.getX().at(1).print("\n Player 2:");
		  auto stat2 = CNG2.getStatistics();

		  results << std::to_string(nodes) + "_" + std::to_string(instance) << ";ZERORegrets;"
					 << std::to_string(stat2.PureNashEquilibrium.get()) << ";"
					 << std::to_string(stat2.Status.get()) << ";"
					 << std::to_string(stat2.WallClockTime.get()) << ";" << CNG2.getSocialWelfare()
					 << ";" << stat2.NumIterations.get() << "\n";

		  // break;
		  results.flush();
		}
	 }
  } catch (ZEROException &e) {
	 throw ZEROException(e);
  }
  results.close();
}