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

#include "games/algorithms/IPG/ipg_zeroregrets.h"



/**
 * @brief This method initializes some fields for the algorithm. Also, it warm starts the
 * initial strategies to pure best responses.
 */
void Algorithms::IPG::ZERORegrets::initialize() {
  // Set the number of iterations to zero.
  this->IPG->Stats.NumIterations.set(0);
  // Reset cuts statistics
  this->Cuts = {std::pair<std::string, int>("EI", 0)};
  this->IPG->Stats.AlgorithmData.Cuts.set(this->Cuts);

  // Initialize the equilibrium MIP
  this->JointProgram = std::move(std::make_unique<GRBModel>(*this->Env));

  // X vars
  this->x = new GRBVar *[this->IPG->NumVariables];
  this->p = new GRBVar[this->IPG->NumPlayers];

  // solutions
  this->xLast = std::vector<arma::vec>(this->IPG->NumPlayers);

  // Assert all players solve pure integer programs and initialize the variables
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 auto nVarPlayer   = this->IPG->PlayersIP.at(i)->getNumVars();
	 this->xLast.at(i) = arma::vec();
	 ZEROAssert(nVarPlayer == this->IPG->PlayersIP.at(i)->getIntegers().size());
	 this->x[i] = new GRBVar[nVarPlayer];
	 this->p[i] = this->JointProgram->addVar(
		  -GRB_INFINITY, +GRB_INFINITY, 0, GRB_CONTINUOUS, "p_" + std::to_string(i));
  }



  // Add the variables and the constraints, for each player
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 auto Bounds = this->IPG->PlayersIP.at(i)->getBounds();

	 auto nVarPlayer = this->IPG->PlayersIP.at(i)->getNumVars();
	 for (int j = 0; j < nVarPlayer; j++)
		this->x[i][j] = this->JointProgram->addVar(Bounds.at(i).first,
																 Bounds.at(i).second,
																 0,
																 GRB_INTEGER,
																 "x" + std::to_string(i) + "_" + std::to_string(j));


	 this->JointProgram->update();
	 Utils::addSparseConstraints(this->IPG->PlayersIP.at(i)->getB(false),
										  this->IPG->PlayersIP.at(i)->getb(false),
										  this->x[i],
										  "Constr_",
										  this->JointProgram.get(),
										  GRB_LESS_EQUAL,
										  nullptr);
  }


  // Add the payoff variable for each player
  GRBQuadExpr payoff = 0;
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 payoff        = 0;
	 auto C        = this->IPG->PlayersIP.at(i)->getC();
	 auto c        = this->IPG->PlayersIP.at(i)->getc();
	 auto d        = this->IPG->PlayersIP.at(i)->getd();
	 int  counterd = 0;
	 for (int j = 0; j < this->IPG->PlayersIP.at(i)->getNumVars(); j++) {
		payoff.addTerm(c.at(j), this->x[i][j]);
		int counter = 0;
		for (unsigned int o = 0; o < this->IPG->NumPlayers; ++o) {
		  if (o != i) {
			 payoff.addTerm(d.at(counterd), this->x[o][j]);
			 counterd += 1;
			 for (int k = 0; k < this->IPG->PlayersIP.at(o)->getNumVars(); k++) {
				payoff.addTerm(C.at(j, counter), this->x[i][j], this->x[o][k]);
				counter += 1;
			 }
		  }
		}
	 }

	 this->JointProgram->addQConstr(this->p[i] == payoff);
  }

  bool        socialcost = true;
  GRBQuadExpr objective  = 0;
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 objective += p[i];
	 if (this->IPG->Stats.AlgorithmData.Objective.get() ==
		  Data::IPG::Objectives::ZERORegrets_PlayerOne) {
		socialcost = false;
		LOG_S(INFO) << "Algorithms::IPG::ZERORegrets::initialize: Setting the default objective to "
							"PlayerOne.";
		break;
	 }
  }
  if (socialcost)
	 LOG_S(INFO) << "Algorithms::IPG::ZERORegrets::initialize: Setting the default objective to "
						 "SocialCost.";
  this->JointProgram->setObjective(objective, GRB_MINIMIZE);
}
/**
 * @brief Checks if there is more time remaining.
 * @param remaining An output filled with the time remaining
 * @return True if there is still time left.
 */
bool Algorithms::IPG::ZERORegrets::checkTime(double &remaining) const {
  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
	 const std::chrono::duration<double> timeElapsed =
		  std::chrono::high_resolution_clock::now() - this->IPG->InitTime;
	 remaining = this->IPG->Stats.AlgorithmData.TimeLimit.get() - timeElapsed.count();
	 if (remaining <= 0) {
		LOG_S(1) << "Algorithms::IPG::CutAndPlay::checkTime: "
						"Time limit hit.";
		this->IPG->Stats.AlgorithmData.Cuts.set(this->Cuts);
		if (this->IPG->Stats.Status.get() == ZEROStatus::Uninitialized)
		  this->IPG->Stats.Status.set(ZEROStatus::TimeLimit);
		return false;
	 } else
		return true;
  } else {
	 remaining = -1;
	 return true;
  }
}

/**
 * @brief Generates an equilibrium inequality starting from xOfI
 * @param player
 * @param xOfI
 * @return
 */
bool Algorithms::IPG::ZERORegrets::addEquilibriumInequality(unsigned int     player,
																				const arma::vec &xOfI) {
  GRBLinExpr payoff = 0;
  payoff            = 0;
  auto C            = this->IPG->PlayersIP.at(player)->getC();
  auto c            = this->IPG->PlayersIP.at(player)->getc();
  auto d            = this->IPG->PlayersIP.at(player)->getd();
  int  counterd     = 0;
  for (int j = 0; j < this->IPG->PlayersIP.at(player)->getNumVars(); j++) {
	 payoff += c.at(j) * xOfI.at(j);
	 int counter = 0;
	 for (unsigned int o = 0; o < this->IPG->NumPlayers; ++o) {
		if (o != player) {
		  payoff += d.at(counterd) * this->x[o][j];
		  counterd += 1;
		  for (int k = 0; k < this->IPG->PlayersIP.at(o)->getNumVars(); k++) {
			 payoff += (C.at(j, counter) * xOfI.at(j)) * this->x[o][k];
			 counter += 1;
		  }
		}
	 }
  }
  this->JointProgram->addQConstr(this->p[player] <= payoff);
  this->JointProgram->update();

  return true;
}



/**
 * @brief Solves the IPG with the Equilibrium CutAndPlay algorithm.
 */
void Algorithms::IPG::ZERORegrets::solve() {


  this->initialize();
  // Main loop condition
  bool solved{false};
  int  Iteration = 0;

  this->JointProgram->set(GRB_IntParam_Threads, this->IPG->Stats.AlgorithmData.Threads.get());
  this->JointProgram->set(GRB_IntParam_OutputFlag, 0);
  while (!solved) {
	 // Increase the number of iterations
	 Iteration++;
	 LOG_S(INFO) << "Algorithms::IPG::ZERORegrets::solve: Iteration ########### " << Iteration
					 << " ###########";
	 this->IPG->Stats.NumIterations.set(Iteration);


	 /* ************************************
	  * Check time and solve the MIP
	  **************************************/
	 int MIPStatus = GRB_UNBOUNDED;

	 if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
		double remaining;
		if (this->checkTime(remaining) && remaining > 0) {
		  this->JointProgram->set(GRB_DoubleParam_TimeLimit, remaining * 0.95);
		} else {
		  this->IPG->Stats.AlgorithmData.Cuts.set(this->Cuts);
		  return;
		}
	 }

	 this->JointProgram->optimize();
	 this->JointProgram->write("Joint.lp");
	 MIPStatus = this->JointProgram->get(GRB_IntAttr_Status);

	 if (MIPStatus != GRB_OPTIMAL and MIPStatus != GRB_TIME_LIMIT) {
		this->IPG->Stats.AlgorithmData.Cuts.set(this->Cuts);
		if (!this->Solved) {
		  this->IPG->Stats.Status.set(ZEROStatus::NashEqNotFound);
		  LOG_S(INFO) << "Algorithms::IPG::ZERORegrets::solve: A Nash Equilibrium has not been "
							  "found.";
		  return;
		} else {
		  LOG_S(INFO) << "Algorithms::IPG::ZERORegrets::solve: Problem closed for infeasibility.";
		  return;
		}
	 }


	 /* ************************************
	  * Support objects
	  **************************************/
	 // How many cuts added in total?
	 unsigned int addedCuts = 0;

	 try {
		if (MIPStatus == GRB_OPTIMAL or MIPStatus == GRB_TIME_LIMIT) {
		  // Reset solution status
		  solved = true;
		  // last objective
		  auto      last_obj       = this->JointProgram->getObjective().getValue();
		  arma::vec last_payoffs   = arma::zeros(this->IPG->NumPlayers);
		  this->IPG->SocialWelfare = 0;
		  // Fetch the solution. This goes outside the main loop
		  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
			 last_payoffs.at(i) = p[i].get(GRB_DoubleAttr_X);
			 this->IPG->SocialWelfare += last_payoffs.at(i);
			 this->xLast.at(i) = arma::zeros(this->IPG->PlayersIP.at(i)->getNumVars());
			 for (int j = 0; j < this->IPG->PlayersIP.at(i)->getNumVars(); j++) {
				this->xLast.at(i).at(j) = x[i][j].get(GRB_DoubleAttr_X);
				// std::cout << "x_" << i << "," << j
				//			 << "=" << this->xLast.at(i).at(j) << std::endl;
			 }
		  }


		  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {

			 arma::vec xMinusI = arma::zeros(0);
			 for (unsigned int o = 0; o < this->IPG->NumPlayers; ++o)
				if (i != o)
				  xMinusI.insert_rows(xMinusI.n_rows, this->xLast.at(o));


			 auto      BR_Program = this->IPG->PlayersIP.at(i)->solveFixed(xMinusI, true);
			 arma::vec BR         = arma::zeros(this->IPG->PlayersIP.at(i)->getNumVars());
			 auto      BR_payoff  = BR_Program->getObjective().getValue();

			 // std::cout << "BR: ";
			 for (int j = 0; j < this->IPG->PlayersIP.at(i)->getNumVars(); j++) {
				BR.at(j) = BR_Program->getVarByName("y_" + std::to_string(j)).get(GRB_DoubleAttr_X);
				// std::cout << BR.at(j) << " ";
			 }
			 // std::cout << std::endl;

			 // std::cout << last_payoffs.at(i) << " versus br of " << BR_payoff << "\n";
			 if ((last_payoffs.at(i) - BR_payoff) >
				  this->IPG->Stats.AlgorithmData.DeviationTolerance.get()) {
				// There is a deviation, not solved!
				solved = false;
				this->addEquilibriumInequality(i, BR);
				addedCuts += 1;
			 }
		  }

		  // cutoff
		  // this->JointProgram->set(GRB_DoubleParam_Cutoff, last_obj);
		}

		this->Cuts.at(0).second += addedCuts;
		LOG_S(INFO) << "Algorithms::IPG::ZERORegrets::initialize: Added " << std::to_string(addedCuts)
						<< " cuts.";
		this->IPG->Stats.AlgorithmData.Cuts.set(this->Cuts);
	 } catch (GRBException &e) {
		throw ZEROException(e);
	 }
  }


  if (solved) {
	 LOG_S(INFO) << "Algorithms::IPG::ZERORegrets::solve: A Nash Equilibrium has been found (PNE).";
	 this->IPG->Stats.Status.set(ZEROStatus::NashEqFound);
	 this->Solved = true;
	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i)
		this->IPG->Solution.at(i) = this->xLast.at(i);

	 this->IPG->Stats.AlgorithmData.Cuts.set(this->Cuts);
  }



  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
	 double remaining;
	 if (!this->checkTime(remaining) || remaining <= 0) {
		this->IPG->Stats.AlgorithmData.Cuts.set(this->Cuts);
		return;
	 }
  }
}
