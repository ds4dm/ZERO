
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

#include "games/algorithms/IPG/ipg_cutandplay.h"
#include "coin/CglGMI.hpp"
#include "coin/CglKnapsackCover.hpp"
#include "coin/CoinPackedMatrix.hpp"
#include "coin/OsiSolverInterface.hpp"
#include <coin/CglGomory.hpp>
#include <coin/CglMixedIntegerRounding2.hpp>
#include <memory>

/**
 * @brief @brief Given @p vertex, it adds a vertex to the field V. If @p checkDuplicate is true,
 * it will check whether the vertex is already contained in the pool.
 * @param vertex The input vertex
 * @param checkDuplicate True if the method needs to check for duplicates
 * @return True if the vertex is added.
 */
bool Algorithms::IPG::IPG_Player::addVertex(const arma::vec &vertex, const bool checkDuplicate) {
  bool go = false;
  if (checkDuplicate)
	 go = Utils::containsRow(this->V, vertex, this->Tolerance);


  if (!go) {
	 int nCols = this->V.n_cols < 1 ? vertex.size() : this->V.n_cols;
	 this->V.resize(this->V.n_rows + 1, nCols);
	 this->V.row(this->V.n_rows - 1) = vertex.t();
	 return true;
  }
  return false;
}

/**
 * @brief @brief Given @p ray, it adds a ray to the field R. If @p checkDuplicate is true,
 * it will check whether the ray is already contained in the pool.
 * @param ray The input ray
 * @param checkDuplicate True if the method needs to check for duplicates
 * @return True if the ray is added.
 */
bool Algorithms::IPG::IPG_Player::addRay(const arma::vec &ray, const bool checkDuplicate) {
  /**
	* @brief Given @p ray, it adds a ray to the field R. If @p checkDuplicate is true, it
	* will check whether the ray is already contained in the bool.
	* @return true if the ray is added.
	*/
  bool go = false;
  if (checkDuplicate)
	 go = Utils::containsRow(this->R, ray, this->Tolerance);

  if (!go) {
	 int nCols = this->R.n_cols < 1 ? ray.size() : this->R.n_cols;
	 this->R.resize(this->R.n_rows + 1, nCols);
	 this->R.row(this->R.n_rows - 1) = ray.t();
	 return true;
  }
  return false;
}


/**
 * @brief Given @p LHS, @p RHS, it adds the inequalities to the field CutPool_A and b, the
 * CoinModel, and the working IP_Param.
 * @param LHS The LHS matrix
 * @param RHS The RHS vector
 * @return true if the inequality is added.
 */
bool Algorithms::IPG::IPG_Player::addCuts(const arma::sp_mat &LHS, const arma::vec &RHS) {

  arma::sp_mat LHSClean = Utils::clearMatrix(LHS, 1e-6, 1 - 1e-6);
  arma::vec    RHSClean = Utils::clearVector(RHS, 1e-6, 1 - 1e-6);
  unsigned int newCuts  = LHS.n_rows;
  unsigned int cuts     = this->CutPool_A.n_rows;
  int          nCols    = this->CutPool_A.n_cols < 1 ? LHSClean.size() : this->CutPool_A.n_cols;

  // Add the constraints to the cut pool
  this->CutPool_A.resize(cuts + newCuts, nCols);
  this->CutPool_b.resize(cuts + newCuts);
  this->CutPool_A.submat(cuts, 0, cuts + newCuts - 1, nCols - 1) = LHSClean;
  this->CutPool_b.subvec(cuts, cuts + newCuts - 1)               = RHSClean;
  // Add the constraints to the Coin model


  /******RECURSIVE CUT GENERATION********
  //Uncomment this to enable recursive cut generation
 auto convertedCuts = Utils::armaToCoinPackedVector(LHS);
 try {
	for (unsigned int i = 0; i < newCuts; ++i)
	  this->CoinModel->addRow(convertedCuts.at(i), 'L', RHS.at(i), 0);

 } catch (CoinError &e) {
	throw ZEROException(ZEROErrorCode::SolverError,
							  "Invalid Coin-OR interface response: " + e.message());
 }
  *******RECURSIVE CUT GENERATION*********/

  // Add the constraints to the parametrized IP
  this->ParametrizedIP->addConstraints(LHSClean, RHSClean);

  //******DEBUG********
  // LHS.print_dense("LHS");
  // RHS.print("RHS");
  //******DEBUG********
  return true;
}

/**
 * @brief @brief Given a player @p player, one of its best responses @p xOfIBestResponses, the
 * strategies of the other players @p xMinusI, it adds an inequality of the type \f[ f^i(x^i,\bar
 * x^{-i}) \geq f^i(\hat x^i, \bar x^{-i})\f] to the cut pool of that player.
 * @param player The player id
 * @param RHS The RHS value
 * @param xMinusI The input \f[ x^{-i} \f]
 * @return True if the cut was added
 */
bool Algorithms::IPG::CutAndPlay::addValueCut(unsigned int     player,
														double           RHS,
														const arma::vec &xMinusI) {



  //******DEBUG********
  // this->Players.at(player)->Incumbent.print("Incumbent of I");
  // xMinusI.print("xMinusI");
  // this->IPG->PlayersIP.at(player)->getc().print("cOfI");
  // this->IPG->PlayersIP.at(player)->getC().print_dense("CofI");
  //******DEBUG********


  arma::vec LHS = -(this->IPG->PlayersIP.at(player)->getc() +
						  this->IPG->PlayersIP.at(player)->getC() * xMinusI);

  // Constant!
  if (Utils::isEqual(arma::max(arma::abs(LHS)), 0)) {
	 LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::addValueCut: "
						 "Constant cut. Discarding. ";
	 return false;
  }

  //******DEBUG********
  // LHS.impl_raw_print("LHS with RHS=" + std::to_string(-RHS));
  //******DEBUG********
  Utils::normalizeIneq(LHS, RHS, false);
  //******DEBUG********
  // LHS.impl_raw_print("LHS with RHS=" + std::to_string(-RHS));
  //******DEBUG********

  if (Utils::nonzeroDecimals(RHS, 6) >= 5) {
	 LOG_S(0) << "Algorithms::IPG::CutAndPlay::addValueCut: "
					 "Numerically unstable. Generating another cut.";
	 if (this->IPG->Stats.AlgorithmData.CutAggressiveness.get() !=
		  Data::IPG::CutsAggressiveness::NotEvenTry) {
		if (this->externalCutGenerator(player, 1, false, true) != 0)
		  return true;
	 }

	 // Force normalization, in case it wasn't before.
	 RHS = Utils::round_nplaces(RHS, 6);
	 for (unsigned int i = 0; i < LHS.size(); ++i)
		LHS.at(i) = Utils::round_nplaces(LHS.at(i), 6);
	 Utils::normalizeIneq(LHS, RHS, true);
	 LOG_S(0) << "Algorithms::IPG::CutAndPlay::addValueCut: "
					 "WARNING: Cannot generate another cut. Adding normalized value-cut.";
  }
  //******DEBUG********
  // LHS.print("Value-Cut: LHS with RHS of" + std::to_string(-RHS));
  //******DEBUG********
  this->Cuts.at(0).second += 1;
  // Again, minus sign on RHS (the inequality is >=)
  return this->Players.at(player)->addCuts(arma::sp_mat{LHS.t()}, arma::vec{-RHS});
}

/**
 * @brief Checks if there is more time remaining.
 * @param remaining An output filled with the time remaining
 * @return True if there is still time left.
 */
bool Algorithms::IPG::CutAndPlay::checkTime(double &remaining) const {
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
 * @brief Initialize the LCP Objective with the quadratic and linear terms. These will be later used
 * if necessary
 */
void Algorithms::IPG::CutAndPlay::initLCPObjective() {

  this->LCP_Q.zeros(this->IPG->NumVariables, this->IPG->NumVariables);
  this->LCP_c.zeros(this->IPG->NumVariables);

  unsigned int varCounter = 0;
  for (unsigned int p = 0; p < this->IPG->NumPlayers; ++p) {
	 // Fill the c vector
	 unsigned int playerVars = this->IPG->PlayerVariables.at(p);
	 this->LCP_c.subvec(varCounter, varCounter + playerVars - 1) =
		  this->IPG->PlayersIP.at(p)->getc();


	 unsigned int otherVarsCounter = 0;
	 for (int o = 0; o < this->IPG->NumPlayers; ++o) {
		if (p != o) {
		  int startCol = std::accumulate(
				this->IPG->PlayerVariables.begin(), this->IPG->PlayerVariables.begin() + o, 0);
		  unsigned int otherVars = this->IPG->PlayerVariables.at(o);
		  this->LCP_Q.submat(
				varCounter, startCol, varCounter + playerVars - 1, startCol + otherVars - 1) =
				this->IPG->PlayersIP.at(p)->getC().submat(
					 0, otherVarsCounter, playerVars - 1, otherVarsCounter + otherVars - 1);
		  otherVarsCounter += otherVars;
		}
	 }
	 varCounter += playerVars;
  }

  //******DEBUG********
  // this->LCP_c.print("This is LCP_c");
  // this->LCP_Q.print_dense("This is LCP_Q");
  //******DEBUG********
}

/**
 * @brief Solves the IPG with the Equilibrium CutAndPlay algorithm.
 */
void Algorithms::IPG::CutAndPlay::solve() {
  this->initialize();
  bool MIPCuts  = true;
  auto cutLevel = this->IPG->Stats.AlgorithmData.CutAggressiveness.get();
  if (cutLevel == Data::IPG::CutsAggressiveness::NoThanks ||
		cutLevel == Data::IPG::CutsAggressiveness::NotEvenTry)
	 MIPCuts = false;
  int cutsAggressiveness = 0;
  if (MIPCuts) {
	 if (this->IPG->Stats.AlgorithmData.CutAggressiveness.get() ==
		  Data::IPG::CutsAggressiveness::Truculent)
		cutsAggressiveness = 3;
	 else
		cutsAggressiveness = 1;
  }

  if (this->Infeasible) {
	 this->IPG->Stats.Status.set(ZEROStatus::NashEqNotFound);
	 LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::solve: A Nash Equilibrium has not been "
						 "found. At least one of the players problem is infeasible.";
	 return;
  }

  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 if (MIPCuts) {
		LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::solve: Adding root cuts.";
		this->externalCutGenerator(i, 5, true, false);
	 }
  }

  // Main loop condition
  bool solved{false};
  // When an MNE is found via MIP, we start searching for the best one in the incumbent
  // approximation
  bool requireOptimality{this->IPG->Stats.AlgorithmData.Objective.get() !=
								 Data::IPG::Objectives::Feasibility};

  // Which players are feasible
  std::vector<int> feasible(this->IPG->NumPlayers, 0);
  // Number of mip cuts added
  std::vector<int> addedMIPCuts(this->IPG->NumPlayers, 0);
  int              Iteration = 0;

  while (!solved) {
	 ZEROStatus status;
	 // Increase the number of iterations
	 Iteration++;
	 LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::solve: Iteration ########### " << Iteration
					 << " ###########";
	 this->IPG->Stats.NumIterations.set(Iteration);


	 /* ************************************
	  * Check time and solve the LCP
	  **************************************/
	 if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
		double remaining;
		if (this->checkTime(remaining) && remaining > 0)
		  status = this->equilibriumLCP(remaining * 0.95, true, true);
		else
		  return;
	 } else
		status = this->equilibriumLCP(-1, true, true);


	 /* ************************************
	  * Numerical issues?
	  **************************************/

	 if (status == ZEROStatus::Numerical || this->IPG->Stats.Status.get() == ZEROStatus::Numerical) {
		this->IPG->Stats.AlgorithmData.Cuts.set(this->Cuts);
		LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::solve: Numerical errors.";
		this->IPG->Stats.Status.set(ZEROStatus::Numerical);
		return;
	 }

	 /* ************************************
	  * Support objects
	  **************************************/
	 // Is the player feasible?
	 std::fill(feasible.begin(), feasible.end(), 0);
	 // Did we add other cuts?
	 std::fill(addedMIPCuts.begin(), addedMIPCuts.end(), 0);
	 // How many cuts added in total?
	 unsigned int addedCuts = 0;
	 // Build xMinusIs
	 std::vector<arma::vec> xMinusI_s;
	 xMinusI_s = {};
	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i)
		xMinusI_s.push_back(this->buildXminusI(i));

	 if (status == ZEROStatus::NashEqFound) {
		/* ************************************
		 * Loop until at least one cut is added, or a solution is found.
		 **************************************/
		while (addedCuts <= 0) {
		  solved = true;

		  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
			 // Check only if the player hasn't proven to be feasible
			 if (feasible.at(i) == 0) {
				// Number of added cuts.
				int EO_cut = 0;
				/* ************************************
				 * Call the EO. Put in EO_cut the number of cuts
				 ***************************************/
				int EO = this->preEquilibriumOracle(
					 i, EO_cut, this->Players.at(i)->Incumbent, xMinusI_s.at(i));

				/* ************************************
				 * Numerical errors?
				 ***************************************/
				if (this->IPG->Stats.Status.get() == ZEROStatus::Numerical && EO == -1) {
				  this->IPG->Stats.AlgorithmData.Cuts.set(this->Cuts);
				  LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::solve: Numerical errors.";
				  return;
				}

				/* ************************************
				 * If we didn't have numerical problems
				 ***************************************/
				if (EO != 1) {
				  // Here we have either 0 or -2 (Infeasible or iteration limit)
				  solved = false;


				  if (EO == 0) {
					 /* ************************************
					  * Infeasible
					  ***************************************/
					 addedCuts += EO_cut; //+ this->separateCoinCuts(i, 5);
					 // if (Iteration > 5)
					 if (MIPCuts) {
						addedCuts += this->externalCutGenerator(
							 i,
							 (Iteration > 5 || (Iteration == 1 && cutsAggressiveness > 1))
								  ? cutsAggressiveness
								  : 1,
							 false,
							 false);
					 }
				  } else {
					 /* ************************************
					  * Iteration limit.
					  ***************************************/
					 // If we have cutting planes turned on, add more cuts.
					 if (MIPCuts && addedMIPCuts.at(i) == 0) {
						addedCuts += this->externalCutGenerator(i, cutsAggressiveness, false, false);
						addedMIPCuts.at(i) = 1;
					 }
				  }
				} else {
				  /* ************************************
					* Player is feasible
					***************************************/
				  feasible.at(i) = 1;
				}
			 }
		  }
		  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
			 double remaining;
			 if (!this->checkTime(remaining) || remaining <= 0) {
				return;
			 }
		  }
		  // If not solved but there are cuts, or it is solved and there are not cuts, exit.
		  if ((addedCuts > 0 && !solved) || (addedCuts == 0 && solved)) {
			 if (addedCuts == 0 && solved) {

				// We have solved the problem. Check for a better equilibrium?
				if (requireOptimality) {
				  LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::solve: Trying to improve the equilibrium.";
				  double remaining;
				  if (this->checkTime(remaining) && remaining > 0)
					 this->equilibriumLCP(remaining * 0.95, false, false);
				}

				this->Solved = true;
				bool pure    = true;
				for (unsigned int i = 0, counter = 0; i < this->IPG->NumPlayers; ++i) {
				  this->IPG->Solution.at(i) = this->Players.at(i)->Incumbent;
				  counter += this->IPG->PlayerVariables.at(i);
				  if (!this->Players.at(i)->Pure) {
					 pure = false;
				  }
				}
				this->Pure      = pure;
				arma::vec realX = this->xLast.subvec(0, this->LCP_c.size() - 1);
				this->IPG->SocialWelfare =
					 arma::as_scalar(this->LCP_c.t() * realX + realX.t() * this->LCP_Q * realX);
				this->IPG->Stats.AlgorithmData.Cuts.set(this->Cuts);
				LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::solve: A Nash Equilibrium has been found ("
								<< (pure == 0 ? "MNE" : "PNE") << ").";
				// Return if we do not need optimality
				this->IPG->Stats.Status.set(ZEROStatus::NashEqFound);
				return;

			 } else {
				// addedCuts > 0 && !solved
				// Not solved but cuts. So far nothing, but need to recompute the MNE. Break from here
			 }
		  }
		} // end while cuts
	 } else if (status == ZEROStatus::NashEqNotFound) {
		/* ************************************
		 * What if not solved? So far, it may be a condition triggered by the cutoff of the optimal
		 *equilibrium
		 *************************************/
		//@todo Implement for unbounded IPGs, where an eq may not exist at a given iteration
		this->IPG->Stats.AlgorithmData.Cuts.set(this->Cuts);
		LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::solve: No Nash Equilibrium has been found";
		this->Solved = false;
		this->IPG->Stats.Status.set(ZEROStatus::NashEqNotFound);
		return;
	 }

	 if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
		double remaining;
		if (!this->checkTime(remaining) || remaining <= 0)
		  return;
	 }
  }
}


/**
 * @brief Given the player id @p player, checks whether the current strategy is feasible or
 * not. In order to do so, a more complex separation technique may be called.
 * @param player The player id
 * @param addedCuts Filled with the number of added cuts
 * @param xOfI The strategy of @p player
 * @param xMinusI The strategy of the other players
 * @return 1 if feasible. 0 if infeasible. 2 if iteration limit was hit.
 */
int Algorithms::IPG::CutAndPlay::preEquilibriumOracle(const unsigned int player,
																  int &              addedCuts,
																  arma::vec &        xOfI,
																  arma::vec &        xMinusI) {

  ZEROAssert(player < this->IPG->NumPlayers);
  LOG_S(2) << "Algorithms::IPG::CutAndPlay::preEquilibriumOracle: (P" << player
			  << ") The oracle has been called. Preprocessing.";

  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
	 double remaining;
	 if (!this->checkTime(remaining) || remaining <= 0)
		return 0;
  }

  unsigned int Ny = this->IPG->PlayerVariables.at(player); // Equals to numVars by definition

  // Remember: the standard is minimization!

  // Update working strategies with "educated guesses"
  auto PureIP = this->IPG->PlayersIP.at(player)->getIPModel(xMinusI, false);
  //******DEBUG********
  // PureIP->write("PureIP.lp");
  //******DEBUG********
  PureIP->set(GRB_IntParam_SolutionLimit, 100);
  PureIP->optimize();
  int status = PureIP->get(GRB_IntAttr_Status);
  if (status == GRB_OPTIMAL) {
	 // Then, we have a best response

	 double IP_Objective  = PureIP->getObjective().getValue();
	 double REL_Objective = this->Players.at(player)->Payoff;

	 if (IP_Objective == GRB_INFINITY) {
		LOG_S(1) << "Algorithms::IPG::CutAndPlay::preEquilibriumOracle (P" << player
					<< ") Unbounded deviation.";
		addedCuts = false;
		return 0;
	 }


	 auto diff = REL_Objective - IP_Objective;
	 if (!Utils::isEqual(
				std::abs(diff), 0, 10 * this->IPG->Stats.AlgorithmData.DeviationTolerance.get())) {
		// There exists a difference between the payoffs

		if (diff > this->IPG->Stats.AlgorithmData.DeviationTolerance.get()) {
		  // This cannot happen!
		  this->IPG->Stats.Status.set(ZEROStatus::Numerical);
		  LOG_S(0) << "Algorithms::IPG::CutAndPlay::preEquilibriumOracle: |NUMERICAL WARNING| Invalid "
						  "payoff relation (better best response) of " +
								std::to_string(diff);
		  return -1;
		} else {
		  LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::preEquilibriumOracle:  (P" << player
						  << ") REL: " << REL_Objective << " vs IP: " << IP_Objective
						  << ". Adding a value cut.";
		  // Infeasible strategy. Add a value-cut

		  this->addValueCut(player, IP_Objective, xMinusI);
		  addedCuts = 1;
		  return 0;
		} // end abs(diff)
	 } else {

		// No discrepancy between payoffs

		int numSols = PureIP->get(GRB_IntAttr_SolCount);
		// Will be set to true if any pure-best response correspond to the current strategy
		bool equal = false;
		// Number of best responses found
		int bestResponses = 0;

		for (unsigned int s = 0; s < numSols; ++s) {

		  PureIP->set(GRB_IntParam_SolutionNumber, s);
		  // Check if the strategies are the same!
		  arma::vec bestResponse(Ny, arma::fill::zeros);
		  for (unsigned int k = 0; k < Ny; ++k)
			 bestResponse.at(k) = PureIP->getVarByName("y_" + std::to_string(k)).get(GRB_DoubleAttr_X);

		  // Add the strategy to the best-responses pool as a vertex
		  if (!this->Players.at(player)->addVertex(bestResponse, true))
			 bestResponses++;

		  if (Utils::isZero(xOfI - bestResponse, this->Tolerance)) {
			 this->Players.at(player)->Pure = true;
			 LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::preEquilibriumOracle:  (P" << player
							 << ") Feasible strategy (BR)";
			 equal = true;
		  }
		}

		LOG_S(2) << "Algorithms::IPG::CutAndPlay::preEquilibriumOracle: (P" << player << ") Found "
					<< bestResponses << " best responses vertices.";

		if (equal)
		  return 1;
		else {
		  // In this case, we need to call the proper oracle.
		  unsigned int iterations = this->IPG->PlayerVariables.at(player) * 5;
		  return this->equilibriumOracle(player, iterations, xOfI, xMinusI, addedCuts);
		}
	 }


  } else if (status == GRB_UNBOUNDED) {
	 LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::preEquilibriumOracle:  (P" << player
					 << ") The problem is unbounded.";
	 throw ZEROException(ZEROErrorCode::Numeric, "Unbounded best response.");
  }
  return 0;
}

/**
 * @brief Returns true if all players are playing a pure strategy in a Nash Equilibrium
 * @return True if the Equilibrium is pure
 */
bool Algorithms::IPG::CutAndPlay::isPureStrategy() const {
  if (this->Solved)
	 return this->Pure;
  else
	 return false;
}

/**
 * @brief The main Equilibrium CutAndPlay loop. Given a player, a maximum number of iterations, a
 * strategy and the other players strategy, it tries to determine if @p xOfI is feasible for @p
 * player.
 * @param player The player id
 * @param iterations The bound on iterations
 * @param xOfI The strategy of @p player
 * @param xMinusI The strategies of other players
 * @param addedCuts The number of added cuts
 * @return 1 if feasible. 0 if infeasible. 2 if iteration limit was hit.
 */
int Algorithms::IPG::CutAndPlay::equilibriumOracle(const unsigned int player,
															  const unsigned int iterations,
															  const arma::vec &  xOfI,
															  const arma::vec &  xMinusI,
															  int &              addedCuts) {

  ZEROAssert(player < this->IPG->NumPlayers);

  LOG_S(2) << "Algorithms::IPG::CutAndPlay::equilibriumOracle: (P" << player << ") Starting separator";



  // Store Membership LP outside the loop
  this->updateMembership(player, xOfI);
  auto   dualMembership = this->Players.at(player)->MembershipLP.get();
  GRBVar y[xOfI.size()]; // Dual membership variables
  for (unsigned int i = 0; i < xOfI.size(); i++)
	 y[i] = dualMembership->getVarByName("alpha_" + std::to_string(i));
  GRBVar betaVar = dualMembership->getVarByName("beta");

  // Update the objective for the membership. Avoid doing it every time
  // Update normalization
  GRBLinExpr expr = -betaVar;
  for (int j = 0; j < xOfI.size(); ++j)
	 expr += xOfI.at(j) * y[j];
  dualMembership->setObjective(expr, GRB_MAXIMIZE);


  //  Store the playerModel outside the loop
  std::unique_ptr<GRBModel> playerModel =
		std::unique_ptr<GRBModel>(this->IPG->PlayersIP.at(player)->getIPModel(xMinusI, false));
  GRBVar l[xOfI.size()]; // Dual membership variables
  for (unsigned int i = 0; i < xOfI.size(); i++)
	 l[i] = playerModel->getVarByName("y_" + std::to_string(i));
  playerModel->set(GRB_IntParam_OutputFlag, 0);
  playerModel->set(GRB_IntParam_SolutionLimit, 100);
  playerModel->set(GRB_IntParam_InfUnbdInfo, 1);
  // Max number of iterations
  for (int k = 0; k < iterations; ++k) {

	 // Check time at each iteration
	 if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
		double remaining;
		if (!this->checkTime(remaining) || remaining <= 0)
		  return 0;
	 }


	 // First, we check whether the point is a convex combination of feasible
	 // KNOWN points

	 // First iteration is out, since we have done it before.
	 if (k > 0)
		this->updateMembership(player, xOfI);

	 //******DEBUG********
	 // dualMembership->write("Convex_" + std::to_string(player) + ".lp");
	 //******DEBUG********

	 // Optimize the PRLP
	 dualMembership->optimize();
	 int membershipStatus = dualMembership->get(GRB_IntAttr_Status);


	 LOG_S(2) << "Algorithms::IPG::CutAndPlay::equilibriumOracle: (P" << player
				 << ") Membership status is " << membershipStatus;



	 // We should check the inequality on the player's model
	 ZEROAssert(membershipStatus == GRB_OPTIMAL);
	 double dualObj = dualMembership->getObjective().getValue();

	 // First: is the point redundant in the description?
	 // Namely, does xOfI belongs to the V-polyhedral approximation?
	 if (Utils::isEqual(dualObj, 0, this->Tolerance)) {
		LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::equilibriumOracle: (P" << player
						<< ") Feasible point. ";
		this->Players.at(player)->Feasible = true;


		arma::vec support;
		support.zeros(this->Players.at(player)->VertexCounter);
		for (unsigned int v = 0; v < this->Players.at(player)->VertexCounter; ++v) {
		  // abs to avoid misunderstanding with sign conventions
		  support.at(v) =
				dualMembership->getConstrByName("V_" + std::to_string(v)).get(GRB_DoubleAttr_Pi);
		}


		if (Utils::isEqual(support.max(), 1, this->Tolerance)) {
		  this->Players.at(player)->Pure = true;
		}
		//******DEBUG********
		// support.print("MNE Support: "+std::to_string(arma::sum(support)));
		//******DEBUG********


		return 1;
	 }


	 // Get the separating hyperplane

	 arma::vec alphaVal(xOfI.size(), arma::fill::zeros);
	 for (unsigned int i = 0; i < xOfI.size(); i++)
		alphaVal.at(i) = y[i].get(GRB_DoubleAttr_X);


	 //  Optimize the resulting inequality over the original feasible set
	 expr = 0;
	 for (unsigned int i = 0; i < xOfI.size(); ++i)
		expr += alphaVal.at(i) * l[i];


	 LOG_S(2) << "Algorithms::IPG::CutAndPlay::equilibriumOracle: (P" << player
				 << ") DualMembership has " << dualMembership->get(GRB_IntAttr_SolCount)
				 << " solution(s) with objective =" << dualObj << ".";


	 playerModel->setObjective(expr, GRB_MAXIMIZE);
	 playerModel->update();
	 playerModel->optimize();

	 int leaderStatus = playerModel->get(GRB_IntAttr_Status);
	 int numSols      = playerModel->get(GRB_IntAttr_SolCount);


	 ZEROAssert((leaderStatus == GRB_OPTIMAL) || (leaderStatus == GRB_SUBOPTIMAL) ||
					(leaderStatus == GRB_UNBOUNDED));
	 if (leaderStatus == GRB_OPTIMAL || (leaderStatus == GRB_SUBOPTIMAL && numSols > 0)) {

		//@todo <numSols or 1?
		for (int s = 0; s < numSols; ++s) {
		  playerModel->set(GRB_IntParam_SolutionNumber, s);

		  // The separating hyperplane plane evaluated at xOfI
		  double betaX = arma::as_scalar(alphaVal.t() * xOfI);
		  // The separating hyperplane evaluated at the optimal player's model
		  double betaPlayer = playerModel->getObjective().getValue();


		  // If the violation is negative: new vertex. Namely, the leader go further than xOfI
		  // If the violation is positive: cutting plane. Namely, xOfI is too far away
		  auto violation = betaX - betaPlayer;


		  //******DEBUG********
		  // xOfI.print("xOfI");
		  // alphaVal.print("alpha");
		  // dualMembership->write("Convex_" + std::to_string(player) + ".sol");
		  //******DEBUG********


		  LOG_S(2) << "Algorithms::IPG::CutAndPlay::equilibriumOracle: (P" << player << ") Violation of "
					  << violation;

		  if (violation >= this->Tolerance) {
			 // We have a cut.
			 // Ciao Moni

			 Utils::normalizeIneq(alphaVal, betaPlayer, true);

			 //******DEBUG********
			 // alphaVal.print("alphaVal with RHS of" + std::to_string(betaPlayer));
			 //******DEBUG********

			 this->Cuts.at(1).second += 1;
			 this->Players.at(player)->addCuts(arma::sp_mat{alphaVal.t()}, arma::vec{betaPlayer});


			 LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::equilibriumOracle: (P" << player
							 << ") Adding a cut.";
			 addedCuts = 1;
			 return 0;
		  } else {

			 // We found a new vertex
			 arma::vec v(this->Players.at(player)->V.n_cols, arma::fill::zeros);
			 for (unsigned int i = 0; i < v.size(); ++i)
				v[i] = l[i].get(GRB_DoubleAttr_X);
			 // v.print("vertex");
			 auto add = this->Players.at(player)->addVertex(v, numSols > 1);
			 LOG_S(2) << "Algorithms::IPG::CutAndPlay::equilibriumOracle: (P" << player
						 << ") adding vertex for Player (" << std::to_string(add) << "). "
						 << (iterations - k - 1) << "/" << iterations << " iterations left";
		  }
		}

	 } // status optimal for playerModel
	 else if (leaderStatus == GRB_UNBOUNDED) {
		// Get the extreme ray
		auto relaxed = playerModel->relax();
		relaxed.optimize();
		for (unsigned int i = 0; i < alphaVal.size(); ++i)
		  alphaVal[i] = relaxed.getVarByName("y_" + std::to_string(i)).get(GRB_DoubleAttr_UnbdRay);
		alphaVal = Utils::normalizeVec(alphaVal);
		LOG_S(2) << "Algorithms::IPG::CutAndPlay::equilibriumOracle: (P" << player << ") new ray";
		this->Players.at(player)->addRay(alphaVal);

	 } // status unbounded for playerModel
  }
  // Iteration limit
  return 2;
}


/**
 * @brief Update the Membership Linear Program for the given player and the verter @p vertex
 * @param player The player id
 * @param vertex  The input point to be checked
 */
void Algorithms::IPG::CutAndPlay::updateMembership(const unsigned int &player,
															  const arma::vec &   vertex) {
  ZEROAssert(player < this->IPG->NumPlayers);
  ZEROAssert(vertex.size() == this->IPG->PlayerVariables.at(player));
  MathOpt::getDualMembershipLP(this->Players.at(player)->MembershipLP,
										 this->Players.at(player)->VertexCounter,
										 this->Players.at(player)->V,
										 this->Players.at(player)->RayCounter,
										 this->Players.at(player)->R,
										 vertex,
										 this->Players.at(player)->containsOrigin);
}


/**
 * @brief Creates and solves the equilibrium LCP wrt the current game approximation
 * @param localTimeLimit A time limit for the computation
 * @param build If true, a new LCP will be built. Otherwise, the last one will be used.
 * @param firstSolution If true, the method will just seek for one solution.
 * @return The ZEROStatus for the equilibrium LCP
 */
ZEROStatus
Algorithms::IPG::CutAndPlay::equilibriumLCP(double localTimeLimit, bool build, bool firstSolution) {


  if (build) {
	 // Empty objects for market clearing
	 arma::sp_mat MC(0, this->IPG->NumVariables), dumA(0, this->IPG->NumVariables);
	 arma::vec    MCRHS(0, arma::fill::zeros), dumB(0, arma::fill::zeros);
	 // Downcast the IP_Param to MP_Param
	 std::vector<std::shared_ptr<MathOpt::MP_Param>> MPCasted;
	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
		auto m = std::dynamic_pointer_cast<MathOpt::MP_Param>(this->Players.at(i)->ParametrizedIP);
		MPCasted.push_back(m);
	 }
	 // Build the Nash Game
	 this->NashGame =
		  std::make_unique<Game::NashGame>(this->Env, MPCasted, MC, MCRHS, 0, dumA, dumB);
	 LOG_S(2) << "Algorithms::IPG::CutAndPlay::equilibriumLCP: Formulated the NashGame";
	 // Build the LCP from the Nash Game
	 this->LCP = std::make_unique<MathOpt::LCP>(this->Env, *this->NashGame);

	 // Record some statistics
	 this->IPG->Stats.NumVar         = this->LCP->getNumCols();
	 this->IPG->Stats.NumConstraints = this->LCP->getNumRows();
  }


  // Discriminate between Solver type and Objective type
  auto Solver        = this->IPG->getStatistics().AlgorithmData.LCPSolver.get();
  auto ObjectiveType = this->IPG->Stats.AlgorithmData.Objective.get();
  if (Solver == Data::LCP::Algorithms::PATH &&
		ObjectiveType != Data::IPG::Objectives::Feasibility) {
	 LOG_S(WARNING) << "Algorithms::IPG::CutAndPlay::equilibriumLCP: The LCP's objective is used only "
							 "for computing the payoff. "
							 "PATH does not support this optimization.";
  }
  switch (ObjectiveType) {
  case Data::IPG::Objectives::Linear: {
	 this->LCP->setMIPLinearObjective(this->LCP_c);
  } break;
  case Data::IPG::Objectives::Quadratic:
	 this->LCP->setMIPQuadraticObjective(this->LCP_c, this->LCP_Q);
	 break;
  default:
	 this->LCP->setMIPFeasibilityObjective();
  }


  arma::vec x, z;
  // Try to warm-start with the previous solution.
  x = this->xLast;
  z = this->zLast;
  // Set the value to the incumbent MNE payoff.
  double obj      = -GRB_INFINITY;
  int    solLimit = firstSolution ? 1 : GRB_MAXINT;

  int workers = 1;
  if (this->IPG->Stats.AlgorithmData.Threads.get() >= 8)
	 workers = (int)floor(this->IPG->Stats.AlgorithmData.Threads.get() / 3);
  auto LCPSolver = this->LCP->solve(Solver, x, z, localTimeLimit, workers, obj, solLimit);
  if (LCPSolver == ZEROStatus::NashEqFound) {
	 LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::equilibriumLCP: tentative Equilibrium found";

	 // Record the primal-dual solution and reset feasibility and pure-flag
	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
		this->Players.at(i)->Incumbent =
			 x.subvec(this->NashGame->getPrimalLoc(i), this->NashGame->getPrimalLoc(i + 1) - 1);
		this->Players.at(i)->DualIncumbent =
			 x.subvec(this->NashGame->getDualLoc(i), this->NashGame->getDualLoc(i + 1) - 1);
		this->Players.at(i)->Feasible = false;
		this->Players.at(i)->Pure     = false;
	 }

	 // Compute the payoffs
	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
		this->Players.at(i)->Payoff = this->IPG->PlayersIP.at(i)->computeObjective(
			 this->Players.at(i)->Incumbent, this->buildXminusI(i), false);
	 }

	 // Record the last responses
	 this->xLast   = x;
	 this->zLast   = z;
	 this->objLast = obj;
	 return ZEROStatus::NashEqFound;

  } else if (LCPSolver == ZEROStatus::Numerical) {
	 LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::equilibriumLCP: Numerical errors.";
	 return ZEROStatus::Numerical;
  } else if (LCPSolver == ZEROStatus::TimeLimit) {
	 if (this->IPG->Stats.Status.get() == ZEROStatus::Uninitialized)
		this->IPG->Stats.Status.set(ZEROStatus::TimeLimit);
	 return ZEROStatus::TimeLimit;
  } else {
	 LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::equilibriumLCP: No tentative equilibrium found";
	 return ZEROStatus::NashEqNotFound;
  }
}

/**
 * @brief This method initializes some fields for the algorithm. Also, it warm starts the
 * initial strategies to pure best responses.
 */
void Algorithms::IPG::CutAndPlay::initialize() {
  // Set the number of iterations to zero.
  this->IPG->Stats.NumIterations.set(0);

  // Create the IPG_Player-s objects
  this->Players = std::vector<std::unique_ptr<IPG_Player>>(this->IPG->NumPlayers);
  // Initialize the working objects with respective values
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 // Initialize the IPG_Player
	 this->Players.at(i) =
		  std::make_unique<IPG_Player>(this->IPG->PlayersIP.at(i)->getNumVars(), this->Tolerance);
	 //@todo be aware of variables' changes in presolve
	 if (this->IPG->Stats.AlgorithmData.Presolve.get())
	   this->IPG->PlayersIP.at(i)->presolve();
	 // Add the working IP
	 auto WorkingIP                      = new MathOpt::IP_Param(*this->IPG->PlayersIP.at(i).get());
	 this->Players.at(i)->ParametrizedIP = std::make_shared<MathOpt::IP_Param>(*WorkingIP);
	 // Add the working MembershipLP
	 this->Players.at(i)->MembershipLP = std::move(std::make_unique<GRBModel>(*this->Env));
	 this->initializeCoinModel(i);
  }


  // Initialize some "educated" guesses as best responses.
  this->initializeEducatedGuesses();
  // Reset cuts statistics
  this->Cuts = {std::pair<std::string, int>("Value", 0),
					 std::pair<std::string, int>("VPoly", 0),
					 std::pair<std::string, int>("MIPCuts", 0)};
  this->IPG->Stats.AlgorithmData.Cuts.set(this->Cuts);
  // Initialize the LCP objectives for further use.
  this->initLCPObjective();
}


/**
 * @brief Given the player id @p i, builds the vector x^{-i} from the current working
 * strategies.
 * @param i The player id
 * @return The other players strategies (except @p i)
 */
arma::vec Algorithms::IPG::CutAndPlay::buildXminusI(const unsigned int i) {
  ZEROAssert(i < this->IPG->NumPlayers);
  arma::vec xMinusI(this->IPG->NumVariables - this->IPG->PlayerVariables.at(i), arma::fill::zeros);
  unsigned int counter = 0;
  for (unsigned int j = 0; j < this->IPG->NumPlayers; ++j) {
	 if (i != j) {
		xMinusI.subvec(counter, counter + this->IPG->PlayerVariables.at(j) - 1) =
			 this->Players.at(j)->Incumbent;
		counter += this->IPG->PlayerVariables.at(j);
	 }
  }
  return xMinusI;
}
/**
 * @brief Given a player @p player, a number of maximum cuts to generate @p maxcuts and a bool
 * @p rootNode, this method generates some valid inequalities for the @p player 's integer program.
 * This method uses Coin-OR CGL. So far, Knapsack covers, GMI and MIR inequalities are used. Also,
 * note that there is no recursive cut generation (meaning, we do not generate cuts from a
 * previous cut pool) as to better manage numerical stability. @p cutOff requires to cut off the
 * current solution for @p player.
 * @param player The current player id
 * @param maxCuts The maximum number of cuts
 * @param rootNode True if the cut generation happens at the root node
 * @param cutOff True if the cuts are required to cutoff the current solution
 * @return The number of added cuts
 */
unsigned int Algorithms::IPG::CutAndPlay::externalCutGenerator(unsigned int player,
																			  int          maxCuts,
																			  bool         rootNode,
																			  bool         cutOff) {

  ZEROAssert(player < this->IPG->NumPlayers);
  auto      xOfI     = this->Players.at(player)->Incumbent;
  auto      xOfIDual = this->Players.at(player)->DualIncumbent;
  auto      xMinusI  = this->buildXminusI(player);
  arma::vec objective;

  objective = (this->Players.at(player)->ParametrizedIP->getC() * xMinusI) +
				  this->Players.at(player)->ParametrizedIP->getc();


  auto CoinModel = this->Players.at(player)->CoinModel;

  auto numVars = this->Players.at(player)->ParametrizedIP->getB(false).n_cols;
  /******RECURSIVE CUT GENERATION********
  Uncomment this to enable recursive cut generation. In general, not a good idea for numerical
  stability auto numConstrs = this->Players.at(player)->ParametrizedIP->getB(false).n_rows;
  *******RECURSIVE CUT GENERATION********/
  auto numConstrs = this->IPG->PlayersIP.at(player)->getB(false).n_rows;


  // Empty objects, for root node cuts
  if (rootNode) {
	 xOfI.zeros(numVars);
	 xOfIDual.zeros(numConstrs);
  }

  auto primal = new double[numVars];
  auto dual   = new double[numConstrs];
  auto c      = new double[numVars];

  for (unsigned int i = 0; i < numVars; ++i) {
	 c[i] = objective.at(i);
	 if (!rootNode)
		primal[i] = xOfI.at(i);
  }
  for (unsigned int i = 0; i < numConstrs; ++i) {
	 if (!rootNode)
		dual[i] = xOfIDual.at(i);
  }


  CoinModel->setObjective(c);
  if (!rootNode) {
	 CoinModel->setColSolution(primal);
	 CoinModel->setRowPrice(dual);
  }


  //******DEBUG********
  // CoinModel->writeLp("CoinModel");
  // this->IPG->PlayersIP.at(player)->getIPModel(xMinusI)->write("GurobiModel.lp");
  //******DEBUG********

  try {
	 CoinModel->solveFromSol();

	 if (!rootNode) {
		auto solCheck = CoinModel->getColSolution();
		for (unsigned int i = 0; i < numVars; ++i) {
		  if (!Utils::isEqual(primal[i], solCheck[i]))
			 throw;
		}


		solCheck = CoinModel->getRowPrice();
		for (unsigned int i = 0; i < numConstrs; ++i) {
		  if (!Utils::isEqual(dual[i], solCheck[i]))
			 throw;
		}
	 }

	 auto        candidateCuts = new OsiCuts;
	 CglTreeInfo info          = CglTreeInfo();
	 info.inTree               = false;
	 info.options              = 4;
	 info.pass                 = 0;



	 CglKnapsackCover kpGen;
	 auto             KPs = new OsiCuts;
	 kpGen.setGlobalCuts(true);
	 kpGen.setAggressiveness(100);
	 kpGen.switchOnExpensive();
	 kpGen.setMaxInKnapsack(xOfI.size() + 10);
	 kpGen.generateCuts(*CoinModel, *KPs, info);


	 for (int(i) = 0; (i) < KPs->sizeCuts(); ++(i))
		if (KPs->rowCut(i).globallyValid())
		  candidateCuts->insert(KPs->rowCut(i));


	 CglMixedIntegerRounding2 MIRGen;
	 auto                     MIRs = new OsiCuts;
	 MIRGen.setGlobalCuts(true);
	 MIRGen.setAggressiveness(100);
	 MIRGen.setDoPreproc(1);
	 MIRGen.setMAXAGGR_(100);
	 MIRGen.generateCuts(*CoinModel, *MIRs, info);


	 for (int(i) = 0; (i) < MIRs->sizeCuts(); ++(i))
		if (MIRs->rowCut(i).globallyValid())
		  candidateCuts->insert(MIRs->rowCut(i));

	 CglGMI GMIGen;
	 auto   GMIs = new OsiCuts;
	 if (!cutOff) {
		GMIGen.getParam().setMAX_SUPPORT(numVars);
		GMIGen.getParam().setMAX_SUPPORT_REL(0.8);
		GMIGen.getParam().setMAXDYN(CoinModel->getInfinity());
	 }
	 GMIGen.getParam().setENFORCE_SCALING(true);
	 GMIGen.setGlobalCuts(true);
	 GMIGen.setAggressiveness(100);
	 if (cutOff) {
		GMIGen.getParam().setAway(1e-4);
		GMIGen.getParam().setMINVIOL(1e-3);
	 }
	 GMIGen.generateCuts(*CoinModel, *GMIs, info);

	 for (int(i) = 0; (i) < GMIs->sizeCuts(); ++(i))
		if (GMIs->rowCut(i).globallyValid())
		  candidateCuts->insert(GMIs->rowCut(i));


	 //******DEBUG********
	 // KPs->printCuts();
	 // GMIs->printCuts();
	 // MIRs->printCuts();
	 //******DEBUG********

	 // Sort by effectiveness
	 candidateCuts->sort();


	 if (cutOff) {
		// Try also simple gomory, and verify the cutoff

		CglGomory gomoryGen;
		auto      GMs = new OsiCuts;
		gomoryGen.setAggressiveness(100);
		gomoryGen.setLimit(xOfI.size() + 10);
		gomoryGen.setLimitAtRoot(xOfI.size() + 10);
		gomoryGen.setAway(1e-3);
		gomoryGen.generateCuts(*CoinModel, *GMs, info);


		for (int(i) = 0; (i) < GMs->sizeCuts(); ++(i))
		  if (GMs->rowCut(i).globallyValid())
			 candidateCuts->insert(GMs->rowCut(i));

		candidateCuts->sort();
		// We need to be sure to get only the cuts that are actively cutting off the solution.
		for (int i = 0; i < candidateCuts->sizeCuts(); ++i) {
		  // Check if we have a violation. Otherwise, erase the cut.
		  double violation = candidateCuts->rowCut(i).violated(primal);
		  if (violation <= 0) {
			 LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::externalCutGenerator: (P" << player
							 << ") Discarding a cut (no violation).";
			 candidateCuts->eraseRowCut(i);
		  }
		}
	 }


	 // Minimum among the candidate cuts and the maximum number of cuts
	 auto numCuts = std::min(candidateCuts->sizeCuts(), maxCuts);

	 if (numCuts > 0) {
		// The final cut matrix and RHS.
		arma::sp_mat LHS;
		LHS.zeros(numCuts, numVars);
		arma::vec RHS(numCuts, arma::fill::zeros);

		// Iterate over the candidate cuts (in max numCuts iterations), sorted by efficacy
		for (int i = 0; i < numCuts; ++i) {
		  // Get the cut
		  auto cut = candidateCuts->rowCut(i);

		  // Get the row from the cut
		  auto row = cut.row();
		  // Get the indices and the elements
		  auto indices  = row.getIndices();
		  auto elements = row.getElements();

		  // Iterate over the elements to get the index-value pair
		  for (int j = 0; j < row.getNumElements(); j++)
			 LHS.at(i, indices[j]) = elements[j];


		  // Discern among the cut senses
		  auto sense = cut.sense();
		  switch (cut.sense()) {
		  case 'E': {
			 // Equality. Two cuts
			 LHS.resize(LHS.n_rows + 1, LHS.n_cols);
			 RHS.resize(RHS.n_rows + 1);
			 LHS.row(LHS.n_rows - 1) = -LHS.row(i);
			 RHS.at(RHS.n_rows - 1)  = -cut.rhs();
			 RHS.at(i)               = cut.rhs();
		  } break;
		  case 'G': {
			 // Switch signs
			 RHS.at(i) = -RHS.at(i);
			 RHS.at(i) = -cut.rhs();
		  } break;
		  case 'L': {
			 // Just complete the RHS
			 RHS.at(i) = cut.rhs();
		  } break;
		  case 'R': {
			 // Ranged inequality. We have both upper and lower bound
			 RHS.at(i) = cut.rhs();
		  } break;
		  default: {
			 candidateCuts->printCuts();
			 throw ZEROException(ZEROErrorCode::InvalidData, "Unknown cut sense.");
		  }
		  }
		}

		// If we got some new cuts, add them to the working integer program
		// Add the cuts
		this->Players.at(player)->addCuts(LHS, RHS);
		// Update the statistics
		this->Cuts.at(2).second += numCuts;
		LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::externalCutGenerator: (P" << player << ") Added "
						<< numCuts << " and generated: " << MIRs->sizeCuts() << "  MIRs  - "
						<< KPs->sizeCuts() << "  KPs - " << GMIs->sizeCuts() << "  GMIs.";
	 } else {
		LOG_S(INFO) << "Algorithms::IPG::CutAndPlay::externalCutGenerator: (P" << player
						<< ") No cuts added.";
	 }
	 return numCuts;
  } catch (CoinError &e) {
	 throw ZEROException(ZEROErrorCode::SolverError,
								"Invalid Coin-OR interface response: " + e.message());
  }
}

/**
 * @brief Initializes some pure-strategies for each player.
 */
void Algorithms::IPG::CutAndPlay::initializeEducatedGuesses() {


  // Reset the working strategies to a pure strategy given by the IP
  // Push back the IP_Param copies in WorkingIPs

  // For any player
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {

	 // Equals to NumVarsI by definition
	 unsigned int NumVarsI = this->IPG->PlayerVariables.at(i);
	 // xMinusI
	 arma::vec xMinusI = this->buildXminusI(i);
	 // Update working strategies with "educated guesses"
	 auto PureIP = this->IPG->PlayersIP.at(i)->getIPModel(xMinusI, false);
	 PureIP->optimize();
	 int status = PureIP->get(GRB_IntAttr_Status);
	 if (status == GRB_INFEASIBLE) {
		// Game ended, player is infeasible
		this->IPG->Stats.Status.set(ZEROStatus::NashEqNotFound);
		this->Infeasible = true;
		return;
	 } else {
		// Model is not infeasible. We can have either a ray or a vertex
		if (status == GRB_OPTIMAL) {
		  // We have a vertex
		  for (unsigned int k = 0; k < NumVarsI; ++k)
			 this->Players.at(i)->Incumbent.at(k) =
				  PureIP->getVarByName("y_" + std::to_string(k)).get(GRB_DoubleAttr_X);

		  this->Players.at(i)->addVertex(this->Players.at(i)->Incumbent, false);
		  LOG_S(2) << "Algorithms::IPG::CutAndPlay::initializeEducatedGuesses: (P" << i
					  << ") Adding a vertex.";
		}

		else if (status == GRB_UNBOUNDED) {
		  // Relax and get the ray
		  GRBModel relaxed = PureIP->relax();
		  relaxed.set(GRB_IntParam_InfUnbdInfo, 1);
		  relaxed.set(GRB_IntParam_DualReductions, 0);
		  relaxed.optimize();
		  arma::vec ray(NumVarsI, arma::fill::zeros);
		  for (unsigned int k = 0; k < NumVarsI; ++k) {
			 ray.at(k) = relaxed.getVarByName("y_" + std::to_string(k)).get(GRB_DoubleAttr_UnbdRay);

			 this->Players.at(i)->addRay(ray, false);
			 LOG_S(2) << "Algorithms::IPG::CutAndPlay::initializeEducatedGuesses():  (P" << i
						 << ") Adding a ray";
		  }
		}

		// Check if the origin is feasible
		arma::vec zeros(this->IPG->PlayersIP.at(i)->getNumVars(), arma::fill::zeros);
		double    zero = this->IPG->PlayersIP.at(i)->computeObjective(zeros, xMinusI, true);
		if (zero != GRB_INFINITY)
		  this->Players.at(i)->containsOrigin = true;
	 }
  }
}

/**
 * @brief This method builds the Coin-OR model used in CutAndPlay::externalCutGenerator for the given
 * player
 * @param player The player's id
 */
void Algorithms::IPG::CutAndPlay::initializeCoinModel(const unsigned int player) {

  ZEROAssert(player < this->IPG->NumPlayers);
  // Source the main ingredients. Avoid getting B with bounds, since we already have the raw bounds.
  auto IP_B          = this->Players.at(player)->ParametrizedIP->getB(false);
  auto IP_Bounds     = this->Players.at(player)->ParametrizedIP->getBounds();
  auto IP_b          = this->Players.at(player)->ParametrizedIP->getb(false);
  auto IP_Integers   = this->Players.at(player)->ParametrizedIP->getIntegers();
  auto IP_numVars    = this->IPG->PlayersIP.at(player)->getNumVars();
  auto IP_numConstrs = IP_B.n_rows;


  // Convert B
  auto B = Utils::armaToCoinSparse(IP_B);
  // Double objects
  auto c  = new double[IP_numVars];
  auto b  = new double[IP_numConstrs];
  auto lb = new double[IP_numVars];
  auto ub = new double[IP_numVars];

  // Filling stage
  for (unsigned int i = 0; i < IP_numVars; ++i) {
	 c[i]  = 0;
	 lb[i] = IP_Bounds.at(i).first > 0 ? IP_Bounds.at(i).first : 0;
	 ub[i] = IP_Bounds.at(i).second >= 0 ? IP_Bounds.at(i).second : GRB_INFINITY;
  }
  for (unsigned int i = 0; i < IP_numConstrs; ++i)
	 b[i] = IP_b.at(i);

  // Solver interface
  auto CoinModel = new OsiGrbSolverInterface();
  // Remember to use the same number of threads...
  GRBsetintparam(CoinModel->getEnvironmentPtr(), "Threads", this->Env->get(GRB_IntParam_Threads));
  // Load the problem from the given objects
  CoinModel->loadProblem(B, lb, ub, c, nullptr, b);
  // Set the integer variables
  for (unsigned int i = 0; i < IP_Integers.size(); ++i)
	 CoinModel->setInteger(IP_Integers.at(i));


  // Reset the log level to zero
  CoinModel->messageHandler()->setLogLevel(0);
  // Move to the target Player
  this->Players.at(player)->CoinModel = std::make_shared<OsiGrbSolverInterface>(*CoinModel);
}