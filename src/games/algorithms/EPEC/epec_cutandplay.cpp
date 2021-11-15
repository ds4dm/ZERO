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


#include "games/algorithms/EPEC/epec_cutandplay.h"

#include <memory>


/**
 * @brief Overrides Algorithms::EPEC::PolyBase::isSolved with a custom method.
 * @param tol Numerical tolerance. Currently not useful
 * @return True if the current outer approximation solution is feasible (then, it is solved)
 */
bool Algorithms::EPEC::CutAndPlay::isSolved(double tol) {
  if (this->Feasible)
	 return true;
  else {
	 bool cuts;
	 return this->isFeasible(cuts);
  }
}


/**
 * @brief Checks whether the current outer approximation equilibrium is feasible and solves the
 * problem. Otherwise, it adds cuts or generate useful information for the next iterations
 * @param addedCuts [out] is true if at least a cut has been added
 * @return
 */
bool Algorithms::EPEC::CutAndPlay::isFeasible(bool &addedCuts) {


  // First, we have a NE from Games::computeNashEq
  if (!this->EPECObject->NashEquilibrium)
	 return false;


  // The returned result
  bool result = true;

  // The best response object. filled for every players
  arma::vec bestResponse;
  // Compute payoffs in the current solution
  arma::vec incumbentPayoffs =
		this->EPECObject->TheNashGame->computeQPObjectiveValues(this->EPECObject->SolutionX, true);
  // For any player
  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i) {
	 // Reset the feasibility
	 this->Trees.at(i)->resetFeasibility();
	 // Compute the payoff of the best response, as well as the best response
	 auto bestModel = this->EPECObject->bestResponseProgram(
		  i, this->EPECObject->SolutionX, this->Trees.at(i)->OriginalLCP.get());
	 const int status = bestModel->get(GRB_IntAttr_Status);

	 if (status == GRB_UNBOUNDED) {
		LOG_S(1) << "Algorithms::EPEC::CutAndPlay::isFeasible (P" << i
					<< ") Unbounded deviation";
		addedCuts = false;
		result    = false;
	 }

	 // Get the best response
	 unsigned int Nx = this->EPECObject->PlayersLCP.at(i)->getNumCols();
	 bestResponse.zeros(Nx);
	 for (unsigned int j = 0; j < Nx; ++j)
		bestResponse.at(j) = bestModel->getVarByName("x_" + std::to_string(j)).get(GRB_DoubleAttr_X);
	 double bestPayoff = bestModel->getObjective().getValue();


	 LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::isFeasible (P" << i << ") Payoff of "
					 << incumbentPayoffs.at(i) << " vs bestResponse of " << bestPayoff;
	 // Since it's minimization, difference between the incumbent and best response payoff
	 double absdiff = incumbentPayoffs.at(i) - bestPayoff;



	 if (!Utils::isEqual(incumbentPayoffs.at(i), bestPayoff, this->Tolerance, 1 - this->Tolerance)) {
		// Discrepancy between payoffs! Need to investigate.


		if (absdiff > 10 * this->Tolerance) {
		  // It means the current payoff is more than then optimal response. Then
		  // this is not a best response. Theoretically, this cannot happen from
		  // an outer approximation. This can however happen for numerical reasons

		  LOG_S(WARNING) << "Algorithms::EPEC::CutAndPlay::isFeasible (P" << i << ")"
							  << " No best response (" << incumbentPayoffs.at(i) << " vs " << bestPayoff
							  << " with absdiff=" << incumbentPayoffs.at(i) - bestPayoff << ")";
		  this->EPECObject->Stats.Status.set(ZEROStatus::Numerical);
		  ZEROException(ZEROErrorCode::Numeric, "Invalid payoffs relation (better best response).");
		  return 0;
		} else
		// In any other case, we are good to go.
		{
		  // It means the current payoff is less than the optimal response. The
		  // approximation is not good, and this point is infeasible. Then, we can
		  // generate a value-cut
		  arma::vec xMinusI;
		  this->EPECObject->getXMinusI(this->EPECObject->SolutionX, i, xMinusI);
		  this->addValueCut(i, bestPayoff, xMinusI);
		  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::isFeasible (P" << i << ") Value cut";
		  result = false;
		}
	 } else {

		// Here we have a best response whose payoff coincides with the one of the
		// equilibrium. The strategy might not be feasible, though.

		// Get xOfI
		arma::vec xOfI;
		this->EPECObject->getXofI(this->EPECObject->SolutionX, i, xOfI, false);

		// Check if we need to add the point to the vertex storage.
		arma::vec vertex = bestResponse.subvec(0, xOfI.size() - 1);
		// Add the best response to the storage
		if (this->Trees.at(i)->addVertex(vertex, true)) {
		  LOG_S(1) << "Algorithms::EPEC::CutAndPlay::isFeasible (P" << i
					  << ") Adding vertex (Best Response)";
		} else {
		  LOG_S(1) << "Algorithms::EPEC::CutAndPlay::isFeasible (P" << i
					  << ") Already known vertex (Best Response)";
		}

		if (!Utils::isZero(xOfI - vertex, this->Tolerance)) {
		  // Then, if the answers aren't equal, we need to refine the
		  // approximation or determine if this strategy is anyhow feasible.
		  // We search for a convex combination of best responses so that we can
		  // certify the answer is inside the convex-hull (or not).

		  int budget = this->EPECObject->PlayersQP.at(i)->getNumVars();
		  if (!this->equilibriumOracle(xOfI, this->EPECObject->SolutionX, i, budget, addedCuts)) {
			 LOG_S(1) << "Algorithms::EPEC::CutAndPlay::isFeasible (P" << i
						 << ") CutAndPlay says NO.";
			 result = false;
		  }
		  // Otherwise, the result is true...

		} else {
		  // Feasible point
		  this->Trees.at(i)->setFeasible();
		  this->Trees.at(i)->setPure();
		  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::isFeasible (P" << i
						  << ") Feasible strategy (Best Response)";
		  return true;
		}
	 }
  }
  return result;
}


/**
 * @brief Updates the membership linear-program in the relative
 * Algorithms::EPEC::CutAndPlay::Trees for the player @p player
 * @param player The player index
 * @param xOfI  The point for which the membership LP should be updated for
 */
void Algorithms::EPEC::CutAndPlay::updateMembership(const unsigned int &player,
																				const arma::vec &   xOfI) {


  auto PlayerTree = Trees.at(player);
  MathOpt::getDualMembershipLP(PlayerTree->MembershipLP,
										 PlayerTree->VertexCounter,
										 PlayerTree->V,
										 PlayerTree->RayCounter,
										 PlayerTree->R,
										 xOfI,
										 PlayerTree->containsOrigin);
}



/**
 * @brief The main Equilibrium CutAndPlay loop. Given a player, a maximum number of iterations, a
 * strategy and the other players strategy, it tries to determine if @p xOfI is feasible for @p
 * player.
 * @param xOfI The incumbent strategy for @p player
 * @param x The full solution vector
 * @param player The player id
 * @param budget A bound on the number of iteration
 * @param addedCuts The number of added cuts
 * @return 1 if feasible. 0 if infeasible. 2 if iteration limit was hit.
 */

bool Algorithms::EPEC::CutAndPlay::equilibriumOracle(
	 arma::vec &xOfI, arma::vec &x, unsigned int player, int budget, bool &addedCuts) {


  // Store the leaderModel outside the loop
  auto leaderModel =
		this->EPECObject->bestResponseProgram(player, x, this->Trees.at(player)->OriginalLCP.get());
  GRBVar l[xOfI.size()]; // Dual membership variables
  for (unsigned int i = 0; i < xOfI.size(); i++)
	 l[i] = leaderModel->getVarByName("x_" + std::to_string(i));
  leaderModel->set(GRB_IntParam_InfUnbdInfo, 1);
  leaderModel->set(GRB_IntParam_DualReductions, 0);
  leaderModel->set(GRB_IntParam_OutputFlag, 0);
  leaderModel->set(GRB_IntParam_SolutionLimit, 100);


  // Store Membership LP outside the loop
  this->updateMembership(player, xOfI);
  auto   dualMembershipModel = this->Trees.at(player)->MembershipLP.get();
  GRBVar y[xOfI.size()]; // Dual membership variables
  for (unsigned int i = 0; i < xOfI.size(); i++)
	 y[i] = dualMembershipModel->getVarByName("alpha_" + std::to_string(i));
  GRBVar betaVar = dualMembershipModel->getVarByName("beta");
  // Update the objective for the membership. Avoid doing it every time
  // Update normalization
  GRBLinExpr expr = -betaVar;
  for (int j = 0; j < xOfI.size(); ++j)
	 expr += xOfI.at(j) * y[j];
  dualMembershipModel->setObjective(expr, GRB_MAXIMIZE);



  for (int k = 0; k < budget; ++k) {
	 // First, we check whether the point is a convex combination of feasible
	 // KNOWN points
	 auto V = this->Trees.at(player)->V;

	 // Avoid for the firs time. We already have it
	 if (k > 0)
		this->updateMembership(player, xOfI);


	 dualMembershipModel->optimize();

	 int status = dualMembershipModel->get(GRB_IntAttr_Status);

	 LOG_S(1) << "Algorithms::EPEC::CutAndPlay::equilibriumOracle: (P" << player << ")"
				 << " MembershipLP status is " << status
				 << " (Objective=" << dualMembershipModel->getObjective().getValue() << ")";

	 if (status == GRB_OPTIMAL) {
		double dualObj = dualMembershipModel->getObjective().getValue();

		// Maximization of the dual membership gives zero. Then, the point is feasible
		if (Utils::isEqual(dualObj, 0, this->Tolerance)) {

		  arma::vec sol(xOfI.size(), arma::fill::zeros);

		  for (unsigned int i = 0; i < xOfI.size(); i++)
			 sol.at(i) = y[i].get(GRB_DoubleAttr_X);

		  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::equilibriumOracle: (P" << player
						  << ")"
						  << " The point is a convex combination of known points!";

		  this->Trees.at(player)->setFeasible();

		  arma::vec support, rays;
		  support.zeros(this->Trees.at(player)->getVertexCount());
		  rays.zeros(this->Trees.at(player)->getRayCount());
		  for (unsigned int v = 0; v < this->Trees.at(player)->getVertexCount(); ++v) {
			 support.at(v) =
				  dualMembershipModel->getConstrByName("V_" + std::to_string(v)).get(GRB_DoubleAttr_Pi);
		  }
		  bool flag = false;
		  for (unsigned int r = 0; r < this->Trees.at(player)->getRayCount(); ++r) {
			 rays.at(r) =
				  dualMembershipModel->getConstrByName("R_" + std::to_string(r)).get(GRB_DoubleAttr_Pi);
			 //******DEBUG********
			 /*if (rays.at(r) > this->Tolerance) {

				LOG_S(WARNING) << "Algorithms::EPEC::CutAndPlay::equilibriumOracle: (P"
									<< player << ")"
									<< " Ray " << r << " has a positive coefficient.";
				flag = true;
			 }
			  */
			 //******DEBUG********
		  }

		  //******DEBUG********
		  // auto test = arma::accu(support);
		  // support.print("support vertices" + std::to_string(test));
		  // rays.print("support rays");
		  //******DEBUG********

		  // Check whether this is pure or not.
		  //******DEBUG********
		  // dualMembershipModel->write("Membership.lp");
		  // this->Trees.at(player)->V.impl_print_dense("V");
		  // this->Trees.at(player)->R.impl_print_dense("R");
		  // xOfI.print("xOfI");
		  //******DEBUG********

		  // assert(arma::sum(support) > 1 - 1e-3);
		  if (Utils::isEqual(support.max(), 1, this->Tolerance) &&
				Utils::isEqual(arma::abs(rays).max(), this->Tolerance * 10)) {
			 this->Trees.at(player)->setPure();
		  } else {
			 if (this->isFeasiblePure(player, xOfI)) {
				this->Trees.at(player)->setPure();
				LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::equilibriumOracle: (P" << player
								<< ")"
								<< " Pure strategy.";
			 } else {
				LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::equilibriumOracle: (P" << player
								<< ")"
								<< " Mixed strategy.";
			 }
		  }
		  return true;

		} else {

		  // This is not a convex combination! The dual membership objective is not zero


		  // Get the Farkas' in the form of the unbounded ray of the dual of the
		  // dualMembershipLP (the primal)
		  LOG_S(1) << "Algorithms::EPEC::CutAndPlay::equilibriumOracle: (P" << player
					  << ")"
						  " The point is NOT a convex combination of known points!";
		  arma::vec cutLHS(xOfI.size(), arma::fill::zeros);

		  // Load the Farkas' certificate
		  for (unsigned int i = 0; i < xOfI.size(); i++)
			 cutLHS.at(i) = y[i].get(GRB_DoubleAttr_X);

		  // Optimize the resulting inequality over the original feasible set
		  expr = 0;
		  for (unsigned int i = 0; i < xOfI.size(); ++i)
			 expr += cutLHS.at(i) * l[i];


		  leaderModel->setObjective(expr, GRB_MAXIMIZE);
		  leaderModel->update();
		  leaderModel->optimize();
		  status = leaderModel->get(GRB_IntAttr_Status);


		  // If the status is optimal, either a vertex or a cut
		  if (status == GRB_OPTIMAL ||
				(status == GRB_SUBOPTIMAL && leaderModel->get(GRB_IntAttr_SolCount) > 0)) {

			 // For any solution found
			 int numSols = leaderModel->get(GRB_IntAttr_SolCount);

			 for (int s = 0; s < numSols; ++s) {

				leaderModel->set(GRB_IntParam_SolutionNumber, s);
				double betaPlayer = leaderModel->getObjective().getValue();

				LOG_S(1) << "Algorithms::EPEC::CutAndPlay::equilibriumOracle: (P" << player
							<< ")"
								" LeaderModel status = "
							<< std::to_string(status) << " with objective=" << betaPlayer;

				double betaXofI = arma::as_scalar(cutLHS.t() * xOfI); // c^T xOfI

				LOG_S(1) << "Algorithms::EPEC::CutAndPlay::equilibriumOracle: (P" << player
							<< ") c^Tv=" << betaPlayer << " -- c^TxOfI=" << betaXofI;
				double violation = betaXofI - betaPlayer;

				if (violation >= this->Tolerance) {
				  // False, but we have a cut :-)
				  // Ciao Moni
				  Utils::normalizeIneq(cutLHS, betaPlayer, true);
				  // Resize for the extra variables...
				  arma::sp_mat cutL = Utils::resizePatch(
						arma::sp_mat{cutLHS}.t(), 1, this->PolyLCP.at(player)->getNumCols());

				  // Add the cut and return false
				  this->PolyLCP.at(player)->addCustomCuts(cutL, arma::vec{betaPlayer});
				  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::equilibriumOracle: (P" << player
								  << ") adding a cut";
				  addedCuts = true;
				  return false;
				} else {

				  // We found a new vertex
				  arma::vec v(V.n_cols, arma::fill::zeros);
				  for (unsigned int i = 0; i < V.n_cols; ++i)
					 v[i] = l[i].get(GRB_DoubleAttr_X);

				  this->Trees.at(player)->addVertex(v, false);
				  LOG_S(1) << "Algorithms::EPEC::CutAndPlay::equilibriumOracle: (P" << player
							  << ") adding vertex for Player. " << (budget - k - 1) << "/" << budget
							  << " iterations left";
				}
			 }


		  } // status optimal for leaderModel
		  else if (status == GRB_UNBOUNDED) {
			 auto relaxed = leaderModel->relax();
			 relaxed.optimize();
			 // Well, we have a ray. But let's normalize it...
			 try {
				for (unsigned int i = 0; i < cutLHS.size(); ++i)
				  cutLHS[i] =
						relaxed.getVarByName("x_" + std::to_string(i)).get(GRB_DoubleAttr_UnbdRay);
			 } catch (GRBException &e) {
				throw ZEROException(e);
			 }
			 cutLHS = Utils::normalizeVec(cutLHS);
			 LOG_S(1) << "Algorithms::EPEC::CutAndPlay::equilibriumOracle: (P" << player
						 << ") new ray";
			 // Add the ray and repeat
			 this->Trees.at(player)->addRay(cutLHS);
		  } // status unbounded for leaderModel

		  // Otherwise, we don't know what happened
		  else {
			 throw ZEROException(ZEROErrorCode::Assertion,
										"Unknown status (" + std::to_string(status) +
											 ") for leaderModel for player " + std::to_string(player));
		  }
		  // no separation
		}
	 } else {
		throw ZEROException(ZEROErrorCode::Assertion,
								  "Unknown status for dualMembershipModel for player " +
										std::to_string(player));
	 }
  }
  return false;
}


/**
 * @brief Adds a value cut to @p player MathOpt::LCP
 * @param player The index of the player
 * @param RHS The RHS of the value cut
 * @param xMinusI The strategies of the other players
 */
void Algorithms::EPEC::CutAndPlay::addValueCut(const unsigned int player,
																		 const double       RHS,
																		 const arma::vec &  xMinusI) {

  arma::vec LHS = (this->EPECObject->LeaderObjective.at(player)->c +
						 this->EPECObject->LeaderObjective.at(player)->C * xMinusI);

  // Constant!
  if (Utils::isEqual(arma::max(LHS), 0)) {
	 LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::addValueCut: "
						 "Constant cut. Discarding. ";
	 return;
  }

  double trueRHS;
  if (Utils::nonzeroDecimals(RHS, 6) >= 6) {
	 LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::addValueCut: "
						 "Numerically unstable cut. This may cause numerical problems. "
					 << player;
	 // Let's try to save what we can.
	 for (unsigned int i = 0; i < LHS.size(); ++i)
		LHS.at(i) = Utils::round_nplaces(LHS.at(i), 5);
	 trueRHS = Utils::round_nplaces(RHS, 5);
	 Utils::normalizeIneq(LHS, trueRHS, true);
  } else {
	 trueRHS = Utils::round_nplaces(RHS, 6);
	 //  LHS.print("LHS with RHS of " + std::to_string(trueRHS));
	 Utils::normalizeIneq(LHS, trueRHS, false);
  }


  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::addValueCut: "
					  "adding cut for Player "
				  << player;

  // Resize
  arma::sp_mat cutLHS =
		Utils::resizePatch(arma::sp_mat{LHS}.t(), 1, this->PolyLCP.at(player)->getNumCols());

  this->PolyLCP.at(player)->addCustomCuts(-cutLHS, arma::vec{-trueRHS});
}


/**
 * @brief Gets the LCP for player @p player, and tries to see whether the origin is feasible. If
 * the answer is yes, sets the corresponding object in the tree to true.
 * @param player The player's id
 */
void Algorithms::EPEC::CutAndPlay::originFeasibility(unsigned int player) {

  arma::vec zeros(this->EPECObject->LeaderObjective.at(player)->C.n_cols, arma::fill::zeros);
  auto      origin = this->EPECObject->PlayersLCP.at(player).get()->LCPasMILP(
      this->EPECObject->LeaderObjective.at(player)->C,
      this->EPECObject->LeaderObjective.at(player)->c,
      zeros,
      false);

  for (unsigned int j = 0; j < this->EPECObject->LeaderObjective.at(player)->c.size(); j++) {
	 origin->addConstr(
		  origin->getVarByName("x_" + std::to_string(j)), GRB_EQUAL, 0, "Fix_x_" + std::to_string(j));
  }
  origin->update();
  origin->optimize();
  if (origin->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
	 this->Trees.at(player)->containsOrigin = true;
	 LOG_S(1) << "Algorithms::EPEC::CutAndPlay::originFeasibility: "
					 "Feasible origin for player "
				 << player;
  } else {
	 LOG_S(1) << "Algorithms::EPEC::CutAndPlay::originFeasibility: "
					 "Infeasible origin for player "
				 << player;
  }
}
/**
 * @brief Given the Game::EPEC instance, this method solves it through the outer approximation
 * scheme.
 */
void Algorithms::EPEC::CutAndPlay::solve() {
  // Set the initial point for all countries as 0 and solve the respective LCPs?
  this->EPECObject->SolutionX.zeros(this->EPECObject->NumVariables);
  bool solved = {false};
  if (this->EPECObject->Stats.AlgorithmData.TimeLimit.get() > 0)
	 this->EPECObject->InitTime = std::chrono::high_resolution_clock::now();

  this->EPECObject->Stats.NumIterations.set(0);


  // Initialize Trees
  // We actually do not use the complex tree structure for this vanilla-version, but we nevertheless
  // give the user the capability of doing so.
  this->Trees     = std::vector<OuterTree *>(this->EPECObject->NumPlayers, nullptr);
  this->Incumbent = std::vector<OuterTree::Node *>(this->EPECObject->NumPlayers, nullptr);
  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; i++) {
	 Trees.at(i) = new OuterTree(this->PolyLCP.at(i)->getNumRows(), this->Env);
	 this->Trees.at(i)->OriginalLCP =
		  std::make_unique<MathOpt::PolyLCP>(this->Env, *EPECObject->PlayersLowerLevels.at(i).get());
	 Incumbent.at(i) = Trees.at(i)->getRoot();
	 this->originFeasibility(i);
  }

  bool branch = true;
  // In this case, branchingLocations is a vector of locations with the length
  // of this->EPECObject->NumPlayers
  std::vector<int>      branchingLocations;
  std::vector<int>      branchingCandidatesNumber;
  int                   cumulativeBranchingCandidates = 0;
  unsigned int          branchingChoices              = 0;
  std::vector<long int> branches;
  while (!solved) {
	 branchingLocations.clear();
	 this->EPECObject->Stats.NumIterations.set(this->EPECObject->Stats.NumIterations.get() + 1);
	 LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::solve: Iteration "
					 << std::to_string(this->EPECObject->Stats.NumIterations.get());

	 branchingLocations            = std::vector<int>(this->EPECObject->NumPlayers, -1);
	 branchingCandidatesNumber     = std::vector<int>(this->EPECObject->NumPlayers, 0);
	 cumulativeBranchingCandidates = 0;

	 for (int j = 0; j < this->EPECObject->NumPlayers; ++j) {
		branchingCandidatesNumber.at(j) =
			 Trees.at(j)->getEncodingSize() - Incumbent.at(j)->getCumulativeBranches();
		cumulativeBranchingCandidates += branchingCandidatesNumber.at(j);
	 }

	 bool infeasibilityDetection = false;
	 branchingChoices            = 0;
	 if (branch) {
		for (int j = 0; j < this->EPECObject->NumPlayers && !infeasibilityDetection; ++j) {
		  // Check if we can branch
		  if (branchingCandidatesNumber.at(j) != 0) {
			 // In the first iteration, no complex branching rule.
			 if (this->EPECObject->Stats.NumIterations.get() == 1) {
				branchingLocations.at(j) = this->getFirstBranchLocation(j, Incumbent.at(j));
				// Check if we detected infeasibility
				if (branchingLocations.at(j) < 0) {
				  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::solve: "
									  "firstBranching proves infeasibility for player  "
								  << j;
				  infeasibilityDetection = true;
				  break;
				}
			 } else {
				if (this->EPECObject->Stats.AlgorithmData.BranchingStrategy.get() ==
					 Data::EPEC::BranchingStrategy::HybridBranching)
				  branchingLocations.at(j) = this->hybridBranching(j, Incumbent.at(j));
				else if (this->EPECObject->Stats.AlgorithmData.BranchingStrategy.get() ==
							Data::EPEC::BranchingStrategy::DeviationBranching) {
				  branchingLocations.at(j) = this->deviationBranching(j, Incumbent.at(j));
				  if (branchingLocations.at(j) == -1)
					 branchingLocations.at(j) = this->getFirstBranchLocation(j, Incumbent.at(j));
				}
				// Check if we detected infeasibility
				if (branchingLocations.at(j) == -2) {
				  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::solve: "
									  "hybridBranching proves infeasibility for player "
								  << j;
				  infeasibilityDetection = true;
				  break;
				}
			 }
		  }
		}

		if (infeasibilityDetection) {
		  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::solve: "
							  "Solved without any equilibrium. Proven infeasibility";
		  this->EPECObject->Stats.Status.set(ZEROStatus::NashEqNotFound);
		  solved = true;
		  break;
		}


		// Check that there is at least a player has a branching selection with
		// hybrid branching
		if (*std::max_element(branchingLocations.begin(), branchingLocations.end()) < 0) {

		  // No branching candidates.
		  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::solve: "
							  "No more hybrid branching candidates for "
							  "any player. Checking if "
							  "any complementarities are left.";
		  this->printCurrentApprox();
		  for (int j = 0; j < this->EPECObject->NumPlayers; ++j)
			 branchingLocations.at(j) = this->getFirstBranchLocation(j, Incumbent.at(j));

		  if (*std::max_element(branchingLocations.begin(), branchingLocations.end()) < 0) {
			 LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::solve: "
								 "No more branching candidates.";
			 this->EPECObject->Stats.Status.set(ZEROStatus::NashEqNotFound);
			 break;
		  }
		}
	 }

	 for (int j = 0; j < this->EPECObject->NumPlayers; ++j) {
		if (branchingLocations.at(j) > -1) {
		  branchingChoices   = branchingChoices + 1;
		  branches           = Trees.at(j)->singleBranch(branchingLocations.at(j), *Incumbent.at(j));
		  auto childEncoding = this->Trees.at(j)->getNodes()->at(branches.at(0)).getEncoding();
		  this->PolyLCP.at(j)->outerApproximate(childEncoding, true);
		  // By definition of hybrid branching, the node should be feasible
		  Incumbent.at(j) = &(this->Trees.at(j)->getNodes()->at(branches.at(0)));
		  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::solve: "
							  "branching candidate for player "
						  << j << " is " << branchingLocations.at(j);
		} else if (!branch) {
		  // if we don't branch.
		  this->PolyLCP.at(j)->outerApproximate(Incumbent.at(j)->getEncoding(), true);
		  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::solve: "
							  "No branching for player "
						  << j;
		}
	 }

	 this->printCurrentApprox();
	 this->EPECObject->makePlayersQPs();
	 this->Feasible   = false;
	 int branchesLeft = cumulativeBranchingCandidates - branchingChoices;
	 if (this->EPECObject->Stats.AlgorithmData.TimeLimit.get() > 0) {
		// Then we should take care of time. Also, let's use an heuristic to compute the time for the
		// current outer approximation.
		const std::chrono::duration<double> timeElapsed =
			 std::chrono::high_resolution_clock::now() - this->EPECObject->InitTime;
		const double timeRemaining =
			 this->EPECObject->Stats.AlgorithmData.TimeLimit.get() - timeElapsed.count();


		double timeForNextIteration = timeRemaining * 0.98;

		//@todo Currently not used.
		/*
		if (branchesLeft > 0 && false)
		  timeForNextIteration = (timeRemaining * 0.2) / (cumulativeBranchingCandidates - 1);
		  */

		LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::solve: Allocating "
						<< timeForNextIteration << "s for the next iteration (" << branchesLeft
						<< " complementarities left).";
		this->EPECObject->computeNashEq(
			 this->EPECObject->Stats.AlgorithmData.PureNashEquilibrium.get(),
			 timeForNextIteration,
			 false,
			 true,
			 false);
	 } else {
		this->EPECObject->computeNashEq(
			 this->EPECObject->Stats.AlgorithmData.PureNashEquilibrium.get(), 0, false, true, false);
	 }

	 if (!this->EPECObject->NashEquilibrium && branchesLeft == 0) {
		if (this->EPECObject->Stats.Status.get() == ZEROStatus::TimeLimit) {
		  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::solve: "
							  "Time Limit Hit.";
		  solved = true;
		  break;
		}
		if (this->EPECObject->Stats.Status.get() == ZEROStatus::NashEqNotFound) {
		  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::solve: "
							  "Solved without any equilibrium.";
		  solved = true;
		  break;
		}
	 }

	 this->Feasible = false;
	 if (this->EPECObject->NashEquilibrium) {
		bool addedCuts{false};
		if (this->isFeasible(addedCuts)) {
		  this->Feasible = true;
		  this->EPECObject->Stats.Status.set(ZEROStatus::NashEqFound);
		  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::solve: "
							  "Solved. ";
		  this->after();
		  return;
		} else {
		  if (addedCuts) {
			 branch = false;
			 LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::solve: "
								 "Cuts were added. Skipping next branching phase. ";
		  } else {
			 branch = true;
		  }
		}
	 } else {
		branch = true;
	 }
	 if (this->EPECObject->Stats.AlgorithmData.TimeLimit.get() > 0) {
		const std::chrono::duration<double> timeElapsed =
			 std::chrono::high_resolution_clock::now() - this->EPECObject->InitTime;
		const double timeRemaining =
			 this->EPECObject->Stats.AlgorithmData.TimeLimit.get() - timeElapsed.count();
		if (timeRemaining <= 0) {
		  this->EPECObject->Stats.Status.set(ZEROStatus::TimeLimit);
		  this->after();
		  return;
		}
	 }
  }
  this->after();
}


/**
 * @brief Given the player index @p player, gets a feasibility quadratic problem enforcing @p x to
 * be in the feasible (approximated) region of the Game::EPEC::PlayersQP
 * @param player  The player index
 * @param x The strategy for the player
 * @return A Gurobi pointer to the model
 */
std::unique_ptr<GRBModel>
Algorithms::EPEC::CutAndPlay::getFeasibilityQP(const unsigned int player,
																		 const arma::vec &  x) {


  arma::vec xMinusI;
  this->EPECObject->getXMinusI(this->EPECObject->SolutionX, player, xMinusI);
  xMinusI.resize(this->EPECObject->PlayersQP.at(player)->getNumParams());
  auto model = this->EPECObject->PlayersQP.at(player)->solveFixed(xMinusI, false);
  // Enforce QP::y to be x, namely the point to belong to the feasible region
  for (unsigned int j = 0; j < x.size(); j++)
	 model->addConstr(model->getVarByName("y_" + std::to_string(j)),
							GRB_EQUAL,
							x.at(j),
							"Fix_y_" + std::to_string(j));
  // Reset the objective
  model->setObjective(GRBLinExpr{0}, GRB_MINIMIZE);
  return model;
}


/**
 * @brief Given the player index @p player, gets a feasibility quadratic problem enforcing @p x to
 * be in the feasible region of the given player
 * @param player  The player index
 * @param x The strategy for the player
 * @return A Gurobi pointer to the model
 */
bool Algorithms::EPEC::CutAndPlay::isFeasiblePure(const unsigned int player,
																			 const arma::vec &  x) {


  auto model = this->PolyLCP.at(player)->LCPasMIP(false, -1, 1, 1);
  for (unsigned int j = 0; j < x.size(); j++)
	 model->addConstr(model->getVarByName("x_" + std::to_string(j)),
							GRB_EQUAL,
							x.at(j),
							"Fix_x_" + std::to_string(j));
  // Reset the objective
  model->setObjective(GRBLinExpr{0}, GRB_MINIMIZE);
  model->optimize();
  const int status = model->get(GRB_IntAttr_Status);
  if (status == GRB_INFEASIBLE)
	 return false;
  else
	 return true;
}


/**
 * @brief Given @p player -- containing the id of the player, returns the branching
 * decision for that node given by a hybrid branching rule. In
 * particular, the method return the complementarity id maximizing a
 * combination of constraint violations and number of violated constraints.
 * @p node contains the tree's node. It isn't const since a branching candidate
 * can be pruned if infeasibility is detected. Note that if the problem is infeasible, namely one
 * complementarity branching candidate results in an infeasible relaxation, then all branching
 * candidates are removed from the list of branching candidates.
 * @param player The player id
 * @param node The pointer to the incumbent OuterTree::Node
 * @return The branching candidate. -1 if none. -2 if infeasible.
 **/
int Algorithms::EPEC::CutAndPlay::hybridBranching(const unsigned int player,
																			 OuterTree::Node *  node) {

  LOG_S(INFO) << "CutAndPlay::hybridBranching: Player " << player;

  int bestId = -1;
  if (this->EPECObject->NashEquilibrium) {
	 arma::vec zeros, x;

	 this->EPECObject->getXofI(this->EPECObject->SolutionX, player, x);
	 if (x.size() != this->EPECObject->LeaderObjective.at(player)->c.n_rows)
		throw ZEROException(ZEROErrorCode::Assertion, "wrong dimensioned x^i");

	 auto              currentEncoding = node->getEncoding();
	 std::vector<bool> incumbentApproximation;
	 double            bestScore = -1.0;

	 for (unsigned int i = 0; i < currentEncoding.size(); i++) {
		// For each complementarity


		if (node->getAllowedBranchings().at(i)) {
		  // Consider it if it is a good candidate for branching (namely, we
		  // didn't branch on it, or it wasn't proven to be infeasible)
		  incumbentApproximation = currentEncoding;
		  // Include this complementarity in the approximation
		  incumbentApproximation.at(i) = true;
		  // Build the approximation
		  this->PolyLCP.at(player)->outerApproximate(incumbentApproximation, true);
		  // If the approximation is infeasible, the
		  if (!this->PolyLCP.at(player)->getFeasOuterApp()) {
			 // The problem is infeasible!
			 LOG_S(INFO) << "CutAndPlay::hybridBranching: Player " << player
							 << " has an infeasible problem (outer relaxation induction)";
			 for (unsigned int j = 0; j < currentEncoding.size(); j++) {
				Trees.at(player)->denyBranchingLocation(*node, j);
			 }
			 return -2;
		  } else {
			 // In this case, we can check if the solution belongs to the outer
			 // approximation
			 this->EPECObject->makePlayerQP(player);
			 // Get the QP model with other players decision QP::x fixed to zero
			 // (since they only appear in the objective);
			 auto model = this->getFeasibilityQP(player, x);
			 model->optimize();
			 const int status = model->get(GRB_IntAttr_Status);
			 if (status == GRB_INFEASIBLE) {
				// If the status is infeasible, bingo! We want to get a measure of
				// the constraint violations given by the current x
				model->feasRelax(1, false, false, true);
				model->optimize();
				if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
				  if (model->getObjective().getValue() > bestScore) {
					 bestId    = i;
					 bestScore = model->getObjective().getValue();
					 LOG_S(INFO) << "CutAndPlay::hybridBranching: Player " << player
									 << " has violation of " << bestScore << " with complementarity " << i;
				  } else {
					 LOG_S(INFO) << "CutAndPlay::hybridBranching: Player " << player
									 << " has violation of " << model->getObjective().getValue()
									 << " with complementarity " << i;
				  }
				} else {
				  LOG_S(WARNING)
						<< "CutAndPlay::hybridBranching: Numerical difficulties in evaluating "
						<< std::to_string(i);
				}
			 } else {
				LOG_S(INFO) << "CutAndPlay::hybridBranching: Player " << player
								<< " has no violation with complementarity " << i;
			 }
		  }
		}
	 }
  }
  return bestId;
}


/**
 * @brief Given @p player -- containing the id of the player, returns the branching
 * decision for that node, where the complementarity is the most (possibly)
 * infeasible one (with both x and z positive). In particular, the method
 * return the (positive) id of the complementarity equation if there is a
 * feasible branching decision at @p node, and a negative value otherwise.
 * @param player The player id
 * @param node The pointer to the incumbent OuterTree::Node
 * @return The branching candidate. Negative if none
 */
int Algorithms::EPEC::CutAndPlay::infeasibleBranching(const unsigned int     player,
																				  const OuterTree::Node *node) {

  int result = -1;
  if (this->EPECObject->NashEquilibrium) {
	 // There exists a Nash Equilibrium for the outer approximation, which is not
	 // a Nash Equilibrium for the game
	 arma::vec x, z;
	 this->EPECObject->getXWithoutHull(this->EPECObject->SolutionX, x);
	 z                                      = this->PolyLCP.at(player)->zFromX(x);
	 std::vector<short int> currentSolution = this->PolyLCP.at(player)->solEncode(x);

	 double maxInfeas = 0;

	 //"The most infeasible" branching
	 for (unsigned int i = 0; i < currentSolution.size(); i++) {
		unsigned int varPos = i >= this->PolyLCP.at(player)->getLStart()
										  ? i + this->PolyLCP.at(player)->getNumberLeader()
										  : i;
		if (x(varPos) > 0 && z(i) > 0 && node->getAllowedBranchings().at(i) &&
			 currentSolution.at(i) == 0) {
		  if ((x(varPos) + z(i)) > maxInfeas) {
			 maxInfeas = x(varPos) + z(i);
			 result    = i;
		  }
		}
	 }
  }
  return result;
}

/**
 * @brief Given @p player -- containing the id of the player, returns the branching
 * decision for that node, where the complementarity helps include the deviation. In particular, the
 * method return the (positive) id of the complementarity equation if there is a feasible branching
 * decision at
 * @p node, and a negative value otherwise.
 * @param player The player id
 * @param node The pointer to the incumbent OuterTree::Node
 * @return The branching candidate. Negative if none
 */
int Algorithms::EPEC::CutAndPlay::deviationBranching(const unsigned int     player,
																				 const OuterTree::Node *node) {


  int result = -1;
  if (this->EPECObject->NashEquilibrium) {
	 // There exists a Nash Equilibrium for the outer approximation, which is not
	 // a Nash Equilibrium for the game
	 arma::vec dev;
	 arma::vec x;
	 this->EPECObject->getXWithoutHull(this->EPECObject->SolutionX, x);
	 std::vector<short int> currentSolution = this->PolyLCP.at(player)->solEncode(x);
	 this->EPECObject->bestResponse(
		  dev, player, this->EPECObject->SolutionX, {}, this->Trees.at(player)->OriginalLCP.get());
	 auto encoding = this->PolyLCP.at(player)->solEncode(dev);

	 for (unsigned int i = 0; i < encoding.size(); i++) {
		if (encoding.at(i) > 0 && node->getAllowedBranchings().at(i) && currentSolution.at(i) == 0) {
		  result = i;
		}
	 }
  }
  return result;
}


/**
 * @brief Given @p player -- containing the id of the player, returns the branching
 * decision for that node, with no complementarity condition enforced. In
 * particular, the method return the (positive) id of the complementarity
 * equation if there is a feasible branching decision at @p node, and a
 * negative value otherwise.
 * @param player The player id
 * @param node The pointer to the incumbent OuterTree::Node
 * @return The branching candidate. Negative if none
 */
int Algorithms::EPEC::CutAndPlay::getFirstBranchLocation(const unsigned int player,
																					  OuterTree::Node *  node) {
  /**
	* Given @p player -- containing the id of the player, returns the branching
	* decision for that node. In
	* particular, the method return the (positive) id of the complementarity
	* equation if there is a feasible branching decision at @p node, and a
	* negative value otherwise. Note that if the problem is infeasible, namely one
	* complementarity branching candidate results in an infeasible relaxation, then all branching
	* candidates are removed from the list of branching candidates.
	* @return a positive int with the id of the complementarity to branch on, or
	* a negative value if none exists.
	*/

  if (node->getCumulativeBranches() == Trees.at(player)->getEncodingSize())
	 return -1;
  auto         model = this->PolyLCP.at(player)->LCPasMIP(true, -1, 1, 1);
  unsigned int nR    = this->PolyLCP.at(player)->getNumRows();
  int          pos   = -nR;
  arma::vec    z, x;
  if (this->PolyLCP.at(player)->extractSols(
			 model.get(), z, x, true)) // If already infeasible, nothing to branch!
  {
	 std::vector<short int> v1 = this->PolyLCP.at(player)->solEncode(z, x);

	 double       maxvalx{-1}, maxvalz{-1};
	 unsigned int maxposx{0}, maxposz{0};
	 for (unsigned int i = 0; i < nR; i++) {
		unsigned int varPos = i >= this->PolyLCP.at(player)->getLStart()
										  ? i + this->PolyLCP.at(player)->getNumberLeader()
										  : i;
		if (x(varPos) > maxvalx && node->getAllowedBranchings().at(i)) {
		  maxvalx = x(varPos);
		  maxposx = i;
		}
		if (z(i) > maxvalz && node->getAllowedBranchings().at(i)) {
		  maxvalz = z(i);
		  maxposz = i;
		}
	 }
	 pos = maxvalz > maxvalx ? maxposz : maxposx;
  } else {
	 // The problem is infeasible!
	 LOG_S(INFO) << "CutAndPlay::getFirstBranchLocation: Player " << player
					 << " has an infeasible problem (outer relaxation induction)";
	 for (unsigned int j = 0; j < node->getEncoding().size(); j++) {
		Trees.at(player)->denyBranchingLocation(*node, j);
	 }
	 return -1;
  }
  return pos;
}



/**
 * @brief Given @p player -- containing the id of the player -- and @p node
 * containing a node, returns the branching decision for that node, with
 * respect to the current node. In particular, the method return the
 * (positive) id of the complementarity equation if there is a feasible
 * branching decision at @p node, and a negative value otherwise.
 * @param player The player id
 * @param node The pointer to the incumbent OuterTree::Node
 * @return A vector of 4 integers with the branching location given by the most
 * Algorithms::EPEC::CutAndPlay::infeasibleBranching,
 * Algorithms::EPEC::CutAndPlay::deviationBranching,
 * Algorithms::EPEC::CutAndPlay::hybridBranching, and
 * Algorithms::EPEC::CutAndPlay::getFirstBranchLocation, respectively. If an int is
 * negative, there is no real candidate.
 */
[[maybe_unused]] std::vector<int>
Algorithms::EPEC::CutAndPlay::getNextBranchLocation(const unsigned int player,
																				OuterTree::Node *  node) {

  std::vector<int> decisions = {-1, -1, -1, -1};
  decisions.at(0)            = this->infeasibleBranching(player, node);
  decisions.at(1)            = this->deviationBranching(player, node);
  decisions.at(2)            = this->hybridBranching(player, node);

  if (decisions.at(0) < 0 && decisions.at(1) < 0 && decisions.at(2) < 0) {
	 LOG_S(INFO) << "Player " << player
					 << ": branching with FirstBranchLocation is the only available choice";
	 decisions.at(3) = this->getFirstBranchLocation(player, node);
  }

  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::getNextBranchinglocation: "
					  "given decisions are: ";
  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::"
					  "getNextBranchinglocation:\t Infeasible="
				  << decisions.at(0);
  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::"
					  "getNextBranchinglocation:\t Deviation="
				  << decisions.at(1);
  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::"
					  "getNextBranchinglocation:\t Hybrid="
				  << decisions.at(2);
  LOG_S(INFO) << "Algorithms::EPEC::CutAndPlay::"
					  "getNextBranchinglocation:\t First="
				  << decisions.at(3);
  return decisions;
}


/**
 * @brief Prints a log message containing the encoding at the current outer
 * approximation iteration
 */
void Algorithms::EPEC::CutAndPlay::printCurrentApprox() {
  LOG_S(INFO) << "Current Node Approximation:";
  for (unsigned int p = 0; p < this->EPECObject->NumPlayers; ++p) {
	 std::stringstream msg;
	 msg << "\tPlayer " << p << ":";
	 for (unsigned int i = 0; i < this->Incumbent.at(p)->getEncoding().size(); i++) {
		msg << "\t" << this->Incumbent.at(p)->getEncoding().at(i);
	 }
	 LOG_S(INFO) << msg.str();
  }
}

/**
 * @brief  Given the vector of branching candidates from
 * Algorithms::EPEC::CutAndPlay::getNextBranchLocation, prints a sum up of them
 * @param vector Output of Algorithms::EPEC::CutAndPlay::getNextBranchLocation
 */
void Algorithms::EPEC::CutAndPlay::printBranchingLog(std::vector<int> vector) {
  LOG_S(INFO) << "Current Branching Log:";
  LOG_S(INFO) << "\tInfeasibleBranching: " << vector.at(0);
  LOG_S(INFO) << "\tDeviationBranching: " << vector.at(1);
  LOG_S(INFO) << "\tHybridBranching: " << vector.at(2);
  LOG_S(INFO) << "\tFirstAvail: " << vector.at(3);
}

/**
 * @brief Checks whether the current solution is a pure-strategy nash equilibrium
 * @param tol A numerical tolerance. Currently not used
 * @return True if the strategy is a pure nash equilibrium
 */
bool Algorithms::EPEC::CutAndPlay::isPureStrategy(double tol) const {
  if (!this->Feasible)
	 return false;
  else {
	 for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i)
		if (!Trees.at(i)->getPure())
		  return false;

	 return true;
  }
}
void Algorithms::EPEC::CutAndPlay::after() {
  bool                      pureStrategy = true;
  std::vector<unsigned int> numComps;
  for (unsigned int i = 0; i < this->EPECObject->getNumPlayers(); ++i) {
	 if (!this->Trees.at(i)->getPure()) {
		pureStrategy = false;
	 }
	 unsigned int counter = 0;
	 for (unsigned int j = 0; j < this->Incumbent.at(i)->getEncoding().size(); j++)
		counter += this->Incumbent.at(i)->getEncoding().at(j);
	 numComps.push_back(counter);
  }
  this->EPECObject->Stats.PureNashEquilibrium.set(pureStrategy);
  this->EPECObject->Stats.AlgorithmData.OuterComplementarities.set(numComps);
  LOG_S(3) << "Algorithms::EPEC::CutAndPlay::after: post-processing results.";
}



/**
 * @brief Given the parent node address @p parent, the @p idComp to branch
	on, and the @p id, creates a new node
 * @param parent The parent node
 * @param idComp The id of the node
 * @param id The The branching candidate
 */
Algorithms::EPEC::OuterTree::Node::Node(Node &parent, unsigned int idComp, unsigned long int id) {
  this->IdComps                      = std::vector<unsigned int>{idComp};
  this->Encoding                     = parent.Encoding;
  this->Encoding.at(idComp)          = true;
  this->AllowedBranchings            = parent.AllowedBranchings;
  this->AllowedBranchings.at(idComp) = false;
  this->Id                           = id;
  this->Parent                       = &parent;
}

/**
 * @brief Constructor for the root node, given the encoding size, namely the number of
 * complementarity equations
 * @param encSize The number of complementarities
 */
Algorithms::EPEC::OuterTree::Node::Node(unsigned int encSize) {
  this->Encoding          = std::vector<bool>(encSize, false);
  this->Id                = 0;
  this->AllowedBranchings = std::vector<bool>(encSize, true);
}


/**
 * @brief If a complementarity equation @p location  has proven to be infeasible
 * or it isn't a candidate for branching, this method prevents any further
 * branching on it for the node @p node.
 * @param node The node pointer
 * @param location The denied branching location
 */
void Algorithms::EPEC::OuterTree::denyBranchingLocation(Algorithms::EPEC::OuterTree::Node &node,
																		  const unsigned int &location) const {
  if (location >= this->EncodingSize)
	 throw ZEROException(ZEROErrorCode::OutOfRange, "idComp is larger than the encoding size");
  if (!node.AllowedBranchings.at(location))
	 LOG_S(WARNING) << "Algorithms::EPEC::OuterTree::denyBranchingLocation: location " << location
						 << "has been already denied.";

  LOG_S(INFO) << "Algorithms::EPEC::OuterTree::denyBranchingLocation: location " << location
				  << "denied.";
  node.AllowedBranchings.at(location) = false;
}


/**
 * @brief Given the @p idComp and the parent node @p t, creates a single
 * child by branching on @p idComp.
 * @param idComp The branching id for the complementarity
 * @param t The pointer to the node
 * @return The node counter stored in a single-element vector
 */
std::vector<long int>
Algorithms::EPEC::OuterTree::singleBranch(const unsigned int                 idComp,
														Algorithms::EPEC::OuterTree::Node &t) {
  if (idComp >= this->EncodingSize)
	 throw ZEROException(ZEROErrorCode::OutOfRange, "idComp is larger than the encoding size");
  if (t.Encoding.at(idComp) != 0) {
	 LOG_S(WARNING) << "OuterTree: cannot branch on this complementary, since it already "
							 "has been processed.";
	 return std::vector<long int>{-1};
  }
  auto child = Node(t, idComp, this->nextIdentifier());

  this->Nodes.push_back(child);
  return std::vector<long int>{this->NodeCounter - 1};
}

/**
 * @brief Adds a vertex to OuterTree::V
 * @param vertex The vector containing the vertex
 * @param checkDuplicates True if the method has to check for duplicates
 * @return True if the vertex was added
 */
bool Algorithms::EPEC::OuterTree::addVertex(const arma::vec &vertex, bool checkDuplicates) {
  if (vertex.size() != this->V.n_cols && this->V.n_rows > 0)
	 throw ZEROException(ZEROErrorCode::OutOfRange, "Ill-dimensioned vertex");

  bool go = false;
  if (checkDuplicates)
	 go = Utils::containsRow(this->V, vertex, 1e-5);


  if (!go) {
	 int nCols = this->V.n_cols < 1 ? vertex.size() : this->V.n_cols;
	 this->V.resize(this->V.n_rows + 1, nCols);
	 this->V.row(this->V.n_rows - 1) = vertex.t();
	 return true;
  }
  return false;
}

/**
 * @brief Adds a ray to OuterTree::R
 * @param ray The vector containing the ray
 */
void Algorithms::EPEC::OuterTree::addRay(const arma::vec &ray) {
  if (ray.size() != this->R.n_cols && this->R.n_rows > 0)
	 throw ZEROException(ZEROErrorCode::OutOfRange, "Ill-dimensioned ray");
  int nCols = this->R.n_cols < 1 ? ray.size() : this->R.n_cols;
  this->R.resize(this->R.n_rows + 1, nCols);
  this->R.row(this->R.n_rows - 1) = ray.t();
}

/**
 * @brief Given the parent node address @p parent, the @p idsComp to branch
 * on (containing all the complementarities ids), and the @p id, creates a
 * new node
 * @param parent The parent node pointer
 * @param idsComp  The vector of branching locations
 * @param id  The node id for the children
 */
Algorithms::EPEC::OuterTree::Node::Node(Node &            parent,
													 std::vector<int>  idsComp,
													 unsigned long int id) {
  this->IdComps           = std::vector<unsigned int>();
  this->Encoding          = parent.Encoding;
  this->AllowedBranchings = parent.AllowedBranchings;
  for (auto &idComp : idsComp) {
	 if (idComp < 0)
		throw ZEROException(ZEROErrorCode::Assertion, "idComp is negative");
	 this->Encoding.at(idComp)          = true;
	 this->AllowedBranchings.at(idComp) = false;
	 this->IdComps.push_back(idComp);
  }
  this->Id     = id;
  this->Parent = &parent;
}
