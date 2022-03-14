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

#include "games/epec.h"

#include "games/algorithms/EPEC/epec_polybase.h"
#include <memory>


/**
 * @brief Finalizes the creation of a Game::EPEC object.
 * @details Performs a bunch of job after all data for a Game::EPEC object are
 * given, namely.
 * Models::EPEC::computeLeaderLocations -	Adds the required dummy
 * variables to each leader's problem so that a game among the leaders can be
 * defined. Calls Game::EPEC::addDummyLead
 * 	-	Makes the market clearing constraint in each country. Calls
 */
void Game::EPEC::finalize() {
  if (this->Finalized)
	 std::cerr << "Warning in Game::EPEC::finalize: Model already Finalized\n";

  this->NumPlayers = static_cast<int>(this->PlayersLowerLevels.size());
  /// Game::EPEC::preFinalize() can be overridden, and that code will run before
  /// calling Game::EPEC::finalize()
  this->preFinalize();

  // Vector with the number of variables of the convex hull
  this->ConvexHullVariables = std::vector<unsigned int>(this->NumPlayers, 0);
  // Reset the statistic of feasible polyhedra for any player
  this->Stats.AlgorithmData.FeasiblePolyhedra.set(std::vector<unsigned int>(this->NumPlayers, 0));

  // Assign locations to the variables
  this->computeLeaderLocations(this->numMCVariables);
  // Initialize leader objective and PlayersQP
  this->LeaderObjective           = std::vector<std::shared_ptr<MathOpt::QP_Objective>>(NumPlayers);
  this->LeaderObjectiveConvexHull = std::vector<std::shared_ptr<MathOpt::QP_Objective>>(NumPlayers);
  this->PlayersQP                 = std::vector<std::shared_ptr<MathOpt::QP_Param>>(NumPlayers);
  this->PlayersLCP                = std::vector<std::shared_ptr<MathOpt::LCP>>(NumPlayers);
  this->SizesWithoutHull          = std::vector<unsigned int>(NumPlayers, 0);

  // For each player
  for (unsigned int i = 0; i < this->NumPlayers; i++) {
	 // Add a trivial variables to the player's lower level
	 this->addDummyLead(i);
	 // Leader objective and its dummy-enlarged version with convex hull variables
	 this->LeaderObjective.at(i)           = std::make_shared<MathOpt::QP_Objective>();
	 this->LeaderObjectiveConvexHull.at(i) = std::make_shared<MathOpt::QP_Objective>();
	 // Fill the original objectives
	 this->makeObjectivePlayer(i, *this->LeaderObjective.at(i).get());
	 // Compute sizes
	 this->SizesWithoutHull.at(i) = *this->LocEnds.at(i);
  }

  // Finalize
  this->Finalized = true;

  /// Game::EPEC::postFinalize() can be overridden, and that code will run after
  /// calling Game::EPEC::finalize()
  this->postFinalize();
}


/**
 * @brief Adds dummy variables to the leader of an EPEC - useful after computing the convex hull.
 @param i The leader id to whom dummy variables should be added
 */

void Game::EPEC::addDummyLead(const unsigned int i) {
  ZEROAssert(i < this->NumPlayers);
  const unsigned int nEPECvars        = this->NumVariables;
  const unsigned int nThisCountryvars = *this->LocEnds.at(i);


  ZEROAssert(nEPECvars >= nThisCountryvars);

  // Add dummy vars associated with "everything else"
  this->PlayersLowerLevels.at(i).get()->addDummy(nEPECvars - nThisCountryvars);
}

/**
 * @brief Assigns the values to LeaderLocations to each player.
 * @param addSpaceForMC contains the number of Market Clearing variables
 */
void Game::EPEC::computeLeaderLocations(const unsigned int addSpaceForMC) {

  this->LeaderLocations       = std::vector<unsigned int>(this->NumPlayers);
  this->LeaderLocations.at(0) = 0;
  for (unsigned int i = 1; i < this->NumPlayers; i++) {
	 this->LeaderLocations.at(i) = this->LeaderLocations.at(i - 1) + *this->LocEnds.at(i - 1);
  }
  this->NumVariables = this->LeaderLocations.back() + *this->LocEnds.back() + addSpaceForMC;
}


/**
 * @brief Given a solution in @p x, computes the other players' solution x minus @p i in @p xMinusI
 * @param x The incumbent solution
 * @param i The player id
 * @param xMinusI The output strategy
 */
void Game::EPEC::getXMinusI(const arma::vec &x, const unsigned int &i, arma::vec &xMinusI) const {
  ZEROAssert(i < this->NumPlayers);
  const unsigned int nEPECvars            = this->NumVariables;
  const unsigned int nThisCountryvars     = *this->LocEnds.at(i);
  const unsigned int nThisCountryHullVars = this->ConvexHullVariables.at(i);
  const auto         nConvexHullVars      = static_cast<const unsigned int>(
      std::accumulate(this->ConvexHullVariables.rbegin(), this->ConvexHullVariables.rend(), 0));

  xMinusI.zeros(nEPECvars -            // All variables in EPEC
					 nThisCountryvars -     // Subtracting this country's variables
					 nConvexHullVars +      // We don't want any convex hull variables
					 nThisCountryHullVars); // We subtract the hull variables
													// associated to the ith player
													// convex hull vars, since we double subtracted

  for (unsigned int j = 0, count = 0, current = 0; j < this->NumPlayers; ++j) {
	 if (i != j) {
		current = *this->LocEnds.at(j) - this->ConvexHullVariables.at(j);
		xMinusI.subvec(count, count + current - 1) =
			 x.subvec(this->LeaderLocations.at(j), this->LeaderLocations.at(j) + current - 1);
		count += current;
	 }
  }
  // We need to keep track of MC_vars also for this country
  for (unsigned int j = 0; j < this->numMCVariables; j++)
	 xMinusI.at(xMinusI.n_rows - this->numMCVariables + j) =
		  x.at(this->NumVariables - this->numMCVariables + j);
}


/**
 * @brief Given the player id @p i and the solution @p x, the method returns in @p
 * xWithoutHull the x vector for the given player, with the convex-hull's
 * variables in case @p hull is false. Also, no MC variables are included
 * @param x The EPEC's solution
 * @param i The player id
 * @param xOfI Output solution vector for i
 * @param hull True if convex-hull variables need to be included in the output
 */
void Game::EPEC::getXofI(const arma::vec    &x,
								 const unsigned int &i,
								 arma::vec          &xOfI,
								 bool                hull) const {
  ZEROAssert(i < this->NumPlayers);
  const unsigned int nThisCountryvars     = *this->LocEnds.at(i);
  const unsigned int nThisCountryHullVars = this->ConvexHullVariables.at(i);

  unsigned int vars, current = 0;
  if (hull) {
	 vars    = nThisCountryvars;
	 current = *this->LocEnds.at(i);
  } else {
	 vars    = nThisCountryvars - nThisCountryHullVars;
	 current = *this->LocEnds.at(i) - this->ConvexHullVariables.at(i);
  }
  xOfI.zeros(vars);
  xOfI.subvec(0, vars - 1) =
		x.subvec(this->LeaderLocations.at(i), this->LeaderLocations.at(i) + current - 1);
}


/**
 * @brief Given the the solution @p x, the method returns in @p
 * xWithoutHull the x vector without the convex-hull's
 * variables.  Also, no MC variables are included
 * @param x The EPEC's solution
 * @param xWithoutHull The output solution without convex hull variables
 */
void Game::EPEC::getXWithoutHull(const arma::vec &x, arma::vec &xWithoutHull) const {
  const unsigned int nEPECvars       = this->NumVariables;
  const auto         nConvexHullVars = static_cast<const unsigned int>(
      std::accumulate(this->ConvexHullVariables.rbegin(), this->ConvexHullVariables.rend(), 0));

  xWithoutHull.zeros(nEPECvars -       // All variables in EPEC
							nConvexHullVars); // We subtract the hull variables
													// associated to the convex hull
													// convex hull vars

  for (unsigned int j = 0, count = 0, current; j < this->NumPlayers; ++j) {
	 current = *this->LocEnds.at(j) - this->ConvexHullVariables.at(j);
	 xWithoutHull.subvec(count, count + current - 1) =
		  x.subvec(this->LeaderLocations.at(j), this->LeaderLocations.at(j) + current - 1);
	 count += current;
  }

  for (unsigned int j = 0; j < this->numMCVariables; j++)
	 xWithoutHull.at(xWithoutHull.n_rows - this->numMCVariables + j) =
		  x.at(this->NumVariables - this->numMCVariables + j);
}


/**
 * @brief Given a player id @p i, the incumbent solution @p x (and optionally a custom player LCP @p
 * customLCP), this method returns the model corresponding to the best response of @p i given the
 * other players' decisions in @p x.
 * @param i The player id
 * @param x The incumbent solution
 * @param customLCP An optional parameter with a custom LCP
 * @return A pointer to the (solved) Gurobi model for the best response
 */
std::unique_ptr<GRBModel> Game::EPEC::bestResponseProgram(const unsigned int i,
																			 const arma::vec   &x,
																			 MathOpt::PolyLCP  *customLCP) const {
  ZEROAssert(this->Finalized);
  ZEROAssert(i < this->NumPlayers);

  arma::vec solOther;
  this->getXMinusI(x, i, solOther);
  auto LCP = (customLCP == nullptr) ? this->PlayersLCP.at(i).get() : customLCP;
  if (this->LeaderObjective.at(i)->Q.n_nonzero > 0)
	 return LCP->LCPasMIQP(this->LeaderObjective.at(i)->Q,
								  this->LeaderObjective.at(i)->C,
								  this->LeaderObjective.at(i)->c,
								  solOther,
								  true);
  else
	 return LCP->LCPasMILP(
		  this->LeaderObjective.at(i)->C, this->LeaderObjective.at(i)->c, solOther, true);
}
/**
 * @brief Returns the best-response value for the
 * player @p player given the decision @p x of all other players.
 * @details
 * Calls Game::EPEC::respond and obtains the std::unique_ptr to GRBModel of
 * best response by player @p player. Then solves the model and returns the
 * appropriate objective value.
 * @param sol The optimal response for @p player
 * @param player The player id
 * @param x The solution vector
 * @param prevDev An optional previous deviation encountered. This is useful to normalize any
 * unbounded best response
 * @param customLCP An optional pointer to a custom player LCP.
 * @returns The optimal objective value for the player @p player.
 */
double Game::EPEC::bestResponse(arma::vec        &sol,
										  unsigned int      player,
										  const arma::vec  &x,
										  const arma::vec  &prevDev,
										  MathOpt::PolyLCP *customLCP) const {
  ZEROAssert(this->Finalized);
  ZEROAssert(player < this->NumPlayers);

  auto      model  = this->bestResponseProgram(player, x, customLCP);
  const int status = model->get(GRB_IntAttr_Status);

  if (status == GRB_UNBOUNDED || status == GRB_OPTIMAL) {
	 //If unbounded or optimal either solution or extreme ray
	 unsigned int Nx = this->PlayersLCP.at(player)->getNumCols();
	 sol.zeros(Nx);
	 for (unsigned int i = 0; i < Nx; ++i)
		sol.at(i) = model->getVarByName("x_" + std::to_string(i)).get(GRB_DoubleAttr_X);

	 if (status == GRB_UNBOUNDED) {
		LOG_S(WARNING) << "Game::EPEC::bestResponse: deviation is "
								"unbounded.";
		GRBLinExpr obj = 0;
		model->setObjective(obj);
		model->optimize();
		if (!prevDev.empty()) {
		  LOG_S(1) << "Generating an improvement basing on the extreme ray.";
		  // Fetch objective function coefficients
		  GRBQuadExpr QuadObj = model->getObjective();
		  arma::vec   objcoeff;
		  for (int i = 0; i < QuadObj.size(); ++i)
			 objcoeff.at(i) = QuadObj.getCoeff(i);

		  // Create objective function objects
		  arma::vec objvalue = prevDev * objcoeff;
		  arma::vec newobjvalue{0};
		  bool      improved{false};

		  // improve following the unbounded ray
		  while (!improved) {
			 for (unsigned int i = 0; i < Nx; ++i)
				sol.at(i) = sol.at(i) +
								model->getVarByName("x_" + std::to_string(i)).get(GRB_DoubleAttr_UnbdRay);
			 newobjvalue = sol * objcoeff;
			 if (newobjvalue.at(0) < objvalue.at(0))
				improved = true;
		  }
		  return newobjvalue.at(0);

		} else {
		  return model->get(GRB_DoubleAttr_ObjVal);
		}
	 }
	 return model->get(GRB_DoubleAttr_ObjVal);
  } else {
	 return GRB_INFINITY;
  }
}

/**
 * @brief Makes the MathOpt::QP_Param corresponding to the @p i-th country.
 *  - First gets the Game::LCP object from @p Game::EPEC::PlayersLowerLevels and
 * makes a MathOpt::QP_Param with this LCP as the lower level
 *  - This is achieved by calling LCP::makeQP and using the objective value
 * object in @p Game::EPEC::LeaderObjective
 *  - Finally the locations are updated owing to the complete convex hull
 * calculated during the call to LCP::makeQP
 * @param i The player's id
 * @note Overloaded as Models::EPEC::makePlayerQP()
 */
void Game::EPEC::makePlayerQP(const unsigned int i) {

  ZEROAssert(this->Finalized);
  ZEROAssert(i < this->NumPlayers);
  this->PlayersQP.at(i)     = std::make_shared<MathOpt::QP_Param>(this->Env);
  const auto &origLeadObjec = *this->LeaderObjective.at(i).get();

  this->LeaderObjectiveConvexHull.at(i) = std::make_shared<MathOpt::QP_Objective>(
		MathOpt::QP_Objective{origLeadObjec.Q, origLeadObjec.C, origLeadObjec.c});
  this->PlayersLCP.at(i)->makeQP(*this->LeaderObjectiveConvexHull.at(i).get(),
											*this->PlayersQP.at(i).get());
}

/**
 * @brief Makes the MathOpt::QP_Param for all the countries. Calls are made to
 * Models::EPEC::makePlayerQP(const unsigned int i) for each valid player id
 */
void Game::EPEC::makePlayersQPs() {
  for (unsigned int i = 0; i < this->NumPlayers; ++i) {
	 this->Game::EPEC::makePlayerQP(i);
  }
  for (unsigned int i = 0; i < this->NumPlayers; ++i) {
	 // Adjusting "stuff" because we now have new convHull variables
	 unsigned long int originalSizeWithoutHull = this->LeaderObjective.at(i)->Q.n_rows;
	 unsigned long int convHullVarCount =
		  this->LeaderObjectiveConvexHull.at(i)->Q.n_rows - originalSizeWithoutHull;

	 LOG_S(1) << "Game::EPEC::makePlayerQP: Added " << convHullVarCount
				 << " convex hull variables to QP #" << i;

	 // Location details
	 this->ConvexHullVariables.at(i) = convHullVarCount;
	 // All other players' QP
	 if (this->NumPlayers > 1) {
		for (int j = 0; j < this->NumPlayers; j++) {
		  if (i != j) {
			 this->PlayersQP.at(j)->addDummy(
				  convHullVarCount,
				  0,
				  this->PlayersQP.at(j)->getNumParams() -
						this->numMCVariables); // The position to add parameters is
													  // towards the end of all parameters,
													  // giving space only for the
													  // numMCVariables number of market
													  // clearing variables
		  }
		}
	 }
  }
  this->updateLocations();
  this->computeLeaderLocations(this->numMCVariables);
}

/**
 * @brief Formulates the LCP to compute an equilibrium. In this LCP, each player's feasible region
 * is given by the PlayersQP associated entry. If the QP stems from a full enumeration, for
 * instance, the solution will be exact.
 */
void Game::EPEC::makeTheLCP() {
  ZEROAssert(this->PlayersQP.front() != nullptr);


  // Preliminary set up to get the LCP ready
  unsigned long int Nvar =
		this->PlayersQP.front()->getNumParams() + this->PlayersQP.front()->getNumVars();
  arma::sp_mat MC(0, Nvar), dumA(0, Nvar);
  arma::vec    MCRHS, empty;
  MCRHS.zeros(0);
  empty.zeros(0);
  this->makeMCConstraints(MC, MCRHS);
  LOG_S(1) << "Game::EPEC::makeTheLCP(): Market Clearing "
				  "constraints are ready";
  std::vector<std::shared_ptr<MathOpt::MP_Param>> MPCasted;
  for (auto &item : this->PlayersQP) {
	 auto m = std::dynamic_pointer_cast<MathOpt::MP_Param>(item);
	 MPCasted.push_back(m);
  }
  this->TheNashGame =
		std::make_unique<Game::NashGame>(this->Env, MPCasted, MC, MCRHS, 0, dumA, empty);
  LOG_S(1) << "Game::EPEC::makeTheLCP(): NashGame is ready";
  this->TheLCP = std::make_unique<MathOpt::LCP>(this->Env, *TheNashGame);
  LOG_S(1) << "Game::EPEC::makeTheLCP(): LCP is ready";
  LOG_S(2) << *TheNashGame;
}


/**
 * Given the object Game::EPEC::LCPModel, this method sets its objective to the social welfare. In
 * specific, we can decide whether to include or not the linear or quadratic part of the welfare
 * @param linear True if the linear part is to be included, false otherwise.
 * @param quadratic True if the quadratic part is to be included, false otherwise.
 */
void Game::EPEC::setWelfareObjective(bool linear = true, bool quadratic = true) {


  if (!linear && !quadratic) {
	 this->LCPModel->setObjective(GRBLinExpr{0}, GRB_MAXIMIZE);
	 return;
  }

  std::vector<std::vector<unsigned int>> xOfIs; // indexes of variables for each player
  std::vector<std::vector<unsigned int>>
				  xMinusIs; // indexes of variables for each "other" player (xminusi)
  GRBLinExpr  linearWelfare = 0;
  GRBQuadExpr quadrWelfare  = 0;


  // Linear part + initialization
  for (unsigned int p = 0; p < this->getNumPlayers(); ++p) {
	 unsigned int              playerVars = this->LeaderObjective.at(p)->c.size();
	 std::vector<unsigned int> xOfI;
	 for (unsigned int i = this->LeaderLocations.at(p), v = 0;
			i < this->LeaderLocations.at(p) + playerVars;
			++i, ++v) {
		//
		xOfI.push_back(i);
		linearWelfare += this->LCPModel->getVarByName("x_" + std::to_string(i)) *
							  this->LeaderObjective.at(p)->c.at(v);
	 }
	 xOfIs.push_back(xOfI);
  }
  // Account for market clearing variables

  for (unsigned int p = 0; p < this->getNumPlayers(); ++p) {
	 // For each player, we build the xMinusI vector
	 std::vector<unsigned int> xMinusI;
	 for (unsigned int o = 0; o < this->getNumPlayers(); ++o) {
		if (p != o) {
		  for (unsigned int i = 0; i < this->LeaderObjective.at(o)->c.size(); ++i)
			 xMinusI.push_back(xOfIs.at(o).at(i));
		}
	 }

	 // Get the MC variables
	 for (unsigned int j = 0; j < this->numMCVariables; j++)
		xMinusI.push_back(this->NumVariables - this->numMCVariables + j);

	 xMinusIs.push_back(xMinusI);
  }

  // Now we can build the objective
  for (unsigned int p = 0; p < this->getNumPlayers(); ++p) {
	 GRBQuadExpr interact = 0;
	 for (arma::sp_mat::const_iterator it = this->LeaderObjective.at(p)->C.begin();
			it != this->LeaderObjective.at(p)->C.end();
			++it) {
		unsigned int xPlayer = xOfIs.at(p).at(it.row());
		unsigned int xOther  = xMinusIs.at(p).at(it.col());
		interact += *it * this->LCPModel->getVarByName("x_" + std::to_string(xPlayer)) *
						this->LCPModel->getVarByName("x_" + std::to_string(xOther));
	 }
	 quadrWelfare += interact;
  }

  if (quadratic) {
	 if (linear) { // both linear and quadratic
		LOG_S(INFO) << "Game::EPEC::setWelfareObjective: Setting linear+quadratic objective.";
		this->LCPModel->setObjective(linearWelfare + quadrWelfare, GRB_MINIMIZE);
	 } else { // just quadratic
		LOG_S(INFO) << "Game::EPEC::setWelfareObjective: Setting quadratic objective.";
		this->LCPModel->setObjective(quadrWelfare, GRB_MINIMIZE);
	 }
	 this->LCPModel->set(GRB_IntParam_NonConvex, 2);
  } else {
	 // Then just linear
	 LOG_S(INFO) << "Game::EPEC::setWelfareObjective: Setting linear objective.";
	 this->LCPModel->setObjective(linearWelfare, GRB_MINIMIZE);
  }
}
/**
 * Given that Game::EPEC::PlayersQP are all filled with a each country's
 * MathOpt::QP_Param problem (either exact or approximate), computes the Nash
 * equilibrium.
 * @p pureNE checks for pure Nash Equilibrium. It does not work with
 * EPEC::Algorithms::CutAndPlay
 * @p localTimeLimit sets the timelimit for the solver. a negative value is infinite time
 * @p check If true, the method calls the isSolved() method related to the active algorithm
 * EPEC::Algorithms
 * @p linearWelfare If true, the objective of the resulting LCP includes the sum of the linear
 * objectives for the players
 * @p quadraticWelfare If true, the objective of the resulting LCP includes the sum of the quadratic
 * objectives for the players
 * @returns true if a Nash equilibrium is found
 */
bool Game::EPEC::computeNashEq(bool   pureNE,
										 double localTimeLimit   = -1,
										 bool   check            = false,
										 bool   linearWelfare    = false,
										 bool   quadraticWelfare = false) {

  // Make the Nash Game between countries
  this->NashEquilibrium = false;
  LOG_S(1) << "Game::EPEC::computeNashEq: Making the Master LCP";
  this->makeTheLCP();
  LOG_S(1) << "Game::EPEC::computeNashEq: Made the Master LCP";

  auto solver = this->Stats.AlgorithmData.LCPSolver.get();
  if (solver == Data::LCP::Algorithms::PATH) {
	 LOG_S(WARNING) << "Game::EPEC::computeNashEq: Cannot use PATH fallback (EPEC Market Clearings "
							 "cannot be handeled). Using MIP";
	 solver = Data::LCP::Algorithms::MIP;
  }

  /*
	* In these cases, we can only use a MIP solver to get multiple solutions or PNEs
	*/

  unsigned int MIPWorkers = 1;
  if (solver == Data::LCP::Algorithms::MIP || solver == Data::LCP::Algorithms::MINLP) {
	 if (this->Stats.AlgorithmData.Threads.get() >= 8) {
		int wrk    = static_cast<int>(std::round(std::floor(this->Stats.AlgorithmData.Threads.get() / 4)));
		MIPWorkers = std::max(wrk, 1);
		LOG_S(INFO) << "Game::EPEC::computeNashEq: ConcurrentMIP set to " << MIPWorkers << ".";
	 }
  }
  bool multipleNE = check;
  if (check &&
		this->Stats.AlgorithmData.Algorithm.get() == Data::EPEC::Algorithms::OuterApproximation) {
	 LOG_S(WARNING) << "Game::EPEC::computeNashEq: (check flag is "
							 "true) Cannot search fore multiple NE with the CutAndPlay.";
	 multipleNE = false;
  }
  if (pureNE) {
	 LOG_S(INFO) << "Game::EPEC::computeNashEq: (PureNashEquilibrium flag is "
						 "true) Searching for a pure NE.";
	 if (this->Stats.AlgorithmData.Algorithm.get() != Data::EPEC::Algorithms::OuterApproximation)
		this->Algorithm->makeThePureLCP();
	 else
		LOG_S(WARNING) << "Game::EPEC::computeNashEq: (PureNashEquilibrium flag is "
								"true) Cannot search fore pure NE with the CutAndPlay.";
  }

  if (this->TheLCP->getNumRows() > 250000) {
	 LOG_S(WARNING) << "Game::EPEC::computeNashEq: Too many complementarities. Aborting";
	 this->Stats.Status.set(ZEROStatus::Numerical);
	 return false;
  }

  this->LCPModel =
		this->TheLCP->LCPasMIP(false, localTimeLimit, MIPWorkers, multipleNE ? GRB_MAXINT : 1);


  this->setWelfareObjective(linearWelfare, quadraticWelfare);
  try {
	 this->LCPModel->set(GRB_IntParam_OutputFlag, 1);
	 this->LCPModel->set(GRB_IntParam_NumericFocus, 1);
	 this->LCPModel->optimize();
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
  try { // Try finding a Nash equilibrium for the approximation
	 this->NashEquilibrium =
		  this->TheLCP->extractSols(this->LCPModel.get(), SolutionZ, SolutionX, true);
  } catch (GRBException &e) {
	 throw ZEROException(e);
  }
  if (this->NashEquilibrium) { // If a Nash equilibrium is found, then update
										 // appropriately
	 if (multipleNE) {
		int scount = this->LCPModel->get(GRB_IntAttr_SolCount);
		LOG_S(INFO) << "Game::EPEC::computeNashEq: number of equilibria is " << scount;
		for (int k = 0, stop = 0; k < scount && stop == 0; ++k) {
		  this->LCPModel->set(GRB_IntParam_SolutionNumber, k);
		  this->NashEquilibrium =
				this->TheLCP->extractSols(this->LCPModel.get(), this->SolutionZ, this->SolutionX, true);
		  if (this->Algorithm->isSolved()) {
			 LOG_S(INFO) << "Game::EPEC::computeNashEq: an "
								 "Equilibrium has been found";
			 stop = 1;
		  }
		}
	 } else {
		this->NashEquilibrium = true;
		LOG_S(INFO) << "Game::EPEC::computeNashEq: An Equilibrium has been found (Status: "
						<< this->LCPModel->get(GRB_IntAttr_Status) << ")";
	 }

  } else {
	 LOG_S(INFO) << "Game::EPEC::computeNashEq: no equilibrium has been found.";
	 int status = this->LCPModel->get(GRB_IntAttr_Status);
	 if (status == GRB_TIME_LIMIT)
		this->Stats.Status.set(ZEROStatus::TimeLimit);
	 else
		this->Stats.Status.set(ZEROStatus::NashEqNotFound);
  }
  return this->NashEquilibrium;
}

/**
 * @brief Warmstart the solution with @p x
 * @todo Complete this implementation
 * @param x The warmstart solution
 * @return True if the warmstart was successful
 */
bool Game::EPEC::warmstart(const arma::vec &x) {

  ZEROAssert(x.size() >= this->getNumVar());
  ZEROAssert(this->Finalized);
  ZEROAssert(this->PlayersQP.front() == nullptr);

  this->SolutionX                  = x;
  std::vector<arma::vec> devns     = std::vector<arma::vec>(this->NumPlayers);
  std::vector<arma::vec> prevDevns = std::vector<arma::vec>(this->NumPlayers);
  this->makePlayersQPs();

  arma::vec devn;

  if (this->Algorithm->isSolved())
	 LOG_S(WARNING) << "Game::EPEC::warmstart: "
							 "The loaded solution is optimal.";
  else
	 LOG_S(WARNING) << "Game::EPEC::warmstart: "
							 "The loaded solution is NOT optimal. Trying to repair.";
  return true;
}

/**
 * @brief Call the delegated Algorithm method and return true if the equilibrium is pure
 * @param tol A numerical tolerance
 * @return True if the incumbent equilibrium is pure
 */
bool Game::EPEC::isPureStrategy(double tol) const { return this->Algorithm->isPureStrategy(tol); }
/**
 * @brief Call the delegated Algorithm method and return true if there exist a feasible equilibrium
 * @param tol A numerical tolerance
 * @return True if the incumbent solution is an equilibrium
 */
bool Game::EPEC::isSolved(double tol) const { return this->Algorithm->isSolved(); }

/**
 * @brief Computes Nash equilibrium using the Algorithm set in
 * Game::EPEC::Algorithm. Checks the value of Game::EPEC::Algorithm and delegates the task to
 * appropriate Algorithm wrappers.
 */
void Game::EPEC::findNashEq() {

  ZEROAssert(this->Finalized);

  if (this->Stats.Status.get() != ZEROStatus::Uninitialized) {
	 LOG_S(ERROR) << "Game::EPEC::findNashEq: a Nash Eq was "
						  "already found. Calling this findNashEq might lead to errors!";
  }


  std::stringstream final_msg;

  switch (this->Stats.AlgorithmData.Algorithm.get()) {

  case Data::EPEC::Algorithms::InnerApproximation: {
	 final_msg << "Inner approximation: run completed. ";
	 this->Algorithm = std::shared_ptr<Algorithms::EPEC::PolyBase>(
		  new class Algorithms::EPEC::InnerApproximation(this->Env, this));
	 this->Algorithm->solve();
  } break;

  case Data::EPEC::Algorithms::CombinatorialPne: {
	 final_msg << "CombinatorialPNE: run completed. ";
	 this->Algorithm = std::shared_ptr<Algorithms::EPEC::PolyBase>(
		  new class Algorithms::EPEC::CombinatorialPNE(this->Env, this));
	 this->Algorithm->solve();
  } break;

  case Data::EPEC::Algorithms::OuterApproximation: {
	 final_msg << "Cut-and-Play: run completed. ";
	 this->Algorithm = std::shared_ptr<Algorithms::EPEC::PolyBase>(
		  new class Algorithms::EPEC::CutAndPlay(this->Env, this));
	 this->Algorithm->solve();
  } break;

  case Data::EPEC::Algorithms::FullEnumeration: {
	 final_msg << "Full enumeration: run completed. ";
	 this->Algorithm = std::shared_ptr<Algorithms::EPEC::PolyBase>(
		  new class Algorithms::EPEC::FullEnumeration(this->Env, this));
	 this->Algorithm->solve();
  } break;
  }
  const std::chrono::duration<double> timeElapsed =
		std::chrono::high_resolution_clock::now() - this->InitTime;
  this->Stats.WallClockTime.set(timeElapsed.count() * std::chrono::milliseconds::period::num /
										  std::chrono::milliseconds::period::den);

  // Handing EPECStatistics object to track performance of algorithm
  if (this->LCPModel) {
	 this->Stats.NumVar         = this->LCPModel->get(GRB_IntAttr_NumVars);
	 this->Stats.NumConstraints = this->LCPModel->get(GRB_IntAttr_NumConstrs);
	 this->Stats.NumNonZero     = this->LCPModel->get(GRB_IntAttr_NumNZs);
  } // Assigning appropriate Status messages after solving

  switch (this->Stats.Status.get()) {
  case ZEROStatus::NashEqNotFound:
	 final_msg << "No Nash equilibrium exists.";
	 break;
  case ZEROStatus::NashEqFound: {
	 final_msg << "Found a Nash equilibrium ("
				  << (this->Stats.PureNashEquilibrium.get() == 0 ? "MNE" : "PNE") << ").";
  } break;
  case ZEROStatus::TimeLimit:
	 final_msg << "Nash equilibrium not found. The time limit was hit.";
	 break;
  case ZEROStatus::Numerical:
	 final_msg << "Nash equilibrium not found. The Numerical issues flag was triggered.";
	 break;
  default:
	 final_msg << "Nash equilibrium not found. Unknown status.";
	 break;
  }
  LOG_S(INFO) << "Game::EPEC::findNashEq: " << final_msg.str();
}


/**
 * @brief Decides the Algorithm to be used for solving the given instance of the
 * problem. The choice of algorithms are documented in Game::EPECalgorithm
 * @param algorithm The input algorithm
 */
void Game::EPEC::setAlgorithm(Data::EPEC::Algorithms algorithm) {
  this->Stats.AlgorithmData.Algorithm.set(algorithm);
}

/**
 * @brief Decides the Algorithm to be used for recovering a PNE out of the
 * InnerApproximation procedure.
 * @param strategy Input ecovery strategy
 */
void Game::EPEC::setRecoverStrategy(Data::EPEC::RecoverStrategy strategy) {
  this->Stats.AlgorithmData.RecoverStrategy.set(strategy);
}

/**
 * @brief Gets the position of the j-th variable in the i-th leader
 * Querying Game::EPEC::LCPModel for x[return-value] variable gives the
 * appropriate variable
 * @param i Leader number
 * @param j Follower number
 * @return The queried position
 */
unsigned int Game::EPEC::getPositionLeadFoll(const unsigned int i, const unsigned int j) const {
  ZEROAssert(i < this->NumPlayers);
  const auto LeaderStart = this->TheNashGame->getPrimalLoc(i);
  return LeaderStart + j;
}

/**
 * @brief Gets the position of the j-th Follower variable in the i-th leader
 * Querying Game::EPEC::LCPModel for x[return-value] variable gives the
 * appropriate variable
 * @param i Leader number
 * @param j Follower number
 * @return The queried position
 */
unsigned int Game::EPEC::getPositionLeadLead(const unsigned int i, const unsigned int j) const {
  ZEROAssert(i < this->NumPlayers);
  const auto LeaderStart = this->TheNashGame->getPrimalLoc(i);
  return LeaderStart + this->PlayersLCP.at(i)->getLStart() + j;
}

/**
 * @brief Gets the value of the j-th variable in i-th leader
 * @param i Leader number
 * @param j Follower number
 * @return The queried position
 */
double Game::EPEC::getValLeadFoll(const unsigned int i, const unsigned int j) const {
  ZEROAssert(i < this->NumPlayers);
  ZEROAssert(this->LCPModel != nullptr);

  return this->LCPModel->getVarByName("x_" + std::to_string(this->getPositionLeadFoll(i, j)))
		.get(GRB_DoubleAttr_X);
}

/**
 * @brief Returns the indices of polyhedra feasible for the leader, from which strategies are played
 * with probability greater than the tolerance.
 * @param i Leader index
 * @param tol Tolerance
 * @return The queried attribute
 */
std::vector<unsigned int> Game::EPEC::mixedStrategyPoly(unsigned int i, double tol) const {
  ZEROAssert(i < this->NumPlayers);
  ZEROAssert(this->LCPModel != nullptr);
  return this->Algorithm->mixedStrategyPoly(i, tol);
}

/**
 * @brief Returns the probability associated with the k -th polyhedron of the i -th leader.
 * @param i Leader index
 * @param k The polyhedron index
 * @return The queried attribute
 */
double Game::EPEC::getValProbab(unsigned int i, unsigned int k) {
  ZEROAssert(i < this->NumPlayers);
  ZEROAssert(this->LCPModel != nullptr);
  return this->Algorithm->getValProbab(i, k);
}

/**
 * @brief For the @p i -th leader, gets the @p k -th pure strategy at position @p j
 * @param i The leader index
 * @param j The position index
 * @param k The pure strategy index
 * @param tol A numerical tolerance
 * @return The queried attribute
 */
double Game::EPEC::getValLeadFollPoly(const unsigned int i,
												  const unsigned int j,
												  const unsigned int k,
												  const double       tol) const {
  ZEROAssert(i < this->NumPlayers);
  ZEROAssert(this->LCPModel != nullptr);
  return this->Algorithm->getValLeadFollPoly(i, j, k, tol);
}

/**
 * @brief For the @p i -th leader, gets the @p k -th pure strategy at leader position @p j
 * @param i The leader index
 * @param j The position index
 * @param k The pure strategy index
 * @param tol A numerical tolerance
 * @return The queried attribute
 */
double Game::EPEC::getValLeadLeadPoly(const unsigned int i,
												  const unsigned int j,
												  const unsigned int k,
												  const double       tol) const {
  ZEROAssert(i < this->NumPlayers);
  ZEROAssert(this->LCPModel != nullptr);
  return this->Algorithm->getValLeadLeadPoly(i, j, k, tol);
}

/**
 * @brief Gets the value of the j-th non-follower variable in i-th leader
 * @param i Leader number
 * @param j Follower number
 * @return The queried position
 */
double Game::EPEC::getValLeadLead(const unsigned int i, const unsigned int j) const {
  ZEROAssert(i < this->NumPlayers);
  ZEROAssert(this->LCPModel != nullptr);
  return this->LCPModel->getVarByName("x_" + std::to_string(this->getPositionLeadLead(i, j)))
		.get(GRB_DoubleAttr_X);
}
/**
 * @brief Sets the branching strategy for the CutAndPlay
 * @param strategy The input strategy
 */
void Game::EPEC::setBranchingStrategy(Data::EPEC::BranchingStrategy strategy) {
  this->Stats.AlgorithmData.BranchingStrategy.set(strategy);
}

/**
 * @brief Convert a Data::EPEC::Algorithms into a string
 * @param al The input algorithm
 * @return  The string literal
 */
std::string std::to_string(const Data::EPEC::Algorithms al) {
  switch (al) {
  case Data::EPEC::Algorithms::FullEnumeration:
	 return std::string("FullEnumeration");
  case Data::EPEC::Algorithms::InnerApproximation:
	 return std::string("InnerApproximation");
  case Data::EPEC::Algorithms::CombinatorialPne:
	 return std::string("CombinatorialPNE");
  case Data::EPEC::Algorithms::OuterApproximation:
	 return std::string("CutAndPlay");
  }
  return "";
}

/**
 * @brief Convert a Data::EPEC::RecoverStrategy into a string
 * @param strategy  The input strategy
 * @return  The string literal
 */

std::string std::to_string(const Data::EPEC::RecoverStrategy strategy) {
  switch (strategy) {
  case Data::EPEC::RecoverStrategy::IncrementalEnumeration:
	 return std::string("IncrementalEnumeration");
  case Data::EPEC::RecoverStrategy::Combinatorial:
	 return std::string("Combinatorial");
  }
  return "";
}

/**
 * @brief Convert a Data::EPEC::BranchingStrategy into a string
 * @param strategy The input strategy
 * @return The string literal
 */
std::string std::to_string(const Data::EPEC::BranchingStrategy strategy) {
  switch (strategy) {
  case Data::EPEC::BranchingStrategy::HybridBranching:
	 return std::string("HybridBranching");
  case Data::EPEC::BranchingStrategy::DeviationBranching:
	 return std::string("DeviationBranching");
  }
  return "";
}
