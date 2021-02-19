/* #############################################
 *             This file is part of
 *                    ZERO
 *
 *             Copyright (c) 2020
 *     Released under the Creative Commons
 *        Zero v1.0 Universal License
 *
 *              Find out more at
 *        https://github.com/ds4dm/ZERO
 * #############################################*/

#include "games/epec.h"
#include "games/algorithms/EPEC/epec_polybase.h"
#include "zero.h"
#include <armadillo>
#include <boost/log/trivial.hpp>
#include <iostream>
#include <memory>

void Game::EPEC::preFinalize()
/**
  @brief Empty function - optionally re-implementable in derived class
@details This function can be optionally implemented by
 the derived class. Code in this class will be run <i>before</i>
 calling Game::EPEC::finalize().
*/
{}

void Game::EPEC::postFinalize()
/**
  @brief Empty function - optionally reimplementable in derived class
@details This function can be optionally implemented by
 the derived class. Code in this class will be run <i>after</i>
 calling Game::EPEC::finalize().
*/
{}

void Game::EPEC::finalize()
/**
 * @brief Finalizes the creation of a Game::EPEC object.
 * @details Performs a bunch of job after all data for a Game::EPEC object are
 * given, namely.
 * Models::EPEC::computeLeaderLocations -	Adds the required dummy
 * variables to each leader's problem so that a game among the leaders can be
 * defined. Calls Game::EPEC::addDummyLead
 * 	-	Makes the market clearing constraint in each country. Calls
 */
{
  if (this->Finalized)
	 std::cerr << "Warning in Game::EPEC::finalize: Model already Finalized\n";

  this->NumPlayers = static_cast<int>(this->PlayersLowerLevels.size());
  ;
  /// Game::EPEC::preFinalize() can be overridden, and that code will run before
  /// calling Game::EPEC::finalize()
  this->preFinalize();

  this->ConvexHullVariables = std::vector<unsigned int>(this->NumPlayers, 0);
  this->Stats.AlgorithmData.FeasiblePolyhedra.set(std::vector<unsigned int>(this->NumPlayers, 0));
  this->computeLeaderLocations(this->numMCVariables);
  // Initialize leader objective and PlayersQP
  this->LeaderObjective           = std::vector<std::shared_ptr<MathOpt::QP_Objective>>(NumPlayers);
  this->LeaderObjectiveConvexHull = std::vector<std::shared_ptr<MathOpt::QP_Objective>>(NumPlayers);
  this->PlayersQP                 = std::vector<std::shared_ptr<MathOpt::QP_Param>>(NumPlayers);
  this->PlayersLCP                = std::vector<std::shared_ptr<MathOpt::LCP>>(NumPlayers);
  this->SizesWithoutHull          = std::vector<unsigned int>(NumPlayers, 0);

  for (unsigned int i = 0; i < this->NumPlayers; i++) {
	 this->addDummyLead(i);
	 this->LeaderObjective.at(i)           = std::make_shared<MathOpt::QP_Objective>();
	 this->LeaderObjectiveConvexHull.at(i) = std::make_shared<MathOpt::QP_Objective>();
	 this->makeObjectivePlayer(i, *this->LeaderObjective.at(i).get());
	 // this->PlayersLCP.at(i) =std::shared_ptr<MathOpt::PolyLCP>(new
	 // PolyLCP(this->Env,*this->PlayersLowerLevels.at(i).get()));
	 this->SizesWithoutHull.at(i) = *this->LocEnds.at(i);
  }

  this->Finalized = true;

  /// Game::EPEC::postFinalize() can be overridden, and that code will run after
  /// calling Game::EPEC::finalize()
  this->postFinalize();
}

void Game::EPEC::addDummyLead(
	 const unsigned int i ///< The leader to whom dummy variables should be added
) {
  /// Adds dummy variables to the leader of an EPEC - useful after computing the
  /// convex hull.
  const unsigned int nEPECvars        = this->NumVariables;
  const unsigned int nThisCountryvars = *this->LocEnds.at(i);
  // this->Locations.at(i).at(Models::LeaderVars::End);

  if (nEPECvars < nThisCountryvars)
	 throw ZEROException(ZEROErrorCode::Assertion,
								"Mismatch between variable counts " + std::to_string(nEPECvars) + " and " +
									 std::to_string(nThisCountryvars));

  this->PlayersLowerLevels.at(i).get()->addDummy(nEPECvars - nThisCountryvars);
}

void Game::EPEC::computeLeaderLocations(const unsigned int addSpaceForMC) {
  this->LeaderLocations       = std::vector<unsigned int>(this->NumPlayers);
  this->LeaderLocations.at(0) = 0;
  for (unsigned int i = 1; i < this->NumPlayers; i++) {
	 this->LeaderLocations.at(i) = this->LeaderLocations.at(i - 1) + *this->LocEnds.at(i - 1);
  }
  this->NumVariables = this->LeaderLocations.back() + *this->LocEnds.back() + addSpaceForMC;
}

void Game::EPEC::getXMinusI(const arma::vec &x, const unsigned int &i, arma::vec &solOther) const {
  const unsigned int nEPECvars            = this->NumVariables;
  const unsigned int nThisCountryvars     = *this->LocEnds.at(i);
  const unsigned int nThisCountryHullVars = this->ConvexHullVariables.at(i);
  const unsigned int nConvexHullVars      = static_cast<const unsigned int>(
      std::accumulate(this->ConvexHullVariables.rbegin(), this->ConvexHullVariables.rend(), 0));

  solOther.zeros(nEPECvars -        // All variables in EPEC
					  nThisCountryvars - // Subtracting this country's variables,
					  // since we only want others'
					  nConvexHullVars +      // We don't want any convex hull variables
					  nThisCountryHullVars); // We subtract the hull variables
													 // associated to the ith player
													 // convex hull vars, since we double subtracted

  for (unsigned int j = 0, count = 0, current = 0; j < this->NumPlayers; ++j) {
	 if (i != j) {
		current = *this->LocEnds.at(j) - this->ConvexHullVariables.at(j);
		solOther.subvec(count, count + current - 1) =
			 x.subvec(this->LeaderLocations.at(j), this->LeaderLocations.at(j) + current - 1);
		count += current;
	 }
  }
  // We need to keep track of MC_vars also for this country
  for (unsigned int j = 0; j < this->numMCVariables; j++)
	 solOther.at(solOther.n_rows - this->numMCVariables + j) =
		  x.at(this->NumVariables - this->numMCVariables + j);
}

void Game::EPEC::getXofI(const arma::vec &   x,
								 const unsigned int &i,
								 arma::vec &         solI,
								 bool                hull) const {
  /**
	* Given the player id @p i and the solution @p x, the method returns in @p
	* xWithoutHull the x vector for the given player, with the convex-hull's
	* variables in case @p hull is false. Also, no MC variables are included
	*/
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
  solI.zeros(vars);
  solI.subvec(0, vars - 1) =
		x.subvec(this->LeaderLocations.at(i), this->LeaderLocations.at(i) + current - 1);
}

void Game::EPEC::getXWithoutHull(const arma::vec &x, arma::vec &xWithoutHull) const {
  /**
	* Given the the solution @p x, the method returns in @p
	* xWithoutHull the x vector without the convex-hull's
	* variables.  Also, no MC variables are included
	*
	*/
  const unsigned int nEPECvars       = this->NumVariables;
  const unsigned int nConvexHullVars = static_cast<const unsigned int>(
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
}

std::unique_ptr<GRBModel> Game::EPEC::respond(const unsigned int i, const arma::vec &x) const {
  if (!this->Finalized)
	 throw ZEROException(ZEROErrorCode::Assertion, "The model was not Finalized");

  if (i >= this->NumPlayers)
	 throw ZEROException(ZEROErrorCode::OutOfRange, "Country number is invalid");

  arma::vec solOther;
  this->getXMinusI(x, i, solOther);
  if (this->LeaderObjective.at(i)->Q.n_nonzero > 0)
	 return this->PlayersLCP.at(i).get()->MPECasMIQP(this->LeaderObjective.at(i)->Q,
																	 this->LeaderObjective.at(i)->C,
																	 this->LeaderObjective.at(i)->c,
																	 solOther,
																	 true);
  else
	 return this->PlayersLCP.at(i).get()->MPECasMILP(
		  this->LeaderObjective.at(i)->C, this->LeaderObjective.at(i)->c, solOther, true);
}

double
Game::EPEC::respondSol(arma::vec &      sol,    ///< [out] Optimal response
							  unsigned int     player, ///< Player whose optimal response is to be computed
							  const arma::vec &x, ///< A std::vector of pure strategies (either for all
							  ///< players or all other players
							  const arma::vec &prevDev
							  ///< [in] if any, the std::vector of previous deviations.
) const {
  /**
	* @brief Returns the optimal objective value that is obtainable for the
	* player @p player given the decision @p x of all other players.
	* @details
	* Calls Game::EPEC::respond and obtains the std::unique_ptr to GRBModel of
	* best response by player @p player. Then solves the model and returns the
	* appropriate objective value.
	* @returns The optimal objective value for the player @p player.
	*/
  auto model = this->respond(player, x);
  LOG_S(1) << "Game::EPEC::respondSol: Writing dat/RespondSol" + std::to_string(player) +
						".lp to disk";
  // model->write("dat/RespondSol" + std::to_string(player) + ".lp");
  const int status = model->get(GRB_IntAttr_Status);
  if (status == GRB_UNBOUNDED || status == GRB_OPTIMAL) {
	 unsigned int Nx = this->PlayersLCP.at(player)->getNumCols();
	 sol.zeros(Nx);
	 for (unsigned int i = 0; i < Nx; ++i)
		sol.at(i) = model->getVarByName("x_" + std::to_string(i)).get(GRB_DoubleAttr_X);

	 if (status == GRB_UNBOUNDED) {
		LOG_S(WARNING) << "Game::EPEC::respondSol: deviation is "
								"unbounded.";
		GRBLinExpr obj = 0;
		model->setObjective(obj);
		model->optimize();
		if (!prevDev.empty()) {
		  LOG_S(1) << "Generating an improvement basing on the extreme ray.";
		  // Fetch objective function coefficients
		  GRBQuadExpr QuadObj = model->getObjective();
		  arma::vec   objcoeff;
		  for (unsigned int i = 0; i < QuadObj.size(); ++i)
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
	 if (status == GRB_OPTIMAL) {
		return model->get(GRB_DoubleAttr_ObjVal);
	 }
  } else {
	 return GRB_INFINITY;
  }
  return GRB_INFINITY;
}

const void Game::EPEC::makePlayerQP(const unsigned int i)
/**
 * @brief Makes the MathOpt::QP_Param corresponding to the @p i-th country.
 * @details
 *  - First gets the Game::LCP object from @p Game::EPEC::PlayersLowerLevels and
 * makes a MathOpt::QP_Param with this LCP as the lower level
 *  - This is achieved by calling LCP::makeQP and using the objective value
 * object in @p Game::EPEC::LeaderObjective
 *  - Finally the locations are updated owing to the complete convex hull
 * calculated during the call to LCP::makeQP
 * @note Overloaded as Models::EPEC::makePlayerQP()
 */
{
  // LOG_S(INFO) << "Starting Convex hull computation of the country
  // "
  // << this->AllLeadPars[i].name << '\n';
  if (!this->Finalized)
	 throw ZEROException(ZEROErrorCode::Assertion, "The model was not Finalized");
  if (i >= this->NumPlayers)
	 throw ZEROException(ZEROErrorCode::OutOfRange, "The player id is out of range");
  // if (!this->PlayersQP.at(i).get())
  {
	 this->PlayersQP.at(i)     = std::make_shared<MathOpt::QP_Param>(this->Env);
	 const auto &origLeadObjec = *this->LeaderObjective.at(i).get();

	 this->LeaderObjectiveConvexHull.at(i).reset(
		  new MathOpt::QP_Objective{origLeadObjec.Q, origLeadObjec.C, origLeadObjec.c});
	 this->PlayersLCP.at(i)->makeQP(*this->LeaderObjectiveConvexHull.at(i).get(),
											  *this->PlayersQP.at(i).get());
  }
}

void Game::EPEC::makePlayersQPs()
/**
 * @brief Makes the MathOpt::QP_Param for all the countries
 * @details
 * Calls are made to Models::EPEC::makePlayerQP(const unsigned int i) for
 * each valid @p i
 * @note Overloaded as EPEC::makePlayerQP(unsigned int)
 */
{
  for (unsigned int i = 0; i < this->NumPlayers; ++i) {
	 this->Game::EPEC::makePlayerQP(i);
  }
  for (unsigned int i = 0; i < this->NumPlayers; ++i) {
	 // LeadLocs &Loc = this->Locations.at(i);
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
		for (unsigned int j = 0; j < this->NumPlayers; j++) {
		  if (i != j) {
			 this->PlayersQP.at(j)->addDummy(
				  convHullVarCount,
				  0,
				  this->PlayersQP.at(j)->getNx() -
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

void Game::EPEC::makeTheLCP() {
  if (this->PlayersQP.front() == nullptr) {
	 LOG_S(ERROR) << "Exception in Game::EPEC::makeTheLCP : "
						  "no country QP has been "
						  "made."
					  << '\n';
	 throw ZEROException(ZEROErrorCode::Assertion, "No country QP has been made");
  }
  // Preliminary set up to get the LCP ready
  unsigned long int Nvar = this->PlayersQP.front()->getNx() + this->PlayersQP.front()->getNy();
  arma::sp_mat      MC(0, Nvar), dumA(0, Nvar);
  arma::vec         MCRHS, dumb;
  MCRHS.zeros(0);
  dumb.zeros(0);
  this->makeMCConstraints(MC, MCRHS);
  LOG_S(1) << "Game::EPEC::makeTheLCP(): Market Clearing "
				  "constraints are ready";
  std::vector<std::shared_ptr<MathOpt::MP_Param>> MPCasted;
  for (auto &item : this->PlayersQP) {
	 auto m = std::dynamic_pointer_cast<MathOpt::MP_Param>(item);
	 MPCasted.push_back(m);
  }
  this->TheNashGame = std::unique_ptr<Game::NashGame>(
		new Game::NashGame(this->Env, MPCasted, MC, MCRHS, 0, dumA, dumb));
  LOG_S(1) << "Game::EPEC::makeTheLCP(): NashGame is ready";
  this->TheLCP = std::unique_ptr<MathOpt::LCP>(new MathOpt::LCP(this->Env, *TheNashGame));
  LOG_S(1) << "Game::EPEC::makeTheLCP(): LCP is ready";


  LOG_S(1) << *TheNashGame;
}


/**
 * Given the object Game::EPEC::LCPModel, this method sets its objective to the social welfare. In
 * specific, we can decide whether to include or not the linear or quadratic part of the welfare
 * @param linear True if the linear part is to be included, false otherwise.
 * @param quadratic True if the quadratic part is to be included, false otherwise.
 */
void Game::EPEC::setWelfareObjective(bool linear = true, bool quadratic = true) {

  if (!linear && !quadratic) {
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
 * @p pureNE checks for pure Nash Equilibria. It does not work with
 * EPEC::Algorithms::OuterApproximation
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
		int wrk    = std::round(std::floor(this->Stats.AlgorithmData.Threads.get() / 4));
		MIPWorkers = std::max(wrk, 1);
		LOG_S(INFO) << "Game::EPEC::computeNashEq: ConcurrentMIP set to " << MIPWorkers << ".";
	 }
  }
  bool multipleNE = check;
  if (check &&
		this->Stats.AlgorithmData.Algorithm.get() == Data::EPEC::Algorithms::OuterApproximation) {
	 LOG_S(WARNING) << "Game::EPEC::computeNashEq: (check flag is "
							 "true) Cannot search fore multiple NE with the OuterApproximation.";
	 multipleNE = false;
  }
  if (pureNE) {
	 LOG_S(INFO) << "Game::EPEC::computeNashEq: (PureNashEquilibrium flag is "
						 "true) Searching for a pure NE.";
	 if (this->Stats.AlgorithmData.Algorithm.get() != Data::EPEC::Algorithms::OuterApproximation)
		static_cast<Algorithms::EPEC::PolyBase *>(this->Algorithm.get())->makeThePureLCP();
	 else
		LOG_S(WARNING) << "Game::EPEC::computeNashEq: (PureNashEquilibrium flag is "
								"true) Cannot search fore pure NE with the OuterApproximation.";
  }

  this->LCPModel =
		this->TheLCP->LCPasMIP(false, localTimeLimit, MIPWorkers, multipleNE ? GRB_MAXINT : 1);


  this->setWelfareObjective(linearWelfare, quadraticWelfare);
  try {
	 this->LCPModel->set(GRB_IntParam_OutputFlag, 1);
	 this->LCPModel->write("dat/TheLCP.lp");
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
		LOG_S(INFO) << "Game::EPEC::computeNashEq: An Equilibrium has been found";
	 }

  } else {
	 LOG_S(INFO) << "Game::EPEC::computeNashEq: no equilibrium has been found.";
	 int status = this->LCPModel->get(GRB_IntAttr_Status);
	 if (status == GRB_TIME_LIMIT)
		this->Stats.Status = ZEROStatus::TimeLimit;
	 else
		this->Stats.Status = ZEROStatus::NashEqNotFound;
  }
  return this->NashEquilibrium;
}

bool Game::EPEC::warmstart(const arma::vec x) {
  //@todo complete implementation

  if (x.size() < this->getNumVar())
	 throw ZEROException(ZEROErrorCode::OutOfRange,
								"The number of variables does not fit the instance");

  if (!this->Finalized) {
	 throw ZEROException(ZEROErrorCode::Assertion, "The EPEC was not Finalized");
  }
  if (this->PlayersQP.front() == nullptr) {
	 LOG_S(WARNING) << "Game::EPEC::warmstart: Generating QP as of warmstart.";
  }

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
  /// @todo Game::EPEC::warmstart - to complete implementation?
  return true;
}
bool Game::EPEC::isPureStrategy(double tol) const {
  /**
	* @brief Call the delegated Algorithm and return true if the equilibrium is
	* pure
	*/
  return this->Algorithm->isPureStrategy(tol);
}
bool Game::EPEC::isSolved(double tol) const {
  /**
	* @brief Call the delegated Algorithm and return true if the EPEC has been
	* solved.
	*/
  return this->Algorithm->isSolved(tol);
}

const void Game::EPEC::findNashEq() {
  /**
	* @brief Computes Nash equilibrium using the Algorithm set in
	* Game::EPEC::Algorithm
	* @details
	* Checks the value of Game::EPEC::Algorithm and delegates the task to
	* appropriate Algorithm wrappers.
	*/

  std::stringstream final_msg;
  if (!this->Finalized)
	 throw ZEROException(ZEROErrorCode::Assertion, "The EPEC was not Finalized");

  if (this->Stats.Status.get() != ZEROStatus::Uninitialized) {
	 LOG_S(ERROR) << "Game::EPEC::findNashEq: a Nash Eq was "
						  "already found. Calling this findNashEq might lead to errors!";
  }

  // Choosing the appropriate algorithm
  switch (this->Stats.AlgorithmData.Algorithm.get()) {

  case Data::EPEC::Algorithms::InnerApproximation: {
	 final_msg << "Inner approximation Algorithm completed. ";
	 this->Algorithm = std::shared_ptr<Algorithms::EPEC::PolyBase>(
		  new class Algorithms::EPEC::InnerApproximation(this->Env, this));
	 this->Algorithm->solve();
  } break;

  case Data::EPEC::Algorithms::CombinatorialPne: {
	 final_msg << "CombinatorialPNE Algorithm completed. ";
	 this->Algorithm = std::shared_ptr<Algorithms::EPEC::PolyBase>(
		  new class Algorithms::EPEC::CombinatorialPNE(this->Env, this));
	 this->Algorithm->solve();
  } break;

  case Data::EPEC::Algorithms::OuterApproximation: {
	 final_msg << "Outer approximation Algorithm completed. ";
	 this->Algorithm = std::shared_ptr<Algorithms::EPEC::PolyBase>(
		  new class Algorithms::EPEC::OuterApproximation(this->Env, this));
	 this->Algorithm->solve();
  } break;

  case Data::EPEC::Algorithms::FullEnumeration: {
	 final_msg << "Full enumeration Algorithm completed. ";
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
	 final_msg << "Nash equilibrium not found. Time limit attained";
	 break;
  case ZEROStatus::Numerical:
	 final_msg << "Nash equilibrium not found. Numerical issues might affect "
					  "this result.";
	 break;
  default:
	 final_msg << "Nash equilibrium not found. Time limit attained";
	 break;
  }
  LOG_S(INFO) << "Game::EPEC::findNashEq: " << final_msg.str();
}

void Game::EPEC::setAlgorithm(Data::EPEC::Algorithms algorithm)
/**
 * Decides the Algorithm to be used for solving the given instance of the
 * problem. The choice of algorithms are documented in Game::EPECalgorithm
 */
{
  this->Stats.AlgorithmData.Algorithm.set(algorithm);
}

void Game::EPEC::setRecoverStrategy(Data::EPEC::RecoverStrategy strategy)
/**
 * Decides the Algorithm to be used for recovering a PNE out of the
 * InnerApproximation procedure.
 */
{
  this->Stats.AlgorithmData.RecoverStrategy.set(strategy);
}

unsigned int Game::EPEC::getPositionLeadFoll(const unsigned int i, const unsigned int j) const {
  /**
	* Get the position of the j-th variable in the i-th leader
	* Querying Game::EPEC::LCPModel for x[return-value] variable gives the
	* appropriate variable
	*/
  const auto LeaderStart = this->TheNashGame->getPrimalLoc(i);
  return LeaderStart + j;
}

unsigned int Game::EPEC::getPositionLeadLead(const unsigned int i, const unsigned int j) const {
  /**
	* Get the position of the j-th Follower variable in the i-th leader
	* Querying Game::EPEC::LCPModel for x[return-value] variable gives the
	* appropriate variable
	*/
  const auto LeaderStart = this->TheNashGame->getPrimalLoc(i);
  return LeaderStart + this->PlayersLCP.at(i)->getLStart() + j;
}

double Game::EPEC::getValLeadFoll(const unsigned int i, const unsigned int j) const {
  /**
	* Get the value of the j-th variable in i-th leader
	*/
  if (!this->LCPModel)
	 throw ZEROException(ZEROErrorCode::Assertion, "The LCP was not made nor solved");
  return this->LCPModel->getVarByName("x_" + std::to_string(this->getPositionLeadFoll(i, j)))
		.get(GRB_DoubleAttr_X);
}

double Game::EPEC::getValLeadLead(const unsigned int i, const unsigned int j) const {
  /**
	* Get the value of the j-th non-follower variable in i-th leader
	*/
  if (!this->LCPModel)
	 throw ZEROException(ZEROErrorCode::Assertion, "The LCP was not made nor solved");
  return this->LCPModel->getVarByName("x_" + std::to_string(this->getPositionLeadLead(i, j)))
		.get(GRB_DoubleAttr_X);
}


std::string std::to_string(const Data::EPEC::Algorithms al) {
  switch (al) {
  case Data::EPEC::Algorithms::FullEnumeration:
	 return std::string("FullEnumeration");
  case Data::EPEC::Algorithms::InnerApproximation:
	 return std::string("InnerApproximation");
  case Data::EPEC::Algorithms::CombinatorialPne:
	 return std::string("CombinatorialPNE");
  case Data::EPEC::Algorithms::OuterApproximation:
	 return std::string("OuterApproximation");
  }
  return "";
}

std::string std::to_string(const Data::EPEC::RecoverStrategy strategy) {
  switch (strategy) {
  case Data::EPEC::RecoverStrategy::IncrementalEnumeration:
	 return std::string("IncrementalEnumeration");
  case Data::EPEC::RecoverStrategy::Combinatorial:
	 return std::string("Combinatorial");
  }
  return "";
}
