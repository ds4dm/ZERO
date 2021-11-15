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


#include "games/algorithms/EPEC/epec_polybase.h"


/**
 * @brief Checks if Game::EPEC is solved, otherwise it returns a proof.
 * @details
 * Analogous to Game::NashGame::isSolved but checks if the given Game::EPEC is
 * solved. If it is, then returns true. If not, it returns the country
 * which has a profitable deviation in @p player and the profitable
 * deviation in @p profitableDeviation. @p Tolerance is the tolerance for the
 * check. If the <i> improved objective </i> after the deviation is less than @p
 * Tolerance, then it is not considered as a profitable deviation.
 *
 * Thus we check if the given point is an @f$\epsilon@f$-equilibrium. Value of
 * @f$\epsilon @f$ can be chosen sufficiently close to 0.
 *
 * @warning Setting @p Tolerance = 0 might even reject a real solution as not
 * solved. This is due to Numerical issues arising from the LCP solver (Gurobi).
 * @param player The id of the player
 * @param profitableDeviation An output (possibly non-changed) vector containing the profitable
 * deviation for the given player
 * @param tol A numerical tolerance
 * @return True if there is no profitable deviation, namely the player is optimal
 */
bool Algorithms::EPEC::PolyBase::isSolved(unsigned int *player,
														arma::vec *   profitableDeviation,
														double        tol) const

{

  if (!this->EPECObject->NashEquilibrium)
	 return false;
  tol = this->EPECObject->Stats.AlgorithmData.DeviationTolerance.get();
  if (tol < 0)
	 tol = 1e-5;
  this->EPECObject->TheNashGame->isSolved(
		this->EPECObject->SolutionX, *player, *profitableDeviation, tol);
  arma::vec objvals =
		this->EPECObject->TheNashGame->computeQPObjectiveValues(this->EPECObject->SolutionX, true);
  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i) {
	 double val =
		  this->EPECObject->bestResponse(*profitableDeviation, i, this->EPECObject->SolutionX);
	 if (val == GRB_INFINITY)
		return false;
	 if (!Utils::isEqual(val, objvals.at(i), tol, 1 - tol)) {
		*player = i;
		LOG_S(0) << "Algorithms::EPEC::PolyBase::isSolved: deviation for player " << i << " -- of "
					<< std::abs(val - objvals.at(i));
		return false;
	 }
  }
  return true;
}

/**
 * @brief Checks whether the current Game::EPEC instance is solved for any player, up to a numerical
 * tolerance
 * @param tol The numerical tolerance
 * @return True if the game is solved
 */
bool Algorithms::EPEC::PolyBase::isSolved(double tol) {
  unsigned int countryNumber;
  arma::vec    ProfDevn;
  bool         ret = this->isSolved(&countryNumber, &ProfDevn, tol);
  return ret;
}

/**
 * @brief Get the position of the @p k -th follower variable of the @p i -th leader, in the @p j -th
 * feasible polyhedron.
 * @param i The leader index
 * @param j The polyhedron index
 * @param k The follower variable index
 * @return The position for the queried attribute
 */
unsigned int Algorithms::EPEC::PolyBase::getPositionLeadFollPoly(const unsigned int i,
																					  const unsigned int j,
																					  const unsigned int k) const {
  const auto LeaderStart = this->EPECObject->TheNashGame->getPrimalLoc(i);
  const auto FollPoly    = this->PolyLCP.at(i).get()->convPolyPosition(k, true);
  return LeaderStart + FollPoly + j;
}


/**
 * @brief Get the position of the @p k -th leader variable of the @p i  leader, in the @p j-th
 * feasible polyhedron.
 * @param i The leader index
 * @param j The polyhedron index
 * @param k The leader variable index
 * @return The position for the queried attribute
 */
unsigned int Algorithms::EPEC::PolyBase::getPositionLeadLeadPoly(const unsigned int i,
																					  const unsigned int j,
																					  const unsigned int k) const {
  const auto LeaderStart = this->EPECObject->TheNashGame->getPrimalLoc(i);
  const auto FollPoly    = this->PolyLCP.at(i).get()->convPolyPosition(k, true);
  return LeaderStart + FollPoly + this->PolyLCP.at(i)->getLStart() + j;
}


/**
 * @brief Get the number of polyhedra used in the approximation for the @p i leader
 * @param i The leader index
 * @return The queried  number of polyhedra
 */
unsigned long int Algorithms::EPEC::PolyBase::getNumPolyLead(const unsigned int i) const {
  return this->PolyLCP.at(i).get()->convNumPoly(true);
}


/**
 * @brief Get the position of the probability associated with the @p k -th polyhedron
 * (@p k -th pure strategy) of the @p i -th leader. However, if the leader has an
 * inner approximation with exactly 1 polyhedron, it returns 0;
 * @param i The leader index
 * @param k The polyhedron index
 * @return The probability position associated with the queried values
 */
unsigned int Algorithms::EPEC::PolyBase::getPositionProbab(const unsigned int i,
																			  const unsigned int k) const {
  const auto PolyProbab = dynamic_cast<MathOpt::PolyLCP *>(this->EPECObject->PlayersLCP.at(i).get())
										->convPolyWeight(k, true);
  if (PolyProbab == 0)
	 return 0;
  const auto LeaderStart = this->EPECObject->TheNashGame->getPrimalLoc(i);
  return LeaderStart + PolyProbab;
}

/**
 * @brief Checks if the current equilibrium strategy in Game::EPEC is a pure strategy
 * @param tol A numerical tolerance
 * @return True if it is a pure equilibrium strategy
 */
bool Algorithms::EPEC::PolyBase::isPureStrategy(const double tol) const {
  for (unsigned int i = 0; i < this->EPECObject->getNumPlayers(); ++i) {
	 if (!isPureStrategy(i, tol))
		return false;
  }
  return true;
}


/**
 * @brief Checks whether the current strategy for the @p i player is a pure strategy
 * @param i  The index of the player
 * @param tol A numerical tolerance
 * @return True if it is a pure equilibrium strategy
 */
bool Algorithms::EPEC::PolyBase::isPureStrategy(const unsigned int i, const double tol) const {
  const unsigned int nPoly = this->getNumPolyLead(i);
  for (unsigned int j = 0; j < nPoly; j++) {
	 const double probab = this->getValProbab(i, j);
	 if (probab > 1 - tol) // Current Strategy is a pure strategy!
		return true;
  }
  return false;
}

/**
 * @brief Returns the indices of polyhedra feasible for the leader, from which strategies are played
 * with probability greater than the tolerance
 * @param i  The index of the player
 * @param tol A numerical tolerance
 * @return The indices of polyhedra with active probabilities
 */
std::vector<unsigned int> Algorithms::EPEC::PolyBase::mixedStrategyPoly(const unsigned int i,
																								const double tol) const {
  std::vector<unsigned int> polys{};
  const unsigned int        nPoly = this->getNumPolyLead(i);
  for (unsigned int j = 0; j < nPoly; j++) {
	 const double probab = this->getValProbab(i, j);
	 if (probab > tol)
		polys.push_back(j);
  }
  std::cout << "\n";
  return polys;
}


/**
 * @brief The probability associated with the @p k -th polyhedron of the @p i -th leader
 * @param i The index of the player
 * @param k The index of the polyhedron
 * @return The queried attribute
 */
double Algorithms::EPEC::PolyBase::getValProbab(const unsigned int i, const unsigned int k) const {
  const unsigned int varname{this->getPositionProbab(i, k)};
  if (varname == 0)
	 return 1;
  return this->EPECObject->SolutionX.at(varname);
}



/**
 * @brief For the @p i -th leader, gets the @p k -th pure strategy at position @p j
 * @param i The leader index
 * @param j The position index
 * @param k The pure strategy index
 * @param tol A numerical tolerance
 * @return The queried attribute
 */
double Algorithms::EPEC::PolyBase::getValLeadFollPoly(const unsigned int i,
																		const unsigned int j,
																		const unsigned int k,
																		const double       tol) const {
  if (!this->EPECObject->LCPModel)
	 throw ZEROException(ZEROErrorCode::Assertion, "LCPModel not made nor solved");
  const double probab = this->getValProbab(i, k);
  if (probab > 1 - tol)
	 return this->EPECObject->getValLeadFoll(i, j);
  else
	 return this->EPECObject->SolutionX.at(this->getPositionLeadFollPoly(i, j, k)) / probab;
}



/**
 * @brief For the @p i -th leader, gets the @p k -th pure strategy at leader position @p j
 * @param i The leader index
 * @param j The position index
 * @param k The pure strategy index
 * @param tol A numerical tolerance
 * @return The queried attribute
 */
double Algorithms::EPEC::PolyBase::getValLeadLeadPoly(const unsigned int i,
																		const unsigned int j,
																		const unsigned int k,
																		const double       tol) const {
  if (!this->EPECObject->LCPModel)
	 throw ZEROException(ZEROErrorCode::Assertion, "LCPModel not made nor solved");
  const double probab = this->getValProbab(i, k);
  if (probab > 1 - tol)
	 return this->EPECObject->getValLeadLead(i, j);
  else
	 return this->EPECObject->SolutionX.at(this->getPositionLeadLeadPoly(i, j, k)) / probab;
}


/**
 * @brief Creates an LCP for the inner-full approximation schemes of MathOpt::PolyLCP that
 * explicitly searches for pure-strategy equilibria. The original LCP is moved to
 * Game::EPEC::LCPModelBase
 */
void Algorithms::EPEC::PolyBase::makeThePureLCP() {
  try {
	 LOG_S(1) << "Game::EPEC::makeThePureLCP: editing the LCP model.";
	 this->EPECObject->LCPModelBase =
		  std::unique_ptr<GRBModel>(new GRBModel(*this->EPECObject->LCPModel));
	 const unsigned int nPolyLead = [this]() {
		unsigned int ell = 0;
		for (unsigned int i = 0; i < this->EPECObject->getNumPlayers(); ++i)
		  ell += (this->getNumPolyLead(i));
		return ell;
	 }();

	 // Add a binary variable for each polyhedron of each leader
	 GRBVar       pure_bin[nPolyLead];
	 GRBLinExpr   objectiveTerm{0};
	 unsigned int count{0}, i, j;
	 for (i = 0; i < this->EPECObject->getNumPlayers(); i++) {
		for (j = 0; j < this->getNumPolyLead(i); ++j) {
		  pure_bin[count] = this->EPECObject->LCPModel->addVar(
				0, 1, 0, GRB_BINARY, "pureBin_" + std::to_string(i) + "_" + std::to_string(j));
		  this->EPECObject->LCPModel->addGenConstrIndicator(
				pure_bin[count],
				1,
				this->EPECObject->LCPModel->getVarByName("x_" +
																	  std::to_string(this->getPositionProbab(i, j))),
				GRB_EQUAL,
				0,
				"Indicator_PNE_" + std::to_string(count));
		  objectiveTerm += pure_bin[count];
		  count++;
		}
	 }
	 this->EPECObject->LCPModel->setObjective(objectiveTerm, GRB_MAXIMIZE);
  } catch (GRBException &e) {
	 throw ZEROException(ZEROErrorCode::SolverError,
								std::to_string(e.getErrorCode()) + e.getMessage());
  }
}