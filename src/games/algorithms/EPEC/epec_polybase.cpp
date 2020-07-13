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


#include "games/algorithms/EPEC/epec_polybase.h"
#include "zero.h"
#include <boost/log/trivial.hpp>

bool Algorithms::EPEC::PolyBase::isSolved(unsigned int *countryNumber,
														arma::vec *   profitableDeviation,
														double        tol) const
/**
 * @brief Checks if Game::EPEC is solved, otherwise it returns a proof.
 * @details
 * Analogous to Game::NashGame::isSolved but checks if the given Game::EPEC is
 * solved. If it is solved, then retruns true. If not, it returns the country
 * which has a profitable deviation in @p countryNumber and the profitable
 * deviation in @p profitableDeviation. @p Tolerance is the tolerance for the
 * check. If the <i> improved objective </i> after the deviation is less than @p
 * Tolerance, then it is not considered as a profitable deviation.
 *
 * Thus we check if the given point is an @f$\epsilon@f$-equilibrium. Value of
 * @f$\epsilon @f$ can be chosen sufficiently close to 0.
 *
 * @warning Setting @p Tolerance = 0 might even reject a real solution as not
 * solved. This is due to Numerical issues arising from the LCP solver (Gurobi).
 */
{
  if (!this->EPECObject->TheNashGame)
	 return false;
  if (!this->EPECObject->NashEquilibrium)
	 return false;
  if (tol < 0)
	 tol = this->EPECObject->Stats.AlgorithmData.DeviationTolerance.get();
  this->EPECObject->TheNashGame->isSolved(
		this->EPECObject->SolutionX, *countryNumber, *profitableDeviation, tol);
  arma::vec objvals =
		this->EPECObject->TheNashGame->computeQPObjectiveValues(this->EPECObject->SolutionX, true);
  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i) {
	 double val = this->EPECObject->respondSol(*profitableDeviation, i, this->EPECObject->SolutionX);
	 if (val == GRB_INFINITY)
		return false;
	 if (std::abs(val - objvals.at(i)) > tol) {
		*countryNumber = i;
		BOOST_LOG_TRIVIAL(trace) << "Algorithms::EPEC::PolyBase::isSolved: deviation for player " << i
										 << " -- of " << std::abs(val - objvals.at(i));
		return false;
	 }
  }
  return true;
}

bool Algorithms::EPEC::PolyBase::isSolved(double tol) const {
  unsigned int countryNumber;
  arma::vec    ProfDevn;
  bool         ret = this->isSolved(&countryNumber, &ProfDevn, tol);
  return ret;
}

unsigned int Algorithms::EPEC::PolyBase::getPositionLeadFollPoly(const unsigned int i,
																					  const unsigned int j,
																					  const unsigned int k) const {
  /**
	* Get the position of the k-th follower variable of the i-th leader, in the
	* j-th feasible polyhedron.
	*
	* Indeed it should hold that @f$ j < @f$
	* Algorithms::EPEC::PolyBase::getNumPolyLead(i)
	*/
  const auto LeaderStart = this->EPECObject->TheNashGame->getPrimalLoc(i);
  const auto FollPoly =
		static_cast<MathOpt::PolyLCP *>(this->PolyLCP.at(i).get())->convPolyPosition(k);
  return LeaderStart + FollPoly + j;
}

unsigned int Algorithms::EPEC::PolyBase::getPositionLeadLeadPoly(const unsigned int i,
																					  const unsigned int j,
																					  const unsigned int k) const {
  /**
	* Get the position of the k-th leader variable of the i-th leader, in the
	* j-th feasible polyhedron.
	*
	* Indeed it should hold that @f$ j < @f$
	* Algorithms::EPEC::PolyBase::getNumPolyLead(i)
	*/
  const auto LeaderStart = this->EPECObject->TheNashGame->getPrimalLoc(i);
  const auto FollPoly =
		static_cast<MathOpt::PolyLCP *>(this->PolyLCP.at(i).get())->convPolyPosition(k);
  return LeaderStart + FollPoly + this->PolyLCP.at(i)->getLStart() + j;
}

unsigned int Algorithms::EPEC::PolyBase::getNumPolyLead(const unsigned int i) const {
  /**
	* Get the number of polyhedra used in the inner approximation of the
	* feasible region of the i-th leader.*
	*/
  return static_cast<MathOpt::PolyLCP *>(this->PolyLCP.at(i).get())->convNumPoly();
}

unsigned int Algorithms::EPEC::PolyBase::getPositionProbab(const unsigned int i,
																			  const unsigned int k) const {
  /**
	* Get the position of the probability associated with the k-th polyhedron
	* (k-th pure strategy) of the i-th leader. However, if the leader has an
	* inner approximation with exactly 1 polyhedron, it returns 0;
	*/
  const auto PolyProbab =
		static_cast<MathOpt::PolyLCP *>(this->EPECObject->PlayersLCP.at(i).get())->convPolyWeight(k);
  if (PolyProbab == 0)
	 return 0;
  const auto LeaderStart = this->EPECObject->TheNashGame->getPrimalLoc(i);
  return LeaderStart + PolyProbab;
}

bool Algorithms::EPEC::PolyBase::isPureStrategy(const double tol) const {
  /**
	* Checks if the returned strategy leader is a pure strategy for the leader
	* i. The strategy is considered a pure strategy, if it is played with a
	* probability greater than 1 - Tolerance;
	*/
  for (unsigned int i = 0; i < this->EPECObject->getNumPlayers(); ++i) {
	 if (!isPureStrategy(i, tol))
		return false;
  }
  return true;
}

bool Algorithms::EPEC::PolyBase::isPureStrategy(const unsigned int i, const double tol) const {
  /**
	* Checks if the returned strategy leader is a pure strategy for the leader
	* i. The strategy is considered a pure strategy, if it is played with a
	* probability greater than 1 - Tolerance;
	*/
  const unsigned int nPoly = this->getNumPolyLead(i);
  for (unsigned int j = 0; j < nPoly; j++) {
	 const double probab = this->getValProbab(i, j);
	 if (probab > 1 - tol) // Current Strategy is a pure strategy!
		return true;
  }
  return false;
}

std::vector<unsigned int> Algorithms::EPEC::PolyBase::mixedStrategyPoly(const unsigned int i,
																								const double tol) const
/**
 * Returns the indices of polyhedra feasible for the leader, from which
 * strategies are played with probability greater than Tolerance.
 */
{
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

double Algorithms::EPEC::PolyBase::getValProbab(const unsigned int i, const unsigned int k) const {
  /**
	* Get the probability associated with the k-th polyhedron
	* (k-th pure strategy) of the i-th leader.
	*/
  const unsigned int varname{this->getPositionProbab(i, k)};
  if (varname == 0)
	 return 1;
  return this->EPECObject->LCPModel->getVarByName("x_" + std::to_string(varname))
		.get(GRB_DoubleAttr_X);
}

double Algorithms::EPEC::PolyBase::getValLeadFollPoly(const unsigned int i,
																		const unsigned int j,
																		const unsigned int k,
																		const double       tol) const {
  /**
	* For the i-th leader, gets the k-th pure strategy for i-th leader at
	* position j
	*/
  if (!this->EPECObject->LCPModel)
	 throw ZEROException(ZEROErrorCode::Assertion, "LCPModel not made nor solved");
  const double probab = this->getValProbab(i, k);
  if (probab > 1 - tol)
	 return this->EPECObject->getValLeadFoll(i, j);
  else
	 return this->EPECObject->LCPModel
					->getVarByName("x_" + std::to_string(this->getPositionLeadFollPoly(i, j, k)))
					.get(GRB_DoubleAttr_X) /
			  probab;
}

double Algorithms::EPEC::PolyBase::getValLeadLeadPoly(const unsigned int i,
																		const unsigned int j,
																		const unsigned int k,
																		const double       tol) const {
  /**
	* For the i-th leader, gets the k-th pure strategy for i-th leader at
	* non-follower leader position j
	*/
  if (!this->EPECObject->LCPModel)
	 throw ZEROException(ZEROErrorCode::Assertion, "LCPModel not made nor solved");
  const double probab = this->getValProbab(i, k);
  if (probab > 1 - tol)
	 return this->EPECObject->getValLeadLead(i, j);
  else
	 return this->EPECObject->LCPModel
					->getVarByName("x_" + std::to_string(this->getPositionLeadLeadPoly(i, j, k)))
					.get(GRB_DoubleAttr_X) /
			  probab;
}

void Algorithms::EPEC::PolyBase::makeThePureLCP(bool indicators) {
  /**
	* Given that Game::EPEC::LCPModel is filled with the final LCP,
	* directs the search toward a pure nash EQ. If such an equilibrium does not
	* exist, then the model will return anyway a MNE. The original LCP is
	* stored in the field Game::EPEC::LCPModelBase. @p Indicators dictates
	* whether the resulting LCP should use indicator constraints instead of
	* general binaries. In general, there are advantages in using the binary
	* variables instead of such constraints, since there is no BigM involved in
	* the formulation.
	*/
  try {
	 BOOST_LOG_TRIVIAL(trace) << "Game::EPEC::makeThePureLCP: editing the LCP model.";
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
	 //@todo fix bug in else
	 for (i = 0; i < this->EPECObject->getNumPlayers(); i++) {
		for (j = 0; j < this->getNumPolyLead(i); ++j) {
		  pure_bin[count] = this->EPECObject->LCPModel->addVar(
				0, 1, 0, GRB_BINARY, "pureBin_" + std::to_string(i) + "_" + std::to_string(j));
		  if (indicators) {
			 this->EPECObject->LCPModel->addGenConstrIndicator(
				  pure_bin[count],
				  1,
				  this->EPECObject->LCPModel->getVarByName(
						"x_" + std::to_string(this->getPositionProbab(i, j))),
				  GRB_EQUAL,
				  0,
				  "Indicator_PNE_" + std::to_string(count));
		  } else {
			 this->EPECObject->LCPModel->addConstr(
				  this->EPECObject->LCPModel->getVarByName(
						"x_" + std::to_string(this->getPositionProbab(i, j))),
				  GRB_GREATER_EQUAL,
				  pure_bin[count]);
		  }
		  objectiveTerm += pure_bin[count];
		  count++;
		}
	 }
	 this->EPECObject->LCPModel->setObjective(objectiveTerm, GRB_MAXIMIZE);
	 if (indicators) {
		BOOST_LOG_TRIVIAL(trace) << "Algorithms::EPEC::PolyBase::makeThePureLCP: using "
											 "indicator constraints.";
	 } else {
		BOOST_LOG_TRIVIAL(trace) << "Algorithms::EPEC::PolyBase::makeThePureLCP: using "
											 "indicator constraints.";
	 }
  } catch (GRBException &e) {
	 throw ZEROException(ZEROErrorCode::SolverError,
								std::to_string(e.getErrorCode()) + e.getMessage());
  }
}