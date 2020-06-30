#include "algorithms/IPG/ipg_oracle.h"

#include <boost/log/trivial.hpp>
#include <chrono>
#include <gurobi_c++.h>
#include <set>
#include <string>

bool Algorithms::IPG::Oracle::isSolved(double tol) const { return this->Feasible; }

void Algorithms::IPG::Oracle::addValueCut(unsigned int player,
														arma::vec    xOfIBestResponse,
														arma::vec    xMinusI) {

  double cutRHS =
		this->IPG->PlayersIP.at(player)->computeObjective(xOfIBestResponse, xMinusI, false);
  arma::vec LHS =
		this->IPG->PlayersIP.at(player)->getc() + this->IPG->PlayersIP.at(player)->getC() * xMinusI;
  arma::sp_mat cutLHS =
		Utils::resizePatch(arma::sp_mat{LHS}.t(), 1, this->IPG->PlayersIP.at(player)->getC().n_cols);
  BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::addValueCut: "
									  "adding cut for Player "
								  << player;
  //@todo add cut
  // this->outerLCP.at(player)->addCustomCuts(-cutLHS, arma::vec{-cutRHS});
}

bool Algorithms::IPG::Oracle::isFeasible(bool &addedCuts, double tol) {

  /**
	* Check whether the current outer approximation equilibrium is feasible and
	* solves the outer approximation. Otherwise, add cuts, or generate useful
	* points for next iterations.
	*/

  // First, we have a NE from Games::computeNashEq
  if (!this->IPG->NashEquilibrium)
	 return false;

  // Then, the feasibility is implied also by the deviations
  bool      result = {true};
  arma::vec bestResponse;
  arma::vec currentPayoffs =
		this->EPECObject->TheNashGame->computeQPObjectiveValues(this->EPECObject->SolutionX, true);
  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i) {
	 this->Trees.at(i)->resetFeasibility();
	 double val = this->EPECObject->respondSol(bestResponse, i, this->EPECObject->SolutionX);
	 if (val == GRB_INFINITY) {
		BOOST_LOG_TRIVIAL(trace) << "Algorithms::EPEC::OuterApproximation:: Unbounded deviation for "
										 << i;
		addedCuts = false;
		return false;
	 }
	 // minimization standard
	 if (std::abs(currentPayoffs.at(i) - val) > tol) {
		// Discrepancy between payoffs! Need to investigate.
		if ((currentPayoffs.at(i) - val) > tol) {
		  // It means the current payoff is more than then optimal response. Then
		  // this is not a best response. Theoretically, this cannot happen from
		  // an outer approximation. This if case is a warning case then
		  //@todo can this happen?

		  BOOST_LOG_TRIVIAL(warning) << "Algorithms::EPEC::OuterApproximation::"
												  "isFeasible: No best response for Player "
											  << i;
		  BOOST_LOG_TRIVIAL(trace) << "Algorithms::EPEC::OuterApproximation:: "
											<< currentPayoffs.at(i) << " vs " << val;
		  result = false;
		  // throw;
		  // throw;
		} else if ((currentPayoffs.at(i) - val) < tol) {
		  // It means the current payoff is less than the optimal response. The
		  // approximation is not good, and this point is infeasible. Then, we can
		  // generate a value-cut
		  arma::vec xMinusI;
		  this->EPECObject->getXMinusI(this->EPECObject->SolutionX, i, xMinusI);
		  this->addValueCut(i, bestResponse, xMinusI);
		  BOOST_LOG_TRIVIAL(info) << "Algorithms::EPEC::OuterApproximation::isFeasible: "
											  "Value cut at for Player "
										  << i;
		  result = false;
		}
	 } else {
		// Here we have a best response whose payoff coincides with the one of the
		// equilibrium. The strategy might not be feasible, though.
		arma::vec xOfI;
		this->EPECObject->getXofI(this->EPECObject->SolutionX, i, xOfI, false);

		// Check if we need to add the point to the vertex storage.
		arma::vec vertex = bestResponse.subvec(0, xOfI.size() - 1);
		vertex.print("Best Response");
		if (!Utils::containsRow(*this->Trees.at(i)->getV(), vertex, this->Tolerance)) {
		  this->Trees.at(i)->addVertex(vertex);
		  BOOST_LOG_TRIVIAL(info) << "Algorithms::EPEC::OuterApproximation::isFeasible: "
											  "Adding vertex as of best response for Player "
										  << i << " (Best Response)";
		} else {
		  BOOST_LOG_TRIVIAL(info) << "Algorithms::EPEC::OuterApproximation::isFeasible: "
											  "Already known best response for Player "
										  << i << " (Best Response)";
		}

		// Check if best response coincides with the strategy in the equilibrium
		bool same = true;
		for (unsigned int k = 0; k < xOfI.size(); ++k) {
		  if (std::abs(xOfI.at(k) - bestResponse.at(k)) > tol) {
			 same = false;
			 break;
		  }
		}
		if (!same) {
		  // Then, if the answers do not coincide, we need to refine the
		  // approximation or determine if this strategy is anyhow feasible.
		  // We search for a convex combination of best responses so that we can
		  // certify the answer is inside the convex-hull (or not).

		  int budget = 15;
		  if (!this->separationOracle(xOfI, this->EPECObject->SolutionX, i, budget, addedCuts)) {
			 BOOST_LOG_TRIVIAL(trace) << "Algorithms::EPEC::OuterApproximation::isFeasible: "
												  "Oracle gave a negative answer for Player "
											  << i;
			 result = false;
		  }

		} else {
		  this->Trees.at(i)->setFeasible();
		  this->Trees.at(i)->setPure();
		  BOOST_LOG_TRIVIAL(info) << "Algorithms::EPEC::OuterApproximation::isFeasible: "
											  "Feasible strategy for Player "
										  << i << " (Best Response)";
		}
	 }
  }
  return result;
}

bool Algorithms::IPG::Oracle::addConstraintsToPool(
	 const arma::sp_mat A,      ///< The LHS of the constraints to add
	 const arma::vec    b,      ///< The RHS of the constraints to add
	 const unsigned int player, ///< The index of the player
	 bool               check   ///< A boolean controling whether the method check for the
										 ///< uniqueness of the constraints
) {
  /**
	* @brief This method adds to the members CutPool_A and CutPool_b the content
	* of @p A and @p b. If @p check is set to true, then the methods checks that
	* the content of A and b is unique. In case a duplicate is detected, the copy
	* is not added to the class' members.
	* @return true in at least one constraint has been added
	*/

  if (this->CutPool_A.size() < player || this->CutPool_b.size() < player)
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Mismatch between CutPool size and player number");

  if (this->CutPool_A.at(player).n_cols != A.n_cols)
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Mismatch between the CutPool of the player and the input matrix");
  if (b.size() != A.n_rows)
	 throw ZEROException(ZEROErrorCode::InvalidData, "Mismatch between the rows of the inputs");

  bool ret{false};
  if (!check)
	 bool ret = true;
  for (unsigned int i = 0; i < A.n_rows; i++) {
	 arma::sp_mat Ai = A.submat(i, 0, i, A.n_cols);
	 if (!Utils::containsConstraint(this->CutPool_A.at(player),
											  this->CutPool_b.at(player),
											  Ai,
											  b.at(player),
											  this->Tolerance)) {
		// This constraint does not exist
		ret                        = true;
		this->CutPool_A.at(player) = arma::join_cols(this->CutPool_A.at(player), Ai);
		this->CutPool_b.at(player) = arma::join_cols(this->CutPool_b.at(player), arma::vec{b.at(i)});
	 }
  }

  return ret;
}
