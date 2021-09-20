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


#include "games/algorithms/EPEC/epec_combPNE.h"

#include <utility>

/**
 * @brief Solves the Game::EPEC instance with the algorithm by excluding some combinations of
 * polyhedra that may have been already tested.
 * @param excludeList A set containing the combinations that should be excluded
 */
void Algorithms::EPEC::CombinatorialPNE::solveWithExcluded(
	 const std::vector<std::set<unsigned long int>> &excludeList) {
  if (this->EPECObject->Stats.AlgorithmData.TimeLimit.get() > 0) {
	 // Checking the function hasn't been called from InnerApproximation
	 if (this->EPECObject->Stats.NumIterations.get() <= 0) {
		this->EPECObject->InitTime = std::chrono::high_resolution_clock::now();
	 }
  }
  std::vector<long int> start(this->EPECObject->NumPlayers, -1);
  this->combPNE(start, excludeList);
  if (this->EPECObject->Stats.Status.get() == ZEROStatus::Uninitialized)
	 this->EPECObject->Stats.Status.set(ZEROStatus::NashEqNotFound);
  this->after();
}

/**
 * @brief This method initializes the algorithm recursion with @p combination. Each element is the
 * index of a polyhedron for the corresponding player. If the index is -1, then the recursion will
 * generate children for any polyhedron of the given player. Otherwise, if there exist a positive
 * value @f$v@f$ in a location @f$l@f$, player @f$l@f$ will only play strategies in the polyhedron
 * @f$v@f$.
 * @param combination A set of either -1 or positive numbers corresponding to the polyhedron of each
 * player. -1 will recurse
 * @param excludeList A set of combinations of polyhedra that should be excluded from the search.
 */
void Algorithms::EPEC::CombinatorialPNE::combPNE(
	 std::vector<long int>                           combination,
	 const std::vector<std::set<unsigned long int>> &excludeList) {
  if ((this->EPECObject->Stats.Status.get() == ZEROStatus::NashEqFound &&
		 this->EPECObject->Stats.PureNashEquilibrium.get()) ||
		this->EPECObject->Stats.Status.get() == ZEROStatus::TimeLimit)
	 return;

  if (this->EPECObject->Stats.AlgorithmData.TimeLimit.get() > 0) {
	 const std::chrono::duration<double> timeElapsed =
		  std::chrono::high_resolution_clock::now() - this->EPECObject->InitTime;
	 const double timeRemaining =
		  this->EPECObject->Stats.AlgorithmData.TimeLimit.get() - timeElapsed.count();
	 if (timeRemaining <= 0) {
		this->EPECObject->Stats.Status.set(ZEROStatus::TimeLimit);
		return;
	 }
  }

  std::vector<long int> childCombination(std::move(combination));
  bool                  found{false};
  unsigned int          i;
  for (i = 0; i < this->EPECObject->NumPlayers; i++) {
	 if (childCombination.at(i) == -1) {
		found = true;
		break;
	 }
  }
  if (found) {
	 for (unsigned int j = 0; j < this->PolyLCP.at(i)->getNumTheoreticalPoly(); ++j) {
		if (this->PolyLCP.at(i)->checkPolyFeas(j, true)) {
		  childCombination.at(i) = j;
		  this->combPNE(childCombination, excludeList);
		}
	 }
  } else {
	 // Combination is filled and ready!
	 // Check that this combination is not in the excluded list
	 LOG_S(1) << "Algorithms::EPEC::CombinatorialPNE::combPNE: "
					 "considering a FULL combination";
	 bool excluded = false;
	 if (!excludeList.empty()) {
		excluded = true;
		for (unsigned int j = 0; j < this->EPECObject->NumPlayers; ++j) {
		  if (excludeList.at(j).find(static_cast<const unsigned long &>(childCombination.at(j))) ==
				excludeList.at(j).end()) {
			 excluded = false;
		  }
		}
	 }

	 if (!excluded) {
		LOG_S(1) << "Algorithms::EPEC::CombinatorialPNE::combPNE: considering a "
						"FEASIBLE combination of polyhedra.";
		for (unsigned long j = 0; j < this->EPECObject->NumPlayers; ++j) {
		  this->PolyLCP.at(j)->clearPolyhedra(true);
		  this->PolyLCP.at(j)->addThePoly(static_cast<const unsigned long &>(childCombination.at(j)),
													 true);
		}
		this->EPECObject->makePlayersQPs();
		bool res = 0;
		if (this->EPECObject->Stats.AlgorithmData.TimeLimit.get() > 0) {
		  const std::chrono::duration<double> timeElapsed =
				std::chrono::high_resolution_clock::now() - this->EPECObject->InitTime;
		  const double timeRemaining =
				this->EPECObject->Stats.AlgorithmData.TimeLimit.get() - timeElapsed.count();
		  res = this->EPECObject->computeNashEq(false, timeRemaining, true, false, false);
		} else
		  res = this->EPECObject->computeNashEq(false, -1.0, true, false, false);

		if (res) {
		  if (this->isSolved()) {
			 // Check that the equilibrium is a pure strategy
			 if ((this->isPureStrategy())) {
				LOG_S(INFO) << "Algorithms::EPEC::CombinatorialPNE::combPNE: "
									"found a pure strategy.";
				this->EPECObject->Stats.Status.set(ZEROStatus::NashEqFound);
				this->EPECObject->Stats.PureNashEquilibrium = true;
				return;
			 }
		  }
		} else if (this->EPECObject->Stats.Status.get() == ZEROStatus::Numerical)
		  return;
	 } else {
		LOG_S(1) << "Algorithms::EPEC::CombinatorialPNE::combPNE:"
						" configuration pruned.";
		return;
	 }
  }
}