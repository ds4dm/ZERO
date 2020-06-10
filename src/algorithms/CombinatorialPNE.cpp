#include "algorithms/combinatorialpne.h"
#include <boost/log/trivial.hpp>
#include <chrono>
#include <gurobi_c++.h>
#include <set>

void Algorithms::CombinatorialPNE::solveWithExcluded(
    const std::vector<std::set<unsigned long int>> &excludeList) {
  /** @brief Solve the referenced EPEC instance with the Combinatorial
   *pure-equilibrium Algorithm
   * @p exclude-list contains the set of excluded polyhedra combinations.
   */
  if (this->EPECObject->Stats.AlgorithmParam.TimeLimit > 0) {
    // Checking the function hasn't been called from InnerApproximation
    if (this->EPECObject->Stats.NumIterations <= 0) {
      this->EPECObject->InitTime = std::chrono::high_resolution_clock::now();
    }
  }
  std::vector<long int> start = {};
  for (int j = 0; j < this->EPECObject->NumPlayers; ++j)
    start.push_back(-1);
  this->combPNE(start, excludeList);
  if (this->EPECObject->Stats.Status == Game::EPECsolveStatus::Uninitialized)
    this->EPECObject->Stats.Status = Game::EPECsolveStatus::NashEqNotFound;
  this->postSolving();
}

void Algorithms::CombinatorialPNE::combPNE(
    std::vector<long int> combination,
    const std::vector<std::set<unsigned long int>> &excludeList) {
  /** @brief Starting from @p combination, the methods builds the recursion to
   * generate the subproblems associated with all the existing combinations of
   * polyhedra. Then, it solves each subproblem, and if a solution is found, it
   * terminates and  stores the solution  into the referenced EPEC object. @p
   * excludeList contains the excluded combinations of polyhedra.
   */
  if ((this->EPECObject->Stats.Status == Game::EPECsolveStatus::NashEqFound &&
       this->EPECObject->Stats.PureNashEquilibrium) ||
      this->EPECObject->Stats.Status == Game::EPECsolveStatus::TimeLimit)
    return;

  if (this->EPECObject->Stats.AlgorithmParam.TimeLimit > 0) {
    const std::chrono::duration<double> timeElapsed =
        std::chrono::high_resolution_clock::now() - this->EPECObject->InitTime;
    const double timeRemaining =
        this->EPECObject->Stats.AlgorithmParam.TimeLimit - timeElapsed.count();
    if (timeRemaining <= 0) {
      this->EPECObject->Stats.Status = Game::EPECsolveStatus::TimeLimit;
      return;
    }
  }

  std::vector<long int> childCombination(combination);
  bool found{false};
  unsigned int i{0};
  for (i = 0; i < this->EPECObject->NumPlayers; i++) {
    if (childCombination.at(i) == -1) {
      found = true;
      break;
    }
  }
  if (found) {
    for (unsigned int j = 0; j < this->PolyLCP.at(i)->getNumTheoreticalPoly();
         ++j) {
      if (this->PolyLCP.at(i)->checkPolyFeas(j)) {
        childCombination.at(i) = j;
        this->combPNE(childCombination, excludeList);
      }
    }
  } else {
    // Combination is filled and ready!
    // Check that this combination is not in the excluded list
    BOOST_LOG_TRIVIAL(trace) << "Algorithms::CombinatorialPNE::combPNE: "
                                "considering a FULL combination";
    bool excluded = false;
    if (!excludeList.empty()) {
      excluded = true;
      for (unsigned int j = 0; j < this->EPECObject->NumPlayers; ++j) {
        if (excludeList.at(j).find(static_cast<const unsigned long &>(
                childCombination.at(j))) == excludeList.at(j).end()) {
          excluded = false;
        }
      }
    }

    if (!excluded) {
      BOOST_LOG_TRIVIAL(trace)
          << "Algorithms::CombinatorialPNE::combPNE: considering a "
             "FEASIBLE combination of polyhedra.";
      for (unsigned long j = 0; j < this->EPECObject->NumPlayers; ++j) {
        this->PolyLCP.at(j)->clearPolyhedra();
        this->PolyLCP.at(j)->addThePoly(
            static_cast<const unsigned long &>(childCombination.at(j)));
      }
      this->EPECObject->makePlayersQPs();
      bool res = 0;
      if (this->EPECObject->Stats.AlgorithmParam.TimeLimit > 0) {
        const std::chrono::duration<double> timeElapsed =
            std::chrono::high_resolution_clock::now() -
            this->EPECObject->InitTime;
        const double timeRemaining =
            this->EPECObject->Stats.AlgorithmParam.TimeLimit -
            timeElapsed.count();
        res = this->EPECObject->computeNashEq(false, timeRemaining, true);
      } else
        res = this->EPECObject->computeNashEq(false, -1.0, true);

      if (res) {
        if (this->isSolved()) {
          // Check that the equilibrium is a pure strategy
          if ((this->isPureStrategy())) {
            BOOST_LOG_TRIVIAL(info) << "Algorithms::CombinatorialPNE::combPNE: "
                                       "found a pure strategy.";
            this->EPECObject->Stats.Status = Game::EPECsolveStatus::NashEqFound;
            this->EPECObject->Stats.PureNashEquilibrium = true;
            return;
          }
        }
      }
    } else {
      BOOST_LOG_TRIVIAL(trace)
          << "Algorithms::CombinatorialPNE::combPNE: configuration pruned.";
      return;
    }
  }
}