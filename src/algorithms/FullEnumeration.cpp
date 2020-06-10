#include "algorithms/fullenumeration.h"
#include "lcp/polylcp.h"
#include <boost/log/trivial.hpp>

void Algorithms::FullEnumeration::solve() {
  /** @brief Solve the referenced EPEC instance with the full enumeration
   * @p excludelist contains the set of excluded polyhedra combinations.
   */
  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i)
    this->PolyLCP.at(i)->enumerateAll(true);
  this->EPECObject->makePlayersQPs();
  BOOST_LOG_TRIVIAL(trace)
      << "Algorithms::FullEnumeration::solve: Starting FullEnumeration search";
  this->EPECObject->computeNashEq(
      this->EPECObject->Stats.AlgorithmParam.PureNashEquilibrium,
      this->EPECObject->Stats.AlgorithmParam.TimeLimit);
  if (this->isSolved()) {
    this->EPECObject->Stats.Status = Game::EPECsolveStatus::NashEqFound;
    if (this->isPureStrategy())
      this->EPECObject->Stats.PureNashEquilibrium = true;
  }
  // Post Solving
  this->postSolving();
}