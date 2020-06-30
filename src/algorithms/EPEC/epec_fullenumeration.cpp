#include "algorithms/EPEC/epec_fullenumeration.h"
#include "lcp/poly_lcp.h"
#include <boost/log/trivial.hpp>

void Algorithms::EPEC::FullEnumeration::solve() {
  /** @brief Solve the referenced EPEC instance with the full enumeration
	* @p excludelist contains the set of excluded polyhedra combinations.
	*/
  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i)
	 this->PolyLCP.at(i)->enumerateAll(true);
  this->EPECObject->makePlayersQPs();
  BOOST_LOG_TRIVIAL(trace) << "Algorithms::EPEC::FullEnumeration::solve: "
										"Starting FullEnumeration search";
  this->EPECObject->computeNashEq(this->EPECObject->Stats.AlgorithmData.PureNashEquilibrium.get(),
											 this->EPECObject->Stats.AlgorithmData.TimeLimit.get());
  if (this->isSolved()) {
	 this->EPECObject->Stats.Status.set(ZEROStatus::NashEqFound);
	 if (this->isPureStrategy())
		this->EPECObject->Stats.PureNashEquilibrium = true;
  }
  // Post Solving
  this->postSolving();
}