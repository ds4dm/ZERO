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


#include "games/algorithms/EPEC/epec_fullenumeration.h"
#include "mathopt/lcp/poly_lcp.h"
#include <boost/log/trivial.hpp>

/**
 * @brief Solves the Game::EPEC by full enumeration.
 */
void Algorithms::EPEC::FullEnumeration::solve() {
  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i)
	 this->PolyLCP.at(i)->exactFullEnumeration(true);
  this->EPECObject->makePlayersQPs();
  LOG_S(1) << "Algorithms::EPEC::FullEnumeration::solve: "
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