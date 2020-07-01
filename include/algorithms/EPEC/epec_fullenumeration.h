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

#pragma once
#include "epec_algorithms.h"
#include "epec_polybase.h"

namespace Algorithms {
  namespace EPEC {
	 ///@brief This class is responsible for the fully enumerative Algorithm
	 class FullEnumeration : public PolyBase {
	 public:
		FullEnumeration(GRBEnv *env, Game::EPEC *EPECObject) : PolyBase(env, EPECObject){};
		void solve();
	 };
  } // namespace EPEC
} // namespace Algorithms