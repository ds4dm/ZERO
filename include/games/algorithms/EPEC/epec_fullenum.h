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
#include "epec_polybase.h"
#include "zero.h"

namespace Algorithms::EPEC {
	 /**
	  * @brief This class manages the full enumeration algorithm for Game::EPEC objects
	  */
	 class FullEnumeration : public PolyBase {
	 public:
		FullEnumeration(GRBEnv *env, Game::EPEC *EPECObject) : PolyBase(env, EPECObject){};
		void solve();
	 };
  } // namespace Algorithms