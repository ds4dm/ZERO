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

#include "abstractseparator.h"
#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Algorithms::IPG::IPG_Separators {
  class Knapsack : public AbstractSeparator {
  public:
	 void separationMain(const arma::vec &x);
    void minimalCover(const unsigned int player, const arma::vec &xOfi, std::vector< int> &C);
  };
} // namespace Game::IPG_Separators