#pragma once
#include "algorithms/polybase.h"
#include <armadillo>
#include <boost/log/trivial.hpp>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Algorithms {

///@brief This class is responsible for the Combinatorial pure-nash Equilibrium
class CombinatorialPNE : public PolyBase {
public:
  CombinatorialPNE(GRBEnv *env, Game::EPEC *EPECObject, bool poly = true)
      : PolyBase(env, EPECObject){};
  ;
  void solve() {
    this->solveWithExcluded(std::vector<std::set<unsigned long int>>{});
  }
  void solveWithExcluded(
      const std::vector<std::set<unsigned long int>> &excludeList = {});

private:
  // Making the method private
  void combPNE(std::vector<long int> combination,
               const std::vector<std::set<unsigned long int>> &excludeList);
};
} // namespace Algorithms