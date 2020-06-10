#pragma once
#include "algorithms/algorithms.h"
#include "algorithms/polybase.h"

namespace Algorithms {

///@brief This class is responsible for the fully enumerative Algorithm
class FullEnumeration : public PolyBase {
public:
  FullEnumeration(GRBEnv *env, Game::EPEC *EPECObject)
      : PolyBase(env, EPECObject){};
  void solve();
};
} // namespace Algorithms