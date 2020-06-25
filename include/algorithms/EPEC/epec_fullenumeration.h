#pragma once
#include "epec_algorithms.h"
#include "epec_polybase.h"

namespace Algorithms {
namespace EPEC {
///@brief This class is responsible for the fully enumerative Algorithm
class FullEnumeration : public PolyBase {
public:
  FullEnumeration(GRBEnv *env, Game::EPEC *EPECObject)
      : PolyBase(env, EPECObject){};
  void solve();
};
} // namespace EPEC
} // namespace Algorithms