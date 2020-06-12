
.. _program_listing_file_include_algorithms_innerapproximation.h:

Program Listing for File innerapproximation.h
=============================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_algorithms_innerapproximation.h>` (``include/algorithms/innerapproximation.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include "algorithms/algorithms.h"
   #include "algorithms/polybase.h"
   #include <armadillo>
   #include <gurobi_c++.h>
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   
   namespace Algorithms {
   
   class InnerApproximation : public PolyBase {
   
   public:
     InnerApproximation(GRBEnv *env, Game::EPEC *EPECObject)
         : PolyBase(env, EPECObject){};
     void solve();
   
   private:
     void start();
     bool addRandomPoly2All(unsigned int aggressiveLevel = 1,
                            bool stopOnSingleInfeasibility = false);
     bool getAllDeviations(std::vector<arma::vec> &deviations,
                           const arma::vec &guessSol,
                           const std::vector<arma::vec> &prevDev = {}) const;
     unsigned int addDeviatedPolyhedron(const std::vector<arma::vec> &deviations,
                                        bool &infeasCheck) const;
   };
   } // namespace Algorithms
