
.. _program_listing_file_include_algorithms_EPEC_epec_innerapproximation.h:

Program Listing for File epec_innerapproximation.h
==================================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_algorithms_EPEC_epec_innerapproximation.h>` (``include/algorithms/EPEC/epec_innerapproximation.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include "epec_algorithms.h"
   #include "epec_polybase.h"
   #include <armadillo>
   #include <gurobi_c++.h>
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   
   namespace Algorithms {
     namespace EPEC {
        class InnerApproximation : public PolyBase {
   
        public:
           InnerApproximation(GRBEnv *env, Game::EPEC *EPECObject) : PolyBase(env, EPECObject){};
           void solve();
   
        private:
           void         start();
           bool         addRandomPoly2All(unsigned int aggressiveLevel           = 1,
                                                    bool         stopOnSingleInfeasibility = false);
           bool         getAllDeviations(std::vector<arma::vec> &      deviations,
                                                   const arma::vec &             guessSol,
                                                   const std::vector<arma::vec> &prevDev = {}) const;
           unsigned int addDeviatedPolyhedron(const std::vector<arma::vec> &deviations,
                                                         bool &                        infeasCheck) const;
        };
     } // namespace EPEC
   } // namespace Algorithms
