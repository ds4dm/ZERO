
.. _program_listing_file_include_games_algorithms_EPEC_epec_innerapproximation.h:

Program Listing for File epec_innerapproximation.h
==================================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_games_algorithms_EPEC_epec_innerapproximation.h>` (``include/games/algorithms/EPEC/epec_innerapproximation.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

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
