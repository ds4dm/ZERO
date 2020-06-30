
.. _program_listing_file_include_algorithms_EPEC_epec_combinatorialpne.h:

Program Listing for File epec_combinatorialpne.h
================================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_algorithms_EPEC_epec_combinatorialpne.h>` (``include/algorithms/EPEC/epec_combinatorialpne.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include "epec_polybase.h"
   #include <armadillo>
   #include <boost/log/trivial.hpp>
   #include <gurobi_c++.h>
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   
   namespace Algorithms {
     namespace EPEC {
   
        class CombinatorialPNE : public PolyBase {
        public:
           CombinatorialPNE(GRBEnv *env, Game::EPEC *EPECObject, bool poly = true)
                : PolyBase(env, EPECObject){};
           ;
           void solve() { this->solveWithExcluded(std::vector<std::set<unsigned long int>>{}); }
           void solveWithExcluded(const std::vector<std::set<unsigned long int>> &excludeList = {});
   
        private:
           // Making the method private
           void combPNE(std::vector<long int>                           combination,
                            const std::vector<std::set<unsigned long int>> &excludeList);
        };
     } // namespace EPEC
   } // namespace Algorithms
