
.. _program_listing_file_include_games_algorithms_EPEC_epec_polybase.h:

Program Listing for File epec_polybase.h
========================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_games_algorithms_EPEC_epec_polybase.h>` (``include/games/algorithms/EPEC/epec_polybase.h``)

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
   #include "epec_algorithms.h"
   #include "mathopt/lcp/lcp.h"
   #include "zero.h"
   #include <boost/log/trivial.hpp>
   
   namespace Algorithms {
     namespace EPEC {
        class PolyBase : public Algorithm {
        protected:
           std::vector<std::shared_ptr<MathOpt::PolyLCP>> PolyLCP{};
   
           void postSolving() override {
             std::vector<unsigned int> feasible;
             for (unsigned int i = 0; i < this->EPECObject->NumPlayers; i++)
                feasible.push_back(this->PolyLCP.at(i)->getFeasiblePolyhedra());
             this->EPECObject->Stats.AlgorithmData.FeasiblePolyhedra.set(feasible);
             this->EPECObject->Stats.PureNashEquilibrium = this->isPureStrategy();
           }
   
        public:
           PolyBase(GRBEnv *env, Game::EPEC *EPECObject) {
             /*
               *  The method will reassign the LCP fields in the EPEC object to new
               * PolyLCP objects
               */
             this->EPECObject = EPECObject;
             this->Env        = env;
             this->PolyLCP    = std::vector<std::shared_ptr<MathOpt::PolyLCP>>(EPECObject->NumPlayers);
             for (unsigned int i = 0; i < EPECObject->NumPlayers; i++) {
                this->PolyLCP.at(i) = std::shared_ptr<MathOpt::PolyLCP>(
                     new class MathOpt::PolyLCP(this->Env, *EPECObject->PlayersLowerLevels.at(i).get()));
                EPECObject->PlayersLCP.at(i) = this->PolyLCP.at(i);
             }
           }
           bool
                 isSolved(unsigned int *countryNumber, arma::vec *profitableDeviation, double tol = -1) const;
           bool isSolved(double tol = -1) const override;
           void makeThePureLCP(bool indicators);
   
           double
           getValLeadFollPoly(unsigned int i, unsigned int j, unsigned int k, double tol = 1e-5) const;
   
           double
           getValLeadLeadPoly(unsigned int i, unsigned int j, unsigned int k, double tol = 1e-5) const;
   
           double getValProbab(unsigned int i, unsigned int k) const;
   
           bool isPureStrategy(unsigned int i, double tol = 1e-5) const;
   
           bool isPureStrategy(double tol = 1e-5) const override;
   
           std::vector<unsigned int> mixedStrategyPoly(unsigned int i, double tol = 1e-5) const;
           unsigned int getPositionLeadFollPoly(unsigned int i, unsigned int j, unsigned int k) const;
   
           unsigned int getPositionLeadLeadPoly(unsigned int i, unsigned int j, unsigned int k) const;
   
           unsigned int getNumPolyLead(unsigned int i) const;
   
           unsigned int getPositionProbab(unsigned int i, unsigned int k) const;
        };
     } // namespace EPEC
   } // namespace Algorithms
   
   #include "epec_combinatorialpne.h"
   #include "epec_fullenumeration.h"
   #include "epec_innerapproximation.h"
