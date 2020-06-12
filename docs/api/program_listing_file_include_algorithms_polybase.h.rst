
.. _program_listing_file_include_algorithms_polybase.h:

Program Listing for File polybase.h
===================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_algorithms_polybase.h>` (``include/algorithms/polybase.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include "algorithms.h"
   #include "epecsolve.h"
   #include <boost/log/trivial.hpp>
   
   namespace Algorithms {
   class PolyBase : public Algorithm {
   protected:
     std::vector<std::shared_ptr<Game::PolyLCP>> PolyLCP{};
   
     void postSolving() override {
       for (unsigned int i = 0; i < this->EPECObject->NumPlayers; i++)
         this->EPECObject->Stats.FeasiblePolyhedra.at(i) =
             this->PolyLCP.at(i)->getFeasiblePolyhedra();
       this->EPECObject->Stats.PureNashEquilibrium = this->isPureStrategy();
     }
   
   public:
     PolyBase(GRBEnv *env, Game::EPEC *EPECObject) {
       /*
        *  The method will reassign the LCP fields in the EPEC object to new
        * PolyLCP objects
        */
       this->EPECObject = EPECObject;
       this->Env = env;
       this->EPECObject->Stats.AlgorithmParam.PolyLcp = true;
       this->PolyLCP =
           std::vector<std::shared_ptr<Game::PolyLCP>>(EPECObject->NumPlayers);
       for (unsigned int i = 0; i < EPECObject->NumPlayers; i++) {
         this->PolyLCP.at(i) =
             std::shared_ptr<Game::PolyLCP>(new class Game::PolyLCP(
                 this->Env, *EPECObject->PlayersLowerLevels.at(i).get()));
         EPECObject->PlayersLCP.at(i) = this->PolyLCP.at(i);
       }
     }
     bool isSolved(unsigned int *countryNumber, arma::vec *profitableDeviation,
                   double tol = -1) const;
     bool isSolved(double tol = -1) const override;
     void makeThePureLCP(bool indicators);
   
     double getValLeadFollPoly(unsigned int i, unsigned int j, unsigned int k,
                               double tol = 1e-5) const;
   
     double getValLeadLeadPoly(unsigned int i, unsigned int j, unsigned int k,
                               double tol = 1e-5) const;
   
     double getValProbab(unsigned int i, unsigned int k) const;
   
     bool isPureStrategy(unsigned int i, double tol = 1e-5) const;
   
     bool isPureStrategy(double tol = 1e-5) const override;
   
     std::vector<unsigned int> mixedStrategyPoly(unsigned int i,
                                                 double tol = 1e-5) const;
     unsigned int getPositionLeadFollPoly(unsigned int i, unsigned int j,
                                          unsigned int k) const;
   
     unsigned int getPositionLeadLeadPoly(unsigned int i, unsigned int j,
                                          unsigned int k) const;
   
     unsigned int getNumPolyLead(unsigned int i) const;
   
     unsigned int getPositionProbab(unsigned int i, unsigned int k) const;
   };
   } // namespace Algorithms
   
   #include "algorithms/combinatorialpne.h"
   #include "algorithms/fullenumeration.h"
   #include "algorithms/innerapproximation.h"
