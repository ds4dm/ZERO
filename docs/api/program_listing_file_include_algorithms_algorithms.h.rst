
.. _program_listing_file_include_algorithms_algorithms.h:

Program Listing for File algorithms.h
=====================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_algorithms_algorithms.h>` (``include/algorithms/algorithms.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include "epecsolve.h"
   
   namespace Algorithms {
   // the class generic stores some common information for algorithms
   class Algorithm {
   protected:
     GRBEnv *Env;            
     Game::EPEC *EPECObject; 
     virtual void
     postSolving() = 0; 
   
   public:
     virtual void solve() = 0; 
     virtual bool isSolved(double tol = -1)
         const = 0; 
     virtual bool isPureStrategy(double tol = -1)
         const = 0; 
   };
   // Second level inheritor for polyhedral inner approximations or full
   // enumeration
   class PolyBase;
   // The following algorithms are children of polybase
   class FullEnumeration;
   class InnerApproximation;
   class CombinatorialPNE;
   
   // Then, second level inheritor for the outer approximation
   class OuterApproximation;
   
   } // namespace Algorithms
