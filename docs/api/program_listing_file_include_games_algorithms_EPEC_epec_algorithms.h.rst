
.. _program_listing_file_include_games_algorithms_EPEC_epec_algorithms.h:

Program Listing for File epec_algorithms.h
==========================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_games_algorithms_EPEC_epec_algorithms.h>` (``include/games/algorithms/EPEC/epec_algorithms.h``)

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
   #include "zero.h"
   
   namespace Algorithms {
     namespace EPEC {
        // the class generic stores some common information for algorithms
        class Algorithm {
        protected:
           GRBEnv *     Env;               
           Game::EPEC * EPECObject;        
           virtual void postSolving() = 0; 
   
        public:
           virtual void solve() = 0; 
           virtual bool
           isSolved(double tol = -1) const = 0; 
   
   
           virtual bool isPureStrategy(
                double tol = -1) const = 0; 
        };
   
     } // namespace EPEC
   } // namespace Algorithms
   
   #include "epec_outerapproximation.h"
   #include "epec_polybase.h"
