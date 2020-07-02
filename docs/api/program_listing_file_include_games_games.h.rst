
.. _program_listing_file_include_games_games.h:

Program Listing for File games.h
================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_games_games.h>` (``include/games/games.h``)

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
   #include "support/codes.h"
   #include <armadillo>
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   
   namespace Game {
   
     template <typename DataObjectType> class AbstractGame {
     protected:
        std::chrono::high_resolution_clock::time_point InitTime;
        ZEROStatistics<DataObjectType> Stats = ZEROStatistics<DataObjectType>(DataObjectType());
        ;                             
        GRBEnv *     Env;             
        unsigned int NumVariables{0}; 
        unsigned int NumPlayers{0};   
        bool NashEquilibrium{false};  
     public:
        AbstractGame(GRBEnv *env) : Env{env} {};
        AbstractGame()                  = default; // No default constructor
        AbstractGame(AbstractGame &)    = delete;  // Abstract class - no copy constructor
        ~AbstractGame()                 = default; // Destructor to free data
        virtual const void findNashEq() = 0;       
        virtual bool       isSolved(double tol = 1e-5)
             const = 0; 
        virtual bool
        isPureStrategy(double tol = 1e-5) const = 0; 
        ZEROStatistics<DataObjectType> getStatistics() const { return this->Stats; }
        void                           setNumThreads(unsigned int t) {
         this->Stats.AlgorithmData.Threads.set(t);
         this->Env->set(GRB_IntParam_Threads, t);
        }
        void setRandomSeed(unsigned int t) { this->Stats.AlgorithmData.RandomSeed.set(t); }
   
        void setIndicators(bool val) { this->Stats.AlgorithmData.IndicatorConstraints.set(val); }
   
        void setPureNashEquilibrium(bool val) { this->Stats.AlgorithmData.PureNashEquilibrium = val; }
        void setDeviationTolerance(double val) {
           this->Stats.AlgorithmData.DeviationTolerance.set(val);
        }
   
        void setTimeLimit(double val) { this->Stats.AlgorithmData.TimeLimit.set(val); }
        int  getNumVar() const noexcept { return this->NumVariables; }
        int  getNumPlayers() const noexcept { return this->NumPlayers; }
     };
   
     /*
       * class aGame : public AbstractGame<Data::aGame::DataObject>{
       *
       *    public:
       *      //Override AbstractGame methods
       *     const void findNashEq() override;
       *     bool       isSolved(double tol = 1e-5) const override;
       *     bool isPureStrategy(double tol = 1e-5) const override;
       *
       * }
       *
       *
       */
   
   } // namespace Game
   
   #include "games/epec.h"
   #include "games/ipg.h"
   #include "games/nash.h"
