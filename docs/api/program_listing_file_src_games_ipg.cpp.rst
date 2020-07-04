
.. _program_listing_file_src_games_ipg.cpp:

Program Listing for File ipg.cpp
================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_games_ipg.cpp>` (``src/games/ipg.cpp``)

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
   
   
   #include "games/ipg.h"
   #include "games/algorithms/IPG/ipg_algorithms.h"
   #include "zero.h"
   #include <armadillo>
   #include <boost/log/trivial.hpp>
   #include <memory>
   
   Game::IPG::IPG(
        GRBEnv *                                        env,    
        std::vector<std::shared_ptr<MathOpt::IP_Param>> players 
   ) {
     this->Env       = env;
     this->PlayersIP = players;
     this->finalize();
   }
   void Game::IPG::finalize() {
     this->NumPlayers      = this->PlayersIP.size();
     this->PlayerVariables = std::vector<unsigned int>(this->NumPlayers);
     this->Solution        = std::vector<arma::vec>(this->NumPlayers);
     this->NumVariables    = 0;
     for (unsigned int i = 0; i < this->NumPlayers; ++i) {
        PlayerVariables.at(i) = this->PlayersIP.at(i)->getNy();
        this->NumVariables += PlayerVariables.at(i);
        this->PlayersIP.at(i)->finalize();
     }
     this->Finalized = true;
   }
   
   void Game::IPG::getXMinusI(
        const arma::vec &x,         
        const unsigned int &i,      
        arma::vec &         xMinusI 
        ) const {
     if (this->NumVariables != x.size())
        throw ZEROException(ZEROErrorCode::Assertion, "Invalid size of x");
   
     xMinusI.zeros(this->NumVariables - this->PlayerVariables.at(i));
   
     for (unsigned int j = 0, posIn = 0, posOut = 0; j < this->NumPlayers; ++j) {
        if (i != j) {
           xMinusI.subvec(posOut, posOut + this->PlayerVariables.at(j) - 1) =
                x.subvec(posIn, posIn + this->PlayerVariables.at(j) - 1);
           posOut += this->PlayerVariables.at(j);
        }
        posIn += this->PlayerVariables.at(j);
     }
   }
   
   void Game::IPG::getXofI(const arma::vec &x, 
                                   const unsigned int &i,   
                                   arma::vec &         xOfI 
                                   ) const {
     if (this->NumVariables != x.size())
        throw ZEROException(ZEROErrorCode::Assertion, "Invalid size of x");
   
     int count = 0;
     for (unsigned int j = 0; j < i; ++j)
        count += this->PlayerVariables.at(j);
   
     xOfI.zeros(this->PlayerVariables.at(i));
     xOfI = x.subvec(count, count + this->PlayerVariables.at(i) - 1);
   }
   
   
   const void Game::IPG::findNashEq() {
     std::stringstream final_msg;
     if (!this->Finalized)
        this->finalize();
   
     switch (this->Stats.AlgorithmData.Algorithm.get()) {
     case Data::IPG::Algorithms::Oracle: {
        final_msg << "Oracle Algorithm completed. ";
        this->Algorithm = std::shared_ptr<Algorithms::IPG::Oracle>(
             new class Algorithms::IPG::Oracle(this->Env, this));
        this->Algorithm->solve();
     } break;
     }
   }
   
   bool Game::IPG::isPureStrategy(double tol) const { return this->Algorithm->isPureStrategy(); }
   bool Game::IPG::isSolved(double tol) const { return this->Algorithm->isSolved(); }
   
   
   std::string std::to_string(const Data::IPG::Algorithms al) {
     switch (al) {
     case Data::IPG::Algorithms::Oracle:
        return std::string("Oracle");
     default:
        return std::string("UNKNOWN_ALGORITHM_") + std::to_string(static_cast<int>(al));
     }
   }
