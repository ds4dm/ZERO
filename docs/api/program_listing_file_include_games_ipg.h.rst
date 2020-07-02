
.. _program_listing_file_include_games_ipg.h:

Program Listing for File ipg.h
==============================

|exhale_lsh| :ref:`Return to documentation for file <file_include_games_ipg.h>` (``include/games/ipg.h``)

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
   #include <armadillo>
   #include <gurobi_c++.h>
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   
   namespace Data {
     namespace IPG {
   
        enum class Algorithms {
           Oracle 
        };
   
        class DataObject : public ZEROAlgorithmData {
        public:
           Attr<Data::IPG::Algorithms> Algorithm = {
                Data::IPG::Algorithms::Oracle}; 
           DataObject(){};
        };
     } // namespace IPG
   } // namespace Data
   
   
   namespace Game {
   
     class IPG : public AbstractGame<Data::IPG::DataObject> {
   
     protected: // Datafields
        std::vector<std::shared_ptr<MathOpt::IP_Param>>
             PlayersIP{}; 
   
        std::vector<unsigned int> PlayerVariables{}; 
   
        bool Finalized{false};           
        std::vector<arma::vec> Solution; 
   
     private:
        void getXMinusI(const arma::vec &x, const unsigned int &i, arma::vec &xMinusI) const;
        void getXofI(const arma::vec &x, const unsigned int &i, arma::vec &xOfI) const;
   
        bool computeNashEq(double localTimeLimit = -1.0, bool check = false);
        void finalize();
   
     public: // functions
        friend class Algorithms::IPG::Oracle;
        IPG(GRBEnv *env) { this->Env = env; };
        IPG(GRBEnv *env, std::vector<std::shared_ptr<MathOpt::IP_Param>> players);
   
        const void findNashEq() override;
        bool       isSolved(double tol = 1e-5) const override;
        bool isPureStrategy(double tol = 1e-5) const override; 
   
        std::unique_ptr<GRBModel> respondModel(const unsigned int i, const arma::vec &x) const;
   
        const std::vector<arma::vec> getX() const { return this->Solution; }
   
        ZEROStatistics<Data::IPG::DataObject> getStatistics() const { return this->Stats; }
   
        void setAlgorithm(Data::IPG::Algorithms algorithm);
     };
   
   } // namespace Game
   
   #include "algorithms/IPG/ipg_oracle.h"
