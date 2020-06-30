
.. _program_listing_file_include_games_nash.h:

Program Listing for File nash.h
===============================

|exhale_lsh| :ref:`Return to documentation for file <file_include_games_nash.h>` (``include/games/nash.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include "zero.h"
   #include <armadillo>
   #include <gurobi_c++.h>
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   
   namespace Game {
     class NashGame {
     private:
        GRBEnv *     Env = nullptr;
        arma::sp_mat LeaderConstraints;                        
        arma::vec    LeaderConstraintsRHS;                     
        unsigned int NumPlayers;                               
        std::vector<std::shared_ptr<QP_Param>> Players;        
        arma::sp_mat                           MarketClearing; 
        arma::vec                              MCRHS; 
   
        std::vector<unsigned int> PrimalPosition;
        std::vector<unsigned int> DualPosition;
        unsigned int MC_DualPosition;
        unsigned int LeaderPosition;
        unsigned int numLeaderVar;
   
        void setPositions();
   
     public: // Constructors
        explicit NashGame(GRBEnv *e) noexcept : Env{e} {};
   
        explicit NashGame(GRBEnv *                               e,
                                std::vector<std::shared_ptr<QP_Param>> players,
                                arma::sp_mat                           MC,
                                arma::vec                              MCRHS,
                                unsigned int                           nLeadVar = 0,
                                arma::sp_mat                           leadA    = {},
                                arma::vec                              leadRHS  = {});
   
        // Copy constructor
        NashGame(const NashGame &N);
   
        ~NashGame() = default;
   
        // Verbose declaration
        friend std::ostream &operator<<(std::ostream &os, const NashGame &N) {
           os << '\n';
           os << "--------------------------------------------------------------------"
                   "---"
               << '\n';
           os << "Nash Game with " << N.NumPlayers << " players" << '\n';
           os << "--------------------------------------------------------------------"
                   "---"
               << '\n';
           os << "Number of primal variables:\t\t\t " << N.getNprimals() << '\n';
           os << "Number of dual variables:\t\t\t " << N.getNumDualVars() << '\n';
           os << "Number of shadow price dual variables:\t\t " << N.getNumShadow() << '\n';
           os << "Number of leader variables:\t\t\t " << N.getNumLeaderVars() << '\n';
           os << "--------------------------------------------------------------------"
                   "---"
               << '\n';
           return os;
        }
   
        inline unsigned int getNprimals() const {
           /***
            * Number of primal variables is the sum of the "y" variables present in
            * each player's Game::QP_Param
            */
           return this->PrimalPosition.back();
        }
   
        inline unsigned int getNumShadow() const { return this->MCRHS.n_rows; }
   
        inline unsigned int getNumLeaderVars() const { return this->numLeaderVar; }
   
        inline unsigned int getNumDualVars() const {
           return this->DualPosition.back() - this->DualPosition.front() + 0;
        }
   
        // Position of variables
        inline unsigned int getPrimalLoc(unsigned int i = 0) const { return PrimalPosition.at(i); }
   
        inline unsigned int getMCDualLoc() const { return MC_DualPosition; }
   
        inline unsigned int getLeaderLoc() const { return LeaderPosition; }
   
        inline unsigned int getDualLoc(unsigned int i = 0) const { return DualPosition.at(i); }
   
        // Members
        const NashGame &formulateLCP(arma::sp_mat &M,
                                               arma::vec &   q,
                                               perps &       Compl,
                                               bool          writeToFile = false,
                                               std::string   M_name      = "dat/LCP.txt",
                                               std::string   q_name      = "dat/q.txt") const;
   
        arma::sp_mat rewriteLeadCons() const;
   
        inline arma::vec getLeadRHS() const { return this->LeaderConstraintsRHS; }
   
        inline arma::vec getMCLeadRHS() const {
           return arma::join_cols(arma::join_cols(this->LeaderConstraintsRHS, this->MCRHS),
                                         -this->MCRHS);
        }
   
        // Check solution and correctness
        std::unique_ptr<GRBModel>
        respond(unsigned int player, const arma::vec &x, bool fullvec = true) const;
   
        double
        respondSol(arma::vec &sol, unsigned int player, const arma::vec &x, bool fullvec = true) const;
   
        arma::vec computeQPObjectiveValues(const arma::vec &x, bool checkFeas = false) const;
   
        bool isSolved(const arma::vec &sol,
                           unsigned int &   violPlayer,
                           arma::vec &      violSol,
                           double           tol = 1e-4) const;
   
        //  Modify NashGame members
        NashGame &addDummy(unsigned int par = 0, int position = -1);
   
        NashGame &addLeadCons(const arma::vec &a, double b);
   
        // Read/Write Nashgame functions
        void write(const std::string &filename, bool append = true, bool KKT = false) const;
   
        void save(const std::string &filename, bool erase = true) const;
   
        long int  load(const std::string &filename, long int pos = 0);
        arma::vec computeQPObjectiveValuesWithoutOthers(const arma::vec &x) const;
     };
   
     std::ostream &operator<<(std::ostream &os, const QP_Param &Q);
   
     std::ostream &operator<<(std::ostream &ost, const perps &C);
   
     void print(const perps &C) noexcept;
   } // namespace Game
