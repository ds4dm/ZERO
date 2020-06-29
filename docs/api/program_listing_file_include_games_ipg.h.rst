
.. _program_listing_file_include_games_ipg.h:

Program Listing for File ipg.h
==============================

|exhale_lsh| :ref:`Return to documentation for file <file_include_games_ipg.h>` (``include/games/ipg.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include "zero.h"
   #include <algorithms/IPG/ipg_oracle.h>
   #include <armadillo>
   #include <gurobi_c++.h>
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   
   namespace Game {
   class IP_Param : public MP_Param
   // Shape of C is Ny\times Nx
   {
   private:
     // Gurobi environment and model
     GRBEnv *Env;
     GRBModel IPModel;          
     arma::vec bounds;          
     std::vector<int> integers; 
     bool madeModel{false};     
   
     // These methods should be inaccessible to the inheritor, since we have a
     // different structure.
     using MP_Param::set;
   
   public: // Constructors
     explicit IP_Param(GRBEnv *env = nullptr) : Env{env}, IPModel{(*env)} {
       this->size();
     }
   
     explicit IP_Param(arma::sp_mat C, arma::sp_mat B, arma::vec b, arma::vec c,
                       arma::vec bounds, std::vector<int> integers,
                       GRBEnv *env = nullptr)
         : Env{env}, IPModel{(*env)} {
       this->Q.zeros(0);
       this->A.zeros(0);
       this->set(Q, C, A, B, c, b);
       this->bounds = bounds;
       this->integers = integers;
       this->size();
       if (!this->dataCheck())
         throw ZEROException(ZEROErrorCode::InvalidData, "dataCheck() failed");
     }
   
     std::vector<int> getIntegers() const { return this->integers; }
     arma::vec getBounds() const { return this->bounds; }
     void makeModel();
     void addConstraints(const arma::sp_mat A, const arma::vec b);
   
     IP_Param(const IP_Param &ipg)
         : MP_Param(ipg), Env{ipg.Env}, IPModel{ipg.IPModel}, madeModel{
                                                                  ipg.madeModel} {
       this->size();
     };
   
     // Override setters
     IP_Param &set(const arma::sp_mat &C, const arma::sp_mat &B,
                   const arma::vec &b, const arma::vec &c, const arma::vec &bounds,
                   const std::vector<int> &integers); // Copy data into this
     IP_Param &set(arma::sp_mat &C, arma::sp_mat &&B, arma::vec &&b, arma::vec &&c,
                   arma::vec &&bounds,
                   std::vector<int> &&integers); // Copy data into this
   
     IP_Param &set(const QP_Objective &obj, const QP_Constraints &cons,
                   const arma::vec &bounds = {},
                   const std::vector<int> &integers = {});
     IP_Param &set(QP_Objective &&obj, QP_Constraints &&cons,
                   arma::vec &&bounds = {}, std::vector<int> &&integers = {});
   
     bool operator==(const IP_Param &IPG2) const;
   
     std::unique_ptr<GRBModel> solveFixed(arma::vec x, bool solve);
   
     double computeObjective(const arma::vec &y, const arma::vec &x,
                             bool checkFeas = true, double tol = 1e-6) const;
   
     inline bool isPlayable(const IP_Param &P) const
     {
       bool b1, b2, b3;
       b1 = (this->Nx + this->Ny) == (P.getNx() + P.getNy());
       b2 = this->Nx >= P.getNy();
       b3 = this->Ny <= P.getNx();
       return b1 && b2 && b3;
     }
   
     IP_Param &addDummy(unsigned int pars, unsigned int vars = 0,
                        int position = -1) override;
   
     void write(const std::string &filename, bool append) const override;
   
     double computeObjectiveWithoutOthers(const arma::vec &y) const;
     arma::vec getConstraintViolations(const arma::vec y, double tol);
   };
   
   class IPG {
   private:
     unsigned int NumVariables{0};
     unsigned int NumPlayers{0};
   
   protected: // Datafields
     std::vector<std::shared_ptr<Game::IP_Param>>
         PlayersIP{}; 
   
     std::vector<unsigned int>
         PlayerVariables{}; 
   
     GRBEnv *Env;
     bool Finalized{false}; 
     bool NashEquilibrium{
         false}; 
     std::chrono::high_resolution_clock::time_point InitTime;
     Data::EPEC::EPECStatistics Stats; 
     std::vector<arma::vec>
         Solution; 
   
   private:
     void getXMinusI(const arma::vec &x, const unsigned int &i,
                     arma::vec &xMinusI) const;
     void getXofI(const arma::vec &x, const unsigned int &i,
                  arma::vec &xOfI) const;
   
     bool computeNashEq(double localTimeLimit = -1.0, bool check = false);
     void finalize();
   
   public: // functions
     friend class Algorithms::IPG::Oracle;
     IPG() = delete;      // No default constructor
     IPG(IPG &) = delete; // Abstract class - no copy constructor
     ~IPG() = default;    // Destructor to free data
   
     IPG(GRBEnv *env)
         : Env{env} {}; 
     IPG(GRBEnv *env, std::vector<std::shared_ptr<Game::IP_Param>> players);
   
     const void findNashEq();
     bool isSolved(double tol = 1e-5) const;
     bool isPureStrategy(
         double tol = 1e-5) const; 
   
     std::unique_ptr<GRBModel> respondModel(const unsigned int i,
                                            const arma::vec &x) const;
   
     const std::vector<arma::vec> getX() const { return this->Solution; }
   
     const Data::EPEC::EPECStatistics getStatistics() const { return this->Stats; }
   
     void setAlgorithm(Data::EPEC::EPECalgorithm algorithm);
   
     Data::EPEC::EPECalgorithm getAlgorithm() const {
       return this->Stats.AlgorithmParam.Algorithm;
     }
   
     void setNumThreads(unsigned int t) {
       this->Stats.AlgorithmParam.Threads = t;
       this->Env->set(GRB_IntParam_Threads, t);
     }
   
     unsigned int getNumThreads() const {
       return this->Stats.AlgorithmParam.Threads;
     }
   
     void setPureNashEquilibrium(bool val) {
       this->Stats.AlgorithmParam.PureNashEquilibrium = val;
     }
   
     bool getPureNashEquilibrium() const {
       return this->Stats.AlgorithmParam.PureNashEquilibrium;
     }
   
     void setDeviationTolerance(double val) {
       this->Stats.AlgorithmParam.DeviationTolerance = val;
     }
   
     double getDeviationTolerance() const {
       return this->Stats.AlgorithmParam.DeviationTolerance;
     }
   
     void setTimeLimit(double val) { this->Stats.AlgorithmParam.TimeLimit = val; }
   
     double getTimeLimit() const { return this->Stats.AlgorithmParam.TimeLimit; }
   
     // Methods to get positions of variables
     // The below are all const functions which return an unsigned int.
     int getNumVar() const noexcept { return this->NumVariables; }
   
     unsigned int getNumPlayers() const noexcept { return this->NumPlayers; }
   };
   
   } // namespace Game
