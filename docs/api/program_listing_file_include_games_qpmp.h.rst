
.. _program_listing_file_include_games_qpmp.h:

Program Listing for File qpmp.h
===============================

|exhale_lsh| :ref:`Return to documentation for file <file_include_games_qpmp.h>` (``include/games/qpmp.h``)

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
   class MP_Param {
   protected:
     // Data representing the parameterized QP
     arma::sp_mat Q, A, B, C;
     arma::vec c, b;
     // Object for sizes and integrity check
     unsigned int Nx, Ny, Ncons;
   
     const unsigned int size();
   
     bool dataCheck(bool forceSymmetry = true) const;
   
     virtual inline bool finalize() {
       this->size();
       return this->dataCheck();
     } 
   
   public:
     // Default constructors
     MP_Param() = default;
   
     MP_Param(const MP_Param &M) = default;
   
     // Getters and setters
     arma::sp_mat getQ() const {
       return this->Q;
     } 
     arma::sp_mat getC() const {
       return this->C;
     } 
     arma::sp_mat getA() const {
       return this->A;
     } 
     arma::sp_mat getB() const {
       return this->B;
     } 
     arma::vec getc() const {
       return this->c;
     } 
     arma::vec getb() const {
       return this->b;
     } 
     unsigned int getNx() const {
       return this->Nx;
     } 
     unsigned int getNy() const {
       return this->Ny;
     } 
   
     MP_Param &setQ(const arma::sp_mat &Q) {
       this->Q = Q;
       return *this;
     } 
     MP_Param &setC(const arma::sp_mat &C) {
       this->C = C;
       return *this;
     } 
     MP_Param &setA(const arma::sp_mat &A) {
       this->A = A;
       return *this;
     } 
     MP_Param &setB(const arma::sp_mat &B) {
       this->B = B;
       return *this;
     } 
     MP_Param &setc(const arma::vec &c) {
       this->c = c;
       return *this;
     } 
     MP_Param &setb(const arma::vec &b) {
       this->b = b;
       return *this;
     } 
   
     // Setters and advanced constructors
     virtual MP_Param &set(const arma::sp_mat &Q, const arma::sp_mat &C,
                           const arma::sp_mat &A, const arma::sp_mat &B,
                           const arma::vec &c,
                           const arma::vec &b); // Copy data into this
     virtual MP_Param &set(arma::sp_mat &&Q, arma::sp_mat &&C, arma::sp_mat &&A,
                           arma::sp_mat &&B, arma::vec &&c,
                           arma::vec &&b); // Move data into this
     virtual MP_Param &set(const QP_Objective &obj, const QP_Constraints &cons);
   
     virtual MP_Param &set(QP_Objective &&obj, QP_Constraints &&cons);
   
     virtual MP_Param &addDummy(unsigned int pars, unsigned int vars = 0,
                                int position = -1);
   
     virtual void write(const std::string &filename, bool append = true) const;
   
     static bool dataCheck(const QP_Objective &obj, const QP_Constraints &cons,
                           bool checkObj = true, bool checkCons = true);
   };
   
   class QP_Param : public MP_Param
   // Shape of C is Ny\times Nx
   {
   private:
     // Gurobi environment and model
     GRBEnv *Env;
     GRBModel QuadModel;
     bool madeyQy;
   
     int makeyQy();
   
   public: // Constructors
     explicit QP_Param(GRBEnv *env = nullptr)
         : Env{env}, QuadModel{(*env)}, madeyQy{false} {
       this->size();
     }
   
     QP_Param(arma::sp_mat Q, arma::sp_mat C, arma::sp_mat A, arma::sp_mat B,
              arma::vec c, arma::vec b, GRBEnv *env = nullptr)
         : Env{env}, QuadModel{(*env)}, madeyQy{false} {
       this->set(Q, C, A, B, c, b);
       this->size();
       if (!this->dataCheck())
         throw ZEROException(ZEROErrorCode::InvalidData, "dataCheck() failed");
     }
   
     QP_Param(const QP_Param &Qu)
         : MP_Param(Qu), Env{Qu.Env}, QuadModel{Qu.QuadModel}, madeyQy{
                                                                   Qu.madeyQy} {
       this->size();
     };
   
     // Override setters
     QP_Param &set(const arma::sp_mat &Q, const arma::sp_mat &C,
                   const arma::sp_mat &A, const arma::sp_mat &B,
                   const arma::vec &c,
                   const arma::vec &b) final; // Copy data into this
     QP_Param &set(arma::sp_mat &&Q, arma::sp_mat &&C, arma::sp_mat &&A,
                   arma::sp_mat &&B, arma::vec &&c,
                   arma::vec &&b) final; // Move data into this
     QP_Param &set(const QP_Objective &obj, const QP_Constraints &cons) final;
   
     QP_Param &set(QP_Objective &&obj, QP_Constraints &&cons) final;
   
     bool operator==(const QP_Param &Q2) const;
   
     // Other methods
     unsigned int KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const;
   
     std::unique_ptr<GRBModel> solveFixed(arma::vec x, bool solve);
   
     double computeObjective(const arma::vec &y, const arma::vec &x,
                             bool checkFeas = true, double tol = 1e-6) const;
   
     inline bool isPlayable(const QP_Param &P) const
     {
       bool b1, b2, b3;
       b1 = (this->Nx + this->Ny) == (P.getNx() + P.getNy());
       b2 = this->Nx >= P.getNy();
       b3 = this->Ny <= P.getNx();
       return b1 && b2 && b3;
     }
   
     QP_Param &addDummy(unsigned int pars, unsigned int vars = 0,
                        int position = -1) override;
   
     void write(const std::string &filename, bool append) const override;
   
     void save(const std::string &filename, bool erase = true) const;
   
     long int load(const std::string &filename, long int pos = 0);
     double computeObjectiveWithoutOthers(const arma::vec &y) const;
     arma::vec getConstraintViolations(const arma::vec x, const arma::vec y,
                                       double tol);
   };
   } // namespace Game
   
   //#include "ipg.h"
