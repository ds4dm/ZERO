
.. _program_listing_file_include_mathopt_lcp_lcp.h:

Program Listing for File lcp.h
==============================

|exhale_lsh| :ref:`Return to documentation for file <file_include_mathopt_lcp_lcp.h>` (``include/mathopt/lcp/lcp.h``)

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
   
   namespace Data {
     namespace LCP {
   
        enum class PolyhedraStrategy {
           Sequential        = 0, 
           ReverseSequential = 1, 
           Random = 2 
        };
     }
   } // namespace Data
   
   namespace MathOpt {
   
     class LCP {
   
     protected:
        // Essential data ironment for MIP/LP solves
        GRBEnv *     Env;   
        arma::sp_mat M;     
        arma::vec    q;     
        perps        Compl; 
        unsigned int LeadStart{1}, LeadEnd{0}, NumberLeader{0};
        arma::sp_mat _A = {};
        arma::vec    _b = {}; 
        arma::sp_mat _Acut = {};
        arma::vec    _bcut = {};           
        bool         MadeRlxdModel{false}; 
        unsigned int nR, nC;
   
        GRBModel RlxdModel; 
   
        bool errorCheck(bool throwErr = true) const;
   
        void defConst(GRBEnv *env);
   
        void makeRelaxed();
   
        /* Solving relaxations and restrictions */
        std::unique_ptr<GRBModel> LCPasMIP(std::vector<unsigned int> FixEq  = {},
                                                       std::vector<unsigned int> FixVar = {},
                                                       bool                      solve  = false);
   
        template <class T> inline bool isZero(const T val) const { return (val >= -Eps && val <= Eps); }
   
        std::unique_ptr<spmat_Vec> Ai; 
        std::unique_ptr<vec_Vec> bi;   
   
        inline std::vector<short int> solEncode(GRBModel *model) const;
   
        unsigned int convexHull(arma::sp_mat &A, arma::vec &b);
   
     public:
        long double BigM{1e7};    
        double      Eps{1e-6};    
        double      EpsInt{1e-8}; 
        bool UseIndicators{true}; 
   
        LCP() = delete;
   
        explicit LCP(GRBEnv *e)
             : Env{e}, RlxdModel(*e){}; 
   
        LCP(GRBEnv *     env,
             arma::sp_mat M,
             arma::vec    q,
             unsigned int leadStart,
             unsigned     leadEnd,
             arma::sp_mat A = {},
             arma::vec    b = {}); // Constructor with M,q,leader posn
        LCP(GRBEnv *     env,
             arma::sp_mat M,
             arma::vec    q,
             perps        Compl,
             arma::sp_mat A = {},
             arma::vec    b = {}); // Constructor with M, q, compl pairs
        LCP(GRBEnv *env, const Game::NashGame &N);
   
        ~LCP() = default;
   
        inline arma::sp_mat  getM() { return this->M; }        
        inline arma::sp_mat *getMstar() { return &(this->M); } 
        inline arma::vec     getq() { return this->q; }        
        inline unsigned int  getNumberLeader() {
         return this->NumberLeader;
        }                                                           
        inline arma::vec *        getqstar() { return &(this->q); } 
        const inline unsigned int getLStart() {
           return LeadStart;
        }                                                       
        const inline unsigned int getLEnd() { return LeadEnd; } 
        inline perps              getCompl() { return this->Compl; } 
        void                      print(std::string end = "\n");     
        inline unsigned int       getNumCols() { return this->M.n_cols; };
   
        inline unsigned int getNumRows() { return this->M.n_rows; };
   
        bool extractSols(GRBModel *model, arma::vec &z, arma::vec &x, bool extractZ = false) const;
   
        /* Getting single point solutions */
        std::unique_ptr<GRBModel> LCPasQP(bool solve = false);
   
        std::unique_ptr<GRBModel> LCPasMIP(bool solve = false);
   
        std::unique_ptr<GRBModel> MPECasMILP(const arma::sp_mat &C,
                                                         const arma::vec &   c,
                                                         const arma::vec &   x_minus_i,
                                                         bool                solve = false);
   
        std::vector<short int> solEncode(const arma::vec &z, const arma::vec &x) const;
   
        std::unique_ptr<GRBModel> MPECasMIQP(const arma::sp_mat &Q,
                                                         const arma::sp_mat &C,
                                                         const arma::vec &   c,
                                                         const arma::vec &   x_minus_i,
                                                         bool                solve = false);
   
        std::unique_ptr<GRBModel> LCPasMIP(std::vector<short int> Fixes, bool solve);
   
        void write(std::string filename, bool append = true) const;
   
        void save(std::string filename, bool erase = true) const;
   
        long int load(std::string filename, long int pos = 0);
   
        virtual void makeQP(QP_Objective &QP_obj, QP_Param &QP);
   
        void addCustomCuts(const arma::sp_mat A, const arma::vec b);
   
        bool containCut(const arma::vec LHS, const double RHS, double tol = 1e-5);
   
        std::vector<short int> solEncode(const arma::vec &x) const;
   
        arma::vec zFromX(const arma::vec x);
     };
   } // namespace MathOpt
   
   namespace std {
     string to_string(Data::LCP::PolyhedraStrategy add);
   }
   
   #include "outer_lcp.h"
   #include "poly_lcp.h"
