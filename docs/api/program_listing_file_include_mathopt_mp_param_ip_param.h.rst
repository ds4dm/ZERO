
.. _program_listing_file_include_mathopt_mp_param_ip_param.h:

Program Listing for File ip_param.h
===================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_mathopt_mp_param_ip_param.h>` (``include/mathopt/mp_param/ip_param.h``)

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
   #include "zero.h"
   #include <armadillo>
   #include <gurobi_c++.h>
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   
   namespace MathOpt {
   
   
     std::ostream &operator<<(std::ostream &os, const IP_Param &I);
   
     class IP_Param : public MP_Param
     // Shape of C is Ny\times Nx
     {
     private:
        // Gurobi environment and model
        GRBEnv *  Env;
        GRBModel  IPModel;          
        arma::vec bounds;           
        arma::vec integers;         
        bool      madeModel{false}; 
   
        // These methods should be inaccessible to the inheritor, since we have a
        // different structure.
        using MP_Param::set;
   
     public: // Constructors
        explicit IP_Param(GRBEnv *env = nullptr) : Env{env}, IPModel{(*env)} { this->size(); }
   
        explicit IP_Param(arma::sp_mat C,
                                arma::sp_mat B,
                                arma::vec    b,
                                arma::vec    c,
                                arma::vec    bounds,
                                arma::vec    integers,
                                GRBEnv *     env = nullptr)
             : Env{env}, IPModel{(*env)} {
           this->Q.zeros(0);
           this->A.zeros(0);
           this->set(Q, C, A, B, c, b);
           this->bounds   = bounds;
           this->integers = integers;
           this->size();
           this->forceDataCheck();
        };
   
        arma::vec getIntegers() const { return this->integers; }
   
        arma::vec getBounds() const { return this->bounds; }
   
        void makeModel();
   
        void addConstraints(const arma::sp_mat A, const arma::vec b);
   
        IP_Param(const IP_Param &ipg)
             : MP_Param(ipg), Env{ipg.Env}, IPModel{ipg.IPModel}, madeModel{ipg.madeModel} {
           this->size();
        };
   
        // Override setters
        IP_Param &set(const arma::sp_mat &C,
                           const arma::sp_mat &B,
                           const arma::vec &   b,
                           const arma::vec &   c,
                           const arma::vec &   bounds,
                           const arma::vec &   integers); // Copy data into this
        IP_Param &set(arma::sp_mat & C,
                           arma::sp_mat &&B,
                           arma::vec &&   b,
                           arma::vec &&   c,
                           arma::vec &&   bounds,
                           arma::vec &&   integers); // Copy data into this
   
        IP_Param &set(const QP_Objective &  obj,
                           const QP_Constraints &cons,
                           const arma::vec &     bounds   = {},
                           const arma::vec &     integers = {});
   
        IP_Param &set(QP_Objective &&  obj,
                           QP_Constraints &&cons,
                           arma::vec &&     bounds   = {},
                           arma::vec &&     integers = {});
   
        bool operator==(const IP_Param &IPG2) const;
   
        std::unique_ptr<GRBModel> solveFixed(arma::vec x, bool solve);
   
        double computeObjective(const arma::vec &y,
                                        const arma::vec &x,
                                        bool             checkFeas = true,
                                        double           tol       = 1e-6) const;
   
        inline bool isPlayable(const IP_Param &P) const
        {
           bool b1, b2, b3;
           b1 = (this->Nx + this->Ny) == (P.getNx() + P.getNy());
           b2 = this->Nx >= P.getNy();
           b3 = this->Ny <= P.getNx();
           return b1 && b2 && b3;
        }
   
        IP_Param &addDummy(unsigned int pars, unsigned int vars = 0, int position = -1) override;
   
        void write(const std::string &filename, bool append) const override;
        long load(const std::string &filename, long pos);
   
        double computeObjectiveWithoutOthers(const arma::vec &y) const;
   
        arma::vec getConstraintViolations(const arma::vec y, double tol);
   
        void forceDataCheck();
     };
   } // namespace MathOpt
