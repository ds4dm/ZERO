
.. _program_listing_file_include_mathopt_mp_param_mp_param.h:

Program Listing for File mp_param.h
===================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_mathopt_mp_param_mp_param.h>` (``include/mathopt/mp_param/mp_param.h``)

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
   
   namespace MathOpt {
     class MP_Param {
     protected:
        // Data representing the parameterized QP
        arma::sp_mat Q, A, B, C;
        arma::vec    c, b;
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
        arma::sp_mat getQ() const { return this->Q; }   
        arma::sp_mat getC() const { return this->C; }   
        arma::sp_mat getA() const { return this->A; }   
        arma::sp_mat getB() const { return this->B; }   
        arma::vec    getc() const { return this->c; }   
        arma::vec    getb() const { return this->b; }   
        unsigned int getNx() const { return this->Nx; } 
        unsigned int getNy() const { return this->Ny; } 
   
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
        virtual MP_Param &set(const arma::sp_mat &Q,
                                     const arma::sp_mat &C,
                                     const arma::sp_mat &A,
                                     const arma::sp_mat &B,
                                     const arma::vec &   c,
                                     const arma::vec &   b); // Copy data into this
        virtual MP_Param &set(arma::sp_mat &&Q,
                                     arma::sp_mat &&C,
                                     arma::sp_mat &&A,
                                     arma::sp_mat &&B,
                                     arma::vec &&   c,
                                     arma::vec &&   b); // Move data into this
        virtual MP_Param &set(const QP_Objective &obj, const QP_Constraints &cons);
   
        virtual MP_Param &set(QP_Objective &&obj, QP_Constraints &&cons);
   
        virtual MP_Param &addDummy(unsigned int pars, unsigned int vars = 0, int position = -1);
   
        virtual void write(const std::string &filename, bool append = true) const;
   
        static bool dataCheck(const QP_Objective &  obj,
                                     const QP_Constraints &cons,
                                     bool                  checkObj  = true,
                                     bool                  checkCons = true);
     };
   
   } // namespace MathOpt
   
   #include "ip_param.h"
   #include "qp_param.h"
