
.. _program_listing_file_include_mathopt_mathopt.h:

Program Listing for File mathopt.h
==================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_mathopt_mathopt.h>` (``include/mathopt/mathopt.h``)

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
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   
   namespace MathOpt {
   
     typedef struct QP_Objective {
        arma::sp_mat Q;
        arma::sp_mat C;
        arma::vec    c;
     } QP_objective;
     typedef struct QP_Constraints {
        arma::sp_mat A, B;
        arma::vec    b;
     } QP_constraints;
   
   
     arma::vec LPSolve(const arma::sp_mat &A,
                             const arma::vec &   b,
                             const arma::vec &   c,
                             int &               status,
                             bool                positivity = false);
   
     unsigned int convexHull(const std::vector<arma::sp_mat *> *Ai,
                                     const std::vector<arma::vec *> *   bi,
                                     arma::sp_mat &                     A,
                                     arma::vec &                        b,
                                     arma::sp_mat                       Acom = {},
                                     arma::vec                          bcom = {});
   
     void compConvSize(arma::sp_mat &                     A,
                             unsigned int                       nFinCons,
                             unsigned int                       nFinVar,
                             const std::vector<arma::sp_mat *> *Ai,
                             const std::vector<arma::vec *> *   bi,
                             const arma::sp_mat &               Acom,
                             const arma::vec &                  bcom);
   
     void print(const perps &C) noexcept;
   } // namespace MathOpt
   
   #include "lcp/lcp.h"
   #include "mathopt/mp_param/mp_param.h"
