
.. _program_listing_file_src_mathopt_mp_param_mp_param.cpp:

Program Listing for File mp_param.cpp
=====================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_mathopt_mp_param_mp_param.cpp>` (``src/mathopt/mp_param/mp_param.cpp``)

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
   
   
   #include "mathopt/mp_param/mp_param.h"
   #include <armadillo>
   #include <iostream>
   #include <memory>
   
   void MathOpt::MP_Param::write(const std::string &filename, bool) const {
     this->getQ().save(filename + "_Q.txt", arma::file_type::arma_ascii);
     this->getC().save(filename + "_C.txt", arma::file_type::arma_ascii);
     this->getA().save(filename + "_A.txt", arma::file_type::arma_ascii);
     this->getB().save(filename + "_B.txt", arma::file_type::arma_ascii);
     this->getc().save(filename + "_c.txt", arma::file_type::arma_ascii);
     this->getb().save(filename + "_b.txt", arma::file_type::arma_ascii);
   }
   
   MathOpt::MP_Param &MathOpt::MP_Param::addDummy(unsigned int pars, unsigned int vars, int position)
   {
     this->Nx += pars;
     this->Ny += vars;
     if (vars) {
        Q = Utils::resizePatch(Q, this->Ny, this->Ny);
        B = Utils::resizePatch(B, this->Ncons, this->Ny);
        c = Utils::resizePatch(c, this->Ny);
     }
     switch (position) {
     case -1:
        if (pars)
           A = Utils::resizePatch(A, this->Ncons, this->Nx);
        if (vars || pars)
           C = Utils::resizePatch(C, this->Ny, this->Nx);
        break;
     case 0:
        if (pars)
           A = arma::join_rows(arma::zeros<arma::sp_mat>(this->Ncons, pars), A);
        if (vars || pars) {
           C = Utils::resizePatch(C, this->Ny, C.n_cols);
           C = arma::join_rows(arma::zeros<arma::sp_mat>(this->Ny, pars), C);
        }
        break;
     default:
        if (pars) {
           arma::sp_mat A_temp =
                arma::join_rows(A.cols(0, position - 1), arma::zeros<arma::sp_mat>(this->Ncons, pars));
           if (static_cast<unsigned int>(position) < A.n_cols) {
             A = arma::join_rows(A_temp, A.cols(position, A.n_cols - 1));
           } else {
             A = A_temp;
           }
        }
        if (vars || pars) {
           C = Utils::resizePatch(C, this->Ny, C.n_cols);
           arma::sp_mat C_temp =
                arma::join_rows(C.cols(0, position - 1), arma::zeros<arma::sp_mat>(this->Ny, pars));
           if (static_cast<unsigned int>(position) < C.n_cols) {
             C = arma::join_rows(C_temp, C.cols(position, C.n_cols - 1));
           } else {
             C = C_temp;
           }
        }
        break;
     };
     return *this;
   }
   
   const unsigned int MathOpt::MP_Param::size()
   {
     if (Q.n_rows < 1)
        this->Ny = this->c.size();
     else
        this->Ny = this->Q.n_rows;
     this->Nx    = this->C.n_cols;
     this->Ncons = this->b.size();
     return this->Ny;
   }
   
   MathOpt::MP_Param &MathOpt::MP_Param::set(const arma::sp_mat &Q,
                                                           const arma::sp_mat &C,
                                                           const arma::sp_mat &A,
                                                           const arma::sp_mat &B,
                                                           const arma::vec &   c,
                                                           const arma::vec &   b)
   {
     this->Q = (Q);
     this->C = (C);
     this->A = (A);
     this->B = (B);
     this->c = (c);
     this->b = (b);
     if (!finalize())
        throw ZEROException(ZEROErrorCode::InvalidData, "finalize() failed");
     return *this;
   }
   
   MathOpt::MP_Param &MathOpt::MP_Param::set(arma::sp_mat &&Q,
                                                           arma::sp_mat &&C,
                                                           arma::sp_mat &&A,
                                                           arma::sp_mat &&B,
                                                           arma::vec &&   c,
                                                           arma::vec &&   b)
   {
     this->Q = std::move(Q);
     this->C = std::move(C);
     this->A = std::move(A);
     this->B = std::move(B);
     this->c = std::move(c);
     this->b = std::move(b);
     if (!finalize())
        throw ZEROException(ZEROErrorCode::InvalidData, "finalize() failed");
     return *this;
   }
   
   MathOpt::MP_Param &MathOpt::MP_Param::set(const QP_Objective &obj, const QP_Constraints &cons) {
     return this->set(obj.Q, obj.C, cons.A, cons.B, obj.c, cons.b);
   }
   
   MathOpt::MP_Param &MathOpt::MP_Param::set(QP_Objective &&obj, QP_Constraints &&cons) {
     return this->set(obj.Q, obj.C, cons.A, cons.B, obj.c, cons.b);
   }
   
   bool MathOpt::MP_Param::dataCheck(bool forceSymmetry) const
   {
     if (forceSymmetry) {
        if (!this->Q.is_symmetric())
           return false;
     }
     if (this->Q.n_cols > 0 && this->Q.n_cols != Ny) {
        return false;
     }
     if (this->A.n_cols > 0 && this->A.n_cols != Nx) {
        return false;
     }
     if (this->B.n_cols != Ny) {
        return false;
     }
     if (this->C.n_rows != Ny) {
        return false;
     }
     if (this->c.size() != Ny) {
        return false;
     }
     if (this->A.n_rows > 0 && this->A.n_rows != Ncons) {
        return false;
     }
     if (this->B.n_rows != Ncons) {
        return false;
     }
     return true;
   }
   
   bool MathOpt::MP_Param::dataCheck(const QP_Objective &  obj,
                                                const QP_Constraints &cons,
                                                bool                  checkobj,
                                                bool                  checkcons) {
     unsigned int Ny    = obj.Q.n_rows;
     unsigned int Nx    = obj.C.n_cols;
     unsigned int Ncons = cons.b.size();
     if (checkobj && obj.Q.n_cols != Ny) {
        return false;
     }
     if (checkobj && obj.C.n_rows != Ny) {
        return false;
     }
     if (checkobj && obj.c.size() != Ny) {
        return false;
     }
     if (checkcons && cons.A.n_cols != Nx) {
        return false;
     } // Rest are matrix size compatibility checks
     if (checkcons && cons.B.n_cols != Ny) {
        return false;
     }
     if (checkcons && cons.A.n_rows != Ncons) {
        return false;
     }
     if (checkcons && cons.B.n_rows != Ncons) {
        return false;
     }
     return true;
   }
