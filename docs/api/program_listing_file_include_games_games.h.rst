
.. _program_listing_file_include_games_games.h:

Program Listing for File games.h
================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_games_games.h>` (``include/games/games.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #include "zero.h"
   #include <armadillo>
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   
   namespace Game {
   
   class PolyLCP; // Forward declaration
   
   arma::vec LPSolve(const arma::sp_mat &A, const arma::vec &b, const arma::vec &c,
                     int &status, bool positivity = false);
   
   unsigned int convexHull(const std::vector<arma::sp_mat *> *Ai,
                           const std::vector<arma::vec *> *bi, arma::sp_mat &A,
                           arma::vec &b, arma::sp_mat Acom = {},
                           arma::vec bcom = {});
   
   void compConvSize(arma::sp_mat &A, unsigned int nFinCons, unsigned int nFinVar,
                     const std::vector<arma::sp_mat *> *Ai,
                     const std::vector<arma::vec *> *bi, const arma::sp_mat &Acom,
                     const arma::vec &bcom);
   
   bool isZero(arma::mat M, double tol = 1e-6) noexcept;
   
   bool isZero(arma::sp_mat M, double tol = 1e-6) noexcept;
   
   // bool isZero(arma::vec M, double Tolerance = 1e-6);
   typedef struct QP_Objective {
     arma::sp_mat Q;
     arma::sp_mat C;
     arma::vec c;
   } QP_objective;
   typedef struct QP_Constraints {
     arma::sp_mat A, B;
     arma::vec b;
   } QP_constraints;
   
   } // namespace Game
