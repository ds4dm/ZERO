
.. _program_listing_file_include_zero.h:

Program Listing for File zero.h
===============================

|exhale_lsh| :ref:`Return to documentation for file <file_include_zero.h>` (``include/zero.h``)

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
   
   #include <armadillo>
   #include <iostream>
   #include <map>
   #include <memory>
   #include <string>
   #include <vector>
   
   #include "support/codes.h"
   #include "support/utils.h"
   #include "support/version.h"
   
   using perps = std::vector<std::pair<unsigned int, unsigned int>>;
   std::ostream &                            operator<<(std::ostream &ost, perps C);
   template <class T> std::ostream &         operator<<(std::ostream &ost, std::vector<T> v);
   template <class T, class S> std::ostream &operator<<(std::ostream &ost, std::pair<T, S> p);
   using spmat_Vec = std::vector<std::unique_ptr<arma::sp_mat>>;
   using vec_Vec   = std::vector<std::unique_ptr<arma::vec>>;
   
   // Forward declarations
   class ZEROException;
   enum class ZEROErrorCode;
   
   namespace MathOpt {
     struct QP_Objective;
     struct QP_Constraints;
     class MP_Param;
     class QP_Param;
     class IP_Param;
     class LCP;
   } // namespace MathOpt
   
   namespace Game {
     class NashGame;
     class PolyLCP;
     class OuterLCP;
     class EPEC;
     class IPG;
   } // namespace Game
   namespace Algorithms {
     class AbstractAlgorithm; 
     namespace EPEC {
        class Algorithm;
        class PolyBase;
        class FullEnumeration;
        class InnerApproximation;
        class CombinatorialPNE;
        class OuterApproximation;
     } // namespace EPEC
     namespace IPG {
        class Oracle;
     } // namespace IPG
   } // namespace Algorithms
   
   class ZEROAlgorithmData;
   namespace Data {
     namespace EPEC {
        class DataObject;
     }
     namespace IPG {
        class DataObject;
     }
     namespace LCP {
        enum class PolyhedraStrategy;
     }
   
   } // namespace Data
   #include "games/games.h"
   #include "mathopt/mathopt.h"
