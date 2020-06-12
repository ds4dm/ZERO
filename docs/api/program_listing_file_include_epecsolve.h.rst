
.. _program_listing_file_include_epecsolve.h:

Program Listing for File epecsolve.h
====================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_epecsolve.h>` (``include/epecsolve.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   #include <armadillo>
   #include <iostream>
   #include <map>
   #include <memory>
   #include <vector>
   
   using perps = std::vector<std::pair<unsigned int, unsigned int>>;
   std::ostream &operator<<(std::ostream &ost, perps C);
   template <class T>
   std::ostream &operator<<(std::ostream &ost, std::vector<T> v);
   template <class T, class S>
   std::ostream &operator<<(std::ostream &ost, std::pair<T, S> p);
   using spmat_Vec = std::vector<std::unique_ptr<arma::sp_mat>>;
   using vec_Vec = std::vector<std::unique_ptr<arma::vec>>;
   
   // Forward declarations
   namespace Game {
   struct QP_Objective;
   struct QP_Constraints;
   class MP_Param;
   class QP_Param;
   class NashGame;
   class LCP;
   class PolyLCP;
   class OuterLCP;
   class EPEC;
   enum class EPECAddPolyMethod {
     Sequential,        
     ReverseSequential, 
     Random 
   };
   } // namespace Game
   namespace Algorithms {
   // Forward declarations
   class Algorithm;
   class PolyBase;
   class FullEnumeration;
   class InnerApproximation;
   class CombinatorialPNE;
   class OuterApproximation;
   } // namespace Algorithms
   #include "games.h"
   #include "lcp/lcp.h"
   #include "utils.h"
   #include "version.h"
