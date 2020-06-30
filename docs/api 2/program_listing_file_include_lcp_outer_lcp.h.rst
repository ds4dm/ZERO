
.. _program_listing_file_include_lcp_outer_lcp.h:

Program Listing for File outer_lcp.h
====================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_lcp_outer_lcp.h>` (``include/lcp/outer_lcp.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   #include <armadillo>
   #include <gurobi_c++.h>
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   #include <zero.h>
   
   namespace Game {
     class OuterLCP : public LCP {
        // using LCP::LCP;
     public:
        OuterLCP(GRBEnv *env, const NashGame &N) : LCP(env, N) {
           this->Ai = std::unique_ptr<spmat_Vec>(new spmat_Vec());
           this->bi = std::unique_ptr<vec_Vec>(new vec_Vec());
           this->clearApproximation();
        };
   
        void clearApproximation() {
           this->Ai->clear();
           this->bi->clear();
           this->Approximation.clear();
           this->feasApprox = false;
        }
   
        bool checkComponentFeas(const std::vector<short int> &encoding);
   
        void outerApproximate(std::vector<bool> encoding, bool clear = true);
   
        bool              addComponent(std::vector<short int> encoding,
                                                 bool                   checkFeas,
                                                 bool                   custom = false,
                                                 spmat_Vec *            custAi = {},
                                                 vec_Vec *              custbi = {});
        inline const bool getFeasApprox() { return this->feasApprox; }
   
     private:
        std::set<unsigned long int> Approximation =
             {}; 
        std::set<unsigned long int> FeasibleComponents =
             {}; 
        std::set<unsigned long int> InfeasibleComponents =
             {}; 
        bool isParent(const std::vector<short> &father, const std::vector<short> &child);
   
        void addChildComponents(const std::vector<short> encoding);
   
        bool feasApprox = false;
     };
   } // namespace Game
