
.. _program_listing_file_include_algorithms_EPEC_epec_fullenumeration.h:

Program Listing for File epec_fullenumeration.h
===============================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_algorithms_EPEC_epec_fullenumeration.h>` (``include/algorithms/EPEC/epec_fullenumeration.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include "epec_algorithms.h"
   #include "epec_polybase.h"
   
   namespace Algorithms {
     namespace EPEC {
        class FullEnumeration : public PolyBase {
        public:
           FullEnumeration(GRBEnv *env, Game::EPEC *EPECObject) : PolyBase(env, EPECObject){};
           void solve();
        };
     } // namespace EPEC
   } // namespace Algorithms
