
.. _program_listing_file_include_algorithms_fullenumeration.h:

Program Listing for File fullenumeration.h
==========================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_algorithms_fullenumeration.h>` (``include/algorithms/fullenumeration.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include "algorithms/algorithms.h"
   #include "algorithms/polybase.h"
   
   namespace Algorithms {
   
   class FullEnumeration : public PolyBase {
   public:
     FullEnumeration(GRBEnv *env, Game::EPEC *EPECObject)
         : PolyBase(env, EPECObject){};
     void solve();
   };
   } // namespace Algorithms
