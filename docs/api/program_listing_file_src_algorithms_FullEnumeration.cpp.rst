
.. _program_listing_file_src_algorithms_FullEnumeration.cpp:

Program Listing for File FullEnumeration.cpp
============================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_algorithms_FullEnumeration.cpp>` (``src/algorithms/FullEnumeration.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #include "algorithms/fullenumeration.h"
   #include "lcp/polylcp.h"
   #include <boost/log/trivial.hpp>
   
   void Algorithms::FullEnumeration::solve() {
     for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i)
       this->PolyLCP.at(i)->enumerateAll(true);
     this->EPECObject->makePlayersQPs();
     BOOST_LOG_TRIVIAL(trace)
         << "Algorithms::FullEnumeration::solve: Starting FullEnumeration search";
     this->EPECObject->computeNashEq(
         this->EPECObject->Stats.AlgorithmParam.PureNashEquilibrium,
         this->EPECObject->Stats.AlgorithmParam.TimeLimit);
     if (this->isSolved()) {
       this->EPECObject->Stats.Status = Game::EPECsolveStatus::NashEqFound;
       if (this->isPureStrategy())
         this->EPECObject->Stats.PureNashEquilibrium = true;
     }
     // Post Solving
     this->postSolving();
   }
