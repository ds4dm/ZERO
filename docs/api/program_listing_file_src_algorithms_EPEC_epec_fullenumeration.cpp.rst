
.. _program_listing_file_src_algorithms_EPEC_epec_fullenumeration.cpp:

Program Listing for File epec_fullenumeration.cpp
=================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_algorithms_EPEC_epec_fullenumeration.cpp>` (``src/algorithms/EPEC/epec_fullenumeration.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #include "algorithms/EPEC/epec_fullenumeration.h"
   #include "lcp/poly_lcp.h"
   #include <boost/log/trivial.hpp>
   
   void Algorithms::EPEC::FullEnumeration::solve() {
     for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i)
        this->PolyLCP.at(i)->enumerateAll(true);
     this->EPECObject->makePlayersQPs();
     BOOST_LOG_TRIVIAL(trace) << "Algorithms::EPEC::FullEnumeration::solve: "
                                           "Starting FullEnumeration search";
     this->EPECObject->computeNashEq(this->EPECObject->Stats.AlgorithmData.PureNashEquilibrium.get(),
                                                this->EPECObject->Stats.AlgorithmData.TimeLimit.get());
     if (this->isSolved()) {
        this->EPECObject->Stats.Status.set(ZEROStatus::NashEqFound);
        if (this->isPureStrategy())
           this->EPECObject->Stats.PureNashEquilibrium = true;
     }
     // Post Solving
     this->postSolving();
   }
