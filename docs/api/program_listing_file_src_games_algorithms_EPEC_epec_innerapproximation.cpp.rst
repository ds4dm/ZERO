
.. _program_listing_file_src_games_algorithms_EPEC_epec_innerapproximation.cpp:

Program Listing for File epec_innerapproximation.cpp
====================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_games_algorithms_EPEC_epec_innerapproximation.cpp>` (``src/games/algorithms/EPEC/epec_innerapproximation.cpp``)

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
   
   
   #include "games/algorithms/EPEC/epec_innerapproximation.h"
   #include "games/algorithms/EPEC/epec_combinatorialpne.h"
   #include "zero.h"
   #include <armadillo>
   #include <boost/log/trivial.hpp>
   #include <chrono>
   #include <gurobi_c++.h>
   #include <set>
   #include <string>
   
   void Algorithms::EPEC::InnerApproximation::solve() {
     this->start();
     this->postSolving();
   }
   
   void Algorithms::EPEC::InnerApproximation::start() {
     // Set the initial point for all countries as 0 and solve the respective LCPs?
     this->EPECObject->SolutionX.zeros(this->EPECObject->NumVariables);
     bool solved = {false};
     bool addRandPoly{false};
     bool infeasCheck{false};
     // When true, a MNE has been found. The algorithm now tries to find a PNE, at
     // the cost of incrementally enumerating the remaining polyhedra, up to the
     // TimeLimit (if any).
     //
     bool incrementalEnumeration{false};
   
     std::vector<arma::vec> prevDevns(this->EPECObject->NumPlayers);
     this->EPECObject->Stats.NumIterations = 0;
     if (this->EPECObject->Stats.AlgorithmData.PolyhedraStrategy.get() ==
           Data::LCP::PolyhedraStrategy::Random) {
        for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i) {
           // 42 is the answer, we all know
           long int seed = this->EPECObject->Stats.AlgorithmData.RandomSeed.get() < 0
                                     ? std::chrono::high_resolution_clock::now().time_since_epoch().count() +
                                             42 + PolyLCP.at(i)->getNumRows()
                                     : this->EPECObject->Stats.AlgorithmData.RandomSeed.get();
           PolyLCP.at(i)->AddPolyMethodSeed = seed;
        }
     }
     if (this->EPECObject->Stats.AlgorithmData.TimeLimit.get() > 0)
        this->EPECObject->InitTime = std::chrono::high_resolution_clock::now();
   
     // Stay in this loop, till you find a Nash equilibrium or prove that there
     // does not exist a Nash equilibrium or you run out of time.
     while (!solved) {
        this->EPECObject->Stats.NumIterations.set(this->EPECObject->Stats.NumIterations.get() + 1);
        BOOST_LOG_TRIVIAL(info) << "Algorithms::EPEC::InnerApproximation::solve: Iteration "
                                        << std::to_string(this->EPECObject->Stats.NumIterations.get());
   
        if (addRandPoly) {
           BOOST_LOG_TRIVIAL(info) << "Algorithms::EPEC::InnerApproximation::solve: using "
                                               "heuristical polyhedra selection";
           bool success =
                this->addRandomPoly2All(this->EPECObject->Stats.AlgorithmData.Aggressiveness.get(),
                                                this->EPECObject->Stats.NumIterations.get() == 1);
           if (!success) {
             this->EPECObject->Stats.Status.set(ZEROStatus::NashEqNotFound);
             return;
           }
        } else { // else we are in the case of finding deviations.
           unsigned int deviatedCountry{0};
           arma::vec    countryDeviation{};
           if (this->isSolved(&deviatedCountry, &countryDeviation)) {
             this->EPECObject->Stats.Status.set(ZEROStatus::NashEqFound);
             this->EPECObject->Stats.PureNashEquilibrium = this->isPureStrategy();
             if ((this->EPECObject->Stats.AlgorithmData.PureNashEquilibrium.get() &&
                    !this->EPECObject->Stats.PureNashEquilibrium.get())) {
                // We are seeking for a pure strategy. Then, here we switch between an
                // incremental
                // enumeration or combinations of pure strategies.
                if (this->EPECObject->Stats.AlgorithmData.RecoverStrategy.get() ==
                     Data::EPEC::RecoverStrategy::IncrementalEnumeration) {
                   BOOST_LOG_TRIVIAL(info) << "Algorithms::EPEC::InnerApproximation::solve: "
                                                       "triggering recover strategy "
                                                       "(IncrementalEnumeration)";
                   incrementalEnumeration = true;
                } else if (this->EPECObject->Stats.AlgorithmData.RecoverStrategy.get() ==
                               Data::EPEC::RecoverStrategy::Combinatorial) {
                   BOOST_LOG_TRIVIAL(info) << "Algorithms::EPEC::InnerApproximation::solve: triggering "
                                                       "recover strategy (Combinatorial)";
                   // In this case, we want to try all the combinations of pure
                   // strategies, except the ones between polyhedra we already tested.
                   std::vector<std::set<unsigned long int>> excludeList;
                   for (unsigned long j = 0; j < this->EPECObject->NumPlayers; ++j) {
                     excludeList.push_back(PolyLCP.at(j)->getAllPolyhedra());
                   }
                   Algorithms::EPEC::CombinatorialPNE combPNE(this->Env, this->EPECObject);
                   combPNE.solveWithExcluded(excludeList);
                   return;
                }
   
             } else
                return;
           }
           // Vector of deviations for the countries
           std::vector<arma::vec> devns = std::vector<arma::vec>(this->EPECObject->NumPlayers);
           this->getAllDeviations(devns, this->EPECObject->SolutionX, prevDevns);
           prevDevns              = devns;
           unsigned int addedPoly = this->addDeviatedPolyhedron(devns, infeasCheck);
           if (addedPoly == 0 && this->EPECObject->Stats.NumIterations.get() > 1 &&
                !incrementalEnumeration) {
             BOOST_LOG_TRIVIAL(error) << " In Algorithms::EPEC::InnerApproximation::solve: Not "
                                                   "Solved, but no deviation? Error!\n This might be due to "
                                                   "Numerical issues (tolerances)";
             this->EPECObject->Stats.Status.set(ZEROStatus::Numerical);
             solved = true;
           }
           if (infeasCheck && this->EPECObject->Stats.NumIterations.get() == 1) {
             BOOST_LOG_TRIVIAL(warning) << " In Algorithms::EPEC::InnerApproximation::solve: Problem is "
                                                     "infeasible";
             this->EPECObject->Stats.Status.set(ZEROStatus::NashEqNotFound);
             return;
           }
        }
   
        this->EPECObject->makePlayersQPs();
   
        // TimeLimit
        if (this->EPECObject->Stats.AlgorithmData.TimeLimit.get() > 0) {
           const std::chrono::duration<double> timeElapsed =
                std::chrono::high_resolution_clock::now() - this->EPECObject->InitTime;
           const double timeRemaining =
                this->EPECObject->Stats.AlgorithmData.TimeLimit.get() - timeElapsed.count();
           addRandPoly =
                !this->EPECObject->computeNashEq(
                     this->EPECObject->Stats.AlgorithmData.PureNashEquilibrium.get(), timeRemaining) &&
                !incrementalEnumeration;
        } else {
           // No Time Limit
           addRandPoly = !this->EPECObject->computeNashEq(
                                   this->EPECObject->Stats.AlgorithmData.PureNashEquilibrium.get()) &&
                             !incrementalEnumeration;
        }
        if (addRandPoly)
           this->EPECObject->Stats.AlgorithmData.LostIntermediateEq.set(
                this->EPECObject->Stats.AlgorithmData.LostIntermediateEq.get() + 1);
        for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i) {
           BOOST_LOG_TRIVIAL(info) << "Country " << i << PolyLCP.at(i)->feasabilityDetailString();
        }
        // This might be reached when a NashEq is found, and need to be verified.
        // Anyway, we are over the TimeLimit and we should stop
        if (this->EPECObject->Stats.AlgorithmData.TimeLimit.get() > 0) {
           const std::chrono::duration<double> timeElapsed =
                std::chrono::high_resolution_clock::now() - this->EPECObject->InitTime;
           const double timeRemaining =
                this->EPECObject->Stats.AlgorithmData.TimeLimit.get() - timeElapsed.count();
           if (timeRemaining <= 0) {
             if (!incrementalEnumeration)
                this->EPECObject->Stats.Status.set(ZEROStatus::TimeLimit);
             return;
           }
        }
     }
   }
   
   bool Algorithms::EPEC::InnerApproximation::addRandomPoly2All(unsigned int aggressiveLevel,
                                                                                    bool         stopOnSingleInfeasibility)
   {
     BOOST_LOG_TRIVIAL(trace) << "Adding Random polyhedra to countries";
     bool infeasible{true};
     for (unsigned int i = 0; i < this->EPECObject->NumPlayers; i++) {
        auto addedPolySet = PolyLCP.at(i)->addAPoly(
             aggressiveLevel, this->EPECObject->Stats.AlgorithmData.PolyhedraStrategy.get());
        if (stopOnSingleInfeasibility && addedPolySet.empty()) {
           BOOST_LOG_TRIVIAL(info) << "Algorithms::EPEC::InnerApproximation::addRandomPoly2All: No Nash "
                                               "equilibrium. due to "
                                               "infeasibility of country "
                                           << i;
           return false;
        }
        if (!addedPolySet.empty())
           infeasible = false;
     }
     return !infeasible;
   }
   
   bool Algorithms::EPEC::InnerApproximation::getAllDeviations(
        std::vector<arma::vec> &      deviations, 
        const arma::vec &             guessSol,   
        const std::vector<arma::vec> &prevDev //<[in] The previous vector of deviations, if any exist.
        ) const
   {
     deviations = std::vector<arma::vec>(this->EPECObject->NumPlayers);
   
     for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i) { // For each country
        // If we cannot compute a deviation, it means model is infeasible!
        if (this->EPECObject->respondSol(deviations.at(i), i, guessSol, prevDev.at(i)) == GRB_INFINITY)
           return false;
        // cout << "Game::EPEC::getAllDeviations: deviations(i): "
        // <<deviations.at(i);
     }
     return true;
   }
   
   unsigned int Algorithms::EPEC::InnerApproximation::addDeviatedPolyhedron(
        const std::vector<arma::vec> &deviations, 
        bool &infeasCheck 
        ) const {
     infeasCheck        = false;
     unsigned int added = 0;
     for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i) { // For each country
        bool ret = false;
        if (!deviations.at(i).empty())
           PolyLCP.at(i)->addPolyFromX(deviations.at(i), ret);
        if (ret) {
           BOOST_LOG_TRIVIAL(trace) << "Algorithms::EPEC::InnerApproximation::"
                                                "addDeviatedPolyhedron: added "
                                                "polyhedron for player "
                                            << i;
           ++added;
        } else {
           infeasCheck = true;
           BOOST_LOG_TRIVIAL(trace) << "Algorithms::EPEC::InnerApproximation::addDeviatedPolyhedron: NO "
                                                "polyhedron added for player "
                                            << i;
        }
     }
     return added;
   }
