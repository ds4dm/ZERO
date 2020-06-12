
.. _program_listing_file_src_lcp_PolyLCP.cpp:

Program Listing for File PolyLCP.cpp
====================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_lcp_PolyLCP.cpp>` (``src/lcp/PolyLCP.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #include "lcp/polylcp.h"
   #include <algorithm>
   #include <armadillo>
   #include <boost/log/trivial.hpp>
   #include <gurobi_c++.h>
   #include <iostream>
   #include <memory>
   #include <random>
   #include <set>
   #include <string>
   
   bool operator==(std::vector<short int> encoding1,
                   std::vector<short int> encoding2)
   {
     if (encoding1.size() != encoding2.size())
       return false;
     for (unsigned int i = 0; i < encoding1.size(); i++) {
       if (encoding1.at(i) != encoding2.at(i))
         return false;
     }
     return true;
   }
   
   bool operator<(std::vector<short int> encoding1,
                  std::vector<short int> encoding2)
   {
     if (encoding1.size() != encoding2.size())
       return false;
     for (unsigned int i = 0; i < encoding1.size(); i++) {
       if (encoding1.at(i) != encoding2.at(i) &&
           encoding1.at(i) * encoding2.at(i) != 0) {
         return false; // encoding1 is not a child of encoding2
       }
     }
     return true; // encoding1 is a child of encoding2
   }
   
   bool operator>(std::vector<int> encoding1, std::vector<int> encoding2) {
     return (encoding2 < encoding1);
   }
   
   Game::PolyLCP &Game::PolyLCP::addPolyFromX(const arma::vec &x, bool &ret)
   {
     const auto numCompl = this->Compl.size();
     auto encoding = this->solEncode(x);
     std::stringstream encStr;
     for (auto vv : encoding)
       encStr << vv << " ";
     BOOST_LOG_TRIVIAL(trace)
         << "Game::PolyLCP::addPolyFromX: Handling deviation with encoding: "
         << encStr.str() << '\n';
     // Check if the encoding polyhedron is already in this->AllPolyhedra
     for (const auto &i : AllPolyhedra) {
       std::vector<short int> bin = Utils::numToVec(i, numCompl);
       if (encoding < bin) {
         BOOST_LOG_TRIVIAL(trace) << "Game::PolyLCP::addPolyFromX: Encoding " << i
                                  << " already in All Polyhedra! ";
         ret = false;
         return *this;
       }
     }
   
     BOOST_LOG_TRIVIAL(trace)
         << "Game::PolyLCP::addPolyFromX: New encoding not in All Polyhedra! ";
     // If it is not in AllPolyhedra
     // First change any zero indices of encoding to 1
     for (short &i : encoding) {
       if (i == 0)
         ++i;
     }
     // And then add the relevant polyhedron
     ret = this->addPolyFromEncoding(encoding, false);
     // ret = true;
     return *this;
   }
   
   bool Game::PolyLCP::addPolyFromEncoding(
       const std::vector<short int>
           encoding, 
       bool checkFeas, 
       bool custom, 
       spmat_Vec *custAi, 
       vec_Vec *custbi 
       )
   {
     unsigned int encodingNumber = Utils::vecToNum(encoding);
     BOOST_LOG_TRIVIAL(trace)
         << "Game::PolyLCP::addPolyFromEncoding: Working on polyhedron #"
         << encodingNumber;
   
     bool eval = false;
     if (checkFeas)
       eval = this->checkPolyFeas(encoding);
     else
       eval = true;
   
     if (eval) {
       if (!custom && !AllPolyhedra.empty()) {
         if (AllPolyhedra.find(encodingNumber) != AllPolyhedra.end()) {
           BOOST_LOG_TRIVIAL(trace) << "Game::PolyLCP::addPolyFromEncoding: "
                                       "Previously added polyhedron #"
                                    << encodingNumber;
           return false;
         }
       }
       std::unique_ptr<arma::sp_mat> Aii =
           std::unique_ptr<arma::sp_mat>(new arma::sp_mat(nR, nC));
       Aii->zeros();
       std::unique_ptr<arma::vec> bii =
           std::unique_ptr<arma::vec>(new arma::vec(nR, arma::fill::zeros));
       for (unsigned int i = 0; i < this->nR; i++) {
         if (encoding.at(i) == 0) {
           throw("Error in Game::PolyLCP::addPolyFromEncoding. 0s not allowed in "
                 "argument vector");
         }
         if (encoding.at(i) == 1) // Equation to be fixed top zero
         {
           for (auto j = this->M.begin_row(i); j != this->M.end_row(i); ++j)
             if (!this->isZero((*j)))
               Aii->at(i, j.col()) =
                   (*j); // Only mess with non-zero elements of a sparse matrix!
           bii->at(i) = -this->q(i);
         } else // Variable to be fixed to zero, i.e. x(j) <= 0 constraint to be
                // added
         {
           unsigned int variablePosition =
               (i >= this->LeadStart) ? i + this->NumberLeader : i;
           Aii->at(i, variablePosition) = 1;
           bii->at(i) = 0;
         }
       }
       if (custom) {
         custAi->push_back(std::move(Aii));
         custbi->push_back(std::move(bii));
       } else {
         AllPolyhedra.insert(encodingNumber);
         this->Ai->push_back(std::move(Aii));
         this->bi->push_back(std::move(bii));
       }
       return true; // Successfully added
     }
     BOOST_LOG_TRIVIAL(trace) << "Game::PolyLCP::addPolyFromEncoding: Checkfeas + "
                                 "Infeasible polyhedron #"
                              << encodingNumber;
     return false;
   }
   
   Game::PolyLCP &Game::PolyLCP::addPoliesFromEncoding(
       const std::vector<short int>
           encoding, 
       bool checkFeas, 
       bool custom, 
       spmat_Vec *custAi, 
       vec_Vec *custbi 
       )
   {
     bool flag = false; // flag that there may be multiple polyhedra, i.e. 0 in
     // some encoding entry
     std::vector<short int> encodingCopy(encoding);
     unsigned int i = 0;
     for (i = 0; i < this->nR; i++) {
       if (encoding.at(i) == 0) {
         flag = true;
         break;
       }
     }
     if (flag) {
       encodingCopy[i] = 1;
       this->addPoliesFromEncoding(encodingCopy, checkFeas, custom, custAi,
                                   custbi);
       encodingCopy[i] = -1;
       this->addPoliesFromEncoding(encodingCopy, checkFeas, custom, custAi,
                                   custbi);
     } else
       this->addPolyFromEncoding(encoding, checkFeas, custom, custAi, custbi);
     return *this;
   }
   
   unsigned long int Game::PolyLCP::getNextPoly(
       Game::EPECAddPolyMethod
           method 
   ) {
     switch (method) {
     case Game::EPECAddPolyMethod::Sequential: {
       while (this->SequentialPolyCounter < this->MaxTheoreticalPoly) {
         const auto isAll =
             AllPolyhedra.find(this->SequentialPolyCounter) != AllPolyhedra.end();
         const auto isInfeas = InfeasiblePoly.find(this->SequentialPolyCounter) !=
                               InfeasiblePoly.end();
         this->SequentialPolyCounter++;
         if (!isAll && !isInfeas) {
           return this->SequentialPolyCounter - 1;
         }
       }
       return this->MaxTheoreticalPoly;
     } break;
     case Game::EPECAddPolyMethod::ReverseSequential: {
       while (this->ReverseSequentialPolyCounter >= 0) {
         const auto isAll =
             AllPolyhedra.find(this->ReverseSequentialPolyCounter) !=
             AllPolyhedra.end();
         const auto isInfeas =
             InfeasiblePoly.find(this->ReverseSequentialPolyCounter) !=
             InfeasiblePoly.end();
         this->ReverseSequentialPolyCounter--;
         if (!isAll && !isInfeas) {
           return this->ReverseSequentialPolyCounter + 1;
         }
       }
       return this->MaxTheoreticalPoly;
     } break;
     case Game::EPECAddPolyMethod::Random: {
       static std::mt19937 engine{this->AddPolyMethodSeed};
       std::uniform_int_distribution<unsigned long int> dist(
           0, this->MaxTheoreticalPoly - 1);
       if ((InfeasiblePoly.size() + AllPolyhedra.size()) ==
           this->MaxTheoreticalPoly)
         return this->MaxTheoreticalPoly;
       while (true) {
         auto randomPolyId = dist(engine);
         const auto isAll = AllPolyhedra.find(randomPolyId) != AllPolyhedra.end();
         const auto isInfeas =
             InfeasiblePoly.find(randomPolyId) != InfeasiblePoly.end();
         if (!isAll && !isInfeas)
           return randomPolyId;
       }
     }
     }
   }
   
   std::set<std::vector<short int>>
   Game::PolyLCP::addAPoly(unsigned long int nPoly, Game::EPECAddPolyMethod method,
                           std::set<std::vector<short int>> polyhedra) {
     // We already have polyhedra AllPolyhedra and in
     // InfeasiblePoly, that are known to be infeasible.
     // Effective maximum of number of polyhedra that can be added
     // at most
     const auto numCompl = this->Compl.size();
   
     if (this->MaxTheoreticalPoly <
         nPoly) {                 // If you cannot add that numVariablesY polyhedra
       BOOST_LOG_TRIVIAL(warning) // Then issue a warning
           << "Warning in Game::PolyLCP::randomPoly: "
           << "Cannot add " << nPoly << " polyhedra. Promising a maximum of "
           << this->MaxTheoreticalPoly;
       nPoly = this->MaxTheoreticalPoly; // and update maximum possibly addable
     }
   
     if (nPoly == 0) // If nothing to be added, then nothing to be done
       return polyhedra;
   
     if (nPoly < 0) // There is no way that this can happen!
     {
       BOOST_LOG_TRIVIAL(error) << "nPoly can't be negative, i.e., " << nPoly;
       throw("Error in Game::PolyLCP::addAPoly: nPoly reached a negative value!");
     }
   
     bool complete{false};
     while (!complete) {
       auto choiceDecimal = this->getNextPoly(method);
       if (choiceDecimal >= this->MaxTheoreticalPoly)
         return polyhedra;
   
       const std::vector<short int> choice =
           Utils::numToVec(choiceDecimal, numCompl);
       auto added = this->addPolyFromEncoding(choice, true);
       if (added) // If choice is added to All Polyhedra
       {
         polyhedra.insert(choice); // Add it to set of added polyhedra
         if (polyhedra.size() == nPoly) {
           return polyhedra;
         }
       }
     }
     return polyhedra;
   }
   bool Game::PolyLCP::addThePoly(const unsigned long int &decimalEncoding) {
     if (this->MaxTheoreticalPoly < decimalEncoding) {
       // This polyhedron does not exist
       BOOST_LOG_TRIVIAL(warning)
           << "Warning in Game::PolyLCP::addThePoly: Cannot add "
           << decimalEncoding << " polyhedra, since it does not exist!";
       return false;
     }
     const unsigned int numCompl = this->Compl.size();
     const std::vector<short int> choice =
         Utils::numToVec(decimalEncoding, numCompl);
     return this->addPolyFromEncoding(choice, true);
   }
   
   Game::PolyLCP &Game::PolyLCP::enumerateAll(
       const bool
           solveLP 
       )
   {
     std::vector<short int> encoding = std::vector<short int>(nR, 0);
     this->Ai->clear();
     this->bi->clear();
     this->addPoliesFromEncoding(encoding, solveLP);
     if (this->Ai->empty()) {
       BOOST_LOG_TRIVIAL(warning)
           << "Empty vector of polyhedra given! Problem might be infeasible."
           << '\n';
       // 0 <= -1 for infeasability
       std::unique_ptr<arma::sp_mat> A(new arma::sp_mat(1, this->M.n_cols));
       std::unique_ptr<arma::vec> b(new arma::vec(1));
       b->at(0) = -1;
       this->Ai->push_back(std::move(A));
       this->bi->push_back(std::move(b));
     }
     return *this;
   }
   
   std::string Game::PolyLCP::feasabilityDetailString() const {
     std::stringstream ss;
     ss << "\tProven feasible: ";
     for (auto vv : this->AllPolyhedra)
       ss << vv << ' ';
     // ss << "\tProven infeasible: ";
     // for (auto vv : this->InfeasiblePoly)
     // ss << vv << ' ';
   
     return ss.str();
   }
   
   unsigned long Game::PolyLCP::convNumPoly() const {
     return this->AllPolyhedra.size();
   }
   
   unsigned int Game::PolyLCP::convPolyPosition(const unsigned long int i) const {
     const unsigned int nPoly = this->convNumPoly();
     if (i > nPoly) {
       BOOST_LOG_TRIVIAL(error) << "Error in Game::PolyLCP::convPolyPosition: "
                                   "Invalid argument. Out of bounds for i";
       throw("Error in Game::PolyLCP::convPolyPosition: Invalid "
             "argument. Out of bounds for i");
     }
     const unsigned int nC = this->M.n_cols;
     return nC + i * nC;
   }
   
   unsigned int Game::PolyLCP::convPolyWeight(const unsigned long int i) const {
     const unsigned int nPoly = this->convNumPoly();
     if (nPoly <= 1) {
       return 0;
     }
     if (i > nPoly) {
       throw("Error in Game::PolyLCP::convPolyWeight: "
             "Invalid argument. Out of bounds for i");
     }
     const unsigned int nC = this->M.n_cols;
   
     return nC + nPoly * nC + i;
   }
   
   bool Game::PolyLCP::checkPolyFeas(
       const unsigned long int
           &decimalEncoding 
   ) {
     return this->checkPolyFeas(
         Utils::numToVec(decimalEncoding, this->Compl.size()));
   }
   
   bool Game::PolyLCP::checkPolyFeas(
       const std::vector<short int>
           &encoding 
   ) {
     unsigned long int encodingNumber = Utils::vecToNum(encoding);
   
     if (InfeasiblePoly.find(encodingNumber) != InfeasiblePoly.end()) {
       BOOST_LOG_TRIVIAL(trace)
           << "Game::PolyLCP::checkPolyFeas: Previously known "
              "infeasible polyhedron. "
           << encodingNumber;
       return false;
     }
   
     if (FeasiblePoly.find(encodingNumber) != FeasiblePoly.end()) {
       BOOST_LOG_TRIVIAL(trace)
           << "Game::PolyLCP::checkPolyFeas: Previously known "
              "feasible polyhedron."
           << encodingNumber;
       return true;
     }
   
     unsigned int count{0};
     try {
       makeRelaxed();
       GRBModel model(this->RlxdModel);
       for (auto i : encoding) {
         if (i > 0)
           model.getVarByName("z_" + std::to_string(count))
               .set(GRB_DoubleAttr_UB, 0);
         if (i < 0)
           model
               .getVarByName("x_" + std::to_string(count >= this->LeadStart
                                                       ? count + NumberLeader
                                                       : count))
               .set(GRB_DoubleAttr_UB, 0);
         count++;
       }
       model.set(GRB_IntParam_OutputFlag, 0);
       model.optimize();
       if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
         FeasiblePoly.insert(encodingNumber);
         return true;
       } else {
         BOOST_LOG_TRIVIAL(trace)
             << "Game::PolyLCP::checkPolyFeas: Detected infeasibility of "
             << encodingNumber << " (GRB_STATUS=" << model.get(GRB_IntAttr_Status)
             << ")";
         InfeasiblePoly.insert(encodingNumber);
         return false;
       }
     } catch (const char *e) {
       std::cerr << "Error in Game::PolyLCP::checkPolyFeas: " << e << '\n';
       throw;
     } catch (std::string e) {
       std::cerr << "String: Error in Game::PolyLCP::checkPolyFeas: " << e << '\n';
       throw;
     } catch (std::exception &e) {
       std::cerr << "Exception: Error in Game::PolyLCP::checkPolyFeas: "
                 << e.what() << '\n';
       throw;
     } catch (GRBException &e) {
       std::cerr << "GRBException: Error in Game::PolyLCP::checkPolyFeas: "
                 << e.getErrorCode() << ": " << e.getMessage() << '\n';
       throw;
     }
     return false;
   }
