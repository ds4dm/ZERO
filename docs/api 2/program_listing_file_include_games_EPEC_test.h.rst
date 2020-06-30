
.. _program_listing_file_include_games_EPEC_test.h:

Program Listing for File EPEC_test.h
====================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_games_EPEC_test.h>` (``include/games/EPEC_test.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   #include "lcp/lcp.h"
   #include "zero.h"
   #include <armadillo>
   #include <gurobi_c++.h>
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   #include <support/codes.h>
   
   namespace Data {
     namespace EPEC {
   
        enum class Algorithms {
           FullEnumeration, 
           InnerApproximation, 
           CombinatorialPne, 
           OuterApproximation 
        };
   
        enum class RecoverStrategy {
           IncrementalEnumeration, 
           Combinatorial           
        };
   
        class DataObject : public ZEROAlgorithmData {
        public:
           Attr<Data::EPEC::Algorithms> Algorithm = {
                Data::EPEC::Algorithms::FullEnumeration}; 
           Attr<Data::EPEC::RecoverStrategy> RecoverStrategy = {
                Data::EPEC::RecoverStrategy::IncrementalEnumeration}; 
           Attr<Data::LCP::PolyhedraStrategy>
                                    PolyhedraStrategy; 
           Attr<unsigned int> Aggressiveness{
                1}; 
           Attr<std::vector<unsigned int>> FeasiblePolyhedra =
                std::vector<unsigned int>(); 
           Attr<bool> BoundPrimals{false};  
           Attr<double> BoundBigM{1e5};           
           Attr<int>    LostIntermediateEq = {0}; 
           DataObject() : PolyhedraStrategy{static_cast<LCP::PolyhedraStrategy>(0)} {};
        };
   
     } // namespace EPEC
   } // namespace Data
   
   namespace Game {
   
     class EPEC {
     private:
        std::vector<unsigned int>  SizesWithoutHull{};
        std::unique_ptr<Game::LCP> TheLCP;       
        std::unique_ptr<GRBModel>  LCPModel;     
        std::unique_ptr<GRBModel>  LCPModelBase; 
        unsigned int                                 NumVariables{0};
        unsigned int                                 NumPlayers{0};
        std::shared_ptr<Algorithms::EPEC::Algorithm> Algorithm{};
   
     protected: // Datafields
        std::vector<std::shared_ptr<Game::NashGame>> PlayersLowerLevels{};
        std::vector<std::shared_ptr<Game::LCP>>      PlayersLCP{};
   
        std::vector<std::shared_ptr<Game::QP_Param>>
             PlayersQP{}; 
        std::vector<std::shared_ptr<Game::QP_Objective>>
             LeaderObjective{}; 
        std::vector<std::shared_ptr<Game::QP_Objective>>
             LeaderObjectiveConvexHull{}; 
   
        std::unique_ptr<Game::NashGame> TheNashGame; 
   
        std::vector<unsigned int> LeaderLocations{}; 
        std::vector<const unsigned int *> LocEnds{};
        std::vector<unsigned int>         ConvexHullVariables{};
        unsigned int                      numMCVariables{0};
   
        GRBEnv *Env;
        bool    Finalized{false};
        bool NashEquilibrium{false}; 
        std::chrono::high_resolution_clock::time_point InitTime;
        ZEROStatistics<Data::EPEC::DataObject>         Stats = ZEROStatistics<Data::EPEC::DataObject>(
           Data::EPEC::DataObject()); 
        arma::vec SolutionZ,                   
             SolutionX;                         
        bool warmstart(arma::vec x);           
   
     private:
        void       addDummyLead(unsigned int i); 
        const void makePlayerQP(unsigned int i);
   
        void makePlayersQPs();
   
        void makeTheLCP();
   
        void computeLeaderLocations(unsigned int addSpaceForMC = 0);
   
        void getXMinusI(const arma::vec &x, const unsigned int &i, arma::vec &solOther) const;
   
        bool computeNashEq(bool pureNE = false, double localTimeLimit = -1.0, bool check = false);
   
     protected:                                  // functions
        explicit EPEC(GRBEnv *env) : Env{env} {}; 
   
        // virtual function to be implemented by the inheritor.
        virtual void makeObjectivePlayer(const unsigned int i, Game::QP_Objective &QP_obj) = 0;
   
        // virtual function to be optionally implemented by the inheritor.
        virtual void preFinalize();
   
        virtual void postFinalize();
   
        virtual void updateLocations() = 0; // If any location tracking system is implemented, that
        // can be called from in here.
        virtual void makeMCConstraints(arma::sp_mat &MC, arma::vec &RHS) const {
           MC.zeros();
           RHS.zeros();
        };
   
     public: // functions
        // Friends algorithmic classes
        friend class Algorithms::EPEC::PolyBase;
   
        friend class Algorithms::EPEC::InnerApproximation;
   
        friend class Algorithms::EPEC::OuterApproximation;
   
        friend class Algorithms::EPEC::CombinatorialPNE;
   
        friend class Algorithms::EPEC::FullEnumeration;
   
        EPEC()       = delete;  // No default constructor
        EPEC(EPEC &) = delete;  // Abstract class - no copy constructor
        ~EPEC()      = default; // Destructor to free data
   
        void finalize();
   
        const void findNashEq();
        bool       isSolved(double tol = 1e-5) const;
   
        std::unique_ptr<GRBModel> respond(const unsigned int i, const arma::vec &x) const;
   
        double respondSol(arma::vec &      sol,
                                unsigned int     player,
                                const arma::vec &x,
                                const arma::vec &prevDev = {}) const;
   
        const arma::vec getX() const { return this->SolutionX; }
   
        void reset() { this->SolutionX.ones(); }
   
        const arma::vec getZ() const { return this->SolutionZ; }
   
        bool isPureStrategy(double tol = 1e-5) const; 
   
        ZEROStatistics<Data::EPEC::DataObject> getStatistics() const { return this->Stats; }
   
        void setAlgorithm(Data::EPEC::Algorithms algorithm);
   
        void setRecoverStrategy(Data::EPEC::RecoverStrategy strategy);
   
        void setAggressiveness(unsigned int a) { this->Stats.AlgorithmData.Aggressiveness = a; }
   
        void setNumThreads(unsigned int t) {
           this->Stats.AlgorithmData.Threads.set(t);
           this->Env->set(GRB_IntParam_Threads, t);
        }
   
        void setRandomSeed(unsigned int t) { this->Stats.AlgorithmData.RandomSeed.set(t); }
   
        void setIndicators(bool val) { this->Stats.AlgorithmData.IndicatorConstraints.set(val); }
   
        void setPureNashEquilibrium(bool val) { this->Stats.AlgorithmData.PureNashEquilibrium = val; }
   
        void setBoundPrimals(bool val) { this->Stats.AlgorithmData.BoundPrimals.set(val); }
   
        void setBoundBigM(double val) { this->Stats.AlgorithmData.BoundBigM.set(val); }
   
        void setDeviationTolerance(double val) {
           this->Stats.AlgorithmData.DeviationTolerance.set(val);
        }
   
        void setTimeLimit(double val) { this->Stats.AlgorithmData.TimeLimit.set(val); }
   
        void setAddPolyMethod(Data::LCP::PolyhedraStrategy add) {
           this->Stats.AlgorithmData.PolyhedraStrategy.set(add);
        }
        // Methods to get positions of variables
        // The below are all const functions which return an unsigned int.
        int getNumVar() const noexcept { return this->NumVariables; }
   
        unsigned int getNumLeaders() const noexcept {
           return static_cast<int>(this->PlayersLowerLevels.size());
        }
   
        unsigned int getPositionLeadFoll(unsigned int i, unsigned int j) const;
   
        unsigned int getPositionLeadLead(unsigned int i, unsigned int j) const;
   
        // The following obtain the variable values
        double getValLeadFoll(unsigned int i, unsigned int j) const;
   
        double getValLeadLead(unsigned int i, unsigned int j) const;
   
        const LCP &getLCPDescription() const { return *this->TheLCP.get(); }
   
        const GRBModel &getLCPModel() const { return *this->LCPModel.get(); }
   
        void writeLCPModel(const std::string &filename) const { this->LCPModel->write(filename); }
   
        void getXWithoutHull(const arma::vec &x, arma::vec &xWithoutHull) const;
        void
        getXofI(const arma::vec &x, const unsigned int &i, arma::vec &solI, bool hull = false) const;
     };
   }; // namespace Game
   
   namespace std {
   
     string to_string(Data::EPEC::Algorithms al);
   
     string to_string(Data::EPEC::RecoverStrategy st);
   
   }; // namespace std
