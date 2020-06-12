
.. _program_listing_file_include_games.h:

Program Listing for File games.h
================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_games.h>` (``include/games.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   #include "epecsolve.h"
   #include <armadillo>
   #include <gurobi_c++.h>
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   
   template <class T>
   std::ostream &operator<<(std::ostream &ost, std::vector<T> v) {
     for (auto elem : v)
       ost << elem << " ";
     ost << '\n';
     return ost;
   }
   
   template <class T, class S>
   std::ostream &operator<<(std::ostream &ost, std::pair<T, S> p) {
     ost << "<" << p.first << ", " << p.second << ">";
     return ost;
   }
   
   namespace Game {
   class PolyLCP; // Forward declaration
   
   arma::vec LPSolve(const arma::sp_mat &A, const arma::vec &b, const arma::vec &c,
                     int &status, bool positivity = false);
   
   unsigned int convexHull(const std::vector<arma::sp_mat *> *Ai,
                           const std::vector<arma::vec *> *bi, arma::sp_mat &A,
                           arma::vec &b, arma::sp_mat Acom = {},
                           arma::vec bcom = {});
   
   void compConvSize(arma::sp_mat &A, unsigned int nFinCons, unsigned int nFinVar,
                     const std::vector<arma::sp_mat *> *Ai,
                     const std::vector<arma::vec *> *bi, const arma::sp_mat &Acom,
                     const arma::vec &bcom);
   
   bool isZero(arma::mat M, double tol = 1e-6) noexcept;
   
   bool isZero(arma::sp_mat M, double tol = 1e-6) noexcept;
   
   // bool isZero(arma::vec M, double Tolerance = 1e-6);
   typedef struct QP_Objective {
     arma::sp_mat Q;
     arma::sp_mat C;
     arma::vec c;
   } QP_objective;
   typedef struct QP_Constraints {
     arma::sp_mat A, B;
     arma::vec b;
   } QP_constraints;
   
   class MP_Param {
   protected:
     // Data representing the parameterized QP
     arma::sp_mat Q, A, B, C;
     arma::vec c, b;
     // Object for sizes and integrity check
     unsigned int Nx, Ny, Ncons;
   
     const unsigned int size();
   
     bool dataCheck(bool forceSymmetry = true) const;
   
     virtual inline bool finalize() {
       this->size();
       return this->dataCheck();
     } 
   
   public:
     // Default constructors
     MP_Param() = default;
   
     MP_Param(const MP_Param &M) = default;
   
     // Getters and setters
     arma::sp_mat getQ() const {
       return this->Q;
     } 
     arma::sp_mat getC() const {
       return this->C;
     } 
     arma::sp_mat getA() const {
       return this->A;
     } 
     arma::sp_mat getB() const {
       return this->B;
     } 
     arma::vec getc() const {
       return this->c;
     } 
     arma::vec getb() const {
       return this->b;
     } 
     unsigned int getNx() const {
       return this->Nx;
     } 
     unsigned int getNy() const {
       return this->Ny;
     } 
   
     MP_Param &setQ(const arma::sp_mat &Q) {
       this->Q = Q;
       return *this;
     } 
     MP_Param &setC(const arma::sp_mat &C) {
       this->C = C;
       return *this;
     } 
     MP_Param &setA(const arma::sp_mat &A) {
       this->A = A;
       return *this;
     } 
     MP_Param &setB(const arma::sp_mat &B) {
       this->B = B;
       return *this;
     } 
     MP_Param &setc(const arma::vec &c) {
       this->c = c;
       return *this;
     } 
     MP_Param &setb(const arma::vec &b) {
       this->b = b;
       return *this;
     } 
   
     // Setters and advanced constructors
     virtual MP_Param &set(const arma::sp_mat &Q, const arma::sp_mat &C,
                           const arma::sp_mat &A, const arma::sp_mat &B,
                           const arma::vec &c,
                           const arma::vec &b); // Copy data into this
     virtual MP_Param &set(arma::sp_mat &&Q, arma::sp_mat &&C, arma::sp_mat &&A,
                           arma::sp_mat &&B, arma::vec &&c,
                           arma::vec &&b); // Move data into this
     virtual MP_Param &set(const QP_Objective &obj, const QP_Constraints &cons);
   
     virtual MP_Param &set(QP_Objective &&obj, QP_Constraints &&cons);
   
     virtual MP_Param &addDummy(unsigned int pars, unsigned int vars = 0,
                                int position = -1);
   
     virtual void write(const std::string &filename, bool append = true) const;
   
     static bool dataCheck(const QP_Objective &obj, const QP_Constraints &cons,
                           bool checkObj = true, bool checkCons = true);
   };
   
   class QP_Param : public MP_Param
   // Shape of C is Ny\times Nx
   {
   private:
     // Gurobi environment and model
     GRBEnv *Env;
     GRBModel QuadModel;
     bool madeyQy;
   
     int makeyQy();
   
   public: // Constructors
     explicit QP_Param(GRBEnv *env = nullptr)
         : Env{env}, QuadModel{(*env)}, madeyQy{false} {
       this->size();
     }
   
     QP_Param(arma::sp_mat Q, arma::sp_mat C, arma::sp_mat A, arma::sp_mat B,
              arma::vec c, arma::vec b, GRBEnv *env = nullptr)
         : Env{env}, QuadModel{(*env)}, madeyQy{false} {
       this->set(Q, C, A, B, c, b);
       this->size();
       if (!this->dataCheck())
         throw("Error in QP_Param::QP_Param: Invalid data for constructor");
     }
   
     QP_Param(const QP_Param &Qu)
         : MP_Param(Qu), Env{Qu.Env}, QuadModel{Qu.QuadModel}, madeyQy{
                                                                   Qu.madeyQy} {
       this->size();
     };
   
     // Override setters
     QP_Param &set(const arma::sp_mat &Q, const arma::sp_mat &C,
                   const arma::sp_mat &A, const arma::sp_mat &B,
                   const arma::vec &c,
                   const arma::vec &b) final; // Copy data into this
     QP_Param &set(arma::sp_mat &&Q, arma::sp_mat &&C, arma::sp_mat &&A,
                   arma::sp_mat &&B, arma::vec &&c,
                   arma::vec &&b) final; // Move data into this
     QP_Param &set(const QP_Objective &obj, const QP_Constraints &cons) final;
   
     QP_Param &set(QP_Objective &&obj, QP_Constraints &&cons) final;
   
     bool operator==(const QP_Param &Q2) const;
   
     // Other methods
     unsigned int KKT(arma::sp_mat &M, arma::sp_mat &N, arma::vec &q) const;
   
     std::unique_ptr<GRBModel> solveFixed(arma::vec x, bool solve);
   
     double computeObjective(const arma::vec &y, const arma::vec &x,
                             bool checkFeas = true, double tol = 1e-6) const;
   
     inline bool isPlayable(const QP_Param &P) const
     {
       bool b1, b2, b3;
       b1 = (this->Nx + this->Ny) == (P.getNx() + P.getNy());
       b2 = this->Nx >= P.getNy();
       b3 = this->Ny <= P.getNx();
       return b1 && b2 && b3;
     }
   
     QP_Param &addDummy(unsigned int pars, unsigned int vars = 0,
                        int position = -1) override;
   
     void write(const std::string &filename, bool append) const override;
   
     void save(const std::string &filename, bool erase = true) const;
   
     long int load(const std::string &filename, long int pos = 0);
     double computeObjectiveWithoutOthers(const arma::vec &y) const;
     arma::vec getConstraintViolations(arma::vec x, arma::vec y, double tol);
   };
   
   class NashGame {
   private:
     GRBEnv *Env = nullptr;
     arma::sp_mat LeaderConstraints; 
     arma::vec LeaderConstraintsRHS; 
     unsigned int NumPlayers;        
     std::vector<std::shared_ptr<QP_Param>>
         Players;                 
     arma::sp_mat MarketClearing; 
     arma::vec MCRHS;             
   
     std::vector<unsigned int> PrimalPosition;
     std::vector<unsigned int> DualPosition;
     unsigned int MC_DualPosition;
     unsigned int LeaderPosition;
     unsigned int numLeaderVar;
   
     void setPositions();
   
   public: // Constructors
     explicit NashGame(GRBEnv *e) noexcept : Env{e} {};
   
     explicit NashGame(GRBEnv *e, std::vector<std::shared_ptr<QP_Param>> players,
                       arma::sp_mat MC, arma::vec MCRHS, unsigned int nLeadVar = 0,
                       arma::sp_mat leadA = {}, arma::vec leadRHS = {});
   
     // Copy constructor
     NashGame(const NashGame &N);
   
     ~NashGame() = default;
   
     // Verbose declaration
     friend std::ostream &operator<<(std::ostream &os, const NashGame &N) {
       os << '\n';
       os << "--------------------------------------------------------------------"
             "---"
          << '\n';
       os << "Nash Game with " << N.NumPlayers << " players" << '\n';
       os << "--------------------------------------------------------------------"
             "---"
          << '\n';
       os << "Number of primal variables:\t\t\t " << N.getNprimals() << '\n';
       os << "Number of dual variables:\t\t\t " << N.getNumDualVars() << '\n';
       os << "Number of shadow price dual variables:\t\t " << N.getNumShadow()
          << '\n';
       os << "Number of leader variables:\t\t\t " << N.getNumLeaderVars() << '\n';
       os << "--------------------------------------------------------------------"
             "---"
          << '\n';
       return os;
     }
   
     inline unsigned int getNprimals() const {
       /***
        * Number of primal variables is the sum of the "y" variables present in
        * each player's Game::QP_Param
        */
       return this->PrimalPosition.back();
     }
   
     inline unsigned int getNumShadow() const { return this->MCRHS.n_rows; }
   
     inline unsigned int getNumLeaderVars() const { return this->numLeaderVar; }
   
     inline unsigned int getNumDualVars() const {
       return this->DualPosition.back() - this->DualPosition.front() + 0;
     }
   
     // Position of variables
     inline unsigned int getPrimalLoc(unsigned int i = 0) const {
       return PrimalPosition.at(i);
     }
   
     inline unsigned int getMCDualLoc() const { return MC_DualPosition; }
   
     inline unsigned int getLeaderLoc() const { return LeaderPosition; }
   
     inline unsigned int getDualLoc(unsigned int i = 0) const {
       return DualPosition.at(i);
     }
   
     // Members
     const NashGame &formulateLCP(arma::sp_mat &M, arma::vec &q, perps &Compl,
                                  bool writeToFile = false,
                                  std::string M_name = "dat/LCP.txt",
                                  std::string q_name = "dat/q.txt") const;
   
     arma::sp_mat rewriteLeadCons() const;
   
     inline arma::vec getLeadRHS() const { return this->LeaderConstraintsRHS; }
   
     inline arma::vec getMCLeadRHS() const {
       return arma::join_cols(
           arma::join_cols(this->LeaderConstraintsRHS, this->MCRHS), -this->MCRHS);
     }
   
     // Check solution and correctness
     std::unique_ptr<GRBModel> respond(unsigned int player, const arma::vec &x,
                                       bool fullvec = true) const;
   
     double respondSol(arma::vec &sol, unsigned int player, const arma::vec &x,
                       bool fullvec = true) const;
   
     arma::vec computeQPObjectiveValues(const arma::vec &x,
                                        bool checkFeas = false) const;
   
     bool isSolved(const arma::vec &sol, unsigned int &violPlayer,
                   arma::vec &violSol, double tol = 1e-4) const;
   
     //  Modify NashGame members
     NashGame &addDummy(unsigned int par = 0, int position = -1);
   
     NashGame &addLeadCons(const arma::vec &a, double b);
   
     // Read/Write Nashgame functions
     void write(const std::string &filename, bool append = true,
                bool KKT = false) const;
   
     void save(const std::string &filename, bool erase = true) const;
   
     long int load(const std::string &filename, long int pos = 0);
     arma::vec computeQPObjectiveValuesWithoutOthers(const arma::vec &x) const;
   };
   
   std::ostream &operator<<(std::ostream &os, const QP_Param &Q);
   
   std::ostream &operator<<(std::ostream &ost, const perps &C);
   
   void print(const perps &C) noexcept;
   } // namespace Game
   
   // The EPEC stuff
   namespace Game {
   
   enum class EPECsolveStatus {
     NashEqNotFound, 
     NashEqFound,    
     TimeLimit,      
     Numerical,      
     Uninitialized   
   };
   
   enum class EPECalgorithm {
     FullEnumeration, 
     InnerApproximation, 
     CombinatorialPne, 
     OuterApproximation 
   };
   
   enum class EPECRecoverStrategy {
     IncrementalEnumeration, 
     Combinatorial 
   };
   
   struct EPECAlgorithmParams {
     Game::EPECalgorithm Algorithm = Game::EPECalgorithm::FullEnumeration;
     Game::EPECRecoverStrategy RecoverStrategy =
         EPECRecoverStrategy::IncrementalEnumeration;
     bool PolyLcp{
         true}; 
     Game::EPECAddPolyMethod AddPolyMethod = Game::EPECAddPolyMethod::Sequential;
     bool BoundPrimals{false}; 
     double BoundBigM{1e5}; 
     double DeviationTolerance{
         51e-4}; 
     long int AddPolyMethodSeed{
         -1}; 
     bool Indicators{true}; 
     double TimeLimit{
         -1}; 
     unsigned int Threads{
         0}; 
     unsigned int Aggressiveness{
         1}; 
     bool PureNashEquilibrium{
         false}; 
   };
   
   struct EPECStatistics {
     Game::EPECsolveStatus Status = Game::EPECsolveStatus::Uninitialized;
     int NumVar = {-1};        
     int NumIterations = {-1}; 
     int NumConstraints = {-1}; 
     int NumNonZero = {-1}; 
     int LostIntermediateEq = {0}; 
     bool NumericalIssues = {
         false}; 
     std::vector<unsigned int> FeasiblePolyhedra =
         {}; 
     double WallClockTime = {0};
     bool PureNashEquilibrium{false}; 
     EPECAlgorithmParams AlgorithmParam =
         {}; 
   };
   
   class EPEC {
   private:
     std::vector<unsigned int> SizesWithoutHull{};
     std::unique_ptr<Game::LCP> TheLCP; 
     std::unique_ptr<GRBModel>
         LCPModel; 
     std::unique_ptr<GRBModel>
         LCPModelBase; 
     unsigned int NumVariables{0};
     unsigned int NumPlayers{0};
     std::shared_ptr<Algorithms::Algorithm> Algorithm{};
   
   protected: // Datafields
     std::vector<std::shared_ptr<Game::NashGame>> PlayersLowerLevels{};
     std::vector<std::shared_ptr<Game::LCP>> PlayersLCP{};
   
     std::vector<std::shared_ptr<Game::QP_Param>>
         PlayersQP{}; 
     std::vector<std::shared_ptr<Game::QP_Objective>>
         LeaderObjective{}; 
     std::vector<std::shared_ptr<Game::QP_Objective>>
         LeaderObjectiveConvexHull{}; 
   
     std::unique_ptr<Game::NashGame> TheNashGame; 
   
     std::vector<unsigned int> LeaderLocations{}; 
     std::vector<const unsigned int *> LocEnds{};
     std::vector<unsigned int> ConvexHullVariables{};
     unsigned int numMCVariables{0};
   
     GRBEnv *Env;
     bool Finalized{false};
     bool NashEquilibrium{
         false}; 
     std::chrono::high_resolution_clock::time_point InitTime;
     EPECStatistics Stats{};      
     arma::vec SolutionZ,         
         SolutionX;               
     bool warmstart(arma::vec x); 
   
   private:
     void addDummyLead(unsigned int i); 
     const void makePlayerQP(unsigned int i);
   
     void makePlayersQPs();
   
     void makeTheLCP();
   
     void computeLeaderLocations(unsigned int addSpaceForMC = 0);
   
     void getXMinusI(const arma::vec &x, const unsigned int &i,
                     arma::vec &solOther) const;
   
     bool computeNashEq(bool pureNE = false, double localTimeLimit = -1.0,
                        bool check = false);
   
   protected: // functions
     EPEC(GRBEnv *env)
         : Env{env} {}; 
   
     // virtual function to be implemented by the inheritor.
     virtual void makeObjectivePlayer(const unsigned int i,
                                      Game::QP_Objective &QP_obj) = 0;
   
     // virtual function to be optionally implemented by the inheritor.
     virtual void preFinalize();
   
     virtual void postFinalize();
   
     virtual void
     updateLocations() = 0; // If any location tracking system is implemented, that
     // can be called from in here.
     virtual void makeMCConstraints(arma::sp_mat &MC, arma::vec &RHS) const {
       MC.zeros();
       RHS.zeros();
     };
   
   public: // functions
     // Friends algorithmic classes
     friend class Algorithms::PolyBase;
   
     friend class Algorithms::InnerApproximation;
   
     friend class Algorithms::OuterApproximation;
   
     friend class Algorithms::CombinatorialPNE;
   
     friend class Algorithms::FullEnumeration;
   
     EPEC() = delete;       // No default constructor
     EPEC(EPEC &) = delete; // Abstract class - no copy constructor
     ~EPEC() = default;     // Destructor to free data
   
     void finalize();
   
     const void findNashEq();
     bool isSolved(double tol = 1e-5) const;
   
     std::unique_ptr<GRBModel> respond(const unsigned int i,
                                       const arma::vec &x) const;
   
     double respondSol(arma::vec &sol, unsigned int player, const arma::vec &x,
                       const arma::vec &prevDev = {}) const;
   
     const arma::vec getX() const { return this->SolutionX; }
   
     void reset() { this->SolutionX.ones(); }
   
     const arma::vec getZ() const { return this->SolutionZ; }
   
     bool isPureStrategy(
         double tol = 1e-5) const; 
   
     const EPECStatistics getStatistics() const { return this->Stats; }
   
     void setAlgorithm(Game::EPECalgorithm algorithm);
   
     Game::EPECalgorithm getAlgorithm() const {
       return this->Stats.AlgorithmParam.Algorithm;
     }
   
     void setRecoverStrategy(Game::EPECRecoverStrategy strategy);
   
     Game::EPECRecoverStrategy getRecoverStrategy() const {
       return this->Stats.AlgorithmParam.RecoverStrategy;
     }
   
     void setAggressiveness(unsigned int a) {
       this->Stats.AlgorithmParam.Aggressiveness = a;
     }
   
     unsigned int getAggressiveness() const {
       return this->Stats.AlgorithmParam.Aggressiveness;
     }
   
     void setNumThreads(unsigned int t) {
       this->Stats.AlgorithmParam.Threads = t;
       this->Env->set(GRB_IntParam_Threads, t);
     }
   
     unsigned int getNumThreads() const {
       return this->Stats.AlgorithmParam.Threads;
     }
   
     void setAddPolyMethodSeed(unsigned int t) {
       this->Stats.AlgorithmParam.AddPolyMethodSeed = t;
     }
   
     unsigned long getAddPolyMethodSeed() const {
       return this->Stats.AlgorithmParam.AddPolyMethodSeed;
     }
   
     void setIndicators(bool val) { this->Stats.AlgorithmParam.Indicators = val; }
   
     bool getIndicators() const { return this->Stats.AlgorithmParam.Indicators; }
   
     void setPureNashEquilibrium(bool val) {
       this->Stats.AlgorithmParam.PureNashEquilibrium = val;
     }
   
     bool getPureNashEquilibrium() const {
       return this->Stats.AlgorithmParam.PureNashEquilibrium;
     }
   
     void setBoundPrimals(bool val) {
       this->Stats.AlgorithmParam.BoundPrimals = val;
     }
   
     bool getBoundPrimals() const {
       return this->Stats.AlgorithmParam.BoundPrimals;
     }
   
     void setBoundBigM(double val) { this->Stats.AlgorithmParam.BoundBigM = val; }
   
     double getBoundBigM() const { return this->Stats.AlgorithmParam.BoundBigM; }
   
     void setDeviationTolerance(double val) {
       this->Stats.AlgorithmParam.DeviationTolerance = val;
     }
   
     double getDeviationTolerance() const {
       return this->Stats.AlgorithmParam.DeviationTolerance;
     }
   
     void setTimeLimit(double val) { this->Stats.AlgorithmParam.TimeLimit = val; }
   
     double getTimeLimit() const { return this->Stats.AlgorithmParam.TimeLimit; }
   
     void setAddPolyMethod(Game::EPECAddPolyMethod add) {
       this->Stats.AlgorithmParam.AddPolyMethod = add;
     }
   
     Game::EPECAddPolyMethod getAddPolyMethod() const {
       return this->Stats.AlgorithmParam.AddPolyMethod;
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
   
     void writeLCPModel(const std::string &filename) const {
       this->LCPModel->write(filename);
     }
   
     void getXWithoutHull(const arma::vec &x, arma::vec &xWithoutHull) const;
     void getXofI(const arma::vec &x, const unsigned int &i, arma::vec &solI,
                  bool hull = false) const;
   };
   } // namespace Game
   
   namespace std {
   string to_string(Game::EPECsolveStatus st);
   
   string to_string(Game::EPECalgorithm al);
   
   string to_string(Game::EPECRecoverStrategy st);
   
   string to_string(Game::EPECAlgorithmParams al);
   
   string to_string(Game::EPECAddPolyMethod add);
   }; // namespace std
   
   /* Example for QP_Param */
   
   /* Example of NashGame */
