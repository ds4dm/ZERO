
.. _program_listing_file_include_algorithms_outerapproximation.h:

Program Listing for File outerapproximation.h
=============================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_algorithms_outerapproximation.h>` (``include/algorithms/outerapproximation.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   #include "algorithms.h"
   #include "epecsolve.h"
   #include <armadillo>
   #include <gurobi_c++.h>
   #include <iostream>
   #include <memory>
   #include <set>
   #include <string>
   
   class OuterTree {
   public:
     struct Node {
     public:
       friend class OuterTree;
   
       Node(unsigned int encSize);
   
       Node(Node &parent, unsigned int idComp, unsigned long int id);
       Node(Node &parent, std::vector<int> idComps, unsigned long int id);
   
       inline std::vector<unsigned int> getIdComps() const {
         return this->IdComps;
       } 
       inline unsigned long int getId() const {
         return this->Id;
       } 
       inline unsigned long int getCumulativeBranches() const {
         return std::count(this->AllowedBranchings.begin(),
                           this->AllowedBranchings.end(), false);
       } 
       inline std::vector<bool> getEncoding() const {
         return this->Encoding;
       } 
   
       inline std::vector<bool> getAllowedBranchings() const {
         return this->AllowedBranchings;
       } 
   
       inline Node *getParent() const {
         return this->Parent;
       } 
   
     private:
       std::vector<unsigned int>
           IdComps; 
       std::vector<bool>
           Encoding; 
       std::vector<bool>
           AllowedBranchings; 
       unsigned long int
           Id;       
       Node *Parent; 
     };
   
   private:
     Node Root = Node(0);           
     unsigned int EncodingSize = 0; 
     unsigned int NodeCounter = 1;  
     std::vector<Node> Nodes{};     
     arma::sp_mat V{}; 
     arma::sp_mat
         R{}; 
     unsigned int VertexCounter = 0; 
     unsigned int RayCounter = 0;    
     GRBModel *MembershipLP; 
     bool MembershipInit{false};
     bool isPure{false};
     bool isFeasible{false};
   
     unsigned int nextIdentifier() {
       this->NodeCounter++;
       return (this->NodeCounter - 1);
     } 
   
   public:
     OuterTree(unsigned int encSize, GRBEnv *env)
         : MembershipLP(new GRBModel(*env)) {
       this->Root = Node(encSize);
       this->EncodingSize = encSize;
       this->Nodes.push_back(this->Root);
     } 
   
     GRBModel *getMembershipLP() { return this->MembershipLP; }
   
     const bool getMembershipInit() { return this->MembershipInit; }
   
     inline void setMembershipInit() { this->MembershipInit = true; }
   
     inline void resetFeasibility() {
       this->isPure = false;
       this->isFeasible = false;
     }
   
     inline bool getPure() const { return this->isPure; }
   
     inline void setFeasible() { this->isFeasible = true; }
   
     inline void setPure() { this->isPure = true; }
   
     const inline unsigned int getEncodingSize() {
       return this->EncodingSize;
     } 
   
     inline const arma::sp_mat *getV() { return &this->V; }
   
     inline const arma::sp_mat *getR() { return &this->R; }
   
     void incrementVertices(unsigned int increment) {
       this->VertexCounter += increment;
     }
     void incrementRays(unsigned int increment) { this->RayCounter += increment; }
   
     inline const unsigned int getVertexCount() { return this->VertexCounter; }
     inline const unsigned int getRayCount() { return this->RayCounter; }
   
     inline const unsigned int getNodeCount() { return this->NodeCounter; }
   
     inline void addVertex(arma::vec vertex) {
       this->V = arma::join_cols(this->V, arma::sp_mat{vertex.t()});
     }
   
     inline void addRay(arma::vec ray) {
       this->R = arma::join_cols(this->R, arma::sp_mat{ray.t()});
     }
   
     inline bool containsRay(arma::vec ray, double tol) {
       if (ray.size() != this->R.n_cols)
         return false;
       for (int i = 0; i < this->R.n_rows; ++i) {
         bool res = true;
         for (int j = 0; j < this->R.n_cols; ++j) {
           if (std::abs(ray.at(j) - this->R.row(i).at(j)) > 1e-5) {
             res = false;
             break;
           }
         }
         if (res)
           return true;
       }
       return false;
     }
   
     inline bool containsVertex(arma::vec vertex, double tol) {
       if (vertex.size() != this->V.n_cols)
         return false;
       for (int i = 0; i < this->V.n_rows; ++i) {
         bool res = true;
         for (int j = 0; j < this->V.n_cols; ++j) {
           if (std::abs(vertex.at(j) - this->V.row(i).at(j)) > tol) {
             res = false;
             break;
           }
         }
         if (res)
           return true;
       }
       return false;
     }
   
     inline Node *const getRoot() {
       return &this->Root;
     } 
   
     inline std::vector<Node> *getNodes() { return &this->Nodes; };
   
     void denyBranchingLocation(Node &node, const unsigned int &location);
     void denyBranchingLocations(Node &node, const std::vector<int> &locations);
   
     std::vector<long int> singleBranch(const unsigned int idComp, Node &t);
   
     std::vector<long int> multipleBranch(const std::vector<int> idsComp, Node &t);
   };
   
   namespace Algorithms {
   class OuterApproximation : public Algorithm {
   private:
     std::vector<std::shared_ptr<Game::OuterLCP>> outerLCP{};
     std::vector<OuterTree *> Trees;
     std::vector<OuterTree::Node *> Incumbent;
     bool Feasible{false};
     double Tolerance = 1e-6;
   
   public:
     double getTol() const { return Tolerance; }
     void setTol(double tol) { this->Tolerance = tol; }
   
   private:
     std::vector<int> getNextBranchLocation(const unsigned int player,
                                            OuterTree::Node *node);
     int getFirstBranchLocation(const unsigned int player,
                                const OuterTree::Node *node);
   
   protected:
     void postSolving() override{
         //@todo implement
     };
   
   public:
     friend class EPEC;
   
     OuterApproximation(GRBEnv *env, Game::EPEC *EpecObj) {
       this->EPECObject = EpecObj;
       this->Env = env;
       /*
        *  The constructor re-builds the LCP fields in the EPEC object as new
        * OuterLCP objects
        */
       this->EPECObject->Stats.AlgorithmParam.PolyLcp = false;
       this->outerLCP =
           std::vector<std::shared_ptr<Game::OuterLCP>>(EPECObject->NumPlayers);
       for (unsigned int i = 0; i < EPECObject->NumPlayers; i++) {
         this->outerLCP.at(i) = std::shared_ptr<Game::OuterLCP>(new Game::OuterLCP(
             this->Env, *EPECObject->PlayersLowerLevels.at(i).get()));
         EPECObject->PlayersLCP.at(i) = this->outerLCP.at(i);
       }
   
     }; 
     void solve() override;
   
     //@todo define these for the outer approximation
     bool isSolved(double tol = 1e-4) const override;
     bool isFeasible(bool &addedCuts, double tol = 1e-4);
     bool isPureStrategy(double tol = 1e-4) const override;
   
     void printCurrentApprox();
     int hybridBranching(const unsigned int player, OuterTree::Node *node);
     int infeasibleBranching(const unsigned int player,
                             const OuterTree::Node *node);
     int deviationBranching(const unsigned int player,
                            const OuterTree::Node *node);
     void printBranchingLog(std::vector<int> vector);
     std::unique_ptr<GRBModel> getFeasQP(const unsigned int player, arma::vec x);
     void addValueCut(unsigned int player, arma::vec xOfIBestResponse,
                      arma::vec xMinusI);
     bool separationOracle(arma::vec &xOfI, arma::vec &x, unsigned int player,
                           int budget, bool &addedCuts);
     GRBModel *getDualMembershipLP(unsigned int player, arma::vec vertex,
                                   bool normalization = true);
     arma::vec normalizeRay(const arma::vec ray);
   };
   } // namespace Algorithms
