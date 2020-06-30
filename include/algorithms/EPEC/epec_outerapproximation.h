#pragma once

#include "epec_algorithms.h"
#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Algorithms {
  namespace EPEC {

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
		  }                                                           ///< Getter method for idComp
		  inline unsigned long int getId() const { return this->Id; } ///< Getter method for id
		  inline unsigned long int getCumulativeBranches() const {
			 return std::count(this->AllowedBranchings.begin(), this->AllowedBranchings.end(), false);
		  } ///< Returns the number of variables that cannot be candidate for the
		  ///< branching decisions, namely the ones on which a branching decision
		  ///< has already been taken, or for which the resulting child node is
		  ///< infeasible.
		  inline std::vector<bool> getEncoding() const {
			 return this->Encoding;
		  } ///< Getter method for the encoding.

		  inline std::vector<bool> getAllowedBranchings() const {
			 return this->AllowedBranchings;
		  } ///< Getter method for the allowed branchings

		  inline Node *getParent() const {
			 return this->Parent;
		  } ///< Getter method for the parent node

		private:
		  std::vector<unsigned int> IdComps; ///< Contains the branching decisions taken at the node
		  std::vector<bool> Encoding;        ///< An encoding of bool. True if the complementarity
		  ///< condition is included in the current node outer
		  ///< approximation, false otherwise.
		  std::vector<bool> AllowedBranchings; ///< A vector where true means that the corresponding
		  ///< complementarity is a candidate for banching at
		  ///< the current node
		  unsigned long int Id; ///< A long int giving the numerical identifier for the node
		  Node *Parent;         ///< A pointer to the parent node.
		};

	 private:
		Node Root                 = Node(0); ///< The root node of the tree
		unsigned int EncodingSize = 0;       ///< The size of the encoding, namely the
		///< number of complementarity equations
		unsigned int NodeCounter = 1; ///< The counter for node ids
		std::vector<Node> Nodes{};    ///< Storage of nodes in the tree
		arma::sp_mat V{};             ///< This object stores points that are inside the feasible
		///< region of the respective leader. Thesee are used to derive
		///< valid cuts, or certify that an equilibrium is inside
		///< (outside) the convex-hull of the feasible region.
		arma::sp_mat R{}; ///< As in V, but instead of vertices, this object contains rays
		unsigned int VertexCounter = 0; ///< The counter for node ids
		unsigned int RayCounter    = 0; ///< The counter for node ids
		GRBModel *MembershipLP;         ///< This member stores the membership LP associated
		///< with the vertices in V
		bool MembershipInit{false};
		bool isPure{false};
		bool isFeasible{false};

		unsigned int nextIdentifier() {
		  this->NodeCounter++;
		  return (this->NodeCounter - 1);
		} ///< Increments the node counter and get the id of the new node.

	 public:
		OuterTree(unsigned int encSize, GRBEnv *env) : MembershipLP(new GRBModel(*env)) {
		  this->Root         = Node(encSize);
		  this->EncodingSize = encSize;
		  this->Nodes.push_back(this->Root);
		} ///< Constructor of the Tree given the encoding size

		GRBModel *getMembershipLP() { return this->MembershipLP; }

		const bool getMembershipInit() { return this->MembershipInit; }

		inline void setMembershipInit() { this->MembershipInit = true; }

		inline void resetFeasibility() {
		  this->isPure     = false;
		  this->isFeasible = false;
		}

		inline bool getPure() const { return this->isPure; }

		inline void setFeasible() { this->isFeasible = true; }

		inline void setPure() { this->isPure = true; }

		const inline unsigned int getEncodingSize() {
		  return this->EncodingSize;
		} ///< Getter for the encoding size

		inline const arma::sp_mat *getV() { return &this->V; }

		inline const arma::sp_mat *getR() { return &this->R; }

		void incrementVertices(unsigned int increment) { this->VertexCounter += increment; }
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

		inline Node *const getRoot() { return &this->Root; } ///< Getter for the root node

		inline std::vector<Node> *getNodes() { return &this->Nodes; };

		void denyBranchingLocation(Node &node, const unsigned int &location);
		void denyBranchingLocations(Node &node, const std::vector<int> &locations);

		std::vector<long int> singleBranch(const unsigned int idComp, Node &t);

		std::vector<long int> multipleBranch(const std::vector<int> idsComp, Node &t);
	 };

	 ///@brief This class is responsible for the outer approximation Algorithm
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
		std::vector<int> getNextBranchLocation(const unsigned int player, OuterTree::Node *node);
		int getFirstBranchLocation(const unsigned int player, const OuterTree::Node *node);

	 protected:
		void postSolving() override{
		    //@todo implement
		};

	 public:
		friend class EPEC;

		OuterApproximation(GRBEnv *env, Game::EPEC *EpecObj) {
		  this->EPECObject = EpecObj;
		  this->Env        = env;
		  /*
		   *  The constructor re-builds the LCP fields in the EPEC object as new
		   * OuterLCP objects
		   */
		  this->outerLCP = std::vector<std::shared_ptr<Game::OuterLCP>>(EPECObject->NumPlayers);
		  for (unsigned int i = 0; i < EPECObject->NumPlayers; i++) {
			 this->outerLCP.at(i) = std::shared_ptr<Game::OuterLCP>(
			     new Game::OuterLCP(this->Env, *EPECObject->PlayersLowerLevels.at(i).get()));
			 EPECObject->PlayersLCP.at(i) = this->outerLCP.at(i);
		  }

		}; ///< Constructor requires a pointer to the Gurobi
		///< Environment and the calling EPEC object
		void solve() override;

		//@todo define these for the outer approximation
		bool isSolved(double tol = 1e-4) const override;
		bool isFeasible(bool &addedCuts, double tol = 1e-4);
		bool isPureStrategy(double tol = 1e-4) const override;

		void printCurrentApprox();
		int hybridBranching(const unsigned int player, OuterTree::Node *node);
		int infeasibleBranching(const unsigned int player, const OuterTree::Node *node);
		int deviationBranching(const unsigned int player, const OuterTree::Node *node);
		void printBranchingLog(std::vector<int> vector);
		std::unique_ptr<GRBModel> getFeasQP(const unsigned int player, arma::vec x);
		void addValueCut(unsigned int player, arma::vec xOfIBestResponse, arma::vec xMinusI);
		bool separationOracle(arma::vec &xOfI, arma::vec &x, unsigned int player, int budget,
		                      bool &addedCuts);
		GRBModel *getDualMembershipLP(unsigned int player, arma::vec vertex,
		                              bool normalization = true);
	 };
  } // namespace EPEC

} // namespace Algorithms