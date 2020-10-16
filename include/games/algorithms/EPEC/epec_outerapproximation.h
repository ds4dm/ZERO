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


#pragma once

#include "epec_polybase.h"
#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Algorithms {
  namespace EPEC {

	 /**
	  * @brief This class manages the outer-approximation tree
	  */
	 class OuterTree {
	 public:
		friend Algorithms::EPEC::OuterApproximation;
		struct Node {
		public:
		  friend class OuterTree;

		  Node(unsigned int encSize);

		  Node(Node &parent, unsigned int idComp, unsigned long int id);
		  Node(Node &parent, std::vector<int> idComps, unsigned long int id);

		  /**
			* @brief Returns the number of variables that cannot be candidate for the branching
			* decisions, namely the ones on which a branching decision has already been taken, or for
			* which the resulting child node is infeasible.
			* @return The number of unsuitable branching candidates
			*/
		  inline unsigned long int getCumulativeBranches() const {
			 return std::count(this->AllowedBranchings.begin(), this->AllowedBranchings.end(), false);
		  }

		  inline std::vector<bool> getEncoding() const {
			 return this->Encoding;
		  } ///< Getter method for the encoding.

		  inline std::vector<bool> getAllowedBranchings() const {
			 return this->AllowedBranchings;
		  } ///< Getter method for the allowed branchings


		private:
		  std::vector<unsigned int> IdComps; ///< Contains the branching decisions taken at the node
		  /**
			* An encoding of bool. True if the complementarity condition is included in the current
			* node outer approximation, false otherwise.
			*/
		  std::vector<bool> Encoding;

		  /**
			* A vector where true means that the corresponding complementarity is a candidate for
			* banching at the current node
			*/
		  std::vector<bool> AllowedBranchings;
		  unsigned long int Id;     ///< A long int giving the numerical identifier for the node
		  Node *            Parent; ///< A pointer to the parent node.
		};

	 private:
		Node         Root         = Node(0); ///< The root node of the tree
		unsigned int EncodingSize = 0;       ///< The size of the encoding, namely the
		///< number of complementarity equations
		unsigned int NodeCounter = 1; ///< The counter for node ids

		/**
		 * Storage of nodes in the tree with the vertices in V
		 */
		std::vector<Node> Nodes{};
		bool isPure{false}; ///< True if the strategy at the current node is a pure-strategy
		bool isFeasible{
			 false}; ///< True if the strategy at the current node is feasible for the original game

		unsigned int nextIdentifier() {
		  this->NodeCounter++;
		  return (this->NodeCounter - 1);
		} ///< Increments the node counter and get the id of the new node.

	 protected:
		/**
		 * Stores the pointer to the MembershipLP associated to the tree.
		 */
		std::unique_ptr<GRBModel> MembershipLP;

		/**
		 * Stores the known extreme vertices of the player's feasible region. These are used to derive
		 * valid cuts, or certify that an equilibrium is inside (outside) the convex-hull of the
		 * feasible region.
		 */
		arma::sp_mat V{};
		arma::sp_mat R{}; ///< As in V, but instead of vertices, this object contains rays
		unsigned int VertexCounter = 0; ///< The counter for node ids
		unsigned int RayCounter    = 0; ///< The counter for node ids

	 public:
		OuterTree(unsigned int encSize, GRBEnv *env) : MembershipLP(new GRBModel(*env)) {
		  this->Root         = Node(encSize);
		  this->EncodingSize = encSize;
		  this->Nodes.push_back(this->Root);
		} ///< Constructor of the Tree given the encoding size


		inline void resetFeasibility() {
		  this->isPure     = false;
		  this->isFeasible = false;
		} ///< Reset the feasibility parameters for the tree

		inline bool getPure() const {
		  return this->isPure;
		} ///< Read-only getter for OuterTree:isPure

		inline void setFeasible() {
		  this->isFeasible = true;
		} ///< Read-only getter for OuterTree:isFeasible

		inline void setPure() { this->isPure = true; } ///< Setter for OuterTree:isPure

		const inline unsigned int getEncodingSize() {
		  return this->EncodingSize;
		} ///< Getter for the encoding size

		inline const arma::sp_mat *getV() { return &this->V; } ///< Getter for OuterTree:V

		inline const arma::sp_mat *getR() { return &this->R; } ///< Getter for OuterTree:R

		inline const unsigned int getVertexCount() {
		  return this->VertexCounter;
		} ///< Getter for OuterTree:VertexCounter
		inline const unsigned int getRayCount() {
		  return this->RayCounter;
		} ///< Getter for OuterTree:RayCounter

		inline const unsigned int getNodeCount() {
		  return this->NodeCounter;
		} ///< Getter for OuterTree:NodeCounter


		inline void addVertex(arma::vec vertex);

		inline void addRay(arma::vec ray);

		inline Node *const getRoot() { return &this->Root; } ///< Getter for the root node

		inline std::vector<Node> *getNodes() { return &this->Nodes; }; ///< Getter for all the nodes

		void denyBranchingLocation(Node &node, const unsigned int &location);

		std::vector<long int> singleBranch(const unsigned int idComp, Node &t);
	 };

	 ///@brief This class is responsible for the outer approximation Algorithm
	 class OuterApproximation : public PolyBase {

	 public:
		OuterApproximation(GRBEnv *env, Game::EPEC *EPECObject) : PolyBase(env, EPECObject){};
		double getTol() const {
		  return this->Tolerance;
		} ///< Read-Only getter for OuterApproximation::Tolerance
		void setTol(double tol) {
		  this->Tolerance = tol;
		} ///< Setter for OuterApproximation::Tolerance

		void solve() override;
		void printCurrentApprox();
		void printBranchingLog(std::vector<int> vector);

		//@todo define these for the outer approximation
		bool isSolved(double tol = 1e-4) const;
		bool isFeasible(bool &addedCuts);
		bool isPureStrategy(double tol = 1e-4) const;


	 private:
		std::vector<OuterTree *>       Trees; ///< The vector of pointer to OuterTree for each player
		std::vector<OuterTree::Node *> Incumbent; ///< The incumbent nodes for each player
		bool   Feasible{false};                   ///< True if a feasible solution has been found
		double Tolerance = 1e-6;                  ///< A numberical tolerance

		std::vector<int> getNextBranchLocation(const unsigned int player, OuterTree::Node *node);
		int getFirstBranchLocation(const unsigned int player, const OuterTree::Node *node);

	 protected:
		void postSolving(){
			 //@todo implement
		};

		void updateMembership(const unsigned int &player,
									 const arma::vec &   xOfI,
									 bool                normalization = true);
		int  hybridBranching(const unsigned int player, OuterTree::Node *node);
		int  infeasibleBranching(const unsigned int player, const OuterTree::Node *node);
		int  deviationBranching(const unsigned int player, const OuterTree::Node *node);
		std::unique_ptr<GRBModel> getFeasQP(const unsigned int player, arma::vec x);
		void addValueCut(const unsigned int player, const double RHS, const arma::vec xMinusI);
		bool separationOracle(
			 arma::vec &xOfI, arma::vec &x, unsigned int player, int budget, bool &addedCuts);
	 };
  } // namespace EPEC

} // namespace Algorithms