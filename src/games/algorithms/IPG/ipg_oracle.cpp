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

#include "games/algorithms/IPG/ipg_oracle.h"
#include "zero.h"
#include <boost/log/trivial.hpp>
#include <chrono>
#include <gurobi_c++.h>
#include <set>
#include <string>

bool Algorithms::IPG::Oracle::addVertex(const unsigned int player,
													 const arma::vec    vertex,
													 const bool         checkDuplicate) {
  if ((checkDuplicate &&
		 !Utils::containsRow(this->Players.at(player).V, vertex, this->Tolerance)) ||
		!checkDuplicate) {
	 this->Players.at(player).V =
		  arma::join_cols(this->Players.at(player).V, arma::sp_mat{vertex.t()});
	 return true;
  }
  return false;
}

bool Algorithms::IPG::Oracle::addRay(const unsigned int player,
												 const arma::vec    ray,
												 const bool         checkDuplicate) {
  if ((checkDuplicate && !Utils::containsRow(this->Players.at(player).R, ray, this->Tolerance)) ||
		!checkDuplicate) {
	 this->Players.at(player).V = arma::join_cols(this->Players.at(player).R, arma::sp_mat{ray.t()});
	 return true;
  }
  return false;
}


bool Algorithms::IPG::Oracle::isSolved(double tol) const { return this->Feasible; }

void Algorithms::IPG::Oracle::addValueCut(unsigned int player,
														arma::vec    xOfIBestResponse,
														arma::vec    xMinusI) {

  double cutRHS =
		this->IPG->PlayersIP.at(player)->computeObjective(xOfIBestResponse, xMinusI, false);
  arma::vec LHS =
		this->IPG->PlayersIP.at(player)->getc() + this->IPG->PlayersIP.at(player)->getC() * xMinusI;
  arma::sp_mat cutLHS =
		Utils::resizePatch(arma::sp_mat{LHS}.t(), 1, this->IPG->PlayersIP.at(player)->getC().n_cols);
  BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::addValueCut: "
									  "adding cut for Player "
								  << player;
  //@todo add cut
  this->IPG->PlayersIP.at(player)->addConstraints(-cutLHS, arma::vec{-cutRHS});
}

void Algorithms::IPG::Oracle::solve() {

  bool solved = {false};
  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0)
	 this->IPG->InitTime = std::chrono::high_resolution_clock::now();
  this->IPG->Stats.NumIterations.set(0);

  // Initialize the working objects
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i)
	 this->Players.push_back(IPG_Player(
		  this->Env, this->IPG->PlayersIP.at(i)->getNy(), this->IPG->PlayersIP.at(i)->getIPModel()));


  // Reset the working strategies to a pure strategy given by the IP
  // Push back the IP_Param copies in WorkingIPs
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 unsigned int Ny      = this->IPG->PlayerVariables.at(i); // Equals to Ny by definition
	 unsigned int Nx      = this->IPG->PlayersIP.at(i)->getNx();
	 arma::vec    xMinusI = this->buildXminusI(i);

	 // Update working strategies with "educated guesses"
	 auto PureIP = this->IPG->PlayersIP.at(i)->solveFixed(xMinusI, true);
	 int  status = PureIP->get(GRB_IntAttr_Status);
	 if (status == GRB_INFEASIBLE) {
		// Game ended, player is infeasible
		this->IPG->Stats.Status.set(ZEROStatus::NashEqNotFound);
		return;
	 } else if (status == GRB_OPTIMAL) {
		// We have a vertex
		for (unsigned int k = 0; k < Ny; ++k) {
		  this->Players.at(i).Incumbent.at(k) =
				PureIP->getVarByName("y_" + std::to_string(k)).get(GRB_DoubleAttr_X);
		  // This is also a free best response
		  this->addVertex(i, this->Players.at(i).Incumbent, true);
		}
	 } else if (status == GRB_UNBOUNDED) {
		PureIP->relax();
		PureIP->set(GRB_IntParam_InfUnbdInfo, 1);
		PureIP->set(GRB_IntParam_DualReductions, 0);
		PureIP->optimize();
		arma::vec ray;
		for (unsigned int k = 0; k < Ny; ++k) {
		  ray.at(k) = PureIP->getVarByName("y_" + std::to_string(k)).get(GRB_DoubleAttr_UnbdRay);
		  // This is also a free ray
		  this->addRay(i, ray, true);
		}
	 }
  }
}

arma::vec Algorithms::IPG::Oracle::buildXminusI(const unsigned int i) {
  /**
	* @brief Given the player id @p i, builds the vector x^{-i} from the current working strategies.
	*/
  arma::vec xMinusI;
  xMinusI.zeros(this->IPG->NumVariables);
  unsigned int counter = 0;
  for (unsigned int j = 0; j < this->IPG->NumPlayers; ++j) {
	 if (i != j) {
		xMinusI.subvec(counter, counter + this->Players.at(j).Incumbent.size()) =
			 this->Players.at(j).Incumbent;
		counter += this->Players.at(j).Incumbent.size();
	 }
  }
  return xMinusI;
}
