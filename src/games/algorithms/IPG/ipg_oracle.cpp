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
#include <string>

bool Algorithms::IPG::IPG_Player::addVertex(const arma::vec vertex, const bool checkDuplicate) {

  bool go{true};
  if (checkDuplicate)
	 go = Utils::containsRow(this->V, vertex, this->Tolerance);

  if (go) {
	 this->V = arma::join_cols(this->V, arma::sp_mat{vertex.t()});
	 return true;
  }
  return false;
}

bool Algorithms::IPG::IPG_Player::addRay(const arma::vec ray, const bool checkDuplicate) {
  bool go{true};
  if (checkDuplicate)
	 go = Utils::containsRow(this->R, ray, this->Tolerance);

  if (go) {
	 this->R = arma::join_cols(this->R, arma::sp_mat{ray.t()});
	 return true;
  }
  return false;
}


bool Algorithms::IPG::Oracle::isSolved(double tol) const { return this->Solved; }

bool Algorithms::IPG::Oracle::addValueCut(unsigned int player,
														arma::vec    xOfIBestResponse,
														arma::vec    xMinusI,
														bool         checkDuplicate) {

  double cutRHS =
		this->IPG->PlayersIP.at(player)->computeObjective(xOfIBestResponse, xMinusI, false);
  arma::vec LHS =
		this->IPG->PlayersIP.at(player)->getc() + this->IPG->PlayersIP.at(player)->getC() * xMinusI;
  arma::sp_mat cutLHS = arma::sp_mat{LHS}.t();

  bool go{true};
  if (checkDuplicate)
	 go = Utils::containsConstraint(this->Players.at(player).getCutPoolA(),
											  this->Players.at(player).getCutPoolb(),
											  LHS,
											  cutRHS,
											  this->Tolerance);

  if (go) {
	 BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::addValueCut: "
										 "adding cut for Player "
									 << player;

	 this->IPG->PlayersIP.at(player)->addConstraints(-cutLHS, arma::vec{-cutRHS});
	 return true;
  }
  return false;
}

void Algorithms::IPG::Oracle::solve() {


  this->initialize();

  while (!this->Solved) {

	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 }
  }
}

void Algorithms::IPG::Oracle::initialize() {
  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0)
	 this->IPG->InitTime = std::chrono::high_resolution_clock::now();
  this->IPG->Stats.NumIterations.set(0);

  // Initialize the working objects
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 IPG_Player player = IPG_Player(this->Env,
											  this->IPG->PlayersIP.at(i)->getNy(),
											  this->IPG->PlayersIP.at(i)->getIPModel(),
											  this->Tolerance);
	 this->Players.push_back(player);
	 // Set GRBCallback
	 this->Players.at(i).getIPModel()->setCallback(&player);
  }


  // Reset the working strategies to a pure strategy given by the IP
  // Push back the IP_Param copies in WorkingIPs
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 unsigned int Ny      = this->IPG->PlayerVariables.at(i); // Equals to Ny by definition
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
		  this->Players.at(i).addVertex(this->Players.at(i).Incumbent, true);
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
		  this->Players.at(i).addRay(ray, true);
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
