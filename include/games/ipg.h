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
#include "zero.h"
#include <armadillo>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Data {
  namespace IPG {

	 enum class Algorithms {
		Oracle ///< Solves the IPG via the separation oracle algorithm
	 };

	 class DataObject : public ZEROAlgorithmData {
	 public:
		Attr<Data::IPG::Algorithms> Algorithm = {
			 Data::IPG::Algorithms::Oracle}; ///< The selected algorithm
		DataObject(){};
	 };
  } // namespace IPG
} // namespace Data


namespace Game {

  class IPG : public AbstractGame<Data::IPG::DataObject> {
	 ///<@brief This class handles Integer Programming Games (IPG), namely multiple agents solving an
	 ///< integer programming game.

  protected: // Datafields
	 std::vector<std::shared_ptr<MathOpt::IP_Param>>
		  PlayersIP{}; ///< The Integer Programs associated to each player

	 std::vector<unsigned int> PlayerVariables{}; ///< The number of variables for each player

	 bool Finalized{false};           ///< When the object is finalized, the solving process
												 ///< can start. No players can be added.
	 std::vector<arma::vec> Solution; ///< Solution variable values, for each player

  private:
	 void getXMinusI(const arma::vec &x, const unsigned int &i, arma::vec &xMinusI) const;
	 void getXofI(const arma::vec &x, const unsigned int &i, arma::vec &xOfI) const;

	 bool computeNashEq(double localTimeLimit = -1.0, bool check = false);
	 void finalize();

  public: // functions
	 friend class Algorithms::IPG::Oracle;
	 IPG(GRBEnv *env) { this->Env = env; };
	 IPG(GRBEnv *env, std::vector<std::shared_ptr<MathOpt::IP_Param>> players);

	 const void findNashEq() override;
	 bool       isSolved(double tol = 1e-5) const override;
	 bool isPureStrategy(double tol = 1e-5) const override; ///< Return a bool indicating whether the
	 ///< equilibrium is a pure strategy

	 std::unique_ptr<GRBModel> respondModel(const unsigned int i, const arma::vec &x) const;

	 const std::vector<arma::vec> getX() const { return this->Solution; }

	 ///@brief Get the EPECStatistics object for the current instance
	 ZEROStatistics<Data::IPG::DataObject> getStatistics() const { return this->Stats; }

	 void setAlgorithm(Data::IPG::Algorithms algorithm);
  };

} // namespace Game

#include "algorithms/IPG/ipg_oracle.h"