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
	 enum class CutsAggressiveness {
		NoThanks,   ///< Do not add "standard" Integer Programming cuts to the game
		KeepItCool, ///< At most one cut per-player
		Truculent   ///< A storm of cuts at each infeasibility detection
	 };
	 enum class Objectives {
		/**
		 * @brief Objective types for the MIP reformulation of the LCPs
		 */
		Feasibility, ///< The LCP objective is feasibility
		Quadratic,   ///< The LCP objective is the quadratic sum of the players objectives
		Linear,      ///< The LCP objective is the linear sum of the players objectives
	 };

	 class DataObject : public ZEROAlgorithmData {
	 public:
		Attr<Data::IPG::CutsAggressiveness> CutAggressiveness = {
			 Data::IPG::CutsAggressiveness::KeepItCool};
		Attr<Data::IPG::Algorithms> Algorithm = {
			 Data::IPG::Algorithms::Oracle};    ///< The selected algorithm
		Attr<Data::LCP::Algorithms> LCPSolver; ///< The preferred LCP Solver
		Attr<Data::IPG::Objectives> Objective = {
			 Data::IPG::Objectives::Linear}; ///< The preferred objective type for the MIP LCP
														///< reformulation
		DataObject() : LCPSolver{static_cast<Data::LCP::Algorithms>(0)} {};
	 };
  } // namespace IPG
} // namespace Data


namespace Game {

  /**
	* @brief This class handles Integer Programming Games (IPG), namely multiple agents solving an
	* integer programming game.
	*/
  class IPG : public AbstractGame<Data::IPG::DataObject> {
  protected: // Datafields
	 std::vector<std::shared_ptr<MathOpt::IP_Param>>
		  PlayersIP{}; ///< The Integer Programs associated to each player

	 std::vector<unsigned int> PlayerVariables{}; ///< The number of variables for each player


	 /**
	  * @brief When the object is Finalized, the solving process can start. No players can be added.
	  */
	 bool                   Finalized{false};
	 std::vector<arma::vec> Solution; ///< Solution variable values, for each player

  private:
	 std::shared_ptr<Algorithms::IPG::Algorithm> Algorithm{};
	 void getXMinusI(const arma::vec &x, const unsigned int &i, arma::vec &xMinusI) const;
	 void getXofI(const arma::vec &x, const unsigned int &i, arma::vec &xOfI) const;


  protected:
	 virtual void preFinalize()  = 0;
	 virtual void postFinalize() = 0;

  public: // functions
	 friend class Algorithms::IPG::Algorithm;
	 friend class Algorithms::IPG::Oracle;
	 void finalize();
	 IPG(GRBEnv *env) { this->Env = env; };
	 IPG(GRBEnv *env, std::vector<std::shared_ptr<MathOpt::IP_Param>> players);

	 const void findNashEq() override;
	 bool       isSolved(double tol = 1e-5) const override;
	 bool isPureStrategy(double tol = 1e-5) const override; ///< Return a bool indicating whether the
	 ///< equilibrium is a pure strategy


	 const std::vector<arma::vec> getX() const { return this->Solution; }

	 ///@brief Get the EPECStatistics object for the current instance
	 ZEROStatistics<Data::IPG::DataObject> getStatistics() const { return this->Stats; }

	 void setAlgorithm(Data::IPG::Algorithms algorithm) {
		this->Stats.AlgorithmData.Algorithm = algorithm;
	 }
	 void setLCPAlgorithm(const Data::LCP::Algorithms algo) {
		this->Stats.AlgorithmData.LCPSolver.set(algo);
	 }
	 void setGameObjective(const Data::IPG::Objectives obj) {
		this->Stats.AlgorithmData.Objective.set(obj);
	 }
	 void setCutsAggressiveness(const Data::IPG::CutsAggressiveness aggressiveness) {
		this->Stats.AlgorithmData.CutAggressiveness.set(aggressiveness);
	 }
  };


} // namespace Game
namespace std {
  string to_string(Data::IPG::Algorithms al);
  string to_string(Data::IPG::Objectives ob);
  string to_string(Data::IPG::CutsAggressiveness ct);

}; // namespace std

#include "algorithms/IPG/ipg_oracle.h"