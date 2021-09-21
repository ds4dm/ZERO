/* #############################################
 *             This file is part of
 *                    ZERO
 *
 *             Copyright (c) 2020
 *     Released under the Creative Commons
 *         CC BY-NC-SA 4.0 License
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

namespace Data::IPG {

  /**
	* @brief IPG Algorithms enum
	*/
  enum class Algorithms {
	 CutAndPlay ///< Solves the IPG via the separation CutAndPlay algorithm
  };

  /**
	* @brief Cuts aggressiveness for Algorithms::IPG::CutAndPlay
	*/
  enum class CutsAggressiveness {
	 NotEvenTry, ///< Do not add "standard" IPG cuts to the game, nor tries to replace value cuts
	 NoThanks,   ///< Do not add "standard" Integer Programming cuts to the game
	 KeepItCool, ///< At most one cut per-player
	 Truculent   ///< A storm of cuts at each infeasibility detection
  };
  /**
	* @brief Objective types for the MIP reformulation of the LCPs
	*/
  enum class Objectives {
	 Feasibility, ///< The LCP objective is feasibility
	 Quadratic,   ///< The LCP objective is the quadratic sum of the players objectives
	 Linear,      ///< The LCP objective is the linear sum of the players objectives
  };

  class DataObject : public ZEROAlgorithmData {
  public:
	 Attr<Data::IPG::CutsAggressiveness> CutAggressiveness = {
		  Data::IPG::CutsAggressiveness::KeepItCool};
	 Attr<Data::IPG::Algorithms> Algorithm = {
		  Data::IPG::Algorithms::CutAndPlay};    ///< The selected algorithm
	 Attr<Data::LCP::Algorithms> LCPSolver; ///< The preferred LCP Solver
	 Attr<Data::IPG::Objectives> Objective = {
		  Data::IPG::Objectives::Linear}; ///< The preferred objective type for the MIP LCP
	 ///< reformulation
	 Attr<std::vector<std::pair<std::string, int>>>
		  Cuts; ///< Statistics about the added cuts. Refer to the indices in
	 ///< IPG::Algorithms::CutAndPlay
	 /**
	  * @brief Standard initializer constructor.
	  */
	 DataObject() : LCPSolver{static_cast<Data::LCP::Algorithms>(0)} {};
  };
} // namespace Data::IPG


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


	 bool                   Finalized{false};
	 std::vector<arma::vec> Solution; ///< Solution variable values, for each player
	 double                 SocialWelfare;   ///< SocialWelfare associated to the incumbent solution

  private:
	 std::shared_ptr<Algorithms::IPG::Algorithm> Algorithm{}; ///< The Algorithm's instance
	 void getXMinusI(const arma::vec &x, const unsigned int &i, arma::vec &xMinusI) const;
	 void getXofI(const arma::vec &x, const unsigned int &i, arma::vec &xOfI) const;


  protected:
	 /**
	  * @brief Virtual (empty) method. Can be implemented by a derived class
	  */
	 virtual void preFinalize() = 0;
	 /**
	  * @brief Virtual (empty) method. Can be implemented by a derived class
	  */
	 virtual void postFinalize() = 0;

  public: // functions
	 friend class Algorithms::IPG::Algorithm;
	 friend class Algorithms::IPG::CutAndPlay;
	 void finalize();
	 /**
	  * @brief Standard initializer
	  * @param env A pointer to the Gurobi Environment
	  */
	 IPG(GRBEnv *env) { this->Env = env; };
	 IPG(GRBEnv *env, std::vector<std::shared_ptr<MathOpt::IP_Param>> players);

	 void findNashEq() override;
	 bool isSolved(double tol = 1e-5) const override;
	 bool isPureStrategy(double tol = 1e-5) const override; ///< Return a bool indicating whether the
	 ///< equilibrium is a pure strategy


	 /**
	  * @brief Gets the X in the incumbent solution
	  * @return A const vector copy of X.
	  */
	 std::vector<arma::vec> getX() const { return this->Solution; }

	 /***
	  * @brief Gets the SocialWelfare associated to the incumbent solution
	  * @return The payoff value
	  */
	 double getSocialWelfare() const { return this->SocialWelfare; }

	 ///@brief Get the EPECStatistics object for the current instance
	 ZEROStatistics<Data::IPG::DataObject> getStatistics() const { return this->Stats; }

	 /**
	  * @brief Sets the Data::IPG::Algorithms for the solution process.
	  * @param algorithm An enum from Data::IPG::Algorithms
	  */
	 void setAlgorithm(Data::IPG::Algorithms algorithm) {
		this->Stats.AlgorithmData.Algorithm = algorithm;
	 }
	 /**
	  * @brief Sets the Data::LCP::Algorithms for the LCP solution process.
	  * @param algo An enum from Data::LCP::Algorithms
	  */
	 void setLCPAlgorithm(const Data::LCP::Algorithms algo) {
		this->Stats.AlgorithmData.LCPSolver.set(algo);
	 }
	 /**
	  * @brief Sets the Data::IPG::Objectives for the LCP objective.
	  * @param obj An enum from Data::IPG::Objectives
	  */
	 void setGameObjective(const Data::IPG::Objectives obj) {
		this->Stats.AlgorithmData.Objective.set(obj);
	 }
	 /**
	  * @brief Sets the Data::IPG::CutsAggressiveness for the cut aggressiveness in
	  * Algorithms::IPG::CutAndPlay
	  * @param aggressiveness An enum from Data::IPG::CutsAggressiveness
	  */
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

#include "algorithms/IPG/ipg_cutandplay.h"
#include "interfaces/ipg_models.h"