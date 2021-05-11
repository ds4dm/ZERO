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
#include <iostream>
#include <memory>
#include <set>
#include <string>

namespace Game {

  /**
	* @brief An abstract class for a game.
	* @tparam DataObjectType The object type for the data
	*  @code
	* class aGame : public AbstractGame<Data::aGame::DataObject>{
	*
	*	 public:
	*	   //Override AbstractGame methods
	*     void findNashEq() override;
	*     bool       isSolved(double tol = 1e-5) const override;
	*     bool isPureStrategy(double tol = 1e-5) const override;
	*
	* }
	* @endcode
	*/
  template <typename DataObjectType> class AbstractGame {
  protected:
	 std::chrono::high_resolution_clock::time_point InitTime;
	 ZEROStatistics<DataObjectType> Stats = ZEROStatistics<DataObjectType>(DataObjectType());
	 ;                             ///< Store run time information
	 GRBEnv *     Env{};           ///< The Gurobi environment
	 unsigned int NumVariables{0}; ///< The number of variables in the game
	 unsigned int NumPlayers{0};   ///< The number of players in the game
	 bool NashEquilibrium{false};  ///< True if computeNashEq returned an equilibrium. Note that this
	 ///< can be the equilibrium of an approximation, and not to the
	 ///< original game. Refer to isSolved() to get a definitive answer
  public:
	 /**
	  * @brief Standard constructor with the Gurobi Environment
	  * @param env A pointer to the Gurobi environment
	  */
	 AbstractGame(GRBEnv *env) : Env{env} {};
	 /**
	  * @brief Default constructor
	  */
	 AbstractGame() = default;
	 /**
	  * @brief Copy constructor. Not active.
	  */
	 AbstractGame(AbstractGame &) = delete;

	 /**
	  * @brief Deconstructor
	  */
	 ~AbstractGame()           = default;
	 virtual void findNashEq() = 0; ///< The main method to start the solving process
	 virtual bool isSolved(double tol = 1e-5)
		  const = 0; ///< Return a bool true if the strategies are all pure, for any player
	 virtual bool
	 isPureStrategy(double tol = 1e-5) const = 0; ///< Return a bool indicating whether the
	 ///< equilibrium is a pure strategy
	 /**
	  * @brief Getter for statistics
	  * @return Returns the appropriate Data Object Type
	  */
	 ZEROStatistics<DataObjectType> getStatistics() const { return this->Stats; }
	 /**
	  * @brief Sets the number of threads for Gurobi
	  * @param t The number of threads
	  */
	 void setNumThreads(unsigned int t) {
		this->Stats.AlgorithmData.Threads.set(t);
		this->Env->set(GRB_IntParam_Threads, t);
	 }
	 /**
	  * @brief Sets the random seed for pseudo-random operations
	  * @param t The seed
	  */
	 void setRandomSeed(unsigned int t) { this->Stats.AlgorithmData.RandomSeed.set(t); }

	 /**
	  * @brief Requires the algorithm to find a pure equilibrium.
	  * @warning This field may not be used by all the inheritor algorithms.
	  * @param val The boolean requirement
	  */
	 void setPureNashEquilibrium(bool val) { this->Stats.AlgorithmData.PureNashEquilibrium = val; }
	 /**
	  * @brief Sets the tolerance for profitable deviations
	  * @param val Deviation tolerance
	  */
	 void setDeviationTolerance(double val) {
		this->Stats.AlgorithmData.DeviationTolerance.set(val);
	 }

	 /**
	  * @brief Sets the timelimit
	  * @param val The timelimit
	  */
	 void setTimeLimit(double val) { this->Stats.AlgorithmData.TimeLimit.set(val); }
	 /**
	  * @brief Gets the number of variables
	  * @return The number of variable
	  */
	 int getNumVar() const noexcept { return this->NumVariables; }
	 /**
	  * @brief Gets the number of players
	  * @return The number of players
	  */
	 int getNumPlayers() const noexcept { return this->NumPlayers; }
  };


} // namespace Game

#include "games/epec.h"
#include "games/ipg.h"
#include "games/nash.h"