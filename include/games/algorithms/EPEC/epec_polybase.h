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
#include "mathopt/lcp/lcp.h"
#include "zero.h"

namespace Algorithms::EPEC {
	 class PolyBase {
		/**
		 *  @brief This is the abstract class manages the algorithms for Game::EPEC. Since they are
		 * all based on MathOpt::PolyLCP, the class keeps a local copy of objects of that class. It
		 * provides a constructor where the Gurobi environment and the EPEC are passed.
		 */
	 protected:
		GRBEnv *    Env;        ///< A pointer to the Gurobi Environment
		Game::EPEC *EPECObject; ///< A pointer to the Game::EPEC instance
		std::vector<std::shared_ptr<MathOpt::PolyLCP>> PolyLCP{}; ///< Local MathOpt::PolyLCP objects


		/**
		 * @brief This method is called after the PolyBase::solve operation. It fills statistics and
		 * can be forcefully overridden by inheritors. The responsibility for calling this method is
		 * left to the inheritor.
		 */
		void after() {
		  std::vector<unsigned int> feasible;
		  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; i++)
			 feasible.push_back(this->PolyLCP.at(i)->getFeasiblePolyhedra());
		  this->EPECObject->Stats.AlgorithmData.FeasiblePolyhedra.set(feasible);
		  if (this->EPECObject->NashEquilibrium)
		    this->EPECObject->Stats.PureNashEquilibrium = this->isPureStrategy();
        LOG_S(3) << "Algorithms::EPEC::PolyBase::after: post-processing results.";
		}

	 public:
		/**
		 * @brief The standard constructor for a PolyBase algorithm. It creates local MathOpt::PolyLCP
		 * objects to work with.
		 * @param env The pointer to the Gurobi environment
		 * @param EPECObject The pointer to the Game::EPEC object
		 */
		PolyBase(GRBEnv *env, Game::EPEC *EPECObject) {
		  this->EPECObject = EPECObject;
		  this->Env        = env;
		  this->PolyLCP    = std::vector<std::shared_ptr<MathOpt::PolyLCP>>(EPECObject->NumPlayers);
		  for (unsigned int i = 0; i < EPECObject->NumPlayers; i++) {
			 this->PolyLCP.at(i) = std::shared_ptr<MathOpt::PolyLCP>(
				  new class MathOpt::PolyLCP(this->Env, *EPECObject->PlayersLowerLevels.at(i).get()));
			 EPECObject->PlayersLCP.at(i) = this->PolyLCP.at(i);
		  }
		}

		virtual void solve() = 0; ///< A general method to solve problems
		bool isSolved(unsigned int *player, arma::vec *profitableDeviation, double tol = -1e-5) const;
		virtual bool isSolved(double tol = 1e-5); ///< A method to check whether the EPEC is solved or
														  ///< not, given a numerical tolerance
		void makeThePureLCP();

		 double
		getValLeadFollPoly(unsigned int i, unsigned int j, unsigned int k, double tol = 1e-5) const;

		 double
		getValLeadLeadPoly(unsigned int i, unsigned int j, unsigned int k, double tol = 1e-5) const;

		 double getValProbab(unsigned int i, unsigned int k) const;

		 bool isPureStrategy(unsigned int i, double tol = 1e-5) const;

		 bool isPureStrategy(double tol = 1e-5) const;

		 std::vector<unsigned int> mixedStrategyPoly(unsigned int i, double tol = 1e-5) const;
		 unsigned int getPositionLeadFollPoly(unsigned int i, unsigned int j, unsigned int k) const;

		 unsigned int getPositionLeadLeadPoly(unsigned int i, unsigned int j, unsigned int k) const;

		 unsigned long int getNumPolyLead(unsigned int i) const;

		 unsigned int getPositionProbab(unsigned int i, unsigned int k) const;
	 };
  } // namespace Algorithms

#include "epec_combPNE.h"
#include "epec_cutandplay.h"
#include "epec_fullenum.h"
#include "epec_innerapp.h"