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
#include "mathopt/lcp/lcp.h"
#include "zero.h"
#include <boost/log/trivial.hpp>

namespace Algorithms {
  namespace EPEC {
	 class PolyBase {
		/**
		 *  @brief This is the abstract class of Algorithms for full enumeration,
		 * inner approximation, and Combinatorial PNE. It provides a constructor where
		 * the Gurobi environment and the EPEC are passed. This is an abstract class.
		 */
	 protected:
		GRBEnv *    Env;        ///< A pointer to the Gurobi Environment
		Game::EPEC *EPECObject; ///< A pointer to the original LCP object
		std::vector<std::shared_ptr<MathOpt::PolyLCP>> PolyLCP{};

		void postSolving() {
		  /**
			* Perform postSolving operations.
			* For instance, it updates the statistics associated with the feasible
			* polyhedra. The responsability for calling this method is left to the
			* inheritor
			*/
		  std::vector<unsigned int> feasible;
		  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; i++)
			 feasible.push_back(this->PolyLCP.at(i)->getFeasiblePolyhedra());
		  this->EPECObject->Stats.AlgorithmData.FeasiblePolyhedra.set(feasible);
		  this->EPECObject->Stats.PureNashEquilibrium = this->isPureStrategy();
		}

	 public:
		PolyBase(GRBEnv *env, Game::EPEC *EPECObject) {
		  /*
			*  The method will reassign the LCP fields in the EPEC object to new
			* PolyLCP objects
			*/
		  this->EPECObject = EPECObject;
		  this->Env        = env;
		  this->PolyLCP    = std::vector<std::shared_ptr<MathOpt::PolyLCP>>(EPECObject->NumPlayers);
		  for (unsigned int i = 0; i < EPECObject->NumPlayers; i++) {
			 this->PolyLCP.at(i) = std::shared_ptr<MathOpt::PolyLCP>(
				  new class MathOpt::PolyLCP(this->Env, *EPECObject->PlayersLowerLevels.at(i).get()));
			 EPECObject->PlayersLCP.at(i) = this->PolyLCP.at(i);
		  }
		}
		virtual void solve() = 0; ///< A method to solve the EPEC
		bool
			  isSolved(unsigned int *countryNumber, arma::vec *profitableDeviation, double tol = -1) const;
		bool isSolved(double tol = -1) const; ///< A method to check whether the EPEC is solved or
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
  } // namespace EPEC
} // namespace Algorithms

#include "epec_combinatorialpne.h"
#include "epec_fullenumeration.h"
#include "epec_innerapproximation.h"
#include "epec_outerapproximation.h"