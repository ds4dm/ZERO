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


#include "games/algorithms/EPEC/epec_innerapp.h"

void Algorithms::EPEC::InnerApproximation::solve() {
  /**
	* @brief This is the high-level method that solves the problem with the inner approximation
	* algorithm.
	*/
  this->start();
  this->after();
}

/**
 * @brief Private main component of the algorithm.  Starting from some profitable deviations from an
 * all-zero strategy vector, the algorithm computes the polyhedra containing such deviations, and
 * add them to the approximation. If an approximate equilibrium is found, then the algorithms keeps
 * adding polyhedra by profitable deviation. Otherwise, it adds a random number of
 * Data::EPEC::DataObject::Aggressiveness polyhedra with the method
 * Data::EPEC::DataObject::PolyhedraStrategy
 */
void Algorithms::EPEC::InnerApproximation::start() {
  // Set the initial point for all countries as 0 and solve the respective LCPs?
  this->EPECObject->SolutionX.zeros(this->EPECObject->NumVariables);
  bool solved = {false};
  bool addRandPoly{false};
  bool infeasCheck{false};
  // When true, a MNE has been found. The algorithm now tries to find a PNE, at
  // the cost of incrementally enumerating the remaining polyhedra, up to the
  // TimeLimit (if any).
  //
  bool incrementalEnumeration{false};

  std::vector<arma::vec> prevDevns(this->EPECObject->NumPlayers);
  this->EPECObject->Stats.NumIterations = 0;
  if (this->EPECObject->Stats.AlgorithmData.PolyhedraStrategy.get() ==
		Data::LCP::PolyhedraStrategy::Random) {
	 for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i) {
		// 42 is the answer, we all know
		long int seed             = this->EPECObject->Stats.AlgorithmData.RandomSeed.get() < 0
												  ? std::chrono::high_resolution_clock::now().time_since_epoch().count() +
                                42 + PolyLCP.at(i)->getNumRows()
												  : this->EPECObject->Stats.AlgorithmData.RandomSeed.get();
		PolyLCP.at(i)->RandomSeed = seed;
	 }
  }
  if (this->EPECObject->Stats.AlgorithmData.TimeLimit.get() > 0)
	 this->EPECObject->InitTime = std::chrono::high_resolution_clock::now();

  // Stay in this loop, till you find a Nash equilibrium or prove that there
  // does not exist a Nash equilibrium or you run out of time.
  while (!solved) {
	 this->EPECObject->Stats.NumIterations.set(this->EPECObject->Stats.NumIterations.get() + 1);
	 LOG_S(INFO) << "Algorithms::EPEC::InnerApproximation::solve: Iteration "
					 << std::to_string(this->EPECObject->Stats.NumIterations.get());

	 if (addRandPoly) {
		LOG_S(INFO) << "Algorithms::EPEC::InnerApproximation::solve: using "
							"heuristic polyhedra selection";
		bool success =
			 this->addRandomPoly2All(this->EPECObject->Stats.AlgorithmData.Aggressiveness.get(),
											 this->EPECObject->Stats.NumIterations.get() == 1);
		if (!success) {
		  this->EPECObject->Stats.Status.set(ZEROStatus::NashEqNotFound);
		  return;
		}
	 } else { // else we are in the case of finding deviations.
		unsigned int deviatedCountry{0};
		arma::vec    countryDeviation{};
		if (this->isSolved(&deviatedCountry, &countryDeviation)) {
		  this->EPECObject->Stats.Status.set(ZEROStatus::NashEqFound);
		  this->EPECObject->Stats.PureNashEquilibrium = this->isPureStrategy();
		  if ((this->EPECObject->Stats.AlgorithmData.PureNashEquilibrium.get() &&
				 !this->EPECObject->Stats.PureNashEquilibrium.get())) {
			 // We are seeking for a pure strategy. Then, here we switch between an
			 // incremental
			 // enumeration or combinations of pure strategies.
			 if (this->EPECObject->Stats.AlgorithmData.RecoverStrategy.get() ==
				  Data::EPEC::RecoverStrategy::IncrementalEnumeration) {
				LOG_S(INFO) << "Algorithms::EPEC::InnerApproximation::solve: "
									"triggering recover strategy "
									"(IncrementalEnumeration)";
				incrementalEnumeration = true;
			 } else if (this->EPECObject->Stats.AlgorithmData.RecoverStrategy.get() ==
							Data::EPEC::RecoverStrategy::Combinatorial) {
				LOG_S(INFO) << "Algorithms::EPEC::InnerApproximation::solve: triggering "
									"recover strategy (Combinatorial)";
				// In this case, we want to try all the combinations of pure
				// strategies, except the ones between polyhedra we already tested.
				std::vector<std::set<unsigned long int>> excludeList;
				for (unsigned long j = 0; j < this->EPECObject->NumPlayers; ++j) {
				  excludeList.push_back(PolyLCP.at(j)->getAllPolyhedra()[0]);
				}
				Algorithms::EPEC::CombinatorialPNE combPNE(this->Env, this->EPECObject);
				combPNE.solveWithExcluded(excludeList);
				return;
			 }

		  } else
			 return;
		}
		// Vector of deviations for the countries
		std::vector<arma::vec> devns = std::vector<arma::vec>(this->EPECObject->NumPlayers);
		this->getAllDeviations(devns, this->EPECObject->SolutionX, prevDevns);
		prevDevns              = devns;
		unsigned int addedPoly = this->addDeviatedPolyhedron(devns, infeasCheck);
		if (addedPoly == 0 && this->EPECObject->Stats.NumIterations.get() > 1 &&
			 !incrementalEnumeration) {
		  LOG_S(ERROR) << " In Algorithms::EPEC::InnerApproximation::solve: Not "
								"Solved, but no deviation? Error!\n This might be due to "
								"Numerical issues (tolerances)";
		  this->EPECObject->Stats.Status.set(ZEROStatus::Numerical);
		  solved = true;
		}
		if (infeasCheck && this->EPECObject->Stats.NumIterations.get() == 1) {
		  LOG_S(WARNING) << " In Algorithms::EPEC::InnerApproximation::solve: Problem is "
								  "infeasible";
		  this->EPECObject->Stats.Status.set(ZEROStatus::NashEqNotFound);
		  return;
		}
	 }

    if (this->EPECObject->Stats.AlgorithmData.TimeLimit.get() > 0) {
	   const std::chrono::duration<double> timeElapsed =
		    std::chrono::high_resolution_clock::now() - this->EPECObject->InitTime;
	   const double timeRemaining =
		    this->EPECObject->Stats.AlgorithmData.TimeLimit.get() - timeElapsed.count();
	   if (timeRemaining <= 0) {
		  if (!incrementalEnumeration)
		    this->EPECObject->Stats.Status.set(ZEROStatus::TimeLimit);
		  return;
	   }
    }
	 this->EPECObject->makePlayersQPs();

	 // TimeLimit
	 if (this->EPECObject->Stats.AlgorithmData.TimeLimit.get() > 0) {
		const std::chrono::duration<double> timeElapsed =
			 std::chrono::high_resolution_clock::now() - this->EPECObject->InitTime;
		const double timeRemaining =
			 this->EPECObject->Stats.AlgorithmData.TimeLimit.get() - timeElapsed.count();
		addRandPoly = !this->EPECObject->computeNashEq(
								this->EPECObject->Stats.AlgorithmData.PureNashEquilibrium.get(),
								timeRemaining,
								false,
								true,
								false) &&
						  !incrementalEnumeration;
	 } else {
		// No Time Limit
		addRandPoly = !this->EPECObject->computeNashEq(
								this->EPECObject->Stats.AlgorithmData.PureNashEquilibrium.get(),
								0,
								false,
								true,
								false) &&
						  !incrementalEnumeration;
	 }

	 //Numerical issues
	 if (addRandPoly && this->EPECObject->Stats.Status.get() == ZEROStatus::Numerical){
		return;
	 }

	 if (addRandPoly)
		this->EPECObject->Stats.AlgorithmData.LostIntermediateEq.set(
			 this->EPECObject->Stats.AlgorithmData.LostIntermediateEq.get() + 1);
	 for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i) {
		LOG_S(INFO) << "Country " << i << PolyLCP.at(i)->feasabilityDetailString();
	 }
	 // This might be reached when a NashEq is found, and need to be verified.
	 // Anyway, we are over the TimeLimit and we should stop
	 if (this->EPECObject->Stats.AlgorithmData.TimeLimit.get() > 0) {
		const std::chrono::duration<double> timeElapsed =
			 std::chrono::high_resolution_clock::now() - this->EPECObject->InitTime;
		const double timeRemaining =
			 this->EPECObject->Stats.AlgorithmData.TimeLimit.get() - timeElapsed.count();
		if (timeRemaining <= 0) {
		  if (!incrementalEnumeration)
			 this->EPECObject->Stats.Status.set(ZEROStatus::TimeLimit);
		  return;
		}
	 }
  }
}


/**
 * @brief  Makes a call to to MathOpt::PolyLCP::addAPoly for each player,  and tries to add a
 * polyhedron to get a better inner approximation for the LCP. @p aggressiveLevel is the maximum
 * number of polyhedra it will try to add to each player. Setting it to an arbitrarily high value
 * will mimic complete enumeration.
 * @param aggressiveLevel The maximum number of polyhedra to be added to each player
 * @param stopOnSingleInfeasibility If set to true, the function will return false if it cannot add
 * a single polyhedron to a country
 * @return True when at least a polyhedron is added
 */
bool Algorithms::EPEC::InnerApproximation::addRandomPoly2All(unsigned int aggressiveLevel,
																				 bool stopOnSingleInfeasibility) {
  LOG_S(1) << "Adding Random polyhedra to countries";
  bool infeasible{true};
  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; i++) {
	 auto addedPoly = PolyLCP.at(i)->addAPoly(
		  aggressiveLevel, this->EPECObject->Stats.AlgorithmData.PolyhedraStrategy.get());
	 if (stopOnSingleInfeasibility && addedPoly == 0) {
		LOG_S(INFO) << "Algorithms::EPEC::InnerApproximation::addRandomPoly2All: No Nash "
							"equilibrium. due to "
							"infeasibility of country "
						<< i;
		return false;
	 }
	 if (addedPoly > 0)
		infeasible = false;
  }
  return !infeasible;
}


/**
 * @brief Given a potential solution vector @p guessSol, it returns the profitable deviations (if
 * any) for all players in @p deviations
 * @param deviations [out] The vector of deviations for all players
 * @param guessSol [in] The guessed solution
 * @param prevDev [in] The previous vector of deviations, if any exist.
 * @return
 */
bool Algorithms::EPEC::InnerApproximation::getAllDeviations(
	 std::vector<arma::vec> &      deviations,
	 const arma::vec &             guessSol,
	 const std::vector<arma::vec> &prevDev) const {
  deviations = std::vector<arma::vec>(this->EPECObject->NumPlayers);

  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i) { // For each country
	 // If we cannot compute a deviation, it means model is infeasible!
	 if (this->EPECObject->bestResponse(deviations.at(i), i, guessSol, prevDev.at(i), nullptr) == GRB_INFINITY)
		return false;
	 // cout << "Game::EPEC::getAllDeviations: deviations(i): "
	 // <<deviations.at(i);
  }
  return true;
}


/**
 * @brief Given a vevtor of profitable deviations for all the players, it adds their corresponding
 * polyhedra to the current approximation
 * @param deviations  A vector of vectors containing the deviations
 * @param infeasCheck  [out] If at least one player cannot add a polyhedron, the method places false
 * in this output parameter
 * @return The number of added polyhedra
 */
unsigned int Algorithms::EPEC::InnerApproximation::addDeviatedPolyhedron(
	 const std::vector<arma::vec> &deviations, bool &infeasCheck) const {

  infeasCheck        = false;
  unsigned int added = 0;
  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i) { // For each country
	 bool ret = false;
	 if (!deviations.at(i).empty())
		ret = PolyLCP.at(i)->addPolyFromX(deviations.at(i), true);
	 if (ret) {
		LOG_S(1) << "Algorithms::EPEC::InnerApproximation::"
						"addDeviatedPolyhedron: added "
						"polyhedron for player "
					<< i;
		++added;
	 } else {
		infeasCheck = true;
		LOG_S(1) << "Algorithms::EPEC::InnerApproximation::addDeviatedPolyhedron: NO "
						"polyhedron added for player "
					<< i;
	 }
  }
  return added;
}