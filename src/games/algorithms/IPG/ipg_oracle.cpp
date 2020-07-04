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


void Algorithms::IPG::IPG_Player::updateIPModel(std::unique_ptr<GRBModel> IPmodel) {
  /**
	* @brief Gets the model in IPModel, and adds the known cuts from the current pool
	*/
  this->Model = std::unique_ptr<GRBModel>(IPmodel.get());
  for (unsigned int i = 0; i < this->CutPool_A.n_rows; i++) {
	 GRBLinExpr LHS{0};
	 for (auto j = this->CutPool_A.begin_row(i); j != this->CutPool_A.end_row(i); ++j)
		LHS += (*j) * this->Model->getVarByName("y_" + std::to_string(j.col()));
	 this->Model->addConstr(LHS, GRB_LESS_EQUAL, this->CutPool_b[i]);
  }
  this->Model->update();
  this->Model->relax();
}

bool Algorithms::IPG::IPG_Player::addVertex(const arma::vec vertex, const bool checkDuplicate) {
  /**
	* @brief Given @p vertex, it adds a vertex to the field R. If @p checkDuplicate is true,
	* it will check whether the vertex is already contained in the bool.
	* @return true if the vertex is added.
	*/
  bool go{true};
  if (checkDuplicate)
	 go = Utils::containsRow(this->V, vertex, this->Tolerance);

  if (go) {
	 this->V = arma::join_cols(this->V, arma::sp_mat{vertex.t()});
	 return true;
  }
  return false;
}

bool Algorithms::IPG::IPG_Player::addCut(const arma::vec LHS,
													  const double    b,
													  const bool      checkDuplicate) {
  /**
	* @brief Given @p LHS, @p b, it adds the inequality to the field CutPool_A and b. If
	* @p checkDuplicate is true, it will check whether the inequality is already contained in the
	* bool.
	* @return true if the inequality is added.
	*/

  bool go{true};
  if (checkDuplicate)
	 go = Utils::containsConstraint(this->CutPool_A, this->CutPool_b, LHS, b, this->Tolerance);

  if (go) {
	 this->CutPool_A = arma::join_cols(this->CutPool_A, arma::sp_mat{LHS.t()});
	 this->CutPool_b = arma::join_cols(this->CutPool_b, arma::vec{b});
	 return true;
  }
  return false;
}


bool Algorithms::IPG::IPG_Player::addRay(const arma::vec ray, const bool checkDuplicate) {
  /**
	* @brief Given @p ray, it adds a ray to the field R. If @p checkDuplicate is true, it
	* will check whether the ray is already contained in the bool.
	* @return true if the ray is added.
	*/
  bool go{true};
  if (checkDuplicate)
	 go = Utils::containsRow(this->R, ray, this->Tolerance);

  if (go) {
	 this->R = arma::join_cols(this->R, arma::sp_mat{ray.t()});
	 return true;
  }
  return false;
}


bool Algorithms::IPG::Oracle::addValueCut(unsigned int player,
														arma::vec    xOfIBestResponse,
														arma::vec    xMinusI,
														bool         checkDuplicate) {
  /**
	* @brief Given a player @p player, one of its best responses @p xOfIBestResponses, the
	* strategies of the other players @p xMinusI, it adds an inequality of the type \f[ f^i(x^i,
	* &\bar& x^{-i}) \geq f^i(\hat x^i, \bar x^{-i})\f] to the cut pool of that player.
	* @p checkDuplicate controls whether the methods search for duplicate inequalities in the
	* pool.
	*/

  double cutRHS =
		this->IPG->PlayersIP.at(player)->computeObjective(xOfIBestResponse, xMinusI, false);
  arma::vec LHS =
		this->IPG->PlayersIP.at(player)->getc() + this->IPG->PlayersIP.at(player)->getC() * xMinusI;
  arma::sp_mat cutLHS = arma::sp_mat{LHS}.t();

  return this->Players.at(player).addCut(LHS, cutRHS, checkDuplicate);
}


void Algorithms::IPG::Oracle::solve() {
  /**
	* @brief Solves the IPG with the Oracle algorithm.
	*/


  this->initialize();
  if (this->Infeasible) {
	 this->IPG->Stats.Status.set(ZEROStatus::NashEqNotFound);
	 BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::solve: A Nash Equilibrium has not been "
										 "found. At least one of the players problem is infeasible.";
	 return;
  }

  bool solved{false};
  while (!solved) {


	 // First, check that the computer strategies are building an equilibrium
	 bool Equilibrium{false};
	 while (!Equilibrium) {
		Equilibrium = true;
		for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
		  arma::vec bestResponse;
		  double    payoff;
		  this->computeStrategy(i, bestResponse, payoff);
		  if (!Utils::isZero(bestResponse - this->Players.at(i).Incumbent, this->Tolerance)) {
			 // The player has a profitable deviation. Update!
			 this->Players.at(i).Incumbent = bestResponse;
			 this->Players.at(i).Payoff    = payoff;
			 this->Players.at(i).Pure      = false;
			 this->Players.at(i).Feasible  = false;
			 Equilibrium                   = false;
		  }
		}
	 }
	 // Now we have an equilibrium, then we need to check whether this is feasible or not

	 solved = true;
	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i)
		if (!this->separationOracle(i))
		  solved = false;
  }

  if (solved) {
	 BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::solve: A Nash Equilibrium has been found";
	 this->Solved = true;
	 bool pure    = true;
	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i)
		if (!this->Players.at(i).Pure) {
		  pure = false;
		  break;
		}
	 this->Pure = pure;
	 this->IPG->Stats.Status.set(ZEROStatus::NashEqFound);
  } else {
	 BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::solve: No Nash Equilibrium has been found";
	 this->Solved = false;
	 this->IPG->Stats.Status.set(ZEROStatus::NashEqNotFound);
  }
}



bool Algorithms::IPG::Oracle::separationOracle(const unsigned int player) {
  /**
	* @brief Given the player id @p player, checks whether the current strategy is feasible or
	* not. In order to do so, a more complex separation technique may be called.
	*/
  BOOST_LOG_TRIVIAL(trace) << "Algorithms::IPG::Oracle::addValueCut: "
										"The Oracle has been called for "
									<< player;

  unsigned int Ny      = this->IPG->PlayerVariables.at(player); // Equals to Ny by definition
  arma::vec    xMinusI = this->buildXminusI(player);
  arma::vec *  xOfI    = &this->Players.at(player).Incumbent;

  // Update working strategies with "educated guesses"
  auto PureIP = this->IPG->PlayersIP.at(player)->solveFixed(xMinusI, true);
  int  status = PureIP->get(GRB_IntAttr_Status);
  if (status == GRB_OPTIMAL) {
	 if (std::abs(PureIP->getObjective().getValue() - this->Players.at(player).Payoff) >
		  this->IPG->Stats.AlgorithmData.DeviationTolerance.get()) {
		// Discrepancy between payoffs! Remember the minimization standard

		// Then the strategy is infeasible! It has a better payoff
		if ((PureIP->getObjective().getValue() - this->Players.at(player).Payoff) >
			 this->IPG->Stats.AlgorithmData.DeviationTolerance.get()) {


		  // Infeasible strategy. Add a value-cut
		  this->addValueCut(player, *xOfI, xMinusI);
		  return false;


		} else if ((this->Players.at(player).Payoff - PureIP->getObjective().getValue()) >
					  this->IPG->Stats.AlgorithmData.DeviationTolerance.get()) {

		  // This cannot happen!
		  throw ZEROException(ZEROErrorCode::Numeric,
									 "Invalid payoff relation (better best response)");
		}
	 } else {
		// Payoffs are the same. Need to stop doing this things


		// Check if the strategies are the same!
		arma::vec bestResponse(Ny, arma::fill::zeros);
		for (unsigned int k = 0; k < Ny; ++k)
		  bestResponse.at(k) = PureIP->getVarByName("y_" + std::to_string(k)).get(GRB_DoubleAttr_X);

		if (Utils::isZero(*xOfI - bestResponse, this->Tolerance)) {
		  this->Players.at(player).Pure = true;
		  BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::separationOracle: "
											  "Feasible strategy for Player "
										  << player << " (Best Response)";
		  return true;
		} else {
		  // In this case, we need to call the proper oracle.
		  int iterations = 15;
		  return this->membershipSeparation(player, iterations, *xOfI, xMinusI);
		}
	 }

  } else
	 return false;
  return false;
}

bool Algorithms::IPG::Oracle::isPureStrategy() const {
  /**
	* @brief Returns true if all players are playing a pure strategy in a Nash Equilibrium
	*/
  if (this->Solved)
	 return this->Pure;
  else
	 return false;
}

bool Algorithms::IPG::Oracle::membershipSeparation(const unsigned int player,
																	const unsigned int iterations,
																	const arma::vec &  xOfI,
																	const arma::vec &  xMinusI) {
  /**
	* @brief Given the player and a bound on the number of iterations, tries to decide whether the
	* given strategy belongs to the feasible region of the player by building the convex-hull with
	* the known rays and vertices. @return true if the point belongs to the feasible region.
	* @p xOfI is the given point to separate.
	*/

  for (int k = 0; k < iterations; ++k) {
	 // First, we check whether the point is a convex combination of feasible
	 // KNOWN points
	 this->updateMembership(player, xOfI, true);
	 xOfI.print("Point to separate: ");

	 this->updateMembership(player, xOfI, true);
	 auto convexModel = *this->Players.at(player).MembershipLP;
	 convexModel.optimize();

	 int status = convexModel.get(GRB_IntAttr_Status);
	 BOOST_LOG_TRIVIAL(trace) << "Algorithms::IPG::Oracle::membershipSeparation: "
										  "MermbershipLP status is "
									  << status;
	 if (status == GRB_OPTIMAL) {
		if (convexModel.getObjective().getValue() == 0 &&
			 convexModel.getConstrByName("Normalization").get(GRB_DoubleAttr_Slack) == 1) {
		  // this->Trees.at(player)->addVertex(xOfI);
		  BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::membershipSeparation: "
											  "The point is a convex combination of known points! Player "
										  << player;
		  this->Players.at(player).Feasible = true;

		  arma::vec support;
		  support.zeros(this->Players.at(player).VertexCounter);
		  for (unsigned int v = 0; v < this->Players.at(player).VertexCounter; ++v) {
			 // abs to avoid misunderstanding with sign conventions
			 support.at(v) =
				  convexModel.getConstrByName("V_" + std::to_string(v)).get(GRB_DoubleAttr_Pi);
		  }
		  support.print("MNE Support: ");
		  if (support.max() == 1)
			 this->Players.at(player).Pure = true;
		  return true;
		}
	 }

	 // Else, the status should be OPTIMAL but without the objective of zero
	 if (status == GRB_OPTIMAL) {
		// Get the Farkas' in the form of the unbounded ray of the dual of the
		// dualMembershipLP (the primal)
		BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::membershipSeparation: "
											"The point is NOT a convex combination of known points! Found "
										<< convexModel.get(GRB_IntAttr_SolCount) << " solutions. Player "
										<< player;
		for (int z = 0; z < convexModel.get(GRB_IntAttr_SolCount); ++z) {
		  convexModel.getEnv().set(GRB_IntParam_SolutionNumber, z);
		  arma::vec cutLHS;
		  cutLHS.zeros(xOfI.size());

		  for (unsigned int i = 0; i < xOfI.size(); i++)
			 cutLHS.at(i) = convexModel.getVarByName("y_" + std::to_string(i)).get(GRB_DoubleAttr_X);
		  cutLHS.print("Separating hyperplane: ");

		  // Optimize the resulting inequality over the original feasible set
		  auto       leaderModel = this->IPG->PlayersIP.at(player)->solveFixed(xMinusI, true);
		  GRBLinExpr expr        = 0;
		  for (unsigned int i = 0; i < xOfI.size(); ++i)
			 expr += cutLHS.at(i) * leaderModel->getVarByName("x_" + std::to_string(i));

		  leaderModel->setObjective(expr, GRB_MAXIMIZE);
		  leaderModel->update();
		  leaderModel->set(GRB_IntParam_InfUnbdInfo, 1);
		  leaderModel->set(GRB_IntParam_DualReductions, 0);
		  leaderModel->set(GRB_IntParam_OutputFlag, 0);
		  leaderModel->write("dat/LeaderModel" + std::to_string(player) + ".lp");
		  leaderModel->optimize();
		  status = leaderModel->get(GRB_IntAttr_Status);

		  if (status == GRB_OPTIMAL) {
			 double cutV = leaderModel->getObjective().getValue();
			 BOOST_LOG_TRIVIAL(trace)
				  << "Algorithms::IPG::Oracle::membershipSeparation: "
					  "LeaderModel status = "
				  << std::to_string(status) << " with objective=" << cutV << " for Player " << player;
			 arma::vec val  = cutLHS.t() * xOfI; // c^T xOfI
			 arma::vec val2 = cutLHS.t() * this->Players.at(player).V.row(0).t();
			 BOOST_LOG_TRIVIAL(trace) << "Algorithms::IPG::Oracle::membershipSeparation: c^Tv=" << cutV
											  << " -- c^TxOfI=" << val.at(0) << " -- c^TV(0)=" << val2.at(0);
			 if (cutV - val.at(0) < -this->Tolerance) {
				// False, but we have a cut :-)
				// Ciao Moni
				cutV              = cutV;
				arma::sp_mat cutL = arma::sp_mat{cutLHS}.t();

				if (!this->Players.at(player).addCut(cutLHS, cutV, true)) {
				  BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::membershipSeparation: "
													  "cut already added for Player "
												  << player;
				  break;

				} else {
				  BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::membershipSeparation: "
													  "adding cut for Player "
												  << player;
				  return false;
				}
			 } else {
				// We found a new vertex
				arma::vec v;
				v.zeros(this->Players.at(player).V.n_cols);
				for (unsigned int i = 0; i < v.size(); ++i) {
				  v[i] = leaderModel->getVarByName("x_" + std::to_string(i)).get(GRB_DoubleAttr_X);
				}

				v.print("Vertex found: ");
				if (Utils::containsRow(this->Players.at(player).V, v, this->Tolerance)) {
				  BOOST_LOG_TRIVIAL(warning) << "Algorithms::IPG::Oracle::membershipSeparation: "
														  "duplicate vertex for  player "
													  << player;
				  //@todo
				  break;
				  // throw;
				} else {
				  this->Players.at(player).addVertex(v);
				  v.print("Vertex");
				  BOOST_LOG_TRIVIAL(info)
						<< "Algorithms::IPG::Oracle::membershipSeparation: "
							"adding vertex for Player. "
						<< (iterations - k - 1) << " iterations left for player " << player;
				  break;
				}
			 }

		  } // status optimal for leaderModel
		  else if (status == GRB_UNBOUNDED) {
			 // Check for a new ray
			 if (!Utils::containsRow(this->Players.at(player).R, cutLHS, this->Tolerance)) {
				BOOST_LOG_TRIVIAL(warning) << "Algorithms::IPG::Oracle::membershipSeparation: "
														"new ray for  player "
													<< player;
				this->Players.at(player).addRay(cutLHS);
				break;
			 } else {
				BOOST_LOG_TRIVIAL(warning) << "Algorithms::IPG::Oracle::membershipSeparation: "
														"duplicate ray for player "
													<< player;
				break;
			 }

		  } // status unbounded for leaderModel

		  else
			 throw ZEROException(ZEROErrorCode::Assertion,
										"Unknown status for leaderModel for player " +
											 std::to_string(player));
		} // end for
		// no separation
	 } else {
		throw ZEROException(ZEROErrorCode::Assertion,
								  "Unknown status for convexModel for player " + std::to_string(player));
	 }
  }
  return false;
}

void Algorithms::IPG::Oracle::updateMembership(const unsigned int &player,
															  const arma::vec &   vertex,
															  bool                normalization) {
  /**
	* @brief Updates the membership LP in the Player vector for the @p player, the point @p
	* xOfI, and
	* @p normalization
	*/
  MathOpt::getDualMembershipLP(this->Players.at(player).MembershipLP,
										 this->Players.at(player).VertexCounter,
										 this->Players.at(player).V,
										 this->Players.at(player).RayCounter,
										 this->Players.at(player).R,
										 vertex,
										 normalization);
}

bool Algorithms::IPG::Oracle::computeStrategy(const unsigned int i, ///< [in] The player id
															 arma::vec &strategy,  ///< [out] The computed strategy
															 double &   payoff     ///< [out] The payoff
) {
  /**
	* @brief Given the @p i as the id of the player, retrieves the relaxation of the integer
	* problem of that player and solves it to find any mixed strategies. The resulting strategy is
	* pushed into @p strategy. In general, the result might not be feasible. In other words, the
	* solution may be outside the mixed-integer convex hull. To check whether this is true or not,
	* ask to the separationOracle. The computed strategy is pushed into an IPG_Player object in the
	* field Players of this class.
	* @returns true if the i-th player problem is feasible
	*/
  unsigned int Ny      = this->IPG->PlayerVariables.at(i); // Equals to Ny by definition
  arma::vec    xMinusI = this->buildXminusI(i);
  strategy.zeros(Ny);

  // Update working strategies
  this->Players.at(i).updateIPModel(this->IPG->PlayersIP.at(i)->solveFixed(xMinusI, false));
  auto Model  = *this->Players.at(i).Model;
  int  status = Model.get(GRB_IntAttr_Status);
  if (status == GRB_INFEASIBLE) {
	 // Game ended, player is infeasible
	 this->IPG->Stats.Status.set(ZEROStatus::NashEqNotFound);
	 return false;
  } else if (status == GRB_OPTIMAL) {
	 BOOST_LOG_TRIVIAL(trace) << "Algorithms::IPG::Oracle::computeStrategy: "
										  "A strategy for player "
									  << i << " has been found.";
	 for (unsigned int k = 0; k < Ny; ++k)
		strategy.at(k) = Model.getVarByName("y_" + std::to_string(k)).get(GRB_DoubleAttr_X);
	 payoff = Model.getObjective().getValue();

  } else if (status == GRB_UNBOUNDED) {
	 BOOST_LOG_TRIVIAL(warning) << "Algorithms::IPG::Oracle::computeStrategy: "
											 "(UNBOUNDED PROBLEM) A strategy for player "
										 << i << " has been found.";
	 for (unsigned int k = 0; k < Ny; ++k)
		strategy.at(k) = Model.getVarByName("y_" + std::to_string(k)).get(GRB_DoubleAttr_UnbdRay);
	 payoff = this->IPG->PlayersIP.at(i)->computeObjective(strategy, xMinusI, false);
  }
  return true;
}

void Algorithms::IPG::Oracle::initialize() {
  /**
	* @brief This method initializes some fields for the algorithm. Also, it warm starts the initial
	* strategies to pure best reponses.
	*/
  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0)
	 this->IPG->InitTime = std::chrono::high_resolution_clock::now();
  this->IPG->Stats.NumIterations.set(0);

  // Initialize the working objects
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i)
	 this->Players.push_back(IPG_Player(this->Env,
													this->IPG->PlayersIP.at(i)->getNy(),
													this->IPG->PlayersIP.at(i)->getIPModel(),
													this->Tolerance));


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
		this->Infeasible = true;
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
