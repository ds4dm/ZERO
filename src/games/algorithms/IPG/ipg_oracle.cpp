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
#include "games/nash.h"
#include "zero.h"
#include <boost/log/trivial.hpp>
#include <chrono>
#include <gurobi_c++.h>
#include <string>


void Algorithms::IPG::IPG_Player::updateIPModel(std::unique_ptr<GRBModel> IPmodel) {
  /**
	* @brief Gets the model in IPModel, and adds the known cuts from the current pool
	*/
  this->Model.release();
  this->Model = std::move(IPmodel);
  for (unsigned int i = 0; i < this->CutPool_A.n_rows; i++) {
	 GRBLinExpr LHS{0};
	 for (auto j = this->CutPool_A.begin_row(i); j != this->CutPool_A.end_row(i); ++j)
		LHS += (*j) * this->Model->getVarByName("y_" + std::to_string(j.col()));
	 this->Model->addConstr(LHS, GRB_LESS_EQUAL, this->CutPool_b[i]);
  }
  this->Model->update();
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

  if (!go) {
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

  if (!go) {
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

  if (!go) {
	 this->R = arma::join_cols(this->R, arma::sp_mat{ray.t()});
	 return true;
  }
  return false;
}


bool Algorithms::IPG::Oracle::addValueCut(unsigned int player,
														double       RHS,
														arma::vec    xMinusI,
														bool         checkDuplicate) {
  /**
	* @brief Given a player @p player, one of its best responses @p xOfIBestResponses, the
	* strategies of the other players @p xMinusI, it adds an inequality of the type \f[ f^i(x^i,
	* &\bar& x^{-i}) \geq f^i(\hat x^i, \bar x^{-i})\f] to the cut pool of that player.
	* @p checkDuplicate controls whether the methods search for duplicate inequalities in the
	* pool.
	*/


  arma::vec LHS =
		this->IPG->PlayersIP.at(player)->getc() + this->IPG->PlayersIP.at(player)->getC() * xMinusI;


  bool go = this->Players.at(player)->addCut(-LHS, -RHS, checkDuplicate);
  if (checkDuplicate)
	 go = !go;
  else
	 go = false;

  if (!go)
	 return this->IPG->PlayersIP.at(player)->addConstraint(
		  -LHS, -RHS, checkDuplicate, this->Tolerance);
  else
	 return false;
}

bool Algorithms::IPG::Oracle::checkTime(double &remaining) const {
  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
	 const std::chrono::duration<double> timeElapsed =
		  std::chrono::high_resolution_clock::now() - this->IPG->InitTime;
	 remaining = this->IPG->Stats.AlgorithmData.TimeLimit.get() - timeElapsed.count();
	 if (remaining <= 0) {
		BOOST_LOG_TRIVIAL(trace) << "Algorithms::IPG::Oracle::checkTime: "
											 "Time limit hit.";
		this->IPG->Stats.Status.set(ZEROStatus::TimeLimit);
		return false;
	 } else
		return true;
  } else {
	 remaining = -1;
	 return true;
  }
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

	 // Increase the number of iterations
	 this->IPG->Stats.NumIterations.set(this->IPG->Stats.NumIterations.get() + 1);
	 // Check the time-limit and form the LCP for the simultaneous game
	 if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
		double remaining;
		if (this->checkTime(remaining) && remaining > 0)
		  bool Eq = this->equilibriumLCP(remaining);
		else
		  return;
	 } else
		bool Eq = this->equilibriumLCP(-1);


	 // Now we have an equilibrium, then we need to check whether this is feasible or not
	 solved = true;
	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
		if (!this->separationOracle(i)) {
		  solved = false;
		  break;
		}
	 }
	 if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
		double remaining;
		if (!this->checkTime(remaining) || remaining <= 0)
		  return;
	 }
  }

  if (solved) {
	 this->Solved = true;
	 bool pure    = true;
	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i)
		if (!this->Players.at(i)->Pure) {
		  pure = false;
		  break;
		}
	 this->Pure = pure;
	 this->IPG->Stats.Status.set(ZEROStatus::NashEqFound);

	 BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::solve: A Nash Equilibrium has been found ("
									 << (pure == 0 ? "MNE" : "PNE") << ").";
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
  BOOST_LOG_TRIVIAL(trace) << "Algorithms::IPG::Oracle::separationOracle: "
										"The Oracle has been called for "
									<< player;

  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
	 double remaining;
	 if (!this->checkTime(remaining) || remaining <= 0)
		return false;
  }

  unsigned int Ny      = this->IPG->PlayerVariables.at(player); // Equals to Ny by definition
  arma::vec    xMinusI = this->buildXminusI(player);
  arma::vec *  xOfI    = &this->Players.at(player)->Incumbent;

  // xOfI->print("xOfI");
  // xMinusI.print("xMinusI");

  // Update working strategies with "educated guesses"
  auto PureIP = this->IPG->PlayersIP.at(player)->getIPModel(xMinusI, false);
  PureIP->optimize();
  PureIP->write("dat/Oracle_" + std::to_string(player) + "+.lp");
  int status = PureIP->get(GRB_IntAttr_Status);
  if (status == GRB_OPTIMAL) {

	 double IPobj  = PureIP->getObjective().getValue();
	 double RELobj = this->IPG->PlayersIP.at(player)->computeObjective(
		  *xOfI, xMinusI, false); // this->Players.at(player)->Payoff;
	 if (std::abs(IPobj - RELobj) > this->IPG->Stats.AlgorithmData.DeviationTolerance.get()) {
		// Discrepancy between payoffs! Remember the minimization standard

		// Then the strategy is infeasible! It has a better payoff
		if ((IPobj - RELobj) > this->IPG->Stats.AlgorithmData.DeviationTolerance.get()) {

		  BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::separationOracle: "
											  "Infeasible strategy. Adding the value-cut.";
		  // Infeasible strategy. Add a value-cut
		  if (this->addValueCut(player, IPobj, xMinusI))
			 return false;
		  else
			 throw ZEROException(ZEROErrorCode::Unknown, "Unknown loop detected");


		} else if ((RELobj - IPobj) > this->IPG->Stats.AlgorithmData.DeviationTolerance.get()) {

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
		  this->Players.at(player)->Pure = true;
		  BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::separationOracle: "
											  "Feasible strategy for Player "
										  << player << " (Best Response)";
		  return true;
		} else {
		  // In this case, we need to call the proper oracle.
		  unsigned int iterations = 100;
		  return this->membershipSeparation(player, iterations, *xOfI, xMinusI);
		}
	 }

  } else if (status == GRB_UNBOUNDED) {
	 BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::separationOracle: "
										 "The problem is unbounded.";
	 return false;
  }
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

  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0) {
	 double remaining;
	 if (!this->checkTime(remaining) || remaining <= 0)
		return false;
  }
  BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::membershipSeparation: "
									  "Starting separator for player "
								  << player;
  for (int k = 0; k < iterations; ++k) {
	 // First, we check whether the point is a convex combination of feasible
	 // KNOWN points
	 xOfI.print("Point to separate: ");

	 this->updateMembership(player, xOfI, true);
	 auto convexModel = *this->Players.at(player)->MembershipLP;
	 convexModel.optimize();
	 convexModel.write("dat/Convex_" + std::to_string(player) + ".lp");

	 int status = convexModel.get(GRB_IntAttr_Status);
	 BOOST_LOG_TRIVIAL(trace) << "Algorithms::IPG::Oracle::membershipSeparation: "
										  "MermbershipLP status is "
									  << status;
	 if (status == GRB_OPTIMAL && convexModel.get(GRB_IntAttr_SolCount) == 1) {
		convexModel.set(GRB_IntParam_SolutionNumber, 0);
		arma::vec sol(xOfI.size(), arma::fill::zeros);
		for (unsigned int i = 0; i < xOfI.size(); i++)
		  sol.at(i) =
				std::abs(convexModel.getVarByName("y_" + std::to_string(i)).get(GRB_DoubleAttr_X));

		if (convexModel.getObjective().getValue() == 0 && sol.max() == 0)
		/*	     ( (convexModel.getConstrByName("Normalization").get(GRB_DoubleAttr_Slack) == 1) ||
			 (convexModel.getConstrByName("Normalization").get(GRB_DoubleAttr_Slack) == 0 &&
			  convexModel.getVarByName("x").get(GRB_DoubleAttr_X) == 1))
			  */
		{
		  // this->Trees.at(player)->addVertex(xOfI);
		  // sol.print("hyper");
		  BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::membershipSeparation: "
											  "The point is a convex combination of known points! Player "
										  << player;
		  this->Players.at(player)->Feasible = true;

		  arma::vec support;
		  support.zeros(this->Players.at(player)->VertexCounter);
		  for (unsigned int v = 0; v < this->Players.at(player)->VertexCounter; ++v) {
			 // abs to avoid misunderstanding with sign conventions
			 support.at(v) =
				  convexModel.getConstrByName("V_" + std::to_string(v)).get(GRB_DoubleAttr_Pi);
		  }
		  support.print("MNE Support: ");
		  if (support.max() == 1)
			 this->Players.at(player)->Pure = true;
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
		  convexModel.set(GRB_IntParam_SolutionNumber, z);
		  arma::vec cutLHS;
		  cutLHS.zeros(xOfI.size());

		  for (unsigned int i = 0; i < xOfI.size(); i++)
			 cutLHS.at(i) = convexModel.getVarByName("y_" + std::to_string(i)).get(GRB_DoubleAttr_X);
		  cutLHS.print("Separating hyperplane: ");

		  if (cutLHS.max() == 0) {
			 this->Players.at(player)->V.save("dat/V.csv", arma::file_type::csv_ascii);
			 xOfI.save("dat/XofI.csv", arma::file_type::csv_ascii);
		  }
		  // Optimize the resulting inequality over the original feasible set
		  std::unique_ptr<GRBModel> leaderModel =
				std::unique_ptr<GRBModel>(this->IPG->PlayersIP.at(player)->getIPModel(xMinusI, false));
		  GRBLinExpr expr = 0;
		  for (unsigned int i = 0; i < xOfI.size(); ++i)
			 expr += cutLHS.at(i) * leaderModel->getVarByName("y_" + std::to_string(i));

		  leaderModel->setObjective(expr, GRB_MAXIMIZE);
		  leaderModel->update();
		  leaderModel->write("dat/LeaderModel" + std::to_string(player) + ".lp");
		  try {
			 // leaderModel->set(GRB_DoubleAttr_MIPGap, this->Tolerance);
		  } catch (GRBException e) {
			 throw ZEROException(e);
		  }
		  leaderModel->optimize();
		  status = leaderModel->get(GRB_IntAttr_Status);

		  if (status == GRB_OPTIMAL) {
			 double cutV = leaderModel->getObjective().getValue();
			 BOOST_LOG_TRIVIAL(trace)
				  << "Algorithms::IPG::Oracle::membershipSeparation: "
					  "LeaderModel status = "
				  << std::to_string(status) << " with objective=" << cutV << " for Player " << player;
			 arma::vec val  = cutLHS.t() * xOfI; // c^T xOfI
			 arma::vec val2 = cutLHS.t() * this->Players.at(player)->V.row(0).t();
			 BOOST_LOG_TRIVIAL(trace) << "Algorithms::IPG::Oracle::membershipSeparation: c^Tv=" << cutV
											  << " -- c^TxOfI=" << val.at(0) << " -- c^TV(0)=" << val2.at(0);
			 if (cutV - val.at(0) < -this->Tolerance) {
				// False, but we have a cut :-)
				// Ciao Moni
				cutV              = cutV;
				arma::sp_mat cutL = arma::sp_mat{cutLHS}.t();

				if (!this->IPG->PlayersIP.at(player)->addConstraint(
						  cutLHS, cutV, true, this->Tolerance)) {
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
				v.zeros(this->Players.at(player)->V.n_cols);
				for (unsigned int i = 0; i < v.size(); ++i)
				  v[i] = leaderModel->getVarByName("y_" + std::to_string(i)).get(GRB_DoubleAttr_X);

				if (v.max() == v.min() && v.min() == 0) {
				  BOOST_LOG_TRIVIAL(warning) << "Algorithms::IPG::Oracle::membershipSeparation: "
														  "The origin! Feasible point"
													  << player;
				  return true;
				}


				v.print("Vertex found: ");
				std::cout << "Objective: " << leaderModel->getObjective();
				if (Utils::containsRow(this->Players.at(player)->V, v, this->Tolerance)) {
				  BOOST_LOG_TRIVIAL(warning) << "Algorithms::IPG::Oracle::membershipSeparation: "
														  "duplicate vertex for  player "
													  << player;
				  //@todo
				  break;
				  // throw;
				} else {
				  this->Players.at(player)->addVertex(v);
				  v.print("Vertex");
				  BOOST_LOG_TRIVIAL(info) << "Algorithms::IPG::Oracle::membershipSeparation: "
													  "adding vertex for Player. ";
				  break;
				}
			 }

		  } // status optimal for leaderModel
		  else if (status == GRB_UNBOUNDED) {
			 // Check for a new ray
			 if (!Utils::containsRow(this->Players.at(player)->R, cutLHS, this->Tolerance)) {
				BOOST_LOG_TRIVIAL(warning) << "Algorithms::IPG::Oracle::membershipSeparation: "
														"new ray for  player "
													<< player;
				this->Players.at(player)->addRay(cutLHS);
				break;
			 } else {
				BOOST_LOG_TRIVIAL(warning) << "Algorithms::IPG::Oracle::membershipSeparation: "
														"duplicate ray for player "
													<< player;
				break;
			 }

		  } // status unbounded for leaderModel

		  else {
			 throw ZEROException(ZEROErrorCode::Assertion,
										"Unknown status for leaderModel for player " +
											 std::to_string(player));
		  }

		  // GRBModel relaxed = leaderModel->relax();
		  // this->Players.at(player)->updateIPModel(std::move(std::unique_ptr<GRBModel>(&relaxed)));
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
  MathOpt::getDualMembershipLP(this->Players.at(player)->MembershipLP,
										 this->Players.at(player)->VertexCounter,
										 this->Players.at(player)->V,
										 this->Players.at(player)->RayCounter,
										 this->Players.at(player)->R,
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
  BOOST_LOG_TRIVIAL(trace) << "Algorithms::IPG::Oracle::computeStrategy: "
										"Computing startegy for for "
									<< std::to_string(i);
  unsigned int Ny      = this->IPG->PlayerVariables.at(i); // Equals to Ny by definition
  arma::vec    xMinusI = this->buildXminusI(i);
  payoff               = 0;
  strategy.zeros(Ny);
  //  Update working strategies
  std::unique_ptr<GRBModel> response =
		std::unique_ptr<GRBModel>(this->IPG->PlayersIP.at(i)->getIPModel(xMinusI, true));
  // this->Players.at(i)->updateIPModel(std::move(response));
  auto Model = *this->Players.at(i)->Model;
  Model.optimize();
  Model.write("dat/Computing_" + std::to_string(i) + ".lp");
  int status = Model.get(GRB_IntAttr_Status);
  if (status == GRB_INFEASIBLE) {
	 // Game ended, player is infeasible
	 BOOST_LOG_TRIVIAL(trace) << "Algorithms::IPG::Oracle::computeStrategy: "
										  "The game is infeasible! Detected infeasiblity for "
									  << std::to_string(i);
	 this->IPG->Stats.Status.set(ZEROStatus::NashEqNotFound);
	 return false;
  } else if (status == GRB_OPTIMAL) {
	 payoff = Model.getObjective().getValue();
	 BOOST_LOG_TRIVIAL(trace) << "Algorithms::IPG::Oracle::computeStrategy: "
										  "A strategy for player "
									  << i << " has been found (payoff=" << payoff << ")";
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


bool Algorithms::IPG::Oracle::equilibriumLCP(double localTimeLimit) {

  arma::sp_mat MC(0, this->IPG->NumVariables), dumA(0, this->IPG->NumVariables);
  arma::vec    MCRHS(0, arma::fill::zeros), dumB(0, arma::fill::zeros);
  std::vector<std::shared_ptr<MathOpt::MP_Param>> MPCasted;
  for (auto &item : this->IPG->PlayersIP) {
	 auto m = std::dynamic_pointer_cast<MathOpt::MP_Param>(item);
	 MPCasted.push_back(m);
  }
  Game::NashGame Nash = Game::NashGame(this->Env, MPCasted, MC, MCRHS, 0, dumA, dumB);
  BOOST_LOG_TRIVIAL(trace) << "@todo: NashGame is ready";
  auto LCP = std::unique_ptr<MathOpt::LCP>(new MathOpt::LCP(this->Env, Nash));

  arma::vec x, z;
  bool      eq = false;
  if (this->IPG->Stats.AlgorithmData.LCPSolver.get() == Data::IPG::LCPAlgorithms::PATH)
	 eq = LCP->solvePATH(localTimeLimit, z, x);
  else if (this->IPG->Stats.AlgorithmData.LCPSolver.get() == Data::IPG::LCPAlgorithms::MIP) {

	 auto LCPModel = LCP->LCPasMIP(false);
	 if (localTimeLimit > 0) {
		LCPModel->set(GRB_DoubleParam_TimeLimit, localTimeLimit);
	 }
	 LCPModel->set(GRB_IntParam_OutputFlag, 1);
	 LCPModel->setObjective(GRBQuadExpr{0}, GRB_MINIMIZE);
	 LCPModel->set(GRB_IntParam_MIPFocus, 1);
	 LCPModel->optimize();
	 if (LCPModel->get(GRB_IntAttr_Status) != (GRB_SOLUTION_LIMIT || GRB_OPTIMAL))
		return false;
	 eq = LCP->extractSols(LCPModel.get(), z, x, true);
  }
  if (eq) {
	 BOOST_LOG_TRIVIAL(info) << "Game::EPEC::computeNashEq: an Equilibrium has been found";
	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
		this->Players.at(i)->Incumbent = x.subvec(Nash.getPrimalLoc(i), Nash.getPrimalLoc(i + 1) - 1);
		this->Players.at(i)->Incumbent.print("Incumbent of " + std::to_string(i));
		this->Players.at(i)->Feasible = false;
		this->Players.at(i)->Pure     = false;
	 }

	 for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i)
		this->Players.at(i)->Payoff = this->IPG->PlayersIP.at(i)->computeObjective(
			 this->Players.at(i)->Incumbent, this->buildXminusI(i), false);

	 return true;

  } else {
	 BOOST_LOG_TRIVIAL(info) << "Game::EPEC::computeNashEq: NO Equilibrium has been found";
	 return false;
  }
}

void Algorithms::IPG::Oracle::initialize() {
  /**
	* @brief This method initializes some fields for the algorithm. Also, it warm starts the initial
	* strategies to pure best reponses.
	*/
  if (this->IPG->Stats.AlgorithmData.TimeLimit.get() > 0)
	 this->IPG->InitTime = std::chrono::high_resolution_clock::now();
  this->IPG->Stats.NumIterations.set(0);

  this->Players = std::vector<std::unique_ptr<IPG_Player>>(this->IPG->NumPlayers);
  // Initialize the working objects
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 std::unique_ptr<GRBModel> Membership = std::unique_ptr<GRBModel>(new GRBModel(*this->Env));
	 this->Players.at(i)                  = std::unique_ptr<IPG_Player>(
        new IPG_Player(this->IPG->PlayersIP.at(i)->getNy(), this->Tolerance));
	 this->Players.at(i)->MembershipLP = std::move(Membership);
  }



  // Reset the working strategies to a pure strategy given by the IP
  // Push back the IP_Param copies in WorkingIPs
  for (unsigned int i = 0; i < this->IPG->NumPlayers; ++i) {
	 unsigned int Ny      = this->IPG->PlayerVariables.at(i); // Equals to Ny by definition
	 arma::vec    xMinusI = this->buildXminusI(i);

	 // Update working strategies with "educated guesses"
	 auto PureIP = this->IPG->PlayersIP.at(i)->getIPModel(xMinusI, false);
	 PureIP->optimize();
	 int status = PureIP->get(GRB_IntAttr_Status);
	 if (status == GRB_INFEASIBLE) {
		// Game ended, player is infeasible
		this->IPG->Stats.Status.set(ZEROStatus::NashEqNotFound);
		this->Infeasible = true;
		return;
	 } else {
		// Model is not infeasible. We can have either a ray or a vertex
		if (status == GRB_OPTIMAL) {
		  // We have a vertex
		  for (unsigned int k = 0; k < Ny; ++k)
			 this->Players.at(i)->Incumbent.at(k) =
				  PureIP->getVarByName("y_" + std::to_string(k)).get(GRB_DoubleAttr_X);
		  // This is also a f

		  if (this->Players.at(i)->addVertex(this->Players.at(i)->Incumbent, true))
			 BOOST_LOG_TRIVIAL(trace) << "Algorithms::IPG::Oracle::initialize(): "
												  "Added vertex for player "
											  << i;
		}

		else if (status == GRB_UNBOUNDED) {
		  GRBModel relaxed = PureIP->relax();
		  relaxed.set(GRB_IntParam_InfUnbdInfo, 1);
		  relaxed.set(GRB_IntParam_DualReductions, 0);
		  relaxed.optimize();
		  arma::vec ray;
		  for (unsigned int k = 0; k < Ny; ++k) {
			 ray.at(k) = relaxed.getVarByName("y_" + std::to_string(k)).get(GRB_DoubleAttr_UnbdRay);
			 // This is also a free ray
			 if (this->Players.at(i)->addRay(ray, true))
				BOOST_LOG_TRIVIAL(trace) << "Algorithms::IPG::Oracle::initialize(): "
													 "Added ray for player "
												 << i;
		  }
		}
		// Give the new IP
		// this->Players.at(i)->updateIPModel(
		//	 std::move(std::unique_ptr<GRBModel>(new GRBModel(PureIP->relax()))));
	 }
  }
}


arma::vec Algorithms::IPG::Oracle::buildXminusI(const unsigned int i) {
  /**
	* @brief Given the player id @p i, builds the vector x^{-i} from the current working strategies.
	*/
  arma::vec xMinusI;
  xMinusI.zeros(this->IPG->NumVariables - this->IPG->PlayerVariables.at(i));
  unsigned int counter = 0;
  for (unsigned int j = 0; j < this->IPG->NumPlayers; ++j) {
	 if (i != j) {
		xMinusI.subvec(counter, counter + this->IPG->PlayerVariables.at(j) - 1) =
			 this->Players.at(j)->Incumbent;
		counter += this->IPG->PlayerVariables.at(j);
	 }
  }
  return xMinusI;
}