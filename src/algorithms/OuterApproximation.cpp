#include "algorithms/outerapproximation.h"

#include <boost/log/trivial.hpp>
#include <chrono>
#include <gurobi_c++.h>
#include <set>
#include <string>

bool Algorithms::OuterApproximation::isSolved(double tol) const {
  return this->Feasible;
}

bool Algorithms::OuterApproximation::isFeasible(bool &addedCuts, double tol) {

  /**
   * Check whether the current outer approximation equilibrium is feasible and
   * solves the outer approximation. Otherwise, add cuts, or generate useful
   * points for next iterations.
   */

  // First, we have a NE from Games::computeNashEq
  if (!this->EPECObject->NashEquilibrium)
    return false;

  // Then, the feasibility is implied also by the deviations
  bool result = {true};
  arma::vec bestResponse;
  arma::vec currentPayoffs =
      this->EPECObject->TheNashGame->computeQPObjectiveValues(
          this->EPECObject->SolutionX, true);
  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i) {
    this->Trees.at(i)->resetFeasibility();
    double val = this->EPECObject->respondSol(bestResponse, i,
                                              this->EPECObject->SolutionX);
    if (val == GRB_INFINITY) {
      BOOST_LOG_TRIVIAL(trace)
          << "Algorithms::OuterApproximation:: Unbounded deviation for " << i;
      addedCuts = false;
      return false;
    }
    // minimization standard
    if (std::abs(currentPayoffs.at(i) - val) > tol) {
      // Discrepancy between payoffs! Need to investigate.
      if ((currentPayoffs.at(i) - val) > tol) {
        // It means the current payoff is more than then optimal response. Then
        // this is not a best response. Theoretically, this cannot happen from
        // an outer approximation. This if case is a warning case then
        //@todo can this happen?

        BOOST_LOG_TRIVIAL(warning) << "Algorithms::OuterApproximation::"
                                      "isFeasible: No best response for Player "
                                   << i;
        BOOST_LOG_TRIVIAL(trace)
            << "Algorithms::OuterApproximation:: " << currentPayoffs.at(i)
            << " vs " << val;
        result = false;
        // throw;
        // throw;
      } else if ((currentPayoffs.at(i) - val) < tol) {
        // It means the current payoff is less than the optimal response. The
        // approximation is not good, and this point is infeasible. Then, we can
        // generate a value-cut
        arma::vec xMinusI;
        this->EPECObject->getXMinusI(this->EPECObject->SolutionX, i, xMinusI);
        this->addValueCut(i, bestResponse, xMinusI);
        // @todo enable cuts
        BOOST_LOG_TRIVIAL(info)
            << "Algorithms::OuterApproximation::isFeasible: "
               "Value cut at for Player "
            << i;
        result = false;
      }
    } else {
      // Here we have a best response whose payoff coincides with the one of the
      // equilibrium. The strategy might not be feasible, though.
      arma::vec xOfI;
      this->EPECObject->getXofI(this->EPECObject->SolutionX, i, xOfI, false);

      // Check if we need to add the point to the vertex storage.
      arma::vec vertex = bestResponse.subvec(0, xOfI.size() - 1);
      vertex.print("Best Response");
      if (!this->Trees.at(i)->containsVertex(vertex, this->Tolerance)) {
        this->Trees.at(i)->addVertex(vertex);
        BOOST_LOG_TRIVIAL(info)
            << "Algorithms::OuterApproximation::isFeasible: "
               "Adding vertex as of best response for Player "
            << i << " (Best Response)";
      } else {
        BOOST_LOG_TRIVIAL(info)
            << "Algorithms::OuterApproximation::isFeasible: "
               "Already known best response for Player "
            << i << " (Best Response)";
      }

      // Check if best response coincides with the strategy in the equilibrium
      bool same = true;
      for (unsigned int k = 0; k < xOfI.size(); ++k) {
        if (std::abs(xOfI.at(k) - bestResponse.at(k)) > tol) {
          same = false;
          break;
        }
      }
      if (!same) {
        // Then, if the answers do not coincide, we need to refine the
        // approximation or determine if this strategy is anyhow feasible.
        // We search for a convex combination of best responses so that we can
        // certify the answer is inside the convex-hull (or not).

        int budget = 15;
        if (!this->separationOracle(xOfI, this->EPECObject->SolutionX, i,
                                    budget, addedCuts)) {
          BOOST_LOG_TRIVIAL(trace)
              << "Algorithms::OuterApproximation::isFeasible: "
                 "Oracle gave a negative answer for Player "
              << i;
          result = false;
        }

      } else {
        this->Trees.at(i)->setFeasible();
        this->Trees.at(i)->setPure();
        BOOST_LOG_TRIVIAL(info)
            << "Algorithms::OuterApproximation::isFeasible: "
               "Feasible strategy for Player "
            << i << " (Best Response)";
      }
    }
  }
  return result;
}

GRBModel *Algorithms::OuterApproximation::getDualMembershipLP(
    unsigned int player, arma::vec vertex, bool normalization) {
  auto convexModel = this->Trees.at(player)->getMembershipLP();
  const arma::sp_mat *V = this->Trees.at(player)->getV();
  const arma::sp_mat *R = this->Trees.at(player)->getR();

  V->print_dense("V");
  R->print_dense("R");

  if (V->n_rows < 1 && R->n_rows < 1) {
    throw "Algorithms::OuterApproximation::getDualMembershipLP: "
          "no points or rays in the membershipLP of " +
        std::to_string(player);
  }
  if (V->n_cols != vertex.size())
    throw "Algorithms::OuterApproximation::getDualMembershipLP: invalid "
          "dimension of vertex.";

  if (!this->Trees.at(player)->getMembershipInit()) {
    // Initialize the model
    GRBVar y[V->n_cols];
    GRBVar z[R->n_cols];
    GRBVar a[V->n_cols + 1];
    GRBVar x;
    GRBLinExpr expr = 0;
    for (unsigned int i = 0; i < vertex.size(); i++) {
      y[i] = convexModel->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS,
                                 "y_" + std::to_string(i));
      a[i] = convexModel->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS,
                                 "abs(y_" + std::to_string(i) + ")");

      convexModel->addConstr(GRBLinExpr{y[i] - a[i]}, GRB_LESS_EQUAL, 0,
                             "Abs_1_y_" + std::to_string(i));
      convexModel->addConstr(GRBLinExpr{-y[i] - a[i]}, GRB_LESS_EQUAL, 0,
                             "Abs_2_y_" + std::to_string(i));
      expr += a[i];
    }

    x = convexModel->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS,
                            "x");
    a[V->n_cols] =
        convexModel->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "abs(x)");
    convexModel->addConstr(GRBLinExpr{x - a[V->n_cols]}, GRB_GREATER_EQUAL, 0,
                           "Abs_1_x");
    convexModel->addConstr(GRBLinExpr{-x - a[V->n_cols]}, GRB_GREATER_EQUAL, 0,
                           "Abs_2_x");
    expr += a[V->n_cols];
    if (normalization)
      convexModel->addConstr(expr, GRB_LESS_EQUAL, 1, "Normalization");

    // Hyperplanes for vertices
    for (unsigned int i = 0; i < V->n_rows; i++) {
      expr = x;
      for (auto j = V->begin_row(i); j != V->end_row(i); ++j)
        expr += (*j) * y[j.col()];
      convexModel->addConstr(expr, GRB_LESS_EQUAL, 0, "V_" + std::to_string(i));
    }
    this->Trees.at(player)->incrementVertices(V->n_rows);

    for (unsigned int i = 0; i < R->n_rows; i++) {
      for (auto j = R->begin_row(i); j != R->end_row(i); ++j)
        expr += (*j) * y[j.col()];
      convexModel->addConstr(expr, GRB_LESS_EQUAL, 0, "R_" + std::to_string(i));
    }

    this->Trees.at(player)->incrementRays(R->n_rows);

    // For the eventual Farkas' proof of infeasibility
    convexModel->set(GRB_IntParam_InfUnbdInfo, 1);
    convexModel->set(GRB_IntParam_DualReductions, 0);
    convexModel->set(GRB_IntParam_OutputFlag, 0);
    convexModel->set(GRB_IntParam_SolutionLimit, 100);
    this->Trees.at(player)->setMembershipInit();
    BOOST_LOG_TRIVIAL(trace)
        << "Algorithms::OuterApproximation::getDualMembershipLP: created model";
  } else {
    // current number of vertices in the model
    if (this->Trees.at(player)->getVertexCount() < V->n_rows) {
      // Then, we need to update the model by adding new constraints
      GRBLinExpr expr = 0;
      for (unsigned int i = this->Trees.at(player)->getVertexCount();
           i < V->n_rows; i++) {
        expr = convexModel->getVarByName("x");
        for (auto j = V->begin_row(i); j != V->end_row(i); ++j)
          expr +=
              (*j) * convexModel->getVarByName("y_" + std::to_string(j.col()));

        convexModel->addConstr(expr, GRB_LESS_EQUAL, 0,
                               "V_" + std::to_string(i));
      }
      this->Trees.at(player)->incrementVertices(
          V->n_rows - this->Trees.at(player)->getVertexCount());
    }

    // current number of rays in the model
    if (this->Trees.at(player)->getRayCount() < R->n_rows) {
      // Then, we need to update the model by adding new constraints
      GRBLinExpr expr = 0;
      for (unsigned int i = this->Trees.at(player)->getRayCount();
           i < R->n_rows; i++) {
        for (auto j = R->begin_row(i); j != R->end_row(i); ++j)
          expr +=
              (*j) * convexModel->getVarByName("y_" + std::to_string(j.col()));

        convexModel->addConstr(expr, GRB_LESS_EQUAL, 0,
                               "R_" + std::to_string(i));
      }

      this->Trees.at(player)->incrementRays(
          R->n_rows - this->Trees.at(player)->getRayCount());
    }

    BOOST_LOG_TRIVIAL(trace)
        << "Algorithms::OuterApproximation::getDualMembershipLP: updated model";
  }
  convexModel->update();
  GRBLinExpr expr = convexModel->getVarByName("x");
  for (int j = 0; j < vertex.size(); ++j)
    expr += vertex.at(j) * convexModel->getVarByName("y_" + std::to_string(j));

  convexModel->setObjective(expr, GRB_MAXIMIZE);
  convexModel->update();
  return convexModel;
}

bool Algorithms::OuterApproximation::separationOracle(arma::vec &xOfI,
                                                      arma::vec &x,
                                                      unsigned int player,
                                                      int budget,
                                                      bool &addedCuts) {

  for (int k = 0; k < budget; ++k) {
    // First, we check whether the point is a convex combination of feasible
    // KNOWN points

    xOfI.print("Point to separate: ");
    const arma::sp_mat *V = this->Trees.at(player)->getV();
    auto convexModel = this->getDualMembershipLP(player, xOfI, true);

    convexModel->write("dat/Convex" + std::to_string(player) + ".lp");
    convexModel->optimize();

    int status = convexModel->get(GRB_IntAttr_Status);
    BOOST_LOG_TRIVIAL(trace)
        << "Algorithms::OuterApproximation::separationOracle: "
           "MermbershipLP status is "
        << status;
    if (status == GRB_OPTIMAL) {
      if (convexModel->getObjective().getValue() == 0 &&
          convexModel->getConstrByName("Normalization")
                  .get(GRB_DoubleAttr_Slack) == 1) {
        // this->Trees.at(player)->addVertex(xOfI);
        BOOST_LOG_TRIVIAL(info)
            << "Algorithms::OuterApproximation::separationOracle: "
               "The point is a convex combination of known points! Player "
            << player;

        this->Trees.at(player)->setFeasible();

        arma::vec support;
        support.zeros(this->Trees.at(player)->getVertexCount());
        auto test = convexModel->getVarByName("x").get(GRB_DoubleAttr_X);
        for (unsigned int v = 0; v < this->Trees.at(player)->getVertexCount();
             ++v) {
          // abs to avoid misunderstanding with sign conventions
          support.at(v) = convexModel->getConstrByName("V_" + std::to_string(v))
                              .get(GRB_DoubleAttr_Pi);
        }
        support.print("MNE Support: ");
        if (support.max() == 1)
          this->Trees.at(player)->setPure();
        return true;
      }
    }

    // Else, the status should be OPTIMAL but without the objective of zero
    if (status == GRB_OPTIMAL) {
      // Get the Farkas' in the form of the unbounded ray of the dual of the
      // dualMembershipLP (the primal)
      BOOST_LOG_TRIVIAL(info)
          << "Algorithms::OuterApproximation::separationOracle: "
             "The point is NOT a convex combination of known points! Found "
          << convexModel->get(GRB_IntAttr_SolCount) << " solutions. Player "
          << player;
      for (int z = 0; z < convexModel->get(GRB_IntAttr_SolCount); ++z) {
        convexModel->getEnv().set(GRB_IntParam_SolutionNumber, z);
        arma::vec cutLHS;
        cutLHS.zeros(xOfI.size());

        for (unsigned int i = 0; i < xOfI.size(); i++)
          cutLHS.at(i) = convexModel->getVarByName("y_" + std::to_string(i))
                             .get(GRB_DoubleAttr_X);
        cutLHS.print("Separating hyperplane: ");

        // Optimize the resulting inequality over the original feasible set
        auto leaderModel = this->EPECObject->respond(player, x);
        GRBLinExpr expr = 0;
        for (unsigned int i = 0; i < xOfI.size(); ++i)
          expr += cutLHS.at(i) *
                  leaderModel->getVarByName("x_" + std::to_string(i));

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
              << "Algorithms::OuterApproximation::separationOracle: "
                 "LeaderModel status = "
              << std::to_string(status) << " with objective=" << cutV
              << " for Player " << player;
          arma::vec val = cutLHS.t() * xOfI; // c^T xOfI
          arma::vec val2 = cutLHS.t() * V->row(0).t();
          BOOST_LOG_TRIVIAL(trace)
              << "Algorithms::OuterApproximation::separationOracle: c^Tv="
              << cutV << " -- c^TxOfI=" << val.at(0)
              << " -- c^TV(0)=" << val2.at(0);
          if (cutV - val.at(0) < -this->Tolerance) {
            // False, but we have a cut :-)
            // Ciao Moni
            cutV = cutV;
            arma::sp_mat cutL =
                Utils::resizePatch(arma::sp_mat{cutLHS}.t(), 1,
                                   this->outerLCP.at(player)->getNumCols());
            if (this->outerLCP.at(player)->containCut(
                    Utils::resizePatch(cutLHS,
                                       this->outerLCP.at(player)->getNumCols()),
                    cutV)) {
              BOOST_LOG_TRIVIAL(info)
                  << "Algorithms::OuterApproximation::separationOracle: "
                     "cut already added for Player "
                  << player;
              // throw;
              break;

            } else {
              this->outerLCP.at(player)->addCustomCuts(cutL, arma::vec{cutV});
              BOOST_LOG_TRIVIAL(info)
                  << "Algorithms::OuterApproximation::separationOracle: "
                     "adding cut for Player "
                  << player;
              addedCuts = true;
              return false;
            }
          } else {
            // We found a new vertex
            arma::vec v;
            v.zeros(V->n_cols);
            for (unsigned int i = 0; i < V->n_cols; ++i) {
              v[i] = leaderModel->getVarByName("x_" + std::to_string(i))
                         .get(GRB_DoubleAttr_X);
            }

            v.print("Vertex found: ");
            if (this->Trees.at(player)->containsVertex(v, this->Tolerance)) {
              BOOST_LOG_TRIVIAL(warning)
                  << "Algorithms::OuterApproximation::separationOracle: "
                     "duplicate vertex for  player "
                  << player;
              //@todo
              break;
              // throw;
            } else {
              this->Trees.at(player)->addVertex(v);
              v.print("Vertex");
              BOOST_LOG_TRIVIAL(info)
                  << "Algorithms::OuterApproximation::separationOracle: "
                     "adding vertex for Player. "
                  << (budget - k - 1) << " iterations left for player "
                  << player;
              break;
            }
          }

        } // status optimal for leaderModel
        else if (status == GRB_UNBOUNDED) {
          // Check for a new ray
          arma::vec normalizedRay = this->normalizeRay(cutLHS);
          if (!this->Trees.at(player)->containsRay(normalizedRay,
                                                   this->Tolerance)) {
            BOOST_LOG_TRIVIAL(warning)
                << "Algorithms::OuterApproximation::separationOracle: "
                   "new ray for  player "
                << player;
            this->Trees.at(player)->addRay(normalizedRay);
            break;
          } else {
            BOOST_LOG_TRIVIAL(warning)
                << "Algorithms::OuterApproximation::separationOracle: "
                   "duplicate ray for player "
                << player;
            break;
          }

        } // status unbounded for leaderModel

        else {
          BOOST_LOG_TRIVIAL(warning)
              << "Algorithms::OuterApproximation::separationOracle: "
                 "unknown status for leaderModel for player "
              << player;
          throw;
        }
      } // end for
      // no separation
    }

    else {
      BOOST_LOG_TRIVIAL(warning)
          << "Algorithms::OuterApproximation::separationOracle: "
             "unknown status for convexrModel for player "
          << player;
      throw;
    }
  }
  return false;
}

arma::vec Algorithms::OuterApproximation::normalizeRay(const arma::vec ray) {
  double max = ray.max();
  double min = std::abs(ray.min());
  double norm = 0;

  if (max > min)
    norm = max;
  else
    norm = min;

  return ray / norm;
}

void Algorithms::OuterApproximation::addValueCut(unsigned int player,
                                                 arma::vec xOfIBestResponse,
                                                 arma::vec xMinusI) {

  double cutRHS = this->EPECObject->PlayersQP.at(player)->computeObjective(
      Utils::resizePatch(xOfIBestResponse,
                         this->EPECObject->PlayersQP.at(player)->getNy(), 1),
      Utils::resizePatch(xMinusI,
                         this->EPECObject->PlayersQP.at(player)->getNx(), 1),
      false);
  arma::vec LHS = this->EPECObject->LeaderObjective.at(player)->c +
                  this->EPECObject->LeaderObjective.at(player)->C * xMinusI;
  arma::sp_mat cutLHS = Utils::resizePatch(
      arma::sp_mat{LHS}.t(), 1, this->outerLCP.at(player)->getNumCols());
  BOOST_LOG_TRIVIAL(info) << "Algorithms::OuterApproximation::addValueCut: "
                             "adding cut for Player "
                          << player;
  this->outerLCP.at(player)->addCustomCuts(-cutLHS, arma::vec{-cutRHS});
}

void Algorithms::OuterApproximation::solve() {
  /**
   * Given the referenced EPEC instance, this method solves it through the outer
   * approximation Algorithm mainly with the hybrid branching rule, with
   * lookahead feature.
   */
  // Set the initial point for all countries as 0 and solve the respective LCPs?
  this->EPECObject->SolutionX.zeros(this->EPECObject->NumVariables);
  bool solved = {false};
  if (this->EPECObject->Stats.AlgorithmParam.TimeLimit > 0)
    this->EPECObject->InitTime = std::chrono::high_resolution_clock::now();

  this->EPECObject->Stats.NumIterations = 0;
  if (this->EPECObject->Stats.AlgorithmParam.TimeLimit > 0)
    this->EPECObject->InitTime = std::chrono::high_resolution_clock::now();

  // Initialize Trees
  this->Trees = std::vector<OuterTree *>(this->EPECObject->NumPlayers, 0);
  this->Incumbent =
      std::vector<OuterTree::Node *>(this->EPECObject->NumPlayers, 0);
  for (unsigned int i = 0; i < this->EPECObject->NumPlayers; i++) {
    Trees.at(i) = new OuterTree(this->outerLCP.at(i)->getNumRows(), this->Env);
    Incumbent.at(i) = Trees.at(i)->getRoot();
  }

  bool branch = true;
  int comp = 0;
  // In this case, branchingLocations is a vector of locations with the length
  // of this->EPECObject->NumPlayers
  std::vector<int> branchingLocations;
  std::vector<long int> branches;
  while (!solved) {
    branchingLocations.clear();
    ++this->EPECObject->Stats.NumIterations;
    BOOST_LOG_TRIVIAL(info)
        << "Algorithms::OuterApproximation::solve: Iteration "
        << std::to_string(this->EPECObject->Stats.NumIterations);

    comp = 0;
    branchingLocations = std::vector<int>(this->EPECObject->NumPlayers, -1);

    if (branch) {
      for (int j = 0; j < this->EPECObject->NumPlayers; ++j) {
        if (Incumbent.at(j)->getCumulativeBranches() ==
            Trees.at(j)->getEncodingSize())
          comp++;
        else {
          if (this->EPECObject->Stats.NumIterations == 1) {
            branchingLocations.at(j) =
                this->getFirstBranchLocation(j, Incumbent.at(j));
          } else {
            branchingLocations.at(j) =
                this->hybridBranching(j, Incumbent.at(j));
          }
        }
      }

      // Check at least a player has at least a branching candidate
      if (comp == this->EPECObject->NumPlayers) {
        BOOST_LOG_TRIVIAL(info) << "Algorithms::OuterApproximation::solve: "
                                   "Solved without any equilibrium.";
        this->EPECObject->Stats.Status = Game::EPECsolveStatus::NashEqNotFound;
        solved = true;
        break;
      }

      // Check that there is at least a player has a branching selection with
      // hybrid branching
      if (*std::max_element(branchingLocations.begin(),
                            branchingLocations.end()) < 0) {

        // No branching candidates.
        BOOST_LOG_TRIVIAL(info) << "Algorithms::OuterApproximation::solve: "
                                   "No more hybrid branching candidates for "
                                   "any player. Checking if "
                                   "any complementarities are left.";
        this->printCurrentApprox();
        for (int j = 0; j < this->EPECObject->NumPlayers; ++j)
          branchingLocations.at(j) =
              this->getFirstBranchLocation(j, Incumbent.at(j));

        if (*std::max_element(branchingLocations.begin(),
                              branchingLocations.end()) < 0) {
          BOOST_LOG_TRIVIAL(info) << "Algorithms::OuterApproximation::solve: "
                                     "No more branching candidates.";
          this->isSolved();
          break;
        }
      }
    }

    for (int j = 0; j < this->EPECObject->NumPlayers; ++j) {
      if (branchingLocations.at(j) > -1) {
        branches = Trees.at(j)->singleBranch(branchingLocations.at(j),
                                             *Incumbent.at(j));
        auto childEncoding =
            this->Trees.at(j)->getNodes()->at(branches.at(0)).getEncoding();
        this->outerLCP.at(j)->outerApproximate(childEncoding, true);
        // By definition of hybrid branching, the node should be feasible
        Incumbent.at(j) = &(this->Trees.at(j)->getNodes()->at(branches.at(0)));
        BOOST_LOG_TRIVIAL(info) << "Algorithms::OuterApproximation::solve: "
                                   "branching candidate for player "
                                << j << " is " << branchingLocations.at(j);
      } else if (!branch) {
        // if we don't branch.
        this->outerLCP.at(j)->outerApproximate(Incumbent.at(j)->getEncoding(),
                                               true);
        BOOST_LOG_TRIVIAL(info) << "Algorithms::OuterApproximation::solve: "
                                   "No branching for player "
                                << j;
      }
    }

    this->printCurrentApprox();
    this->EPECObject->makePlayersQPs();
    // To make computeNashEq skip any feasibility check
    this->Feasible = true;
    if (this->EPECObject->Stats.AlgorithmParam.TimeLimit > 0) {
      const std::chrono::duration<double> timeElapsed =
          std::chrono::high_resolution_clock::now() -
          this->EPECObject->InitTime;
      const double timeRemaining =
          this->EPECObject->Stats.AlgorithmParam.TimeLimit -
          timeElapsed.count();
      this->EPECObject->computeNashEq(
          this->EPECObject->Stats.AlgorithmParam.PureNashEquilibrium,
          timeRemaining);
    } else {
      this->EPECObject->computeNashEq(
          this->EPECObject->Stats.AlgorithmParam.PureNashEquilibrium);
    }

    this->Feasible = false;
    if (this->EPECObject->NashEquilibrium) {
      bool addedCuts{false};
      if (this->isFeasible(addedCuts)) {
        this->Feasible = true;
        this->EPECObject->Stats.Status = Game::EPECsolveStatus::NashEqFound;
        BOOST_LOG_TRIVIAL(info) << "Algorithms::OuterApproximation::solve: "
                                   "Solved. ";
        return;
      } else {
        if (addedCuts) {
          branch = false;
          BOOST_LOG_TRIVIAL(info)
              << "Algorithms::OuterApproximation::solve: "
                 "Cuts were added. Skipping next branching phase. ";
        } else {
          branch = true;
        }
      }
    } else {
      branch = true;
    }
    if (this->EPECObject->Stats.AlgorithmParam.TimeLimit > 0) {
      const std::chrono::duration<double> timeElapsed =
          std::chrono::high_resolution_clock::now() -
          this->EPECObject->InitTime;
      const double timeRemaining =
          this->EPECObject->Stats.AlgorithmParam.TimeLimit -
          timeElapsed.count();
      if (timeRemaining <= 0) {
        this->EPECObject->Stats.Status = Game::EPECsolveStatus::TimeLimit;
        return;
      }
    }
  }
}

std::unique_ptr<GRBModel>
Algorithms::OuterApproximation::getFeasQP(const unsigned int player,
                                          const arma::vec x) {
  /**
   * Given @p player -- containing the id of the player, returns the feasibility
   *QP associated with the current QP made and the point @p x. In other words,
   *it returns the model deciding whether @p x belongs to the feasible region of
   *the QP or not. Note that Game::make_player_QP(player) has to be called
   *before, otherwise the method will return an error.
   *
   **/

  // this->EPECObject->getXMinusI(this->EPECObject->SolutionX, player, xMinusI);
  arma::vec zeros;
  // Dummy vector of zeros associated to x^{-i}
  zeros.zeros(this->EPECObject->PlayersQP.at(player)->getNx());
  auto model = this->EPECObject->PlayersQP.at(player)->solveFixed(zeros, false);
  // Enforce QP::y to be x, namely the point to belong to the feasible region
  for (unsigned int j = 0; j < x.size(); j++)
    model->addConstr(model->getVarByName("y_" + std::to_string(j)), GRB_EQUAL,
                     x.at(j), "Fix_y_" + std::to_string(j));
  // Reset the objective
  model->setObjective(GRBLinExpr{0}, GRB_MINIMIZE);
  // model->write("dat/test.lp");
  return model;
}

int Algorithms::OuterApproximation::hybridBranching(const unsigned int player,
                                                    OuterTree::Node *node) {
  /**
   * Given @p player -- containing the id of the player, returns the branching
   * decision for that node given by a hybrid branching rule. In
   * particular, the method return the complementarity id maximizing a
   *combination of constraint violations and number of violated constraints.
   *@p node contains the tree's node. It isn't const since a branching candidate
   *can be pruned if infeasibility is detected
   *
   **/

  BOOST_LOG_TRIVIAL(info) << "OuterApproximation::hybridBranching: Player "
                          << player;

  int bestId = -1;
  if (this->EPECObject->NashEquilibrium) {
    arma::vec zeros, x;

    this->EPECObject->getXofI(this->EPECObject->SolutionX, player, x);
    if (x.size() != this->EPECObject->LeaderObjective.at(player)->c.n_rows)
      throw "Algorithms::OuterApproximation::hybridBranching: "
            "wrong-dimensioned "
            "x^i";

    auto currentEncoding = node->getEncoding();
    std::vector<bool> incumbentApproximation;
    double bestScore = -1.0;

    for (unsigned int i = 0; i < currentEncoding.size(); i++) {
      // For each complementarity
      if (node->getAllowedBranchings().at(i)) {
        // Consider it if it is a good candidate for branching (namely, we
        // didn't branch on it, or it wasn't proven to be infeasible)
        incumbentApproximation = currentEncoding;
        // Include this complementarity in the approximation
        incumbentApproximation.at(i) = true;
        // Build the approximation
        this->outerLCP.at(player)->outerApproximate(incumbentApproximation,
                                                    true);
        // If the approximation is infeasible, prune this branching location
        // from the candidates
        if (!this->outerLCP.at(player)->getFeasApprox())
          Trees.at(player)->denyBranchingLocation(*node, i);
        else {
          // In this case, we can check if the solution belongs to the outer
          // approximation
          this->EPECObject->makePlayerQP(player);
          // Get the QP model with other players decision QP::x fixed to zero
          // (since they only appear in the objective);
          auto model = this->getFeasQP(player, x);
          model->optimize();
          const int status = model->get(GRB_IntAttr_Status);
          if (status == GRB_INFEASIBLE) {
            // If the status is infeasible, bingo! We want to get a measure of
            // the constraint violations given by the current x
            model->feasRelax(0, false, false, true);
            model->optimize();
            if (model->getObjective().getValue() > bestScore) {
              bestId = i;
              bestScore = model->getObjective().getValue();
              BOOST_LOG_TRIVIAL(debug)
                  << "OuterApproximation::hybridBranching: Player " << player
                  << " has violation of " << bestScore
                  << " with complementarity " << i;
            }
          } else {
            BOOST_LOG_TRIVIAL(debug)
                << "OuterApproximation::hybridBranching: Player " << player
                << " has no violation with complementarity " << i;
          }
        }
      }
    }
  }
  return bestId;
}

int Algorithms::OuterApproximation::infeasibleBranching(
    const unsigned int player, const OuterTree::Node *node) {
  /**
   * Given @p player -- containing the id of the player, returns the branching
   * decision for that node, where the complementarity is the most (possibly)
   * infeasible one (with both x and z positive). In particular, the method
   * return the (positive) id of the complementarity equation if there is a
   * feasible branching decision at @p node, and a negative value otherwise.
   * @return a positive int with the id of the complementarity to branch on, or
   * a negative value if none exists.
   */

  int result = -1;
  if (this->EPECObject->NashEquilibrium) {
    // There exists a Nash Equilibrium for the outer approximation, which is not
    // a Nash Equilibrium for the game
    arma::vec x, z;
    this->EPECObject->getXWithoutHull(this->EPECObject->SolutionX, x);
    z = this->outerLCP.at(player)->zFromX(x);
    std::vector<short int> currentSolution =
        this->outerLCP.at(player)->solEncode(x);

    double maxInfeas = 0;

    //"The most infeasible" branching
    for (unsigned int i = 0; i < currentSolution.size(); i++) {
      unsigned int varPos =
          i >= this->outerLCP.at(player)->getLStart()
              ? i + this->outerLCP.at(player)->getNumberLeader()
              : i;
      if (x(varPos) > 0 && z(i) > 0 && node->getAllowedBranchings().at(i) &&
          currentSolution.at(i) == 0) {
        if ((x(varPos) + z(i)) > maxInfeas) {
          maxInfeas = x(varPos) + z(i);
          result = i;
        }
      }
    }
  }
  return result;
}

int Algorithms::OuterApproximation::deviationBranching(
    const unsigned int player, const OuterTree::Node *node) {
  /**
   * Given @p player -- containing the id of the player, returns the branching
   * decision for that node, where the complementarity helps to cut off a
   * possible deviation . In particular, the method return the (positive) id of
   * the complementarity equation if there is a feasible branching decision at
   * @p node, and a negative value otherwise.
   * @return a positive int with the id of the complementarity to branch on, or
   * a negative value if none exists.
   */

  int result = -1;
  if (this->EPECObject->NashEquilibrium) {
    // There exists a Nash Equilibrium for the outer approximation, which is not
    // a Nash Equilibrium for the game
    arma::vec dev;
    arma::vec x;
    this->EPECObject->getXWithoutHull(this->EPECObject->SolutionX, x);
    std::vector<short int> currentSolution =
        this->outerLCP.at(player)->solEncode(x);
    this->EPECObject->respondSol(dev, player, this->EPECObject->SolutionX);
    auto encoding = this->outerLCP.at(player)->solEncode(dev);

    for (unsigned int i = 0; i < encoding.size(); i++) {
      if (encoding.at(i) > 0 && node->getAllowedBranchings().at(i) &&
          currentSolution.at(i) == 0) {
        result = i;
      }
    }
  }
  return result;
}

int Algorithms::OuterApproximation::getFirstBranchLocation(
    const unsigned int player, const OuterTree::Node *node) {
  /**
   * Given @p player -- containing the id of the player, returns the branching
   * decision for that node, with no complementarity condition enforced. In
   * particular, the method return the (positive) id of the complementarity
   * equation if there is a feasible branching decision at @p node, and a
   * negative value otherwise.
   * @return a positive int with the id of the complementarity to branch on, or
   * a negative value if none exists.
   */

  if (node->getCumulativeBranches() == Trees.at(player)->getEncodingSize())
    return -1;
  auto model = this->outerLCP.at(player)->LCPasMIP(true);
  unsigned int nR = this->outerLCP.at(player)->getNumRows();
  int pos = -nR;
  arma::vec z, x;
  if (this->outerLCP.at(player)->extractSols(
          model.get(), z, x, true)) // If already infeasible, nothing to branch!
  {
    std::vector<short int> v1 = this->outerLCP.at(player)->solEncode(z, x);

    double maxvalx{-1}, maxvalz{-1};
    unsigned int maxposx{0}, maxposz{0};
    for (unsigned int i = 0; i < nR; i++) {
      unsigned int varPos =
          i >= this->outerLCP.at(player)->getLStart()
              ? i + this->outerLCP.at(player)->getNumberLeader()
              : i;
      if (x(varPos) > maxvalx && node->getAllowedBranchings().at(i)) {
        maxvalx = x(varPos);
        maxposx = i;
      }
      if (z(i) > maxvalz && node->getAllowedBranchings().at(i)) {
        maxvalz = z(i);
        maxposz = i;
      }
    }
    pos = maxvalz > maxvalx ? maxposz : maxposx;
  } else {
    BOOST_LOG_TRIVIAL(debug) << "The problem is infeasible";
  }
  return pos;
}

std::vector<int>
Algorithms::OuterApproximation::getNextBranchLocation(const unsigned int player,
                                                      OuterTree::Node *node) {
  /**
   * Given @p player -- containing the id of the player -- and @p node
   * containing a node, returns the branching decision for that node, with
   * respect to the current node. In particular, the method return the
   * (positive) id of the complementarity equation if there is a feasible
   * branching decision at @p node, and a negative value otherwise.
   * @return a vector of 4 ints with the branching location given by the most
   * infeasible branching, deviation branching, hybrid branching, and the
   * firstlocation branching, respectively
   */
  std::vector<int> decisions = {-1, -1, -1, -1};
  decisions.at(0) = this->infeasibleBranching(player, node);
  decisions.at(1) = this->deviationBranching(player, node);
  decisions.at(2) = this->hybridBranching(player, node);

  if (decisions.at(0) < 0 && decisions.at(1) < 0 && decisions.at(2) < 0) {
    BOOST_LOG_TRIVIAL(info)
        << "Player " << player
        << ": branching with FirstBranchLocation is the only available choice";
    decisions.at(3) = this->getFirstBranchLocation(player, node);
  }

  BOOST_LOG_TRIVIAL(debug)
      << "Algorithms::OuterApproximation::getNextBranchinglocation: "
         "given decisions are: ";
  BOOST_LOG_TRIVIAL(debug) << "Algorithms::OuterApproximation::"
                              "getNextBranchinglocation:\t Infeasible="
                           << decisions.at(0);
  BOOST_LOG_TRIVIAL(debug) << "Algorithms::OuterApproximation::"
                              "getNextBranchinglocation:\t Deviation="
                           << decisions.at(1);
  BOOST_LOG_TRIVIAL(debug)
      << "Algorithms::OuterApproximation::getNextBranchinglocation:\t Hybrid="
      << decisions.at(2);
  BOOST_LOG_TRIVIAL(debug)
      << "Algorithms::OuterApproximation::getNextBranchinglocation:\t First="
      << decisions.at(3);
  return decisions;
}

void Algorithms::OuterApproximation::printCurrentApprox() {
  /**
   * Returns a log message containing the encoding at the current outer
   * approximation iteration
   */
  BOOST_LOG_TRIVIAL(info) << "Current Node Approximation:";
  for (unsigned int p = 0; p < this->EPECObject->NumPlayers; ++p) {
    std::stringstream msg;
    msg << "\tPlayer " << p << ":";
    for (unsigned int i = 0; i < this->Incumbent.at(p)->getEncoding().size();
         i++) {
      msg << "\t" << this->Incumbent.at(p)->getEncoding().at(i);
    }
    BOOST_LOG_TRIVIAL(info) << msg.str();
  }
}
void Algorithms::OuterApproximation::printBranchingLog(
    std::vector<int> vector) {
  /**
   * Given the vector of branching logs, prints a sum up of the decision taken
   */
  BOOST_LOG_TRIVIAL(info) << "Current Branching Log:";
  BOOST_LOG_TRIVIAL(info) << "\tInfeasibleBranching: " << vector.at(0);
  BOOST_LOG_TRIVIAL(info) << "\tDeviationBranching: " << vector.at(1);
  BOOST_LOG_TRIVIAL(info) << "\tHybridBranching: " << vector.at(2);
  BOOST_LOG_TRIVIAL(info) << "\tFirstAvail: " << vector.at(3);
}
bool Algorithms::OuterApproximation::isPureStrategy(double tol) const {
  if (!this->Feasible)
    return false;
  else {
    for (unsigned int i = 0; i < this->EPECObject->NumPlayers; ++i)
      if (!Trees.at(i)->getPure())
        return false;

    return true;
  }
}

OuterTree::Node::Node(Node &parent, unsigned int idComp, unsigned long int id) {
  /**
   * Given the parent node address @param parent, the @param idComp to branch
   * on, and the @param id, creates a new node
   */
  this->IdComps = std::vector<unsigned int>{idComp};
  this->Encoding = parent.Encoding;
  this->Encoding.at(idComp) = true;
  this->AllowedBranchings = parent.AllowedBranchings;
  this->AllowedBranchings.at(idComp) = false;
  this->Id = id;
  this->Parent = &parent;
}

OuterTree::Node::Node(unsigned int encSize) {
  /**
   * Constructor for the root node, given the encoding size, namely the number
   * of complementarity equations
   */
  this->Encoding = std::vector<bool>(encSize, 0);
  this->Id = 0;
  this->AllowedBranchings = std::vector<bool>(encSize, true);
}

void OuterTree::denyBranchingLocation(OuterTree::Node &node,
                                      const unsigned int &location) {
  /**
   * If a complementarity equation @param location  has proven to be infeasible
   * or it isn't a candidate for branching, this method prevents any further
   * branching on it for the node @param node.
   */
  if (location >= this->EncodingSize)
    throw "OuteTree::branch idComp is larger than the encoding size.";
  if (!node.AllowedBranchings.at(location))
    BOOST_LOG_TRIVIAL(warning) << "OuterTree::denyBranchingLocation: location "
                                  "has been already denied.";
  node.AllowedBranchings.at(location) = false;
}

void OuterTree::denyBranchingLocations(OuterTree::Node &node,
                                       const std::vector<int> &locations) {
  /**
   * If a complementarity equation @param location  has proven to be infeasible
   * or it isn't a candidate for branching, this method prevents any further
   * branching on it for the node @param node.
   */
  for (auto &location : locations) {
    if (location < 0)
      throw "OuterTree::denyBranchingLocations a location is negative.";
    this->denyBranchingLocation(node, location);
  }
}

std::vector<long int> OuterTree::singleBranch(const unsigned int idComp,
                                              OuterTree::Node &t) {
  /**
   * Given the @param idComp and the parent node @param t, creates a single
   * child by branching on @param idComp.
   */
  if (idComp >= this->EncodingSize)
    throw "OuterTree::branch idComp is larger than the encoding size.";
  if (t.Encoding.at(idComp) != 0) {
    BOOST_LOG_TRIVIAL(warning)
        << "OuterTree: cannot branch on this complementary, since it already "
           "has been processed.";
    return std::vector<long int>{-1};
  }
  auto child = Node(t, idComp, this->nextIdentifier());

  this->Nodes.push_back(child);
  return std::vector<long int>{this->NodeCounter - 1};
}

std::vector<long int> OuterTree::multipleBranch(const std::vector<int> idsComp,
                                                Node &t) {
  /**
   * Given the @param idComp and the parent node @param t, creates a single
   * child by branching on @param idComp.
   */
  for (auto &idComp : idsComp) {
    if (idComp >= this->EncodingSize)
      throw "Tree::branch idComp is larger than the encoding size.";
    if (t.Encoding.at(idComp) != 0) {
      BOOST_LOG_TRIVIAL(warning)
          << "Tree: cannot branch on this complementary, since it already has "
             "been processed.";
      return std::vector<long int>{-1};
    }
  }
  auto child = Node(t, idsComp, this->nextIdentifier());

  this->Nodes.push_back(child);
  return std::vector<long int>{this->NodeCounter - 1};
}

OuterTree::Node::Node(Node &parent, std::vector<int> idsComp,
                      unsigned long int id) {
  /**
   * Given the parent node address @param parent, the @param idsComp to branch
   * on (containing all the complementarities ids), and the @param id, creates a
   * new node
   */
  this->IdComps = std::vector<unsigned int>();
  this->Encoding = parent.Encoding;
  this->AllowedBranchings = parent.AllowedBranchings;
  for (auto &idComp : idsComp) {
    if (idComp < 0)
      throw "OuterTree::Node::Node  idComp is negative.";
    this->Encoding.at(idComp) = true;
    this->AllowedBranchings.at(idComp) = false;
    this->IdComps.push_back(idComp);
  }
  this->Id = id;
  this->Parent = &parent;
}
