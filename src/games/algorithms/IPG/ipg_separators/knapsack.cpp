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

#include "games/algorithms/IPG/ipg_separators/knapsack.h"

/**
 * @brief Given the solution @p x (or an empty vector) this method call the separations routines for
 * the Knapsack problem.
 * @param x The solutions for all players.
 */
void Algorithms::IPG::IPG_Separators::Knapsack::separationMain(const arma::vec &x) {}



/**
 * @brief Computes a minimal cover for player @p player and its solution @p xOfi. Namely, an
 * inequality in the form
 * @f$ \sum_{j\in C}x^p_j \leq |C|-1 @f$ valid for player @p player.
 * The variables in the cover are stored in @p C
 * @param player The player id
 * @param xOfi The solution vector
 * @param C The output indexes of cover variables
 */
 void Algorithms::IPG::IPG_Separators::Knapsack::minimalCover(const unsigned int player,
																				 const arma::vec &  xOfi,
																				 std::vector<int> & C) {

  auto separator = GRBModel(this->Env);
  auto playerVars = this->IPG->PlayerVariables.at(player);
  GRBVar z[playerVars];
}
