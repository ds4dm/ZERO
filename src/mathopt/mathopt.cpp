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

#include "mathopt/mathopt.h"


/**
 * @brief Computes the convex hull of a finite union of polyhedra where
 * each polyhedra @f$P_i@f$ is of the form
 * @f{eqnarray}{
 * A^ix &\leq& b^i\\
 * x &\geq& 0
 * @f}.
 * It uses Balas' approach to compute the convex hull.
 * @param Ai Inequality constraints LHSs that define polyhedra whose convex hull is to be found
 * @param bi Inequality constraints RHSs that define polyhedra whose convex hull is to be found
 * @param A Pointer to store the output of the convex hull LHS
 * @param b Pointer to store the output of the convex hull RHS
 * @param Acom Any common constraints to all the polyhedra - lhs.
 * @param bcom Any common constraints to ALL the polyhedra - RHS.
 * @return
 */
unsigned int MathOpt::convexHull(const std::vector<arma::sp_mat *> *Ai,
											const std::vector<arma::vec *> *   bi,
											arma::sp_mat &                     A,
											arma::vec &                        b,
											const arma::sp_mat &               Acom,
											const arma::vec &                  bcom) {
  // Count number of polyhedra and the space we are in!
  const unsigned int nPoly{static_cast<unsigned int>(Ai->size())};
  // Error check
  ZEROAssert(nPoly > 0);
  // consider
  const unsigned int nC{static_cast<unsigned int>(Ai->front()->n_cols)};
  const unsigned int nComm{static_cast<unsigned int>(Acom.n_rows)};

  ZEROAssert(!(nComm > 0 && Acom.n_cols != nC));
  ZEROAssert(!(nComm > 0 && nComm != bcom.n_rows));
  ZEROAssert(nPoly == bi->size());

  // Count the number of variables in the convex hull.
  unsigned int nFinCons{0}, nFinVar{0};

  for (unsigned int i = 0; i != nPoly; i++) {

	 ZEROAssert(Ai->at(i)->n_cols == nC);
	 ZEROAssert(Ai->at(i)->n_rows == bi->at(i)->n_rows);
	 nFinCons += Ai->at(i)->n_rows;
  }
  // For common constraint copy
  nFinCons += nPoly * nComm;

  const unsigned int FirstCons = nFinCons;

  // 2nd constraint in Eqn 4.31 of Conforti et al. - twice so we have 2 ineq instead of
  // 1 eq constr
  nFinCons += nC * 2;
  // 3rd constr in Eqn 4.31. Again as two ineq constr.
  nFinCons += 2;
  // Common constraints
  // nFinCons += Acom.n_rows;

  nFinVar = nPoly * nC + nPoly + nC; // All x^i variables + delta variables+ original x variables
  A.zeros(nFinCons, nFinVar);
  b.zeros(nFinCons);
  // A.zeros(nFinCons, nFinVar); b.zeros(nFinCons);
  // Implements the first constraint more efficiently using better constructors
  // for sparse matrix
  MathOpt::compConvSize(A, nFinCons, nFinVar, Ai, bi, Acom, bcom);

  bool printEvery = true;
  if (nPoly > 10)
	 printEvery = false;

  // Counting rows completed
  /****************** SLOW LOOP BEWARE *******************/
  for (unsigned int i = 0; i < nPoly; i++) {

	 if (printEvery || i % 10 == 0)
		LOG_S(3) << "MathOpt::convexHull: Handling Polyhedron " << i + 1 << " out of " << nPoly;
	 for (unsigned int j = 0; j < nC; j++) {
		A.at(FirstCons + 2 * j, nC + (i * nC) + j)     = 1;
		A.at(FirstCons + 2 * j + 1, nC + (i * nC) + j) = -1;
	 }
	 // Third constraint in (4.31)
	 A.at(FirstCons + nC * 2, nC + nPoly * nC + i)     = 1;
	 A.at(FirstCons + nC * 2 + 1, nC + nPoly * nC + i) = -1;
  }
  /****************** SLOW LOOP BEWARE *******************/
  // Second Constraint RHS
  for (unsigned int j = 0; j < nC; j++) {
	 A.at(FirstCons + 2 * j, j)     = -1;
	 A.at(FirstCons + 2 * j + 1, j) = 1;
  }
  // Third Constraint RHS
  b.at(FirstCons + nC * 2)     = 1;
  b.at(FirstCons + nC * 2 + 1) = -1;
  LOG_S(1) << "MathOpt::convexHull: Done";
  return nPoly; ///< Perform increasingly better inner approximations in
  ///< iterations
}

/**
 * @brief Generates the matrix "A" in MathOpt::convexHull using batch
 * insertion constructors.
 * Motivation behind this: Response from
 * armadillo:-https://gitlab.com/conradsnicta/armadillo-code/issues/111
 * @param A The output matrix a
 * @param nFinCons Number of rows in final matrix A
 * @param nFinVar Number of columns in the final matrix A
 * @param Ai Input inequality constraints LHSs defining the polyhedra whose convex hull is to be
 * found
 * @param bi Input inequality  RHSs defining the polyhedra whose convex hull is to be found
 * @param Acom LHS of the common constraints for all polyhedra
 * @param bcom RHS of the common constraints for all polyhedra
 */
void MathOpt::compConvSize(arma::sp_mat &                     A,
									const unsigned int                 nFinCons,
									const unsigned int                 nFinVar,
									const std::vector<arma::sp_mat *> *Ai,
									const std::vector<arma::vec *> *   bi,
									const arma::sp_mat &               Acom,
									const arma::vec &                  bcom) {
  const unsigned int nPoly{static_cast<unsigned int>(Ai->size())};
  const unsigned int nC{static_cast<unsigned int>(Ai->front()->n_cols)};
  unsigned int       N{0}; // Total number of nonzero elements in the final matrix
  const unsigned int numCommon{static_cast<unsigned int>(Acom.n_nonzero + bcom.n_rows)};
  for (unsigned int i = 0; i < nPoly; i++) {
	 N += Ai->at(i)->n_nonzero;
	 N += bi->at(i)->n_rows;
  }
  N += numCommon * nPoly; // The common constraints have to be copied for each polyhedron.

  // Now computed N which is the total number of nonzeros.
  arma::umat locations; // location of nonzeros
  arma::vec  val;       // nonzero values
  locations.zeros(2, N);
  val.zeros(N);

  unsigned int count{0}, rowCount{0}, colCount{nC};
  for (unsigned int i = 0; i < nPoly; i++) {
	 for (auto it = Ai->at(i)->begin(); it != Ai->at(i)->end(); ++it) // First constraint
	 {
		locations(0, count) = rowCount + it.row();
		locations(1, count) = colCount + it.col();
		val(count)          = *it;
		++count;
	 }
	 for (unsigned int j = 0; j < bi->at(i)->n_rows; ++j) // RHS of first constraint
	 {
		locations(0, count) = rowCount + j;
		locations(1, count) = nC + nC * nPoly + i;
		val(count)          = -bi->at(i)->at(j);
		++count;
	 }
	 rowCount += Ai->at(i)->n_rows;

	 // For common constraints
	 for (auto it = Acom.begin(); it != Acom.end(); ++it) // First constraint
	 {
		locations(0, count) = rowCount + it.row();
		locations(1, count) = colCount + it.col();
		val(count)          = *it;
		++count;
	 }
	 for (unsigned int j = 0; j < bcom.n_rows; ++j) // RHS of first constraint
	 {
		locations(0, count) = rowCount + j;
		locations(1, count) = nC + nC * nPoly + i;
		val(count)          = -bcom.at(j);
		++count;
	 }
	 rowCount += Acom.n_rows;

	 colCount += nC;
  }
  A = arma::sp_mat(locations, val, nFinCons, nFinVar);
}


/**
 * @brief 	* @brief Given a vector @p R of rays, and @p V or vertices, builds a model in @p
 * ConvexModel that certifies whether @p vertex belongs to the convex-hull generated by @p
 * V and @p R. In case @p numV and/or @p numR are specified, it just updates the model
 * in @p ConvexModel with the missing vertices and rays. The model is always normalized.
 * @param convexModel  The pointer to the model
 * @param numV The number of vertices in the model
 * @param V The matrix containing vertices (as rows)
 * @param numR The number of rays in the model
 * @param R The matrix containing rays (as rows)
 * @param vertex The vertex to separate
 * @param containsOrigin True if the origin is a feasible vertex
 */
void MathOpt::getDualMembershipLP(std::unique_ptr<GRBModel> &convexModel,
											 unsigned int &             numV,
											 const arma::sp_mat &       V,
											 unsigned int &             numR,
											 const arma::sp_mat &       R,
											 const arma::vec &          vertex,
											 bool                       containsOrigin) {


  ZEROAssert(!(V.n_rows < 1 && R.n_rows < 1));
  ZEROAssert(V.n_cols == vertex.size());

  if (numV == 0 && numR == 0) {
	 // Initialize the model
	 convexModel->reset(true);
	 GRBVar     alpha[V.n_cols];
	 GRBVar     a[V.n_cols + 1];
	 GRBVar     beta;
	 GRBLinExpr expr = 0;
	 for (unsigned int i = 0; i < vertex.size(); i++) {
		alpha[i] = convexModel->addVar(
			 -GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "alpha_" + std::to_string(i));


		a[i] = convexModel->addVar(
			 0, GRB_INFINITY, 0, GRB_CONTINUOUS, "abs(alpha_" + std::to_string(i) + ")");

		// Abs: a[i] = abs(alpha[i])
		convexModel->addConstr(a[i], GRB_GREATER_EQUAL, alpha[i], "Abs_1_alpha_" + std::to_string(i));
		convexModel->addConstr(
			 a[i], GRB_GREATER_EQUAL, -alpha[i], "Abs_2_alpha_" + std::to_string(i));
		expr += a[i];
	 }

	 beta        = convexModel->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "beta");
	 a[V.n_cols] = convexModel->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "abs(beta)");
	 convexModel->addConstr(a[V.n_cols], GRB_GREATER_EQUAL, beta, "Abs_1_beta");
	 convexModel->addConstr(a[V.n_cols], GRB_GREATER_EQUAL, -beta, "Abs_2_beta");

	 // Never had a ray to the normalization!
	 // expr += a[V.n_cols];

	 //   Normalization, to be filled later
	 // convexModel->addConstr(0, GRB_LESS_EQUAL, 1, "Normalization");
	 // Absolute normalization
	 convexModel->addConstr(expr, GRB_EQUAL, 1, "Normalization_Abs");

	 if (containsOrigin)
		convexModel->addConstr(beta, GRB_GREATER_EQUAL, 0, "V_Origin");

	 // Hyperplanes for vertices
	 for (unsigned int i = 0; i < V.n_rows; i++) {
		expr = -beta;
		for (auto j = V.begin_row(i); j != V.end_row(i); ++j)
		  expr += (*j) * alpha[j.col()];
		convexModel->addConstr(expr, GRB_LESS_EQUAL, 0, "V_" + std::to_string(i));
	 }
	 numV = V.n_rows;

	 // Without beta-term for vertices
	 for (unsigned int i = 0; i < R.n_rows; i++) {
		expr = 0;
		for (auto j = R.begin_row(i); j != R.end_row(i); ++j)
		  expr += (*j) * alpha[j.col()];
		convexModel->addConstr(expr, GRB_LESS_EQUAL, 0, "R_" + std::to_string(i));
	 }

	 numR = R.n_rows;

	 // For the eventual Farkas' proof of infeasibility
	 convexModel->set(GRB_IntParam_InfUnbdInfo, 1);
	 convexModel->set(GRB_IntParam_DualReductions, 0);
	 convexModel->set(GRB_IntParam_OutputFlag, 0);
	 convexModel->set(GRB_IntParam_SolutionLimit, 1);
	 LOG_S(2) << "MathOpt::getDualMembershipLP: created model";
  } else {
	 // current number of vertices in the model
	 if (numV < V.n_rows) {
		// Then, we need to update the model by adding new constraints
		GRBLinExpr expr = 0;
		GRBVar     beta = convexModel->getVarByName("beta");
		for (unsigned int i = numV; i < V.n_rows; i++) {
		  expr = -beta;
		  for (auto j = V.begin_row(i); j != V.end_row(i); ++j)
			 expr += (*j) * convexModel->getVarByName("alpha_" + std::to_string(j.col()));

		  convexModel->addConstr(expr, GRB_LESS_EQUAL, 0, "V_" + std::to_string(i));
		}
		numV = V.n_rows;
	 }

	 // current number of rays in the model
	 if (numR < R.n_rows) {
		// Then, we need to update the model by adding new constraints
		GRBLinExpr expr = 0;
		for (unsigned int i = numR; i < R.n_rows; i++) {
		  expr = 0;
		  for (auto j = R.begin_row(i); j != R.end_row(i); ++j)
			 expr += (*j) * convexModel->getVarByName("alpha_" + std::to_string(j.col()));

		  convexModel->addConstr(expr, GRB_LESS_EQUAL, 0, "R_" + std::to_string(i));
		}

		numR = R.n_rows;
	 }
  }
  convexModel->update();

  LOG_S(3) << "MathOpt::getDualMembershipLP: updated model";
  convexModel->update();
}

/**
 * @brief Given a vector @p R of rays, and @p V or vertices, builds a model in @p
 * ConvexModel that certifies whether @p vertex belongs to the convex-hull generated by @p
 * V and @p R. In case @p numV and/or @p numR are specified, it just updates the model
 * in @p ConvexModel with the missing vertices and rays. The model is always normalized.
 * From Chvátal, V., Cook, W. and Espinoza, D., 2013. Local cuts for mixed-integer programming.
 * Mathematical Programming Computation, 5(2), pp.171-200.
 * @param convexModel  The pointer to the model
 * @param numV The number of vertices in the model
 * @param V The matrix containing vertices (as rows)
 * @param numR The number of rays in the model
 * @param R The matrix containing rays (as rows)
 * @param vertex The vertex to separate
 * @param containsOrigin True if the origin is a feasible vertex
 */

void MathOpt::getPrimalMembershipLP(std::unique_ptr<GRBModel> &convexModel,
												unsigned int &             numV,
												const arma::sp_mat &       V,
												unsigned int &             numR,
												const arma::sp_mat &       R,
												const arma::vec &          vertex,
												bool                       containsOrigin) {

  ZEROAssert(!(V.n_rows < 1 && R.n_rows < 1));
  ZEROAssert(V.n_cols == vertex.size());

  if (numV == 0 && numR == 0) {
	 // Initialize the model
	 convexModel->reset(true);
	 GRBVar     lambda[V.n_rows], mu[R.n_rows], s, w[vertex.size()];
	 GRBLinExpr expr = 0;

	 // Convex multipliers
	 for (unsigned int i = 0; i < V.n_rows; i++)
		lambda[i] =
			 convexModel->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "lambda_" + std::to_string(i));

	 // Conic multipliers
	 for (unsigned int i = 0; i < R.n_rows; i++)
		mu[i] = convexModel->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "mu_" + std::to_string(i));

	 // Normalization 1
	 s = convexModel->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "s");
	 // Normalization 2
	 for (unsigned int i = 0; i < vertex.size(); i++)
		w[i] = convexModel->addVar(-1, 1, 0, GRB_CONTINUOUS, "w_" + std::to_string(i));

	 convexModel->setObjective(GRBLinExpr{s}, GRB_MAXIMIZE);

	 // For any dimension of the input point
	 for (unsigned int i = 0; i < vertex.size(); i++) {
		expr = w[i] + s * vertex.at(i); // Reset the expression to s(x*[i]) + w[i]

		// Add term for each vertex
		for (unsigned int j = 0; j < V.n_rows; ++j) {
		  GRBVar vars[] = {lambda[j]};
		  // Remember the minus
		  double coeff[] = {-V.at(j, i)};
		  expr.addTerms(coeff, vars, 1);
		}

		// Add term for each ray
		for (unsigned int j = 0; j < R.n_rows; ++j) {
		  GRBVar vars[] = {mu[j]};
		  // Remember the minus
		  double coeff[] = {-R.at(j, i)};
		  expr.addTerms(coeff, vars, 1);
		}
		convexModel->addConstr(expr, GRB_EQUAL, 0, "PR_Combination_" + std::to_string(i));
	 }


	 // Convex Combination
	 expr = -s;
	 double coeff[V.n_rows];
	 for (auto &el : coeff)
		el = 1;
	 expr.addTerms(coeff, lambda, V.n_rows);
	 convexModel->addConstr(expr, GRB_EQUAL, 0, "Convex_Combination");


	 // Update the working sizes
	 numV = V.n_rows;
	 numR = R.n_rows;

	 // For the eventual Farkas' proof of infeasibility
	 convexModel->set(GRB_IntParam_InfUnbdInfo, 1);
	 convexModel->set(GRB_IntParam_DualReductions, 0);
	 convexModel->set(GRB_IntParam_OutputFlag, 0);
	 convexModel->set(GRB_IntParam_SolutionLimit, 1);
	 LOG_S(2) << "MathOpt::getPrimalMembershipLP: created model";

  } else {
	 // Get the constraints for the point-ray combination
	 std::vector<GRBConstr> PR_Combs;
	 for (unsigned int i = 0; i < vertex.size(); ++i)
		PR_Combs.push_back(convexModel->getConstrByName("PR_Combination_" + std::to_string(i)));

	 // Get the s var
	 GRBVar              s = convexModel->getVarByName("s");
	 std::vector<GRBVar> ss(vertex.size(), s);


	 if (numV < V.n_rows) {
		// Then, we need to update the model by adding terms to PR_Combination and Convex_Combination


		// Add the convex combination for points
		std::vector<GRBConstr> Constrs = PR_Combs;
		Constrs.push_back(convexModel->getConstrByName("Convex_Combination"));

		double coeff[vertex.size() + 1];
		// Last coefficient is for the convex combination
		coeff[vertex.size()] = 1;
		for (unsigned int i = numV; i < V.n_rows; i++) {
		  // Change the Coefficient wrt to the point
		  for (unsigned int j = 0; j < vertex.size(); ++j)
			 coeff[j] = V.at(i, j);

		  convexModel->addVar(0,
									 GRB_INFINITY,
									 0,
									 GRB_CONTINUOUS,
									 vertex.size() + 1,
									 &Constrs[0],
									 coeff,
									 "lambda_" + std::to_string(i));
		}
		// Update the working size of the vertices
		numV = V.n_rows;
	 }

	 // current number of rays in the model
	 if (numR < R.n_rows) {
		// Then, we need to update the model by adding new constraints


		double coeff[vertex.size()];
		// coefficients in the point-ray combination
		for (unsigned int i = numR; i < R.n_rows; i++) {
		  for (unsigned int j = 0; j < vertex.size(); ++j)
			 coeff[j] = R.at(i, j);
		  convexModel->addVar(
				0, GRB_INFINITY, 0, GRB_CONTINUOUS, 1, &PR_Combs[0], coeff, "mu_" + std::to_string(i));
		}

		// Update the working size of the rays
		numR = R.n_rows;
	 }
	 convexModel->chgCoeffs(&PR_Combs[0], &ss[0], &vertex[0], static_cast<int>(vertex.size()));
	 LOG_S(3) << "MathOpt::getPrimalMembershipLP: updated model";
  }
  convexModel->reset();
  convexModel->update();
  LOG_S(3) << "MathOpt::getPrimalMembershipLP: ready.";
}