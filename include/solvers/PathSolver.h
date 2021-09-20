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
extern "C" {
#include <path/MCP_Interface.h>
#include <path/Output_Interface.h>
#include <path/Presolve_Interface.h>
}
namespace Solvers {


  /**
	* @brief This class manages the external solver PATH, for mixed-complementarity problems.
	*/
  class PATH {

  private:
	 /**
	  * @brief This struct manages an LCP where each complementarity is given by @f$x \perp
	  * z=(Mx+q)@f$
	  */
	 typedef struct {
		int n;   ///< Number of x variables, namely complementarities
		int nnz; ///< Number of non-zeros in M

		double *x;  ///< Pointer to x values
		double *lb; ///< n-dimensional array of lower bounds
		double *ub; ///< n-dimensional array of upper bounds

		int *   m_start; ///< Fortran like start of M
		int *   m_len; ///< Fortran like length of M
		int *   m_row; ///< Fortran like row of M
		double *m_data; ///< Fortran like data of M

		double *q; ///< q vector
	 } PATHProblem;


	 ZEROStatus    status = ZEROStatus::Numerical; ///< Status for the solver' instance
	 bool          Filled; ///< Boolean controlling instance's filling
	 PATHProblem   Problem;
	 MCP_Interface PATH_Interface = {
		  this,
		  reinterpret_cast<void (*)(void *, int *, int *)>(this->problem_size),
		  reinterpret_cast<void (*)(void *, int, double *, double *, double *)>(this->bounds),
		  reinterpret_cast<int (*)(void *, int, double *, double *)>(&this->function_evaluation),
		  reinterpret_cast<int (*)(
				void *, int, double *, int, double *, int *, int *, int *, int *, double *)>(
				&this->jacobian_evaluation),
		  static_cast<void (*)(void *)>(this->start),
		  NULL,
		  NULL,
		  NULL,
		  NULL}; ///< The MCP_Interface initializer. See PATH documentation for more
	 Presolve_Interface PATH_Presolve = {
		  NULL,
		  NULL,
		  NULL,
		  NULL,
		  NULL,
		  reinterpret_cast<void (*)(void *, int, int *)>(this->mcp_typ),
		  NULL}; ///< The Presolve_Interface initializer. See PATH documentation for more

	 Output_Interface PATH_Output = {
		  NULL,
		  reinterpret_cast<void (*)(void *, int, char *)>(this->messageCB),
		  NULL}; ///< The Output_Interface initializer. See PATH documentation for more

	 int  CreateLMCP(int    n,
						  int    m_nnz,
						  int    m_i[],
						  int    m_j[],
						  double m_ij[],
						  double q[],
						  double lb[],
						  double ub[],
						  double x[],
						  double z[],
						  int    verbose,
						  double timeLimit);
	 void sort(int rows, int cols, int elements, int *row, int *col, double *data);
	 void C_problem_size(int *n, int *nnz);
	 void C_bounds(int n, double *z, double *lb, double *ub);
	 int  C_function_evaluation(int n, double *x, double *f);
	 int  C_jacobian_evaluation(int     n,
										 double *x,
										 int     wantf,
										 double *f,
										 int *   nnz,
										 int *   col_start,
										 int *   col_len,
										 int *   row,
										 double *data);

  public:
	 PATH() = delete;
	 PATH(const arma::sp_mat &  M,
			const arma::vec &     q,
			const perps &         Compl,
			const VariableBounds &Bounds,
			arma::vec &           z,
			arma::vec  &           x,
			double                timeLimit,
			bool verbose);
	 ZEROStatus getStatus() const { return this->status; } ///< Read-only getter for the status


	 static void *messageCB(void *dat, int mode, char *buf);
	 static void *mcp_typ(void *dat, int nnz, int *typ);



	 static int jacobian_evaluation(void *  dat,
											  int     n,
											  double *x,
											  int     wantf,
											  double *f,
											  int *   nnz,
											  int *   col_start,
											  int *   col_len,
											  int *   row,
											  double *data);



	 static int function_evaluation(void *dat, int n, double *x, double *f);


	 static void start(void *dat);

	 static void bounds(void *dat, int n, double *x, double *lb, double *ub);

	 static void problem_size(void *dat, int *n, int *nnz);
  };
} // namespace Solvers