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
extern "C" {
#include <path/MCP_Interface.h>
#include <path/Output_Interface.h>
#include <path/Presolve_Interface.h>
}
namespace Solvers {

  class PATH {

  private:
	 typedef struct {
		int n;
		int nnz;

		double *z;
		double *lb;
		double *ub;

		int *   m_start;
		int *   m_len;
		int *   m_row;
		double *m_data;

		double *q;
	 } PATHProblem;


	 ZEROStatus    status;
	 int           Filled;
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
		  NULL};
	 Presolve_Interface PATH_Presolve = {
		  NULL,
		  NULL,
		  NULL,
		  NULL,
		  NULL,
		  reinterpret_cast<void (*)(void *, int, int *)>(this->mcp_typ),
		  NULL};

	 Output_Interface PATH_Output = {
		  NULL, reinterpret_cast<void (*)(void *, int, char *)>(this->messageCB), NULL};

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
	 int  C_function_evaluation(int n, double *z, double *f);
	 int  C_jacobian_evaluation(int     n,
										 double *z,
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
			arma::vec &           x,
			arma::vec             z,
			double                timeLimit);
	 ZEROStatus getStatus() const { return this->status; }

	 static void *messageCB(void *dat, int mode, char *buf) {
		auto *self = static_cast<Solvers::PATH *>(dat);
		fprintf(stdout, "%s", buf);
		return dat;
	 }
	 static void *mcp_typ(void *dat, int nnz, int *typ);
	 static int   jacobian_evaluation(void *  dat,
												 int     n,
												 double *z,
												 int     wantf,
												 double *f,
												 int *   nnz,
												 int *   col_start,
												 int *   col_len,
												 int *   row,
												 double *data) {
      auto *self = static_cast<Solvers::PATH *>(dat);
      return self->C_jacobian_evaluation(n, z, wantf, f, nnz, col_start, col_len, row, data);
	 }

	 static int function_evaluation(void *dat, int n, double *z, double *f) {
		auto *self = static_cast<Solvers::PATH *>(dat);
		return self->C_function_evaluation(n, z, f);
	 }
	 static void start(void *dat) {
		auto *self   = static_cast<Solvers::PATH *>(dat);
		self->Filled = 0;
	 }

	 static void bounds(void *dat, int n, double *z, double *lb, double *ub);
	 static void problem_size(void *dat, int *n, int *nnz) {
		auto *self = static_cast<Solvers::PATH *>(dat);
		self->C_problem_size(n, nnz);
	 }
  };
} // namespace Solvers