extern "C" {
#include "License.h"
#include "MCP_Interface.h"
#include "Macros.h"
#include "Options.h"
#include "Output.h"
#include "Output_Interface.h"
#include "Path.h"
#include "PathOptions.h"
}
#include "solvers/PathSolver.h"

/**
 * @brief This function is from PATH. Sorts the data in M.
 * @param rows Row count
 * @param cols Column counts
 * @param elements Element counts
 * @param row Row data pointer
 * @param col Column data pointer
 * @param data M-data pointer
 */
void Solvers::PATH::sort(int rows, int cols, int elements, int *row, int *col, double *data) {
  double *m_data;
  int    *m_start;
  int    *m_len;
  int    *m_row;

  int i = 0, cs = 0, ce = 0;

  m_start = new int[cols + 1];
  m_len   = new int[cols + 1];
  m_row   = new int[elements + 1];
  m_data  = new double[elements + 1];

  for (i = 0; i < cols; i++) {
	 m_len[i] = 0;
  }

  for (i = 0; i < elements; i++) {
	 if ((col[i] < 1) || (col[i] > cols)) {
		throw("column incorrect.\n");
	 }

	 if ((row[i] < 1) || (row[i] > rows)) {
		throw("column incorrect.\n");
	 }

	 m_len[col[i] - 1]++;
  }

  m_start[0] = 0;
  for (i = 1; i < cols; i++) {
	 m_start[i]   = m_start[i - 1] + m_len[i - 1];
	 m_len[i - 1] = 0;
  }
  m_len[i - 1] = 0;

  for (i = 0; i < elements; i++) {
	 cs         = col[i] - 1;
	 ce         = m_start[cs] + m_len[cs];
	 m_row[ce]  = row[i];
	 m_data[ce] = data[i];
	 m_len[cs]++;
  }

  elements = 0;
  for (i = 0; i < cols; i++) {
	 cs = m_start[i];
	 ce = cs + m_len[i];

	 while (cs < ce) {
		row[elements]  = m_row[cs];
		col[elements]  = i + 1;
		data[elements] = m_data[cs];
		elements++;
		cs++;
	 }
  }
}

/**
 * @brief Internal method to create the linear mixed-complemetarity problem.
 * @param n Number of x variables
 * @param m_nnz Number of non-zeros in M
 * @param m_i Row indexes for M non-zeros. Fortran style: start is 1 (not 0)
 * @param m_j Column indexes for M non-zeros. Fortran style: start is 1 (not 0)
 * @param m_ij The data in M corresponding to row @p m_i and column @p m_j
 * @param q The q vector
 * @param lb Vector of lower bounds on x
 * @param ub Vector of upper bounds on x
 * @param x Output vector of x variables
 * @param z Output vector of z equation values
 * @param verbose True if PATH will be verbose
 * @param timeLimit A double timelimit
 * @return The PATH MCP_Termination code
 */
int Solvers::PATH::CreateLMCP(int    n,
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
										double timeLimit) {
  MCP_Termination termination;
  Information     information;

  Output_SetLog(stdout);



  auto o = Options_Create();
  Path_AddOptions(o);
  Options_Default(o);



  double *xSol, *zSol;
  double  dnnz;
  int     i;


  if (verbose == 0) {
	 Options_SetBoolean(o, "output", False);
	 Options_SetBoolean(o, "output_errors", False);
	 Options_SetBoolean(o, "output_crash_iterations", False);
	 Options_SetBoolean(o, "output_major_iterations", False);
	 Options_SetBoolean(o, "output_minor_iterations", False);
	 Options_SetBoolean(o, "output_options", False);
	 Options_SetBoolean(o, "output_warnings", False);
  }
  if (timeLimit > 0)
	 Options_SetDouble(o, "time_limit", timeLimit);


  Options_SetInt(o, "cumulative_iteration_limit", 4e6);
  Options_SetInt(o, "major_iteration_limit", 1e6);
  Options_SetInt(o, "minor_iteration_limit", 1e5);
  Options_SetInt(o, "crash_iteration_limit", 5e2);
  Options_SetInt(o, "nms_memory_size", 1e2);
  Options_SetInt(o, "lemke_rank_deficiency_iterations", 1e2);
  // Very important!
  Options_Set(o, "lemke_start_type advanced");
  Options_SetBoolean(o, "crash_perturb", True);
  Options_SetDouble(o, "proximal_perturbation", 1);
  Options_SetInt(o, "crash_nbchange_limit", 50);
  Options_SetDouble(o, "convergence_tolerance", 1e-6);



  if (n == 0) {
	 Options_Destroy(o);
	 this->status = ZEROStatus::Numerical;
	 return MCP_Solved;
  }

  dnnz = MIN(1.0 * m_nnz, 1.0 * n * n);

  if (verbose != 0) {
	 fprintf(stdout,
				"%d row/cols, %f non-zeros, %3.2f%% dense.\n\n",
				n,
				dnnz,
				100.0 * dnnz / (1.0 * n * n));
  }

  this->Problem.n       = n;
  this->Problem.nnz     = m_nnz;
  this->Problem.x       = x;
  this->Problem.q       = q;
  this->Problem.lb      = lb;
  this->Problem.ub      = ub;
  this->Problem.m_start = new int[n];
  this->Problem.m_len   = new int[n];
  this->Problem.m_row   = new int[m_nnz + 1];
  this->Problem.m_data  = new double[m_nnz + 1];

  sort(n, n, m_nnz, m_i, m_j, m_ij);

  int m_index = 0, m_count = 0;
  for (i = 0; i < n; i++) {
	 this->Problem.m_start[i] = m_count + 1;
	 this->Problem.m_len[i]   = 0;

	 while ((m_index < m_nnz) && (m_j[m_index] <= i + 1)) {
		if (m_ij[m_index] != 0) {
		  this->Problem.m_len[i]++;
		  this->Problem.m_row[m_count]  = m_i[m_index];
		  this->Problem.m_data[m_count] = m_ij[m_index];
		  m_count++;
		}
		m_index++;
	 }
  }
  this->Problem.nnz = m_count;



  License_SetString(
		"2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0");
  auto m = MCP_Create(this->Problem.n, this->Problem.nnz + 1);
  MCP_SetInterface(m, &this->PATH_Interface);
  MCP_SetPresolveInterface(m, &this->PATH_Presolve);
  if (verbose != 0)
	 Output_SetInterface(&this->PATH_Output);
  MCP_Jacobian_Structure_Constant(m, static_cast<Boolean>(true));



  Information info;
  if (verbose != 0)
	 info.generate_output = Output_Log | Output_Status | Output_Listing;
  info.use_start  = False;
  info.use_basics = False;

  try {
	 termination = Path_Solve(m, &info);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::SolverError, "PATHWrapper threw an exception");
  }


  switch (termination) {
  case MCP_Solved: {
	 this->status = ZEROStatus::Solved;
	 xSol         = MCP_GetX(m);
	 zSol         = MCP_GetF(m);
	 for (i = 0; i < n; i++) {
		x[i] = xSol[i];

		z[i] = zSol[i];
	 }
  } break;
  case MCP_TimeLimit:
	 this->status = ZEROStatus::TimeLimit;
	 break;
  case MCP_Infeasible:
	 this->status = ZEROStatus::NotSolved;
  default:
	 this->status = ZEROStatus::Numerical;
  }



  MCP_Destroy(m);
  Options_Destroy(o);
  Path_Destroy();

  return termination;
}

/**
 * @brief Main public method to use the solver.
 * @param M The matrix M in the LCP
 * @param q The vector q in the LCP
 * @param Compl Pairs of complementarities <Eqn, Var>
 * @param Bounds Bounds on variables
 * @param z Output vector of x
 * @param x Output vector of z equation values
 * @param timeLimit A double timelimit
 * @param verbose True if verbosity is required
 */
Solvers::PATH::PATH(const arma::sp_mat   &M,
						  const arma::vec      &q,
						  const perps          &Compl,
						  const VariableBounds &Bounds,
						  arma::vec            &z,
						  arma::vec            &x,
						  double                timeLimit,
						  bool                  verbose) {

  // Note that for path z and x are inverted (z are variables, x equations definitions).
  int n   = 0;
  int nnz = 0;

  // M.save("LCP_M.dat", arma::csv_ascii);
  // q.save("LCP_q.dat",arma::csv_ascii);

  this->Filled = false;
  std::vector<double> _q, _Mij, _lb, _ub, _xsol, _zsol;
  std::vector<int>    _Mi, _Mj, _xmap, _zmap;
  unsigned int        row = 1; // Fortran style, we start from 1



  for (const auto p : Compl) {
	 // For each complementarity
	 // z[p.first] \perp x[p.second]


	 {
		/*if (M.row(p.first).n_nonzero >0) {
		  LOG_S(WARNING) << "Solvers::PathLCP: Empty row associated with z_"
											  << std::to_string(p.first);
		}

		else {
			*/
		int lb = Bounds.at(p.second).first;
		int ub = Bounds.at(p.second).second;

		// if (lb != ub)
		{
		  ++n;
		  // Avoid inserting fixed variables

		  _lb.push_back(lb > 0 ? lb : 0);
		  _ub.push_back(ub >= 0 ? ub : 1e20); // PATH will treat 1e20 as infinite. Note that negative UBs means infinity

		  _q.push_back(q.at(p.first));
		  _xmap.push_back(p.second);
		  _zmap.push_back(p.first);

		  _xsol.push_back(0);
		  _zsol.push_back(0);

		  ++row;
		} // end bounds
	 }   // end empty row
  }     // end while


  // Sparse iterator for M
  for (arma::sp_mat::const_iterator it = M.begin(); it != M.end(); ++it) {
	 _Mi.push_back(it.row() + 1);
	 _Mj.push_back(it.col() + 1);
	 _Mij.push_back(*it);
	 ++nnz;
  }

  int stat = 0;
  /*
  _Mi.push_back(0);
  _Mj.push_back(0);
  _Mij.push_back(0);
  _q.push_back(0);
  _lb.push_back(0);
  _ub.push_back(0);
	*/

  try {
	 stat = this->CreateLMCP(n,
									 nnz,
									 &_Mi[0],
									 &_Mj[0],
									 &_Mij[0],
									 &_q[0],
									 &_lb[0],
									 &_ub[0],
									 &_xsol[0],
									 &_zsol[0],
									 verbose,
									 timeLimit);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::SolverError, "PATH threw an exception");
  }


  if (stat == 1) {
	 LOG_S(1) << "Solvers::PathLCP: Found a solution";
	 z.zeros(M.n_cols);
	 x.zeros(M.n_rows);
	 for (unsigned int i = 0; i < _xsol.size(); ++i) {
		z.at(_xmap.at(i)) = _xsol.at(i);
		x.at(_zmap.at(i)) = _zsol.at(i);
	 }
	 // z.print("z");
	 // x.print("x");
	 this->status = ZEROStatus::NashEqFound;
  } else {
	 LOG_S(1) << "Solvers::PathLCP: No solution found (STATUS=" << std::to_string(status) << ")";
  }
}


/**
 * @brief Assigns bounds to variables. See PATH documentation for more
 * @param n Number of variables
 * @param x Pointer to variables
 * @param lb Lower bounds on variables
 * @param ub Upper bounds on variables
 */
void Solvers::PATH::C_bounds(int n, double *x, double *lb, double *ub) {
  int i;

  for (i = 0; i < n; i++) {
	 x[i]  = this->Problem.x[i];
	 lb[i] = this->Problem.lb[i];
	 ub[i] = this->Problem.ub[i];
	 // std::cout <<"\nx_"<<std::to_string(i)<<" in ["<<std::to_string(lb[i])<<",
	 // "<<std::to_string(ub[i])<<"]";
  }
  // std::cout << "\n done\n";
}


/**
 * @brief Static wrapper for PATH::c_function_evaluation. See PATH documentation for more
 * @param n Number of variables
 * @param x Vector of variables for the lcp. Here we call them
 * @param f The output object for the function value
 * @return Unused. This is a callback function
 */
int Solvers::PATH::C_function_evaluation(int n, double *x, double *f) {
  int    col, colStart, colEnd, row;
  double value;

  for (col = 0; col < n; col++) {
	 f[col] = Problem.q[col];
  }

  for (col = 0; col < n; col++) {
	 value = x[col];

	 if (value != 0) {
		colStart = Problem.m_start[col] - 1;
		colEnd   = colStart + Problem.m_len[col];

		while (colStart < colEnd) {
		  row = Problem.m_row[colStart] - 1;
		  f[row] += Problem.m_data[colStart] * value;
		  colStart++;
		}
	 }
  }

  return 0;
}

/**
 * @brief Evaluates the jacobian at a given point. Eventually, returns the value of the
 * complementarity. See PATH documentation for more
 * @param n Number of variables
 * @param x Vector of variables for the lcp. Here we call them
 * @param wantf True (positive) if the function value is needed
 * @param f The output object for the function value
 * @param nnz Number of non-zeros
 * @param col_start Column start vector
 * @param col_len Column length vector
 * @param row Row vector
 * @param data M_Data vector
 * @return Unused. This is a callback function
 */
int Solvers::PATH::C_jacobian_evaluation(int     n,
													  double *x,
													  int     wantf,
													  double *f,
													  int    *nnz,
													  int    *col_start,
													  int    *col_len,
													  int    *row,
													  double *data) {
  int element;

  if (wantf) {
	 this->C_function_evaluation(n, x, f);
  }

  if (!Filled) {
	 for (element = 0; element < Problem.n; element++) {
		col_start[element] = Problem.m_start[element];
		col_len[element]   = Problem.m_len[element];
	 }

	 for (element = 0; element < Problem.nnz; element++) {
		row[element]  = Problem.m_row[element];
		data[element] = Problem.m_data[element];
	 }

	 Filled = true;
  }

  *nnz = Problem.nnz;
  return 0;
}

/**
 * @brief Presolving type for the variable. See PATH documentation for more
 * @param dat  The data passed by the PATH callback. Usually, it is an instance of the class
 * @param nnz Number of non-zeros in M
 * @param typ The output vector for presolve
 * @return Not used. An instance of this
 */
void *Solvers::PATH::mcp_typ(void *dat, int nnz, int *typ) {
  int i;

  for (i = 0; i < nnz; i++) {
	 typ[i] = PRESOLVE_LINEAR;
  }
  return dat;
}

/**
 * @brief Fills the problem size. See PATH documentation for more
 * @param n Number of variables
 * @param nnz Number of non-zeros in M
 */
void Solvers::PATH::C_problem_size(int *n, int *nnz) {
  *n   = this->Problem.n;
  *nnz = this->Problem.nnz + 1;
}
/**
 * @brief Static wrapper for PATH::c_bounds. See PATH documentation for more
 * @param dat  The data passed by the PATH callback. Usually, it is an instance of the class
 * @param n Number of variables
 * @param x Pointer to variables
 * @param lb Lower bounds on variables
 * @param ub Upper bounds on variables
 */
void Solvers::PATH::bounds(void *dat, int n, double *x, double *lb, double *ub) {
  auto *self = static_cast<Solvers::PATH *>(dat);
  self->C_bounds(n, x, lb, ub);
  return;
}
/**
 * @brief Static wrapper for PATH::c_problem_size. See PATH documentation for more
 * @param dat  The data passed by the PATH callback. Usually, it is an instance of the class
 * @param n Number of variables
 * @param nnz Number of non-zeros in M
 */
void Solvers::PATH::problem_size(void *dat, int *n, int *nnz) {
  auto *self = static_cast<Solvers::PATH *>(dat);
  self->C_problem_size(n, nnz);
}

/**
 * @brief Starting function for PATH. This is called whenever the solver is initialized
 * @param dat  The data passed by the PATH callback. Usually, it is an instance of the class
 */
void Solvers::PATH::start(void *dat) {
  auto *self   = static_cast<Solvers::PATH *>(dat);
  self->Filled = 0;
}

/**
 * @brief Static wrapper for PATH::c_function_evaluation. See PATH documentation for more
 * @param dat  The data passed by the PATH callback. Usually, it is an instance of the class
 * @param n Number of variables
 * @param x Vector of variables for the lcp. Here we call them
 * @param f The output object for the function value
 * @return Unused. This is a callback function
 */
int Solvers::PATH::function_evaluation(void *dat, int n, double *x, double *f) {
  auto *self = static_cast<Solvers::PATH *>(dat);
  return self->C_function_evaluation(n, x, f);
}

/**
 * @brief Static wrapper for PATH::c_jacobian_evaluation. See PATH documentation for more
 * @param dat  The data passed by the PATH callback. Usually, it is an instance of the class
 * @param n Number of variables
 * @param x Vector of variables for the lcp. Here we call them
 * @param wantf True (positive) if the function value is needed
 * @param f The output object for the function value
 * @param nnz Number of non-zeros
 * @param col_start Column start vector
 * @param col_len Column length vector
 * @param row Row vector
 * @param data M_Data vector
 * @return Unused. This is a callback function
 */
int Solvers::PATH::jacobian_evaluation(void   *dat,
													int     n,
													double *x,
													int     wantf,
													double *f,
													int    *nnz,
													int    *col_start,
													int    *col_len,
													int    *row,
													double *data) {
  auto *self = static_cast<Solvers::PATH *>(dat);
  return self->C_jacobian_evaluation(n, x, wantf, f, nnz, col_start, col_len, row, data);
}

/**
 * @brief Message callback for path
 * @param dat The data passed by the PATH callback. Usually, it is an instance of the class
 * @param mode Print-mode. Currently unused, see PATH documentation for more
 * @param buf  The char buffer
 * @return Not used. An instance of this
 */
void *Solvers::PATH::messageCB(void *dat, int mode, char *buf) {
  // std::cout << buf;
  return dat;
}
