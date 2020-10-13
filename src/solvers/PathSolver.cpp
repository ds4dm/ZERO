
#include "solvers/PathSolver.h"
#include "path/MCP_Interface.h"
#include "path/Macros.h"
#include "path/Options.h"
#include "path/Output.h"
#include "path/Output_Interface.h"
#include "path/Path.h"
#include "path/PathOptions.h"
#include "path/Types.h"
#include "support/codes.h"
#include <boost/function.hpp>
#include <boost/log/trivial.hpp>



void Solvers::PATH::sort(int rows, int cols, int elements, int *row, int *col, double *data) {
  double *m_data;
  int *   m_start;
  int *   m_len;
  int *   m_row;

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

  return;
}



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

  Output_Printf(
		Output_Log | Output_Status | Output_Listing, "%s: PathWrapper LCP Link\n", Path_Version());


  if (timeLimit > 0)
	 Options_SetDouble(o, "time_limit", timeLimit);


  if (n == 0) {
	 fprintf(stdout, "\n ** EXIT - No variables.\n");
	 Options_Destroy(o);
	 return MCP_Solved;
  }

  dnnz = MIN(1.0 * m_nnz, 1.0 * n * n);
  if (dnnz > INT_MAX) {
	 Output_Printf(Output_Log | Output_Status, "\n ** EXIT - model too large.\n");
	 Options_Destroy(o);
	 return MCP_Error;
  }

  fprintf(stdout,
			 "%d row/cols, %f non-zeros, %3.2f%% dense.\n\n",
			 n,
			 dnnz,
			 100.0 * dnnz / (1.0 * n * n));

  this->Problem.n       = n;
  this->Problem.nnz     = m_nnz;
  this->Problem.z       = x;
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



  auto m = MCP_Create(this->Problem.n, this->Problem.nnz + 1);
  MCP_SetInterface(m, &this->PATH_Interface);
  MCP_SetPresolveInterface(m, &this->PATH_Presolve);
  Output_SetInterface(&this->PATH_Output);
  MCP_Jacobian_Structure_Constant(m, static_cast<Boolean>(true));



  Information info;
  info.generate_output = Output_Log | Output_Status | Output_Listing;
  info.use_start       = True;
  info.use_basics      = True;

  try {
	 termination = Path_Solve(m, &info);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::SolverError, "PATHWrapper threw an exception");
  }


  if (termination == MCP_Solved) {
	 xSol = MCP_GetX(m);
	 zSol = MCP_GetF(m);
	 fprintf(stdout, "%d is  %f \n", n, xSol[i]);
	 for (i = 0; i < n; i++) {
		x[i] = xSol[i];

		z[i] = zSol[i];
	 }
  }


  MCP_Destroy(m);
  Options_Destroy(o);
  Path_Destroy();

  return termination;
}

Solvers::PATH::PATH(const arma::sp_mat &  M,
						  const arma::vec &     q,
						  const perps &         Compl,
						  const VariableBounds &Bounds,
						  arma::vec &           x,
						  arma::vec             z,
						  double                timeLimit) {
  int n   = 0;
  int nnz = 0;


  this->Filled = false;
  std::vector<double> _q, _Mij, _lb, _ub, _xsol, _zsol;
  std::vector<int>    _Mi, _Mj, _xmap, _zmap;
  unsigned int        row = 1; // Fortran style, we start from 1



  for (const auto p : Compl) {
	 // For each complementarity
	 // z[p.first] \perp x[p.second]


	 {
		/*if (M.row(p.first).n_nonzero >0) {
		  BOOST_LOG_TRIVIAL(warning) << "Solvers::PathLCP: Empty row associated with z_"
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
		  _ub.push_back(ub >= 0 ? ub : 1e20); // PATH will treat 1e20 as infinite


		  // std::cout << "Row" << std::to_string(row);
		  for (auto v = M.begin_row(p.first); v != M.end_row(p.first); ++v) {
			 if (*v != 0) {
				_Mi.push_back(row);
				_Mj.push_back(v.col() + 1);
				_Mij.push_back(*v);
				std::cout << "\t" + std::to_string(*v) + "*x_" + std::to_string(v.col());
				++nnz;
			 }
		  }
		  _q.push_back(q.at(p.first));
		  // std::cout << "\t+" + std::to_string(q.at(p.first)) + "\t\tPERP" <<
		  // std::to_string(p.second)
		  //			<< "\n\n";
		  _xmap.push_back(p.second);
		  _zmap.push_back(p.first);

		  _xsol.push_back(0);
		  _zsol.push_back(0);

		  ++row;
		} // end bounds
	 }   // end empty row
  }     // end while


  int stat = 0;
  _Mi.push_back(0);
  _Mj.push_back(0);
  _Mij.push_back(0);
  _q.push_back(0);
  _lb.push_back(0);
  _ub.push_back(0);

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
									 true,
									 timeLimit);
  } catch (...) {
	 throw ZEROException(ZEROErrorCode::SolverError, "PATH threw an exception");
  }


  if (stat == 1) {
	 BOOST_LOG_TRIVIAL(trace) << "Solvers::PathLCP: Found a solution";
	 x.zeros(M.n_cols);
	 z.zeros(M.n_rows);
	 for (unsigned int i = 0; i < _xsol.size(); ++i) {
		x.at(_xmap.at(i)) = _xsol.at(i);
		z.at(_zmap.at(i)) = _zsol.at(i);
	 }
	 // z.print("z");
	 // x.print("x");
	 this->status = ZEROStatus::NashEqFound;
  } else {
	 BOOST_LOG_TRIVIAL(trace) << "Solvers::PathLCP: No solution found (STATUS="
									  << std::to_string(status) << ")";
	 this->status = ZEROStatus::NashEqNotFound;
  }
}
void Solvers::PATH::C_bounds(int n, double *z, double *lb, double *ub) {
  int i;

  for (i = 0; i < n; i++) {
	 z[i]  = this->Problem.z[i];
	 lb[i] = this->Problem.lb[i];
	 ub[i] = this->Problem.ub[i];
	 // std::cout <<"\nx_"<<std::to_string(i)<<" in ["<<std::to_string(lb[i])<<",
	 // "<<std::to_string(ub[i])<<"]";
  }
  // std::cout << "\n done\n";
}
int Solvers::PATH::C_function_evaluation(int n, double *z, double *f) {
  int    col, colStart, colEnd, row;
  double value;

  for (col = 0; col < n; col++) {
	 f[col] = Problem.q[col];
  }

  for (col = 0; col < n; col++) {
	 value = z[col];

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
int Solvers::PATH::C_jacobian_evaluation(int     n,
													  double *z,
													  int     wantf,
													  double *f,
													  int *   nnz,
													  int *   col_start,
													  int *   col_len,
													  int *   row,
													  double *data) {
  int element;

  if (wantf) {
	 this->C_function_evaluation(n, z, f);
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

	 Filled = 1;
  }

  *nnz = Problem.nnz;
  return 0;
}

void *Solvers::PATH::mcp_typ(int nnz, int *typ) {
  int i;

  for (i = 0; i < nnz; i++) {
	 typ[i] = PRESOLVE_LINEAR;
  }
}

void Solvers::PATH::C_problem_size(int *n, int *nnz) {
  *n   = this->Problem.n;
  *nnz = this->Problem.nnz + 1;
}
void Solvers::PATH::bounds(void *dat, int n, double *z, double *lb, double *ub) {
  auto *self = static_cast<Solvers::PATH *>(dat);
  self->C_bounds(n, z, lb, ub);
  return;
}
