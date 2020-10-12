#include "solvers/PathSolver.h"
#include "path/MCP_Interface.h"
#include "path/Macros.h"
#include "path/Options.h"
#include "path/Output.h"
#include "path/Output_Interface.h"
#include "path/Path.h"
#include "path/PathOptions.h"
#include "path/Types.h"
#include <climits>
#include <cstdio>

typedef struct {
  int variables;

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
} Problem_t;

static Problem_t problem;
static int       filled;

extern "C" {


static CB_FUNC(void) problem_size(void *id, int *n, int *nnz) {
  *n   = problem.n;
  *nnz = problem.nnz;
  return;
} // problem_size

static CB_FUNC(void) bounds(void *id, int n, double *z, double *lb, double *ub) {
  int i;

  for (i = 0; i < n; i++) {
	 z[i]  = problem.z[i];
	 lb[i] = problem.lb[i];
	 ub[i] = problem.ub[i];
  }
  return;
} // bounds

static CB_FUNC(int) function_evaluation(void *v, int n, double *z, double *f) {
  int    col, colStart, colEnd, row;
  double value;

  for (col = 0; col < n; col++) {
	 f[col] = problem.q[col];
  }

  for (col = 0; col < n; col++) {
	 value = z[col];

	 if (value != 0) {
		colStart = problem.m_start[col] - 1;
		colEnd   = colStart + problem.m_len[col];

		while (colStart < colEnd) {
		  row = problem.m_row[colStart] - 1;
		  f[row] += problem.m_data[colStart] * value;
		  colStart++;
		}
	 }
  }

  return 0;
}

static CB_FUNC(int) jacobian_evaluation(void *  v,
													 int     n,
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
	 function_evaluation(v, n, z, f);
  }

  if (!filled) {
	 for (element = 0; element < problem.n; element++) {
		col_start[element] = problem.m_start[element];
		col_len[element]   = problem.m_len[element];
	 }

	 for (element = 0; element < problem.nnz; element++) {
		row[element]  = problem.m_row[element];
		data[element] = problem.m_data[element];
	 }

	 filled = 1;
  }

  *nnz = problem.nnz;
  return 0;
}

static CB_FUNC(void) mcp_typ(void *d, int nnz, int *typ) {
  int i;

  for (i = 0; i < nnz; i++) {
	 typ[i] = PRESOLVE_LINEAR;
  }
  return;
}

static CB_FUNC(void) messageCB(void *data, int mode, char *buf) {
  fprintf(stdout, "%s", buf);
} /* messageCB */
}

static Output_Interface outputInterface = {NULL, messageCB, NULL};


static MCP_Interface mcp_interface = {NULL,
												  problem_size,
												  bounds,
												  function_evaluation,
												  jacobian_evaluation,
												  NULL,
												  NULL,
												  NULL,
												  NULL,
												  NULL};

static Presolve_Interface mcp_presolve = {NULL, NULL, NULL, NULL, NULL, mcp_typ, 0};


int PathLCP(int    n,
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
  Options_Interface *o;
  MCP *              m;
  MCP_Termination    termination;
  Information        information;

  Output_SetLog(stdout);


  o = Options_Create();
  Path_AddOptions(o);
  Options_Default(o);


  double *xSol, *zSol;
  double  dnnz;
  int     i;

  Output_SetInterface(&outputInterface);
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

  problem.n       = n;
  problem.z       = x;
  problem.q       = q;
  problem.lb      = lb;
  problem.ub      = ub;
  problem.m_start = new int[n];
  problem.m_len   = new int[n];
  problem.m_row   = new int[m_nnz + 1];
  problem.m_data  = new double[m_nnz + 1];

  int m_index = 0, m_count = 0;
  for (i = 0; i < n; i++) {
	 problem.m_start[i] = m_count + 1;
	 problem.m_len[i]   = 0;

	 while ((m_index < m_nnz) && (m_j[m_index] <= i + 1)) {
		if (m_ij[m_index] != 0) {
		  problem.m_len[i]++;
		  problem.m_row[m_count]  = m_i[m_index];
		  problem.m_data[m_count] = m_ij[m_index];
		  m_count++;
		}
		m_index++;
	 }
  }
  problem.nnz = m_count;

  // Path_Size(problem.n, problem.nnz);
  m = MCP_Create(problem.n, problem.nnz + 1);
  MCP_SetInterface(m, &mcp_interface);
  MCP_SetPresolveInterface(m, &mcp_presolve);
  MCP_Jacobian_Structure_Constant(m, static_cast<Boolean>(true));



  Information info;
  info.generate_output = Output_Log | Output_Status | Output_Listing;
  info.use_start       = True;
  info.use_basics      = True;

  termination = Path_Solve(m, &info);



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

  return termination;
}
