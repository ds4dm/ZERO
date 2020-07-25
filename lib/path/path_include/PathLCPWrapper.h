#pragma once
int SimpleLCP(int     variables,
				  int     m_nnz,
				  int *   m_i,
				  int *   m_j,
				  double *m_ij,
				  double *q,
				  double *lb,
				  double *ub,
				  double *z);
