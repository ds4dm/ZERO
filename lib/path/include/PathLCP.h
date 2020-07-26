#ifndef __PATHLCP_H_
#define __PATHLCP_H_

#ifdef __cplusplus
extern "C" {
#endif

int PathLCP(int     variables,
				int     m_nnz,
				int *   m_i,
				int *   m_j,
				double *m_ij,
				double *q,
				double *lb,
				double *ub,
				double *z,
				int     verbose,
				double  timeLimit);
#ifdef __cplusplus
}
#endif
#endif