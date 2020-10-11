#pragma once

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
				double timeLimit);