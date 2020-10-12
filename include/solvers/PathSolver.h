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