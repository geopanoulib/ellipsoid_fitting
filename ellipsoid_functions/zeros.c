/**
 * \file		zeros.c
 * \brief       Matrix initialization
 */

/**
 *
 * Copyright (c) 2024, Jason Koci
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHA
 * NTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 * 
 * Author:	Jason Koci <iasonaskotsis@hotmail.com>
 */

#include "ellipsoid_functions.h"

/**
 * \brief           Initializes the values of a matrix by resetting the elements of that matrix to zero
 * \param[in]       a: The matrix a 
 * \param[in]       n: The number of rows of the matrix a
 * \param[in]       m: The number of columns of the matrix a
 */
void zeros(type *a, int n, int m)
{
	register int i, j;
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			a[i * m + j] = 0.0L; 
}

