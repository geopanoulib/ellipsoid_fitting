/**
 * \file		multiply.c
 * \brief       Matrix multiplication
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
 * \brief           It multiplies two matrices or vectors
 * \param[in]       a: The matrix a [dimensions (n x m)]
 * \param[in]       b: The matrix b [dimensions (m x k)]
 * \param[in]       c: The matrix c, where c = a * b [dimensions (n x k)]
 * \param[in]       n: The number of rows of the matrix a
 * \param[in]       m: The number of columns of the matrix a
 * \param[in]       k: The number of columns of the matrix b
 */
void multiply(type *a, type *b, type *c, int n, int m, int k)
{
	register int i, j, v;
	type sum;
			
	for (v = 0; v < n; v++)
	{
        	for(j = 0; j < k; j++)
        	{
              		sum = 0.0L;
              		for (i = 0; i < m; i++)
              			sum += a[v * m + i] * b[i * k + j];
        		c[v * k + j] = sum;
          	}
	}		 	   
}

