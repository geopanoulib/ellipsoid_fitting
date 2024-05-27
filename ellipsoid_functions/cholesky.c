/**
 * \file		cholesky.c
 * \brief       Matrix inversion
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
 * \brief           It caclulates the inverse matrix using the cholesky inverse method
 * \param[in]       N: The symmetric matrix to be inversed
 * \param[in]       Ninv: Inverse matrix
 * \param[in]       n: Dimension of the N matrix
 */
void cholesky(type *N, type *Ninv, int n)
{
	register int i, j, k;
     	type sumki, sumkij, sumdc, sumb;
     	type c[n][n], d[n][n];
        
     	/* calculation of the matrix c (upper triangular) */
     	for (i = 0; i < n; i++)
     	{
     		sumki = 0.0L;
     		for (k = 0; k <= i - 1; k++)
     			sumki += c[k][i] * c[k][i];
     		c[i][i] = sqrt(N[i * n + i] - sumki);
     		d[i][i] = 1.0L / c[i][i];
     		for (j = i + 1; j < n; j++)
     		{
     			sumkij = 0.0L;
     			for (k = 0; k <= i - 1; k++)
     				sumkij += c[k][i] * c[k][j];
     			c[i][j] = (N[i * n + j] - sumkij) / c[i][i];
     		}	
     	}
     	/* calculation of the matrix d (upper triangular) */
     	for (i = 0; i < n; i++)
     		for (j = i + 1; j < n; j++)
     		{
     			sumdc = 0.0L;
     			for (k = i; k <= j - 1; k++)
     				sumdc += d[i][k] * c[k][j];
     			d[i][j] = -sumdc / c[j][j];  
     		}
     	/* calculation of the inverse matrix */
     	for (i = 0; i < n; i++)
     		for (j = i; j < n; j++)
     		{
     			sumb = 0.0L;
     			for (k = j; k < n; k++)
     				sumb += d[i][k] * d[j][k];
     			Ninv[i * n + j] = sumb; 
     		}
     	/* calling the function symmetric() to turn the upper triangular matrix Ninv into a symmetric one */		
     	symmetric(Ninv, n); 
}

