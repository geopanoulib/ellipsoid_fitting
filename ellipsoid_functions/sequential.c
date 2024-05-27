/**
 * \file		sequential.c
 * \brief       Sequential adjustments function
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
 * \brief           Given an initial solution, it calculates the
 * 					revised solution by applying the sequential
 * 					adjustment technique
 * \param[in]       file: The data file
 * \param[in]       x1: The structure <solution> which contain
 * 					the previous solution of the sequential adjustment
 * \return			The revised solution of the sequential adjustment
 */
struct solution sequential(FILE *file, struct solution x1)
{
	struct group M2;
	struct solution x;
	type invNbar[9][9], dx1[9], u2Nu2 = 0.0L;
	register int i, j;
	
	/* Calculating the matrix N2 and the vector u2 of the added measurements */
	M2 = direct_calculation(file, &x1.x[0]);
	/* Calculating the final matrix N (N = N1 + N2) */
	for (i = 0; i < 9; i++)
		for(j = 0; j < 9; j++)
			x.Nbar[i][j] = x1.Nbar[i][j] + M2.N_bar[i][j];
	/* Inversion of the matrix N by calling the function cholesky() */
	cholesky(&x.Nbar[0][0], &invNbar[0][0], 9);
	/* Multiplication of the N and u2 matrices by calling the function multiply() */
	multiply(&invNbar[0][0], &M2.U_bar[0], &dx1[0], 9, 9, 1);
	x.r = x1.r + M2.c; /* Degrees of freedom */
	/* Calculation of the revised solution */
	for (i = 0; i < 9; i++)
	{
		x.x[i] = x1.x[i] + dx1[i];
		u2Nu2 += dx1[i] * M2.U_bar[i];
	}
	/* Calculating the revised a-posteriori variance factor */
	x.s02 = (x1.r * x1.s02 - u2Nu2 + M2.sum_piwi2) / x.r;
	return x;
}

