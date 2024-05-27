/**
 * \file		matrix_summary.c
 * \brief       Summation of the matrices N and u
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
 * \brief           Calculates the sum of the matrices N and the
 * 					vectors u of each group of measurements
 * \param[in]       A: A pointer of type <group> which contains
 * 					all the matrices and elements for each group
 * 					of measurements
 * \param[in]       n: The number of rows of the matrix N
 * \return			A structure of type <group> which contains
 * 					the matrices N and u, referring to all groups
 * 					of measurements
 */
struct group summary(struct group *A, int n)
{	
	struct group SUM;
	register int i, j, k;
	
	/* Resetting the elements of the matrix N (of all froups) to zero */
	zeros(&SUM.N_bar[0][0], 9, 9);
	/* Resetting the elements of the vector u (of all groups) to zero */
	zeros(&SUM.U_bar[0], 9, 1);
	SUM.sum_piwi2 = 0.0L;
	/* Summation procedure */
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < 9; j++)
		{
			for (k = 0; k < 9; k++)
				SUM.N_bar[j][k] += A[i].N_bar[j][k];
			SUM.U_bar[j] += A[i].U_bar[j];	
		}
		SUM.sum_piwi2 += A[i].sum_piwi2;			
	}
	return SUM;
}

