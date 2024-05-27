/**
 * \file		display.c
 * \brief       Vector and matrix printer function
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
 * \brief           It prints a given matrix or vector using the minimum space needed to display the elements
 * \param[in]       M: The matrix or vector to be print
 * \param[in]       n: Number of rows
 * \param[in]       m: Number of columns
 * \param[in]       precision: Number of significant digits to be print
 * \param[in]       N: The name of the matrix or the vector
 */

void display(type *M, int n, int m, int precision, char N[])
{
	register int i, j;
	type max_column[m];
	
	max_abs_column(M, &max_column[0], n, m); /* calculating the maximum absolute value of every column by calling the function max_abs_column() */
	printf("\n\n%s = \n\n", N);
	for (i = 0; i < n; i++)
	{
		printf(" | ");
		for (j = 0; j < m; j++)
			printf("%*.*Lf ", digitc(max_column[j]) + precision + 2, precision, M[i * m + j]);
		printf(" | \n");
	}
}

