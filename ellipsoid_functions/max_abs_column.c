/**
 * \file		max_abs_column.c
 * \brief       Maximum absolute value in a column
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
 * \brief           It calculates the maximum absolute value of each column for a given matrix
 * \param[in]       d: The matrix
 * \param[in]       colmn: The vector of maximum values
 * \param[in]       n: number of rows of d
 * \param[in]       m: number of columns of d
 */
void max_abs_column(type *d, type *column, int n, int m)
{
	register int i, j;
	type max;
	
	for (j = 0; j < m; j++){
		max = 0.0L;
		for (i = 0; i < n; i++)
			if(MYABS(d[i * m + j]) > max)
				max = MYABS(d[i * m + j]);
		column[j] = max;
	}
}

