/**
 * \file		symmetric.c
 * \brief       Symmetric matrix
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
 * \brief           Converts an upper triangular matrix into a symmetric one
 * \param[in]       a: The matrix a
 * \param[in]       n: The number of rows or columns
 */
void symmetric(type *a, int n)
{
	register int i, j;
	
		for (i = 1; i < n; i++)
			for (j = 0; j < i; j++)
				a[i * n + j] = a[j * n + i];			
}

