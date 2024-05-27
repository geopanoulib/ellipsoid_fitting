/**
 * \file            alpha_sort.c
 * \brief           Input files sorting
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
 * \brief           Sorts the input files alphabetically
 * \param[in]       s: The names of the files in vector form
 * \param[in]       args: The total number of files
 */
void alpha_sort(char *s[], int args)
{
	int val;
	register int i, j;
	long unsigned int pos1, pos2;
	bool cas1, cas2, cas3;
	char *lower, *bigger;
	
	for (i = 1; i < args; i++)
		for (j = i + 1; j < args; j++){
			pos1 = strlen(s[i]) - strlen(strrchr(s[i], '.'));
			pos2 = strlen(s[j]) - strlen(strrchr(s[j], '.'));
			val = strcmp(s[i], s[j]);
			cas1 = pos1 == pos2 && val < 0;
			cas2 = pos1 < pos2 && val < 0;
			cas3 = pos1 < pos2 && val > 0;
			if(cas1 || cas2 || cas3){
				lower = s[i];
				bigger = s[j];
			}else {
				lower = s[j];
				bigger = s[i];
			} 
			s[i] = lower;
			s[j] = bigger;
		}
}

