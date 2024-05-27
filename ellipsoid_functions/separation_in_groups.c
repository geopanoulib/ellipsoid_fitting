/**
 * \file		separation_in_groups.c
 * \brief       Ellipsoid fitting (technique: separation in groups)
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
 * \brief           Fitting a triaxial ellipsoid by applying the separation in groups (of measurements) technique
 * \param[in]       argc: The number of arguments that main was called (integer)
 * \param[in]       argv: The vector of arguments (string) which
 * 					contains the names of the data files (binary files)
 * \return			An integer value equal to zero
 */
int main(int argc, char *argv[])
{
	register int i, j;
	int c, n, m, r, t = argc - 1, iteration;
	type in_val[9], ds[9], N_inv[9][9], Vx[9][9];
	type uTds, sigma0i, sigma0ip1, stx, sty, stz, sax, say, saz, sthetax, sthetay, sthetaz;
	FILE *files[t];
	struct group mat[t];
	struct group final_mat;
	
	/* Sorting the names of the included data files (binary files) in alphabetical order */
	alpha_sort(argv, argc);
	/* Data Files control */
	for (i = 0; i < t; i++)
		if((files[i] = fopen(argv[i + 1], "rb")) == NULL)
		{
			printf("\nCant open the file %s", argv[i + 1]);
			exit(1);
		}
	/* Calculating the total number of points and the initial values by calling the function initial_values() */
	c = initial_values(files, t, &in_val[0]);
	/* Printing the initial values of the triaxial ellipsoid */
	display(&in_val[0], 9, 1, 4, "Initial Values");
	n = 3 * c; /* Total number of measurements */
	m = 9 + 2 * c; /* Total number of unknowns */
	r = n - m; /* Degrees of freedom */
	iteration = 0;
	sigma0i = 1.0L;
	sigma0ip1 = 2.0L;
	/* Iterative adjustment procedure */
	do {
		sigma0i = sigma0ip1;
		for (i = 0; i < t; i++)
			mat[i] = direct_calculation(files[i], &in_val[0]);
				
		final_mat = summary(mat, t);	
		cholesky(&final_mat.N_bar[0][0], &N_inv[0][0], 9);
		multiply(&N_inv[0][0], &final_mat.U_bar[0], &ds[0], 9, 9, 1);
		for(i = 0; i < 9; i++)
			in_val[i] += ds[i];
			
		multiply(&final_mat.U_bar[0], &ds[0], &uTds, 1, 9, 1);
		sigma0ip1 = sqrt((final_mat.sum_piwi2 - uTds) / r);	
		iteration++;
	} while(MYABS(sigma0i - sigma0ip1) > CONVTOL && iteration < 10);
	
	/* Closing all the data files */
	for (i = 0; i < t; i++)
		fclose(files[i]);
	
	/* Calculating the variance-covariance matrix */
	for(i = 0; i < 9; i++)
		for(j = 0; j < 9; j++)
			Vx[i][j] = sigma0ip1 * sigma0ip1 * N_inv[i][j];
	
	/* Calculating each parameter's std */		
	stx = sqrt(Vx[0][0]);
	sty = sqrt(Vx[1][1]);
	stz = sqrt(Vx[2][2]);
	sax = sqrt(Vx[3][3]);
	say = sqrt(Vx[4][4]);
	saz = sqrt(Vx[5][5]);
	sthetax = sqrt(Vx[6][6]) * RDEG; /* deg */
	sthetay = sqrt(Vx[7][7]) * RDEG; /* deg */
	sthetaz = sqrt(Vx[8][8]) * RDEG; /* deg */
	in_val[6] *= RDEG; /* deg */
	in_val[7] *= RDEG; /* deg */
	in_val[8] *= RDEG; /* deg */
	
	/* Printing results */
	printf("\nNumber of files = %d\nIncluded files :", t);
	for (i = 0; i < t; i++)
		printf("\n%s", argv[i + 1]);
	printf("\n\nc = %d points", c);
	printf("\nn = %d measurements", n);
	printf("\nm = %d unknowns", m);
	printf("\nr = %d degrees of freedom", r);
	printf("\nIterations = %d", iteration);
	printf("\nExecution time = %ld [s]", clock() / CLOCKS_PER_SEC);
	printf("\n\nElipsoid Parameters:");
	printf("\ntx = %-.4Lf +/- %-.5Lf [m]", in_val[0], stx);
	printf("\nty = %-.4Lf +/- %-.5Lf [m]", in_val[1], sty);
	printf("\ntz = %-.4Lf +/- %-.5Lf [m]", in_val[2], stz);
	printf("\nax = %-.4Lf +/- %-.5Lf [m]", in_val[3], sax);
	printf("\nay = %-.4Lf +/- %-.5Lf [m]", in_val[4], say);
	printf("\naz = %-.4Lf +/- %-.5Lf [m]", in_val[5], saz);
	printf("\ntheta_x = %-.4Lf +/- %-.5Lf [deg]", in_val[6], sthetax);
	printf("\ntheta_y = %-.4Lf +/- %-.5Lf [deg]", in_val[7], sthetay);
	printf("\ntheta_z = %-.4Lf +/- %-.5Lf [deg]", in_val[8], sthetaz);
	printf("\ns0_aposteriori = +/- %-.4Lf [m]", sigma0ip1);
	display(&Vx[0][0], 9, 9, 7, "Vx");
	return 0;
}

