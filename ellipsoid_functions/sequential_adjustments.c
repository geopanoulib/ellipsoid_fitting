/**
 * \file		sequential_adjustments.c
 * \brief       Ellipsoid fitting (technique: sequential adjustments)
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
 * \brief           Fitting a triaxial ellipsoid by applying
 * 					the sequential adjustments technique
 * \param[in]       argc: The number of arguments that main was called (integer)
 * \param[in]       argv: The vector of the arguments (string) which
 * 					contains the names of the data files (binary files)
 * \return			An integer value equal to zero
 */
int main(int argc, char *argv[])
{
	int c, n, m, t = argc - 1, iteration;
	register int i, j;
	FILE *files[t];
	type in_val[9], ds[9], N_inv[9][9];
	type Vx[9][9];
	type uTds, sigma0i, sigma0ip1, stx, sty, stz, sax, say, saz, sthetax, sthetay, sthetaz;
	struct solution x, x1;
	struct group mat;
	
	/* Sorting the names of the included data files (binary files) in alphabetical order */
	alpha_sort(argv, argc);
	/* File read control */
	for (i = 0; i < t; i++)
		if((files[i] = fopen(argv[i + 1], "rb")) == NULL)
		{
			printf("\n\tCant open the file %s", argv[i + 1]);
			exit(1);
		}
	/* Calculating the number of points of the first group and the initial (of the first solution) values by calling the function initial_values() */
	c = initial_values(files, 1, &in_val[0]);
	x1.r = c - 9;
	iteration = 0;
	sigma0i = 1.0L;
	sigma0ip1 = 2.0L;
	/* Iterative adjustment procedure (only for the first group) */
	do {
		sigma0i = sigma0ip1;
		mat = direct_calculation(files[0], &in_val[0]);		
		cholesky(&mat.N_bar[0][0], &N_inv[0][0], 9);
		multiply(&N_inv[0][0], &mat.U_bar[0], &ds[0], 9, 9, 1);
		multiply(&mat.U_bar[0], &ds[0], &uTds, 1, 9, 1);
		sigma0ip1 = sqrt((mat.sum_piwi2 - uTds) / x1.r);
		
		for(i = 0; i < 9; i++)
			in_val[i] += ds[i];
		rewind(files[0]);
		iteration++;
	} while(MYABS(sigma0i - sigma0ip1) > CONVTOL && iteration < 10);
	for (i = 0; i < 9; i++)
	{
		for (j = 0; j < 9; j++)
			x1.Nbar[i][j] = mat.N_bar[i][j];
		x1.x[i] = in_val[i];
	}
	
	printf("\nNumber of files = %d\nIncluded files :", t);
	for (i = 0; i < t; i++)
		printf("\n%s", argv[i + 1]);
		
	x1.s02 = mat.sum_piwi2 / x1.r;
	/* Printing the first solution that contains the measurements of the group 1 */
	display(&x1.x[0], 9, 1, 4, "x1");
	printf("\nc1 = %d\nr1 = %d", c, x1.r);
	printf("\ns01 = +/- %-.5Lf\nIterations = %d\n", sigma0ip1, iteration);
	
	/* Sequential adjustments procedure */
	for (i = 1; i < t; i++)
	{
		x = sequential(files[i], x1);
		printf("\n#--------------------------#");
		printf("\nc%d = %d\nr%d = %d", i + 1, x.r + 9, i + 1, x.r);
		printf("\ns0_%d = +/- %-.5Lf", i + 1, (type)sqrt(x.s02));
		x1 = x;
	}
	/* Inversion of the matrix N by calling the function cholesky() */
	cholesky(&x.Nbar[0][0], &N_inv[0][0], 9);
	
	/* Calculating the variance-covariance matrix */
	for (i = 0; i < 9; i++)
		for (j = 0; j < 9; j++)
			Vx[i][j] = x.s02 * N_inv[i][j];
	/* Closing all the data files */
	for (i = 0; i < t; i++)
		fclose(files[i]);
		
	/* Calculating each parameter's std */		
	stx = sqrt(Vx[0][0]);
	sty = sqrt(Vx[1][1]);
	stz = sqrt(Vx[2][2]);
	sax = sqrt(Vx[3][3]);
	say = sqrt(Vx[4][4]);
	saz = sqrt(Vx[5][5]);
	sthetax = sqrt(Vx[6][6]) * RDEG;
	sthetay = sqrt(Vx[7][7]) * RDEG;
	sthetaz = sqrt(Vx[8][8]) * RDEG;
	x.x[6] *= RDEG; 
	x.x[7] *= RDEG;
	x.x[8] *= RDEG;

	/* Printing the results */
	c = x.r + 9;
	n = 3 * c;
	m = 9 + 2 * c;
	printf("\n\nc = %d points", c);
	printf("\nn = %d measurements", n);
	printf("\nm = %d unknowns", m);
	printf("\nr = %d degrees of freedom", x.r);
	printf("\nExecution time = %ld [s]", clock() / CLOCKS_PER_SEC);
	printf("\n\nElipsoid Parameters:");
	printf("\ntx = %-.4Lf +/- %-.5Lf [m]", x.x[0], stx);
	printf("\nty = %-.4Lf +/- %-.5Lf [m]", x.x[1], sty);
	printf("\ntz = %-.4Lf +/- %-.5Lf [m]", x.x[2], stz);
	printf("\nax = %-.4Lf +/- %-.5Lf [m]", x.x[3], sax);
	printf("\nay = %-.4Lf +/- %-.5Lf [m]", x.x[4], say);
	printf("\naz = %-.4Lf +/- %-.5Lf [m]", x.x[5], saz);
	printf("\ntheta_x = %-.4Lf +/- %-.5Lf [deg]", x.x[6], sthetax);
	printf("\ntheta_y = %-.4Lf +/- %-.5Lf [deg]", x.x[7], sthetay);
	printf("\ntheta_z = %-.4Lf +/- %-.5Lf [deg]", x.x[8], sthetaz);
	printf("\ns0_aposteriori = +/- %-.4Lf [m]", (type)sqrt(x.s02));
	display(&Vx[0][0], 9, 9, 7, "Vx");
	return 0;	
}

