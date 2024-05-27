/**
 * \file	ellipsoid_functions.h
 * \brief       This header file contains all the necessary
 *		libraries, functions, structures and constants
 *		for the least-squares fitting of a triaxial ellipsoid
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#define type long double /* Data type */
#define MYABS(x) (((x)>0) ? (x):-(x)) /* A macro function that returns the absolute value of a number (inline function) */
#define RDEG 180.0L / M_PI /* A constant value for the conversion from rad to degrees */
#define CONVTOL 1e-5 /* Convergence tolerance */

int digitc(type);
void max_abs_column(type *, type *, int, int);
void display(type *, int, int, int, char []);
void zeros(type *, int, int);
void symmetric(type *, int);
void cholesky(type *, type *, int);
void multiply(type *, type *, type *, int, int, int);
int initial_values(FILE **, int, type *);
void alpha_sort(char *[], int);

struct solution sequential(FILE *, struct solution);
struct group direct_calculation(FILE *, type *);
struct group summary(struct group *, int);

/* A structure for the Cartesian coordinates and their weights */
struct cart_coord {
 	double x;
	double y;
	double z;
	double w;
};

/* A structure for the elements of each group of measurements */
struct group {
	int c;
	type N_bar[9][9];
	type U_bar[9];
	type sum_piwi2;
};

/* A structure for the elements of a particular solution after an adjustment process */
struct solution {
	int r;
	type x[9];
	type Nbar[9][9];
	type s02;
};

