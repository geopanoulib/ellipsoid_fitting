/**
 * \file		initial_values.c
 * \brief       Initial values for the ellipsoid parameters
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
 * \brief           Calculates the initial values of the parameters of 
 * 					the triaxial ellipsoid by using the data files of the measurements,
 * 					i.e. the points for the ellipsoid fitting
 * \param[in]       fp: Vector of the data files pointers
 * \param[in]       file_num: Number of the data files
 * \param[in]       values: Vector of the initial values of the triaxial ellipsoid
 * \return			The number of points
 */
int initial_values(FILE *fp[], int file_num, type *values)
{
	register int i;
	int cnt = 0;
	type xi, yi, zi, xi2, yi2, zi2, N[9][9], U[9], C[9], Ninv[9][9];
	type cxx, cyy, czz, cxy, cxz, cyz, cx, cy, cz, f1, f2, f3, g2, g3, h3, e;
	type tx, ty, tz, ax, ay, az, theta_x, theta_y, theta_z, qxx, qxy, qxz, qyy, qyz, qzz, d;
	type q1, q2, w, Q;
	type A1, B1, C1, A2, B2, C2, A3, B3, C3, E1, E2, E3;
	type a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11;
	type a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22;
	type a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34;
	struct cart_coord pp;
	
	/* Initializing the values of the matrix N and the vector u */
	a1 = a2 = a3 = a4 = a5 = a6 = a7 = a8 = a9 = a10 = a11 = 0.0L;
	a12 = a13 = a14 = a15 = a16 = a17 = a18 = a19 = a20 = a21 = a22 = 0.0L;
	a23 = a24 = a25 = a26 = a27 = a28 = a29 = a30 = a31 = a32 = a33 = a34 = 0.0L;
	/* Reading the files (binary files) and calculating the elements of the N and u matrices */
	for (i = 0; i < file_num; i++)
		while (fread(&pp, 1, sizeof(pp), fp[i]) > 0) {
			xi = pp.x;
			yi = pp.y;
			zi = pp.z;
			xi2 = xi * xi;
			yi2 = yi * yi;
			zi2 = zi * zi;
			a1 += xi2 * xi2; 
			a2 += xi2 * yi2;
			a3 += xi2 * zi2;
			a4 += xi2 * xi * yi;
			a5 += xi2 * xi * zi;
			a6 += xi2 * yi * zi;
			a7 += xi2 * xi;
			a8 += xi2 * yi;
			a9 += xi2 * zi;
			a10 += yi2 * yi2;
			a11 += yi2 * zi2;
			a12 += xi * yi2 * yi;
			a13 += xi * yi2 * zi;
			a14 += yi2 * yi * zi;
			a15 += xi * yi2;
			a16 += yi2 * yi;
			a17 += yi2 * zi;
			a18 += zi2 * zi2;
			a19 += xi * yi * zi2;
			a20 += xi * zi2 * zi;
			a21 += yi * zi2 * zi;
			a22 += xi * zi2;
			a23 += yi * zi2;
			a24 += zi2 * zi;
			a25 += xi * yi * zi;
			a26 += xi2;
			a27 += xi * yi;
			a28 += xi * zi;
			a29 += yi2;
			a30 += yi * zi;
			a31 += zi2;
			a32 += xi;
			a33 += yi;
			a34 += zi;
			cnt++;
		}
	/* Designing the matrix N */
	N[0][0] = a1;
	N[0][1] = a2;
	N[0][2] = a3;
	N[0][3] = a4;
	N[0][4] = a5;
	N[0][5] = a6;
	N[0][6] = a7;
	N[0][7] = a8;
	N[0][8] = a9;
	N[1][1] = a10;
	N[1][2] = a11;
	N[1][3] = a12;
	N[1][4] = a13;
	N[1][5] = a14;
	N[1][6] = a15;
	N[1][7] = a16;
	N[1][8] = a17;
	N[2][2] = a18;
	N[2][3] = a19;
	N[2][4] = a20;
	N[2][5] = a21;
	N[2][6] = a22;
	N[2][7] = a23;
	N[2][8] = a24;
	N[3][3] = a2;
	N[3][4] = a6;
	N[3][5] = a13;
	N[3][6] = a8;
	N[3][7] = a15;
	N[3][8] = a25;
	N[4][4] = a3;
	N[4][5] = a19;
	N[4][6] = a9;
	N[4][7] = a25;
	N[4][8] = a22;
	N[5][5] = a11;
	N[5][6] = a25;
	N[5][7] = a17;
	N[5][8] = a23;
	N[6][6] = a26;
	N[6][7] = a27;
	N[6][8] = a28;
	N[7][7] = a29;
	N[7][8] = a30;
	N[8][8] = a31;
	/* Designing the vector u */
	U[0] = a26;
	U[1] = a29;
	U[2] = a31;
	U[3] = a27;
	U[4] = a28;
	U[5] = a30;
	U[6] = a32;
	U[7] = a33;
	U[8] = a34;
	/* Converting the upper triangular matrix N into a symmetric one by calling the function symmetric() */
	symmetric(&N[0][0], 9);
	/* Inversion of the matrix N by calling the function cholesky() */
	cholesky(&N[0][0], &Ninv[0][0], 9);
	/* Multiplication of the N and u matrices by calling the function multiply() */
	multiply(&Ninv[0][0], &U[0], &C[0], 9, 9, 1);
	/* Calculation of the initial values of the triaxial ellipsoid */
	cxx = C[0];
	cyy = C[1];
	czz = C[2];
	cxy = C[3];
	cxz = C[4];
	cyz = C[5];
	cx = C[6];
	cy = C[7];
	cz = C[8];
	f1 = 4. * cyy * czz - cyz * cyz;
	f2 = cxz * cyz - 2 * cxy * czz;
	f3 = cxy * cyz - 2 * cxz * cyy;
	g2 = 4. * cxx * czz - cxz * cxz;
	g3 = cxy * cxz - 2 * cxx * cyz;
	h3 = 4. * cxx * cyy - cxy * cxy;
	e = 2. * cxx * f1 + cxy * f2 + cxz * f3;
	if(e == 0.0L)
		printf("\n\t Error in function Initial values , e = 0\n");
	tx = - (f1 * cx + f2 * cy + f3 * cz) / e;
	ty = - (f2 * cx + g2 * cy + g3 * cz) / e;
	tz = - (f3 * cx + g3 * cy + h3 * cz) / e;
	d = 1.0L + cxx * tx * tx + cyy * ty * ty + czz * tz * tz + cxy * tx * ty + cxz * tx * tz + cyz * ty * tz;
	qxx = 2.0L * d * f1 / e;
	qxy = 2.0L * d * f2 / e;
	qxz = 2.0L * d * f3 / e;
	qyy = 2.0L * d * g2 / e;
	qyz = 2.0L * d * g3 / e;
	qzz = 2.0L * d * h3 / e;
	q1 = 1.0L * (qxx + qyy + qzz) / 3.0L;
	q2 = 1.0L * (qyy * qzz + qxx * qzz + qxx * qyy - qyz * qyz - qxz * qxz - qxy * qxy) / 3.0L;
	Q = qxx * (qyy * qzz - qyz * qyz) + qxy * (qxz * qyz - qxy * qzz) + qxz * (qxy * qyz - qxz * qyy);
	w = acos((Q + 2 * q1 * q1 * q1 - 3 * q1 * q2) / (2 * pow(q1 * q1 - q2, 1.5)));
	ax = sqrt(q1 + 2 * sqrt(q1 * q1 - q2) * cos(w / 3));
	ay = sqrt(q1 + 2 * sqrt(q1 * q1 - q2) * cos((w - 2 * M_PI) / 3));
	az = sqrt(q1 + 2 * sqrt(q1 * q1 - q2) * cos((w + 2 * M_PI) / 3));
	A1 = qxy * qxz - qyz * qxx + ax * ax * qyz;
	B1 = qxy * qyz - qxz * qyy + ax * ax * qxz;
	C1 = qxz * qyz - qxy * qzz + ax * ax * qxy;
    	A2 = qxy * qxz - qyz * qxx + ay * ay * qyz;
	B2 = qxy * qyz - qxz * qyy + ay * ay * qxz;		
	C2 = qxz * qyz - qxy * qzz + ay * ay * qxy;
	A3 = qxy * qxz - qyz * qxx + az * az * qyz;
	B3 = qxy * qyz - qxz * qyy + az * az * qxz;
	C3 = qxz * qyz - qxy * qzz + az * az * qxy;
	E1 = sqrt(1.0L / A1 / A1 + 1.0L / B1 / B1 + 1.0L / C1 / C1);
	E2 = sqrt(1.0L / A2 / A2 + 1.0L / B2 / B2 + 1.0L / C2 / C2);
	E3 = sqrt(1.0L / A3 / A3 + 1.0L / B3 / B3 + 1.0L / C3 / C3);
	theta_x = atan(-C3 / B3);
	theta_y = atan((A1 * E1 * A2 * E2) / (A3 * E3 * sqrt(A1 * A1 * E1 * E1 + A2 * A2 * E2 * E2)));
	theta_z = atan(-A1 * E1 / A2 / E2);
	/* Vector of the initial values */
	values[0] = tx;
	values[1] = ty;
	values[2] = tz;
	values[3] = ax;
	values[4] = ay;
	values[5] = az;
	values[6] = (theta_x < 0) ? -theta_x:theta_x;
	values[7] = (theta_y < 0) ? -theta_y:theta_y;
	values[8] = (theta_z < 0) ? -theta_z:theta_z;
	/* Returning the files position indicators to the beginning of the each file */
	for (i = 0; i < file_num; i++)
		rewind(fp[i]);
	return cnt;
}

