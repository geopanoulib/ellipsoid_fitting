/**
 * \file		direct_calculation.c
 * \brief       Direct calculation of the matrix N and the vector u 
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
 * \brief           Calculation of the matrix N and the vector u
 * 					by applying the direct calculation technique
 * \param[in]       fp: Data file pointer
 * \param[in]       values: Vector of the initial values of the triaxial ellipsoid
 * \return			A structure of type <group> which contains matrices
 * 					and elements for a specific group of measurements
 */
struct group direct_calculation(FILE *fp, type *values)
{
	type xi, yi, zi, pi;
	type theta_x, theta_y, theta_z, tx, ty, tz, ax, ay, az; 
	type r11, r12, r13, r21, r22, r23, r31, r32, r33; 
	type pxx, pyy, pzz, pxy, pxz, pyz, A, B, C; 
	type dpxx_dax, dpxx_day, dpxx_daz, dpxx_dthy, dpxx_dthz; 
	type dpyy_dax, dpyy_day, dpyy_daz, dpyy_dthy, dpyy_dthz;
	type dpzz_dax, dpzz_day, dpzz_daz, dpzz_dthy, dpzz_dthz;
	type dpxy_dax, dpxy_day, dpxy_daz, dpxy_dthy, dpxy_dthz;
	type dpxz_dax, dpxz_day, dpxz_daz, dpxz_dthy, dpxz_dthz;
	type dpyz_dax, dpyz_day, dpyz_daz, dpyz_dthy, dpyz_dthz;
	type dFi_dtx, dFi_dty, dFi_dtz, dFi_dax, dFi_day, dFi_daz, dFi_dthx, dFi_dthy, dFi_dthz;
	type Fi, DX, DY, DZ, wi, DX2, DY2, DZ2, DXDY, DXDZ, DYDZ;
	type p_bari, sinx, cosx, siny, cosy, sinz, cosz, Wi;
	type NC1, NC2, NC3, NC4, NC5, NC6, NC7, NC8;
	type N00, N01, N02, N03, N04, N05, N06, N07, N08, U0;
	type N11, N12, N13, N14, N15, N16, N17, N18, U1; 
	type N22, N23, N24, N25, N26, N27, N28, U2;
	type N33, N34, N35, N36, N37, N38, U3;
	type N44, N45, N46, N47, N48 ,U4;
	type N55, N56, N57, N58, U5;
	type N66, N67, N68, U6;
	type N77, N78, U7;
	type N88, U8;
	struct cart_coord pp;	
	struct group matr;
	/* Initializing the values of the matrix N and vector u */
	N00 = N01 = N02 = N03 = N04 = N05 = N06 = N07 = N08 = U0 = N88 = U8 = 0.0L;
	N11 = N12 = N13 = N14 = N15 = N16 = N17 = N18 = U1 = N77 = N78 = U7 = 0.0L;
	N22 = N23 = N24 = N25 = N26 = N27 = N28 = U2 = N66 = N67 = N68 = U6 = 0.0L;
	N33 = N34 = N35 = N36 = N37 = N38 = U3 = N55 = N56 = N57 = N58 = U5 = 0.0L;
	N44 = N45 = N46 = N47 = N48 = U4 = matr.sum_piwi2 = 0.0L;
	/* Initial values */
	tx = values[0]; /* [m] */
	ty = values[1];	/* [m] */
	tz = values[2];	/* [m] */
	ax = values[3]; /* [m] */
	ay = values[4];	/* [m] */
	az = values[5];	/* [m] */
	theta_x = values[6]; /* [rad] */
	theta_y = values[7]; /* [rad] */
	theta_z = values[8]; /* [rad] */
	/* Calculating the elements of the rotation matrix */
	sinx = sin(theta_x);
	siny = sin(theta_y);
	sinz = sin(theta_z);
	cosx = cos(theta_x);
	cosy = cos(theta_y);
	cosz = cos(theta_z);
    r11 = cosy * cosz;
	r12 = cosx * sinz + sinx * siny * cosz;
	r13 = sinx * sinz - cosx * siny * cosz;
	r21 = -cosy * sinz;
	r22 = cosx * cosz - sinx * siny * sinz;
	r23 = sinx * cosz + cosx * siny * sinz;
	r31 = siny;
	r32 = -sinx * cosy;
	r33 = cosx * cosy;
	/* Calculating the p coefficients */
	pxx = r11 * r11 / ax / ax + r21 * r21 / ay / ay + r31 * r31 / az / az;
	pyy = r12 * r12 / ax / ax + r22 * r22 / ay / ay + r32 * r32 / az / az;
	pzz = r13 * r13 / ax / ax + r23 * r23 / ay / ay + r33 * r33 / az / az;
	pxy = r11 * r12 / ax / ax + r21 * r22 / ay / ay + r31 * r32 / az / az;
	pxz = r11 * r13 / ax / ax + r21 * r23 / ay / ay + r31 * r33 / az / az;
	pyz = r12 * r13 / ax / ax + r22 * r23 / ay / ay + r32 * r33 / az / az;
	A = 1.0L / az / az - 1.0L / ax / ax;
	B = 1.0L / ay / ay - 1.0L / az / az;
	C = 1.0L / ax / ax - 1.0L / ay / ay;
	/* Calculating the partial derivatives with respect to ax */
	dpxx_dax = -2.0L * r11 * r11 / ax / ax / ax;
	dpyy_dax = -2.0L * r12 * r12 / ax / ax / ax;
	dpzz_dax = -2.0L * r13 * r13 / ax / ax / ax;
	dpxy_dax = -2.0L * r11 * r12 / ax / ax / ax;
	dpxz_dax = -2.0L * r11 * r13 / ax / ax / ax;
	dpyz_dax = -2.0L * r12 * r13 / ax / ax / ax;
	/* Calculating the partial derivatives with respect to ay */
	dpxx_day = -2.0L * r21 * r21 / ay / ay / ay;
	dpyy_day = -2.0L * r22 * r22 / ay / ay / ay;
	dpzz_day = -2.0L * r23 * r23 / ay / ay / ay;
	dpxy_day = -2.0L * r21 * r22 / ay / ay / ay;
	dpxz_day = -2.0L * r21 * r23 / ay / ay / ay;
	dpyz_day = -2.0L * r22 * r23 / ay / ay / ay;
	/* Calculating the partial derivatives with respect to az */
	dpxx_daz = -2.0L * r31 * r31 / az / az / az;
	dpyy_daz = -2.0L * r32 * r32 / az / az / az;
	dpzz_daz = -2.0L * r33 * r33 / az / az / az;
	dpxy_daz = -2.0L * r31 * r32 / az / az / az;
	dpxz_daz = -2.0L * r31 * r33 / az / az / az;
	dpyz_daz = -2.0L * r32 * r33 / az / az / az;
	/* Calculating the partial derivatives with respect to theta_y */
	dpxx_dthy = 2.0L * (r11 * r31 * cosz * A + r21 * r31 * sinz * B);
	dpyy_dthy = 2.0L * (r12 * r32 * cosz * A + r22 * r32 * sinz * B);
	dpzz_dthy = 2.0L * (r13 * r33 * cosz * A + r23 * r33 * sinz * B);
	dpxy_dthy = (r11 * r32 + r12 * r31) * cosz * A + (r21 * r32 + r22 * r31) * sinz * B;
	dpxz_dthy = (r11 * r33 + r13 * r31) * cosz * A + (r21 * r33 + r23 * r31) * sinz * B;
	dpyz_dthy = (r12 * r33 + r13 * r32) * cosz * A + (r22 * r33 + r23 * r32) * sinz * B;
	/* Calculating the partial derivatives with respect to theta_z */
	dpxx_dthz = 2.0L * r11 * r21 * C;
	dpyy_dthz = 2.0L * r12 * r22 * C;
	dpzz_dthz = 2.0L * r13 * r23 * C;
	dpxy_dthz = (r11 * r22 + r12 * r21) * C;
	dpxz_dthz = (r11 * r23 + r13 * r21) * C;
	dpyz_dthz = (r12 * r23 + r13 * r22) * C;
	matr.c = 0;
	/* Reading the file (binary file) and calculating the elements of the N and u matrices */
	while (fread(&pp, 1, sizeof(pp), fp) > 0){
		xi = pp.x;
		yi = pp.y;
		zi = pp.z;
		pi = pp.w;
		DX = xi - tx;
		DY = yi - ty;
		DZ = zi - tz;
		DX2 = DX * DX;
		DY2 = DY * DY;
		DZ2 = DZ * DZ;
		DXDY = DX * DY;
		DXDZ = DX * DZ;
		DYDZ = DY * DZ;
		/* Calculating the partial derivatives of the function F with respect to tx, ty and tz */
		dFi_dtx = -2.0L * pxx * DX -2.0L * pxy * DY -2.0L * pxz * DZ;
		dFi_dty = -2.0L * pxy * DX -2.0L * pyy * DY -2.0L * pyz * DZ;
		dFi_dtz = -2.0L * pxz * DX -2.0L * pyz * DY -2.0L * pzz * DZ;
		/* Calculating the partial derivatives of the function F with respect to ax, ay and az */
		dFi_dax = dpxx_dax * DX2 + dpyy_dax * DY2 + dpzz_dax * DZ2 + 2 * dpxy_dax * DXDY + 2 * dpxz_dax * DXDZ + 2 * dpyz_dax * DYDZ;
		dFi_day = dpxx_day * DX2 + dpyy_day * DY2 + dpzz_day * DZ2 + 2 * dpxy_day * DXDY + 2 * dpxz_day * DXDZ + 2 * dpyz_day * DYDZ;
		dFi_daz = dpxx_daz * DX2 + dpyy_daz * DY2 + dpzz_daz * DZ2 + 2 * dpxy_daz * DXDY + 2 * dpxz_daz * DXDZ + 2 * dpyz_daz * DYDZ;
		/* Calculating the partial derivatives of the function F with respect to thetax, thetay and thetaz */
		dFi_dthx = -2.0L * pyz * DY2 + 2.0L * pyz * DZ2 - 2.0L * pxz * DXDY + 2 * pxy * DXDZ + 2.0L * pyy * DYDZ - 2.0L * pzz * DYDZ;
		dFi_dthy = dpxx_dthy * DX2 + dpyy_dthy * DY2 + dpzz_dthy * DZ2 + 2 * dpxy_dthy * DXDY + 2 * dpxz_dthy * DXDZ + 2 * dpyz_dthy * DYDZ;
		dFi_dthz = dpxx_dthz * DX2 + dpyy_dthz * DY2 + dpzz_dthz * DZ2 + 2 * dpxy_dthz * DXDY + 2 * dpxz_dthz * DXDZ + 2 * dpyz_dthz * DYDZ;
		wi = dFi_dtx * dFi_dtx + dFi_dty * dFi_dty + dFi_dtz * dFi_dtz;
		wi /= pi;
		p_bari = 1.0L / wi;
		/* Calculating the function F(tx, ty, tz, ax, ay, az, thetax, thetay, thetaz) */
		Fi = pxx * DX2 + pyy * DY2 + pzz * DZ2 + 2 * pxy * DXDY + 2 * pxz * DXDZ + 2 * pyz * DYDZ - 1.0L;
		wi = -Fi;
		Wi = wi * p_bari;
		matr.sum_piwi2 += wi * Wi;
		NC1 = dFi_dtx * p_bari;
		NC2 = dFi_dty * p_bari;
		NC3 = dFi_dtz * p_bari;
		NC4 = dFi_dax * p_bari;
		NC5 = dFi_day * p_bari;
		NC6 = dFi_daz * p_bari;
		NC7 = dFi_dthx * p_bari;
		NC8 = dFi_dthy * p_bari;
		/* Calculating each element of the matrix N */
		N00 += NC1 * dFi_dtx;
		N01 += NC1 * dFi_dty;
		N02 += NC1 * dFi_dtz;
		N03 += NC1 * dFi_dax;
		N04 += NC1 * dFi_day;
		N05 += NC1 * dFi_daz;
		N06 += NC1 * dFi_dthx;
		N07 += NC1 * dFi_dthy;
		N08 += NC1 * dFi_dthz;
		N11 += NC2 * dFi_dty;
		N12 += NC2 * dFi_dtz;
		N13 += NC2 * dFi_dax;
		N14 += NC2 * dFi_day;
		N15 += NC2 * dFi_daz;
		N16 += NC2 * dFi_dthx;
		N17 += NC2 * dFi_dthy;
		N18 += NC2 * dFi_dthz;
		N22 += NC3 * dFi_dtz;
		N23 += NC3 * dFi_dax;
		N24 += NC3 * dFi_day;
		N25 += NC3 * dFi_daz;
		N26 += NC3 * dFi_dthx;
		N27 += NC3 * dFi_dthy;
		N28 += NC3 * dFi_dthz;
		N33 += NC4 * dFi_dax;
		N34 += NC4 * dFi_day;
		N35 += NC4 * dFi_daz;
		N36 += NC4 * dFi_dthx;
		N37 += NC4 * dFi_dthy;
		N38 += NC4 * dFi_dthz;
		N44 += NC5 * dFi_day;
		N45 += NC5 * dFi_daz;
		N46 += NC5 * dFi_dthx;
		N47 += NC5 * dFi_dthy;
		N48 += NC5 * dFi_dthz;
		N55 += NC6 * dFi_daz;
		N56 += NC6 * dFi_dthx;
		N57 += NC6 * dFi_dthy;
		N58 += NC6 * dFi_dthz;
		N66 += NC7 * dFi_dthx;
		N67 += NC7 * dFi_dthy;
		N68 += NC7 * dFi_dthz;
		N77 += NC8 * dFi_dthy;
		N78 += NC8 * dFi_dthz;
		N88 += dFi_dthz * dFi_dthz * p_bari;
		/* Calculating each element of the vector u */
		U0 += dFi_dtx * Wi;
		U1 += dFi_dty * Wi;
		U2 += dFi_dtz * Wi;
		U3 += dFi_dax * Wi;
		U4 += dFi_day * Wi;
		U5 += dFi_daz * Wi;
		U6 += dFi_dthx * Wi;
		U7 += dFi_dthy * Wi;
		U8 += dFi_dthz * Wi;
		matr.c++;
	}
	/* Designing the matrix N */
	matr.N_bar[0][0] = N00;
	matr.N_bar[0][1] = N01;
	matr.N_bar[0][2] = N02;
	matr.N_bar[0][3] = N03;
	matr.N_bar[0][4] = N04;
	matr.N_bar[0][5] = N05;
	matr.N_bar[0][6] = N06;
	matr.N_bar[0][7] = N07;
	matr.N_bar[0][8] = N08;
	matr.N_bar[1][1] = N11;
	matr.N_bar[1][2] = N12;
	matr.N_bar[1][3] = N13;
	matr.N_bar[1][4] = N14;
	matr.N_bar[1][5] = N15;
	matr.N_bar[1][6] = N16;
	matr.N_bar[1][7] = N17;
	matr.N_bar[1][8] = N18;
	matr.N_bar[2][2] = N22;
	matr.N_bar[2][3] = N23;
	matr.N_bar[2][4] = N24;
	matr.N_bar[2][5] = N25;
	matr.N_bar[2][6] = N26;
	matr.N_bar[2][7] = N27;
	matr.N_bar[2][8] = N28;
	matr.N_bar[3][3] = N33;
	matr.N_bar[3][4] = N34;
	matr.N_bar[3][5] = N35;
	matr.N_bar[3][6] = N36;
	matr.N_bar[3][7] = N37;
	matr.N_bar[3][8] = N38;
	matr.N_bar[4][4] = N44;
	matr.N_bar[4][5] = N45;
	matr.N_bar[4][6] = N46;
	matr.N_bar[4][7] = N47;
	matr.N_bar[4][8] = N48;
	matr.N_bar[5][5] = N55;
	matr.N_bar[5][6] = N56;
	matr.N_bar[5][7] = N57;
	matr.N_bar[5][8] = N58;
	matr.N_bar[6][6] = N66;
	matr.N_bar[6][7] = N67;
	matr.N_bar[6][8] = N68;
	matr.N_bar[7][7] = N77;
	matr.N_bar[7][8] = N78;
	matr.N_bar[8][8] = N88;
	/* Designing the vector u */
	matr.U_bar[0] = U0;
	matr.U_bar[1] = U1;
	matr.U_bar[2] = U2;
	matr.U_bar[3] = U3;
	matr.U_bar[4] = U4;
	matr.U_bar[5] = U5;
	matr.U_bar[6] = U6;
	matr.U_bar[7] = U7;
	matr.U_bar[8] = U8;
	/* Converting the upper triangular matrix N into a symmetric one by calling the function symmetric() */
	symmetric(&matr.N_bar[0][0], 9);
	/* Returning the file position indicator to the beginning of the file */
	rewind(fp);
	return matr;
}

