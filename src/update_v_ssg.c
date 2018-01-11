/*------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS3D.
 *
 * IFOS3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * IFOS3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with IFOS3D. See file COPYING and/or
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   updating velocity values by a staggered grid finite difference scheme of
 *   nth order accuracy in space and second order accuracy in time
 *  ----------------------------------------------------------------------*/

#include "fd.h"





extern SLAVE_FUN(update_v_kernel_fusion_slave)(Param_vel *);
extern SLAVE_FUN(update_v_kernel_absorb_salve)(Param_absorb *);

double update_v(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt, float *Fv, float *Fs, float  ***  rho,  float  *** rjp, float  *** rkp, float  *** rip, float **  srcpos_loc, float ** signals, float ** signaly, float ** signalz, int nsrc, float *** absorb_coeff, int back) {


	/*extern FILE *FP;*/
	extern float DT, DX, DY, DZ, ALPHA, BETA;
	extern int NX, NY, NZ;
	double time=0.0; /*, time1=0.0;*/
	/*double time2=0.0;*/
	extern int  FDORDER,  ABS_TYPE, FDCOEFF; extern int MYID;  // ,LOG,*/



	int i, j, k, l;
	float  amp, alpha_rad, beta_rad;
	float b1, b2, b3, b4, b5, b6, dx, dy, dz;
	float sxx_x, sxy_y, sxz_z, syy_y, sxy_x, syz_z;
	float szz_z, sxz_x, syz_y;

	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2; printf("ABS:%d, FDOR:%d\n", ABS_TYPE, FDORDER);}
	 //printf("debugging MYID:%d %s: ------------------------------------------------------------->\n", MYID,  __FUNCTION__);
    /************************************************* part of array fusion  ****************************************/

	//int nrow = NY + 2 * (ll*FDORDER/2) + 1;
	int ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
    //int size = nrow * ncol * ndep;
	float *Frp = transform_3(rip, rjp, rkp, 1, NY, 1, NX, 1, NZ);
	float *Fabsorb = transform3(absorb_coeff, 1, NY, 1, NX, 1, NZ);
	int strip = ndep;
	int slice = ncol * ndep;
	int strip_rp = NZ;
	int slice_rp = NZ * NX;




	int dim_x = nx2 - nx1 + 1;
	int dim_y = ny2 - ny1 + 1;
	int dim_z = nz2 - nz1 + 1;
	//if(MYID == 0)
		//printf("---###-----------------------------------> nz1: %d, nz2:%d dim_z: %d\n", nz1, nz2, dim_z);


	Param_vel param;


	/*unsigned int * test_float;*/



	/*if (LOG)
	if (MYID==0) time1=MPI_Wtime();*/


    switch (FDORDER){

	case 4 :

	    dx=DT/DX;
		dy=DT/DY;
		dz=DT/DZ;
		b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients*/
		if(FDCOEFF==2){
		b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/

		float *vel = Fv + 3 * (ny1 * slice + nx1 * strip + nz1);
		float *stress = Fs + 6 * (ny1 * slice + nx1 * strip + nz1);
		float *rp = Frp + 3 * (ny1 * slice_rp + nx1 * strip_rp + nz1) ;


		param.dim_x = dim_x;
		param.dim_y = dim_y;
		param.dim_z = dim_z;
		param.slice = slice;
		param.strip = strip;
		param.slice_rp = slice_rp;
		param.strip_rp = strip_rp;
		param.dx = dx;
		param.dy = dy;
		param.dz = dz;
		param.FDCOEFF = FDCOEFF;
		param.vel = vel;
		param.stress = stress;
		param.rp = rp;


	  	athread_spawn(update_v_kernel_fusion_slave, &param);
	  	athread_join();


	  	break;



     }
	/* Adding body force components to corresponding particle velocities */



	if(back==0){
	for (l=1;l<=nsrc;l++) {
		i=(int)srcpos_loc[1][l];
		j=(int)srcpos_loc[2][l];
		k=(int)srcpos_loc[3][l];
		int idx = j * slice + i * strip + k;
		int idx_rp = j * slice_rp + i * strip_rp + k;
		amp=DT*signals[l][nt]/(DX*DY*DZ);


		switch ((int)srcpos_loc[7][l]){
		case 2 :
			Fv[idx * 3 + 0]+=amp*Frp[idx_rp * 3 + 0];

			/*/(rip[j][i][k]);  *//* single force in x , DX^3 because of force density*/

			/*test_float = (unsigned int*) (void*) &vx[j][i][k];
			if(((*test_float & 0x7f800000) == 0x0) && ((*test_float & 0x007fffff) != 0x0))fprintf(FP,"Achtung:ILLEGALER FLOAT");*/

			break;
		case 3 :
			Fv[idx * 3 + 2]+=amp*Frp[idx_rp * 3 + 1];  /* single force in z  */
			break;
		case 4 :
			Fv[idx * 3 + 1]+=amp*Frp[idx_rp * 3 + 2];  /* single force in y, vertical direction*/
			break;
		case 5 :
			alpha_rad=ALPHA*PI/180; /* custom force */
			beta_rad=BETA*PI/180;
			Fv[idx * 3 + 0]+=cos(beta_rad)*cos(alpha_rad)*amp;
			Fv[idx * 3 + 1]+=sin(alpha_rad)*amp;
			Fv[idx * 3 + 2]+=sin(beta_rad)*cos(alpha_rad)*amp;
			break;
		}
	}



	}

	if (back==1){

		for (l=1;l<=nsrc;l++) {
			i=(int)srcpos_loc[1][l];
			j=(int)srcpos_loc[2][l];
			k=(int)srcpos_loc[3][l];
			int idx = j * slice + i * strip + k;
			Fv[idx * 3 + 0]+=signals[l][nt];
			Fv[idx * 3 + 2]+=signalz[l][nt];
			Fv[idx * 3 + 1]+=signaly[l][nt];
		}
	}




	/* absorbing boundary condition (exponential damping) */


/*
	Param_absorb param1;

	param1.stress = Fs + 6 * (ny1 * slice + nx1 * strip + nz1);
	param1.vel = Fv + 3 * (ny1 * slice + nx1 * strip + nz1);
	param1.absorb = Fabsorb + (ny1 * slice_rp + nx1 * strip_rp + nz1);
	param1.dimx = dim_x;
	param1.dimy = dim_y;
	param1.dimz = dim_z;
	param1.slice = slice;
	param1.strip = strip;
	param1.slice_a = slice_rp;
	param1.strip_a = strip_rp;

	if (ABS_TYPE==2){
	  	athread_spawn(update_v_kernel_absorb_salve, &param1);
	  	athread_join();
	}
*/

	if (ABS_TYPE==2){
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){
					int idx = j * slice + i * strip + k;
					int idx_absorb = j * slice_rp + i * strip_rp + k;
					Fv[idx * 3 + 0]*=Fabsorb[idx_absorb];
					Fv[idx * 3 + 1]*=Fabsorb[idx_absorb];
					Fv[idx * 3 + 2]*=Fabsorb[idx_absorb];

					Fs[idx * 6 + 0]*=Fabsorb[idx_absorb];
					Fs[idx * 6 + 1]*=Fabsorb[idx_absorb];
					Fs[idx * 6 + 2]*=Fabsorb[idx_absorb];
					Fs[idx * 6 + 3]*=Fabsorb[idx_absorb];
					Fs[idx * 6 + 4]*=Fabsorb[idx_absorb];
					Fs[idx * 6 + 5]*=Fabsorb[idx_absorb];

				}
			}
		}
	}

    free_trans_3(Frp, 1, NY, 1, NX, 1, NZ );
    free_trans(Fabsorb, 1, NY, 1, NX, 1, NZ);


	return time;

}
