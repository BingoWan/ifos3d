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
 *   updating stress values by a staggered grid finite difference scheme of 
 *   nth order accuracy in space and second order accuracy in time
 *   elastic version
 *  ----------------------------------------------------------------------*/

#include "fd.h"




extern SLAVE_FUN(update_s_kernel_fusion_slave)(void *);

double update_s_elastic(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, float *Fv, float *Fs, float ***  pi, float ***  u, float ***  uipjp, float ***  ujpkp, float ***  uipkp) {


	extern float DT, DX, DY, DZ;
	extern int  FDORDER,  FDCOEFF, MYID, ABS_TYPE;//,LOG,
	extern int NX, NY, NZ;

	int i, j, k;
	double time=0.0; double time1=0.0; double time2 = 0.0;
	float vxx,vxy,vxz,vyx,vyy,vyz,vzx,vzy,vzz;
	float vxyyx,vyzzy,vxzzx,vxxyyzz,vyyzz,vxxzz,vxxyy;
	float g,f,fipjp,fjpkp,fipkp;

	float b1, b2, b3, b4, b5, b6;
	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2; printf("ABS:%d, FDOR:%d\n", ABS_TYPE, FDORDER);}
	//printf("debugging MYID:%d %s: ------------------------------------------------------------->\n", MYID,  __FUNCTION__);
	/************************************************* part of array fusion  ****************************************/
	//if (LOG)
	if (MYID==0) time1=MPI_Wtime();

	int nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
	int size = nrow * ncol * ndep;
    float *Fpi = transform3(pi, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    float *Fu = transform3(u, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    float *Fup = transform_3(uipjp, ujpkp, uipkp, 1, NY, 1, NX, 1, NZ);
	int strip = ndep;
	int slice = ncol * ndep;
	int strip_up = NZ;
	int slice_up = NX * NZ;
	int strip_pi = NZ + 1;
	int slice_pi = (NX + 1) * (NZ + 1);



	static int init = 1;



	float *vel = Fv + 3 * (ny1 * slice + nx1 * strip + nz1);
	float *stress = Fs + 6 * (ny1 * slice + nx1 * strip + nz1);
	float *up = Fup + 3 * (ny1 * slice_up + nx1 * strip_up + nz1);
	float *Pi = Fpi + (ny1 * slice_pi + nx1 * strip_pi + nz1);
	float *U = Fu + (ny1 * slice_pi + nx1 * strip_pi + nz1);

	int dim_x = nx2 - nx1 + 1;
	int dim_y = ny2 - ny1 + 1;
	int dim_z = nz2 - nz1 + 1;


	Param_str param;


	switch (FDORDER){


	case 4 :

		b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients */
		if(FDCOEFF==2){
			b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/


		/*
		param.dim_x = dim_x;
		param.dim_y = dim_y;
		param.dim_z = dim_z;
		param.slice = slice;
		param.strip = strip;
		param.slice_up = slice_up;
		param.strip_up = strip_up;
		param.slice_pi = slice_pi;
		param.strip_pi = strip_pi;
		param.DT = DT;
		param.DX = DX;
		param.DY = DY;
		param.DZ = DZ;
		param.FDCOEFF = FDCOEFF;
		param.vel = vel;
		param.stress = stress;
		param.up = up;
		param.pi = Pi;
		param.u = U;

		//signal(SIGUSR1, pc_on_sig);
		athread_spawn(update_s_kernel_fusion_slave, &param);
		athread_join();

		break;
		*/

		for (j = ny1; j <=  ny2; j++){
			for (i = nx1; i <= nx2; i++){
				for (k = nz1; k <= nz2; k++){
					int idx = j * slice + i * strip + k;
					int idx_n = j * slice_up + i * strip_up + k;
					int idx_p = j * slice_pi + i * strip_pi + k;
					vxx = (b1 * (Fv[idx * 3 + 0] - Fv[(idx - strip) * 3 + 0]) + b2 * (Fv[(idx + strip) * 3 + 0] - Fv[(idx - 2 * strip) * 3 + 0])) / DX;
					vxy = (b1 * (Fv[(idx + slice) * 3 + 0] - Fv[idx * 3 + 0]) + b2 * (Fv[(idx + 2 * slice) * 3 + 0] - Fv[(idx - slice) * 3 + 0])) / DY;
					vxz = (b1 * (Fv[(idx + 1) * 3 + 0] - Fv[idx * 3 + 0]) + b2 * (Fv[(idx + 2) * 3 + 0] - Fv[(idx -1) * 3 + 0])) / DZ;
					vyx = (b1 * (Fv[(idx + strip) * 3 + 1] - Fv[idx * 3 + 1]) + b2 * (Fv[(idx + 2 * strip) * 3 + 1] - Fv[(idx - strip) * 3 + 1])) / DX;
					vyy = (b1 * (Fv[idx * 3 + 1] - Fv[(idx - slice) * 3 + 1]) + b2 * (Fv[(idx + slice) * 3 + 1] - Fv[(idx - 2 * slice) * 3 + 1])) / DY;
					vyz = (b1 * (Fv[(idx + 1) * 3 + 1] - Fv[idx * 3 + 1]) + b2 * (Fv[(idx + 2) * 3 + 1] - Fv[(idx - 1) * 3 + 1])) / DZ;
					vzx = (b1 * (Fv[(idx + strip) * 3 + 2] - Fv[idx * 3 + 2]) + b2 * (Fv[(idx + 2 * strip) * 3 + 2] - Fv[(idx - strip) * 3 + 2])) / DX;
					vzy = (b1 * (Fv[(idx + slice) * 3 + 2] - Fv[idx * 3 + 2]) + b2 * (Fv[(idx + 2 * slice) * 3 + 2] - Fv[(idx - slice) * 3 + 2])) / DY;
					vzz = (b1 * (Fv[idx * 3 + 2] - Fv[(idx - 1) * 3 + 2]) + b2 * (Fv[(idx + 1) * 3 + 2] - Fv[(idx - 2) * 3 + 2])) / DZ;


					// updating components of the stress tensor, partially
					fipjp=Fup[idx_n * 3 + 0] * DT;
					fjpkp=Fup[idx_n * 3 + 1] * DT;
					fipkp=Fup[idx_n * 3 + 2] * DT;
					g = Fpi[idx_p];
					f = 2.0*Fu[idx_p];


					vxyyx=vxy+vyx;
					vyzzy=vyz+vzy;
					vxzzx=vxz+vzx;
					vxxyyzz=vxx+vyy+vzz;
					vyyzz=vyy+vzz;
					vxxzz=vxx+vzz;
					vxxyy=vxx+vyy;


					Fs[6 * idx + 0] += DT*((g*vxxyyzz)-(f*vyyzz));
					Fs[6 * idx + 1] += DT*((g*vxxyyzz)-(f*vxxzz));
					Fs[6 * idx + 2] += DT*((g*vxxyyzz)-(f*vxxyy));
					Fs[6 * idx + 3] += (fipjp*vxyyx);
					Fs[6 * idx + 4] += (fjpkp*vyzzy);
					Fs[6 * idx + 5] += (fipkp*vxzzx);
				}
			}
		}

		break;

	}



    free_trans_3(Fup, 1, NY, 1, NX, 1, NZ);
    free_trans(Fpi, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    free_trans(Fu, 1, NY + 1, 1, NX + 1, 1, NZ + 1);



	//if (LOG)
	if (MYID==0){
		time2=MPI_Wtime();
		time=time2-time1;
		//fprintf(FP," Real time for stress tensor update: \t\t %4.2f s.\n",time);
	}

	return time;

}
