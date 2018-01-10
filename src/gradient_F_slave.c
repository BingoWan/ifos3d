/*
 * gradient_F_slave.c
 *
 *  Created on: Dec 19, 2017
 *      Author: bingo
 */



#include <slave.h>
#include <dma.h>
#include <math.h>

#define wx 10
#define wy 1
#define wz 10
#define MX 64
//#define DEBUG

typedef struct {
	int dim_x;
	int dim_y;
	int dim_z;
	int cube;
	int slice;
	int strip;
	int slice_pi;
	int strip_pi;
	int slice_grad;
	int strip_grad;
	float DT;
	float DX;
	float DY;
	float DZ;
	int FDCOEFF;
	int nf;
	float *grad;
	float *fvel;
	float *bvel;
	float *pi;
	float *u;
	float *rho;
	float *finv;
	int MYID;
} Param_grad;


inline  addn(int x, int n, int x1) {
	return x + n < x1 ? x + n : x + n - x1;
}


void gradient_kernel_slave(Param_grad *param) {

	int id = athread_get_id(-1);

	int dim_x = param->dim_x;
	int dim_y = param->dim_y;
	int dim_z = param->dim_z;
	int cube = param->cube;
	int ystep = param->slice;
	int xstep = param->strip;
	int ystep_pi = param->slice_pi;
	int xstep_pi = param->strip_pi;
	int ystep_grad = param->slice_grad;
	int xstep_grad = param->strip_grad;
	float DT = param->DT;
	float DX = param->DX;
	float DY = param->DY;
	float DZ = param->DZ;
	int FDCOEFF = param->FDCOEFF;
	int nf = param->nf;

	float *grad_s0 = param->grad;
	float *fvel_s0 = param->fvel;
	float *bvel_s0 = param->bvel;
	float *pi = param->pi;
	float *u = param->u;
	float *rho = param->rho;
	float *finv = param->finv;
	int MYID = param->MYID;



	int NX = (dim_x + wx * MX - 1) / (wx * MX);
	int NZ = (dim_z + wz - 1) / wz;



	const int x0 = 2;
	const int x1 = 2;
	const int y0 = 2;
	const int y1 = 2;
	const int z0 = 2;
	const int z1 = 2;

	int l, iix, iiz, ix, iy, iz;
	int izbeg, izend, izn, ixbeg, ixend, ixn;

	int xstep_6 = xstep * 6;
	int ystep_6 = ystep * 6;
	int xstep_3 = xstep * 3;
	int ystep_3 = ystep * 3;
	int ystep_grad_3 = ystep_grad * 3;
	int xstep_grad_3 = xstep_grad * 3;
	int xstep_l = wz + z0 + z1;
	int ystep_l = (wx + x0 + x1) * (wz + z0 + z1);
	int xstep_l2 = xstep_l * 2;
	int ystep_l3 = ystep_l * 3;
	int xstep_l3 = xstep_l * 3;
	int ystep_l6 = ystep_l * 6;
	int xstep_l6 = xstep_l * 6;


	float *fvel_s_tmp, *fvel_s1, *bvel_s1, * bvel_s_tmp, *grad_s1, *grad_s1_tmp;
    float *pi_s1, *pi_s1_tmp, *u_s1, *u_s1_tmp, *rho_s1, *rho_s1_tmp;
    float fvel_l[(wy + y0 + y1) * (wx + x0 + x1) * (wz + z0 + z1) * 6];
    float bvel_l[(wy + y0 + y1) * (wx + x0 + x1) * (wz + z0 + z1) * 6];
    float grad_l[wx * wz * 3];
    float u_l[wx * wz], pi_l[wx * wz], rho_l[wx * wz];
    float *fvel_l_tmp, *fvel_l_c, *bvel_l_tmp, *bvel_l_c, *grad_l_tmp;
    float *u_l_tmp, *pi_l_tmp, *rho_l_tmp;
    fvel_l_c = fvel_l + x0 * xstep_l6 + z0 * 6;
    bvel_l_c = bvel_l + x0 * xstep_l6 + z0 * 6;

    volatile unsigned long get_reply, put_reply;


    int plane_comp = 0;

	float finv_l[nf];
	float fvxx=0.0,fvxy=0.0,fvxz=0.0,fvyx=0.0,fvyy=0.0,fvyz=0.0,fvzx=0.0,fvzy=0.0,fvzz=0.0;
	float bvxx=0.0,bvxy=0.0,bvxz=0.0,bvyx=0.0,bvyy=0.0,bvyz=0.0,bvzx=0.0,bvzy=0.0,bvzz=0.0;
	float fivxx=0.0,fivxy=0.0,fivxz=0.0,fivyx=0.0,fivyy=0.0,fivyz=0.0,fivzx=0.0,fivzy=0.0,fivzz=0.0;
	float bivxx=0.0,bivxy=0.0,bivxz=0.0,bivyx=0.0,bivyy=0.0,bivyz=0.0,bivzx=0.0,bivzy=0.0,bivzz=0.0;
	float gradlam, gradmu, gradrho;
	float b1,b2,fdummy;


	//To obtain data of finv
	get_reply = 0;
	athread_get(PE_MODE, finv, finv_l, nf * sizeof(float), &get_reply, 0, 0, 0);
	while(get_reply!=1);





	b1 = 9.0 / 8.0; b2 = -1.0 / 24.0; /* Taylor coefficients (4th order)*/
	if(FDCOEFF==2){
		b1=1.1382; b2=-0.046414;
	} /* Holberg coefficients E=0.1 %*/


#ifdef DEBUG
	int ldm_finv = sizeof(float) * nf;
	int ldm_fvel = sizeof(float) * (wx + x0 + x1) * (wy + y0 + y1) * (wz + z0 + z1) * 6;
	int ldm_bvel = sizeof(float) * (wx + x0 + x1) * (wy + y0 + y1) * (wz + z0 + z1) * 6;
	int ldm_grad = sizeof(float) * wx * wz * 3;
	int ldm_rho = sizeof(float) * wx * wz;
	int ldm_pi = sizeof(float) * wx * wz;
	int ldm_u = sizeof(float) * wx * wz;
	int ldm_total = ldm_fvel + ldm_bvel + ldm_grad + ldm_rho + ldm_pi + ldm_u + ldm_finv;
	if(id == 0) {
		//printf("LDM fvel = %d Bytes\n", ldm_fvel);
		//printf("LDM bvel = %d Bytes\n", ldm_bvel);
		//printf("LDM grad = %d Bytes\n", ldm_grad);
		//printf("LDM rho = %d Bytes\n", ldm_rho);
		//printf("LDM pi = %d Bytes\n", ldm_pi);
		//printf("LDM u = %d Bytes\n", ldm_u);
		printf("LDM total of gradient_F = %d Bytes\n", ldm_total);
  }
#endif

	for (l = 1; l <= nf; l++) {
		fdummy = 0.0;
		fdummy = 2.0 * 2.0 * finv_l[l - 1] * finv_l[l - 1] * M_PI * M_PI;
		fdummy=1/fdummy;
		for (iiz = 0; iiz < NZ; iiz++) {
			izbeg = wz * iiz;
			izend = wz * (iiz + 1);
			izend = izend < dim_z ? izend : dim_z;
			izn = izend - izbeg;

			for (iix = 0; iix < NX; iix++) {
				ixbeg = wx * MX * iix + wx * id;
				ixend = wx * MX * iix + wx * (id + 1);
				ixend = ixend < dim_x ? ixend: dim_x;
				ixn = ixend - ixbeg;

				fvel_s1 = fvel_s0 + (l * cube + ixbeg * xstep + izbeg) * 6;
                bvel_s1 = bvel_s0 + (l * cube + ixbeg * xstep + izbeg) * 6;
                grad_s1 = grad_s0 + (ixbeg * xstep_grad + izbeg) * 3;
                u_s1 = u + (ixbeg * xstep_pi + izbeg);
                pi_s1 = pi + (ixbeg * xstep_pi + izbeg);
                rho_s1 = rho + (ixbeg * xstep_pi + izbeg);


                //To obtain data from main memory(fvel[-y0][-x0][-z0]) to local memory(fvel_l[y0][x0][z0]),  every data size is (wz + x0 + x1) * 3.
                fvel_l_tmp = fvel_l;
                bvel_l_tmp = bvel_l;
				for (iy = -y0; iy < y1 + wy; iy++) {
					fvel_s_tmp = fvel_s1 + iy * ystep_6 - x0 * xstep_6 - z0 * 6;
					bvel_s_tmp = bvel_s1 + iy * ystep_6 - x0 * xstep_6 - z0 * 6;
					for (ix = -x0; ix < x1 + wx; ix++) {
						get_reply = 0;
						athread_get(PE_MODE, fvel_s_tmp, fvel_l_tmp, (wz + z0 + z1) * 6 * sizeof(float), &get_reply, 0, 0, 0);
						athread_get(PE_MODE, bvel_s_tmp, bvel_l_tmp, (wz + z0 + z1) * 6 * sizeof(float), &get_reply, 0, 0, 0);
						while(get_reply!=2);
                        fvel_s_tmp += xstep_6;
                        fvel_l_tmp += xstep_l6;
                        bvel_s_tmp += xstep_6;
                        bvel_l_tmp += xstep_l6;

					}
				}




				plane_comp = 0;
				// compute velocity from iy plane to iy++ plane

				for (iy = 0; iy < dim_y; iy++) {
					//To obtain gradient data from main memory(grad[iy][0]-grad[iy][wz] to local memory (grad_l[iy][0]-grad_l[iy][wz]) and so as data of u, pi, rho.
					grad_s1_tmp = grad_s1 + iy * ystep_grad_3;
					grad_l_tmp = grad_l;
					u_s1_tmp = u_s1 + iy * ystep_pi;
					u_l_tmp = u_l;
					pi_s1_tmp = pi_s1 + iy * ystep_pi;
					pi_l_tmp = pi_l;
					rho_s1_tmp = rho_s1 + iy * ystep_pi;
					rho_l_tmp = rho_l;


					for (ix = 0; ix < ixn; ix ++) {
						get_reply = 0;
						athread_get(PE_MODE, grad_s1_tmp, grad_l_tmp, wz * 3 * sizeof(float), &get_reply, 0, 0, 0);
						athread_get(PE_MODE, rho_s1_tmp, rho_l_tmp, wz * sizeof(float), &get_reply, 0, 0, 0);
						athread_get(PE_MODE, u_s1_tmp, u_l_tmp, wz * sizeof(float), &get_reply, 0, 0, 0);
						athread_get(PE_MODE, pi_s1_tmp, pi_l_tmp, wz * sizeof(float), &get_reply, 0, 0, 0);
						while(get_reply!=4);
						grad_s1_tmp += xstep_grad_3;
						grad_l_tmp += wz * 3;
						rho_s1_tmp += xstep_pi;
						rho_l_tmp += wz;
						u_s1_tmp += xstep_pi;
						u_l_tmp += wz;
						pi_s1_tmp += xstep_pi;
						pi_l_tmp += wz;

					}


					//compute velocity at iy plane.
					for (ix = 0; ix < ixn; ix++) {
						for (iz = 0; iz < izn; iz++) {
							int ym2 = plane_comp;
							int ym1 = addn(plane_comp, 1, wy + y0 + y1);
							int y = addn(plane_comp, y0 , wy + y0 + y1);
							int yp1 = addn(plane_comp, y0 + 1, wy + y0 + y1);
							int yp2 = addn(plane_comp, y0 + 2, wy + y0 + y1);

							/*computing positon*/
							int pos_l = y * ystep_l + ix * xstep_l + iz;

							/* z direction */
							int pos_l_km2 = pos_l - 2;
							int pos_l_km1 = pos_l - 1;
							int pos_l_kp1 = pos_l + 1;
							int pos_l_kp2 = pos_l + 2;

							/* x direction */
							int pos_l_im2 = pos_l - xstep_l2;
							int pos_l_im1 = pos_l - xstep_l;
							int pos_l_ip1 = pos_l + xstep_l;
							int pos_l_ip2 = pos_l + xstep_l2;

							/* y direction */
							int pos_l_jm2 = ym2 * ystep_l + ix * xstep_l + iz;
							int pos_l_jm1 = ym1 * ystep_l + ix * xstep_l + iz;
							int pos_l_jp1 = yp1 * ystep_l + ix * xstep_l + iz;
							int pos_l_jp2 = yp2 * ystep_l + ix * xstep_l + iz;

							int pos_grad = ix * wz + iz;
							int pos_pi = ix * wz + iz;




							fvxx = (b1 * (fvel_l_c[pos_l * 6] - fvel_l_c[pos_l_im1 * 6]) + b2 * (fvel_l_c[pos_l_ip1 * 6] - fvel_l_c[pos_l_im2 * 6])) / DX;
							fvxy = (b1 * (fvel_l_c[pos_l_jp1 * 6] - fvel_l_c[pos_l * 6]) + b2 * (fvel_l_c[pos_l_jp2 * 6] - fvel_l_c[pos_l_jm1 * 6])) / DY;
							fvxz = (b1 * (fvel_l_c[pos_l_kp1 * 6] - fvel_l_c[pos_l * 6]) + b2 * (fvel_l_c[pos_l_kp2 * 6] - fvel_l_c[pos_l_km1 * 6])) / DZ;
							fvyx = (b1 * (fvel_l_c[pos_l_ip1 * 6 + 1] - fvel_l_c[pos_l * 6 + 1]) + b2 * (fvel_l_c[pos_l_ip2 * 6 + 1] - fvel_l_c[pos_l_im1 * 6 + 1])) / DX;
							fvyy = (b1 * (fvel_l_c[pos_l * 6 + 1] - fvel_l_c[pos_l_jm1 * 6 + 1]) + b2 * (fvel_l_c[pos_l_jp1 * 6 + 1] - fvel_l_c[pos_l_jm2 * 6 + 1])) / DY;
							fvyz = (b1 * (fvel_l_c[pos_l_kp1 * 6 + 1] - fvel_l_c[pos_l * 6 + 1]) + b2 * (fvel_l_c[pos_l_kp2 * 6 + 1] - fvel_l_c[pos_l_km1 * 6 + 1])) / DZ;
							fvzx = (b1 * (fvel_l_c[pos_l_ip1 * 6 + 2] - fvel_l_c[pos_l * 6 + 2]) + b2 * (fvel_l_c[pos_l_ip2 * 6 + 2] - fvel_l_c[pos_l_im1 * 6 + 2])) / DX;
							fvzy = (b1 * (fvel_l_c[pos_l_jp1 * 6 + 2] - fvel_l_c[pos_l * 6 + 2]) + b2 * (fvel_l_c[pos_l_jp2 * 6 + 2] - fvel_l_c[pos_l_jm1 * 6 + 2])) / DY;
							fvzz = (b1 * (fvel_l_c[pos_l * 6 + 2] - fvel_l_c[pos_l_km1 * 6 + 2]) + b2 * (fvel_l_c[pos_l_kp1 * 6 + 2] - fvel_l_c[pos_l_km2 * 6 + 2])) / DY;

							fivxx = (b1 * (fvel_l_c[pos_l * 6 + 3] - fvel_l_c[pos_l_im1 * 6 + 3]) + b2 * (fvel_l_c[pos_l_ip1 * 6 + 3] - fvel_l_c[pos_l_im2 * 6 + 3])) / DX;
							fivxy = (b1 * (fvel_l_c[pos_l_jp1 * 6 + 3] - fvel_l_c[pos_l * 6 + 3]) + b2 * (fvel_l_c[pos_l_jp2 * 6 + 3] - fvel_l_c[pos_l_jm1 * 6 + 3])) / DY;
							fivxz = (b1 * (fvel_l_c[pos_l_kp1 * 6 + 3] - fvel_l_c[pos_l * 6 + 3]) + b2 * (fvel_l_c[pos_l_kp2 * 6 + 3] - fvel_l_c[pos_l_km1 * 6 + 3])) / DZ;
							fivyx = (b1 * (fvel_l_c[pos_l_ip1 * 6 + 4] - fvel_l_c[pos_l * 6 + 4]) + b2 * (fvel_l_c[pos_l_ip2 * 6 + 4] - fvel_l_c[pos_l_im1 * 6 + 4])) / DX;
							fivyy = (b1 * (fvel_l_c[pos_l * 6 + 4] - fvel_l_c[pos_l_jm1 * 6 + 4]) + b2 * (fvel_l_c[pos_l_jp1 * 6 + 4] - fvel_l_c[pos_l_jm2 * 6 + 4])) / DY;
							fivyz = (b1 * (fvel_l_c[pos_l_kp1 * 6 + 4] - fvel_l_c[pos_l * 6 + 4]) + b2 * (fvel_l_c[pos_l_kp2 * 6 + 4] - fvel_l_c[pos_l_km1 * 6 + 4])) / DZ;
							fivzx = (b1 * (fvel_l_c[pos_l_ip1 * 6 + 5] - fvel_l_c[pos_l * 6 + 5]) + b2 * (fvel_l_c[pos_l_ip2 * 6 + 5] - fvel_l_c[pos_l_im1 * 6 + 5])) / DX;
							fivzy = (b1 * (fvel_l_c[pos_l_jp1 * 6 + 5] - fvel_l_c[pos_l * 6 + 5]) + b2 * (fvel_l_c[pos_l_jp2 * 6 + 5] - fvel_l_c[pos_l_jm1 * 6 + 5])) / DY;
							fivzz = (b1 * (fvel_l_c[pos_l * 6 + 5] - fvel_l_c[pos_l_km1 * 6 + 5]) + b2 * (fvel_l_c[pos_l_kp1 * 6 + 5] - fvel_l_c[pos_l_km2 * 6 + 5])) / DY;

							bvxx = (b1 * (bvel_l_c[pos_l * 6] - bvel_l_c[pos_l_im1 * 6]) + b2 * (bvel_l_c[pos_l_ip1 * 6] - bvel_l_c[pos_l_im2 * 6])) / DX;
							bvxy = (b1 * (bvel_l_c[pos_l_jp1 * 6] - bvel_l_c[pos_l * 6]) + b2 * (bvel_l_c[pos_l_jp2 * 6] - bvel_l_c[pos_l_jm1 * 6])) / DY;
							bvxz = (b1 * (bvel_l_c[pos_l_kp1 * 6] - bvel_l_c[pos_l * 6]) + b2 * (bvel_l_c[pos_l_kp2 * 6] - bvel_l_c[pos_l_km1 * 6])) / DZ;
							bvyx = (b1 * (bvel_l_c[pos_l_ip1 * 6 + 1] - bvel_l_c[pos_l * 6 + 1]) + b2 * (bvel_l_c[pos_l_ip2 * 6 + 1] - bvel_l_c[pos_l_im1 * 6 + 1])) / DX;
							bvyy = (b1 * (bvel_l_c[pos_l * 6 + 1] - bvel_l_c[pos_l_jm1 * 6 + 1]) + b2 * (bvel_l_c[pos_l_jp1 * 6 + 1] - bvel_l_c[pos_l_jm2 * 6 + 1])) / DY;
							bvyz = (b1 * (bvel_l_c[pos_l_kp1 * 6 + 1] - bvel_l_c[pos_l * 6 + 1]) + b2 * (bvel_l_c[pos_l_kp2 * 6 + 1] - bvel_l_c[pos_l_km1 * 6 + 1])) / DZ;
							bvzx = (b1 * (bvel_l_c[pos_l_ip1 * 6 + 2] - bvel_l_c[pos_l * 6 + 2]) + b2 * (bvel_l_c[pos_l_ip2 * 6 + 2] - bvel_l_c[pos_l_im1 * 6 + 2])) / DX;
							bvzy = (b1 * (bvel_l_c[pos_l_jp1 * 6 + 2] - bvel_l_c[pos_l * 6 + 2]) + b2 * (bvel_l_c[pos_l_jp2 * 6 + 2] - bvel_l_c[pos_l_jm1 * 6 + 2])) / DY;
							bvzz = (b1 * (bvel_l_c[pos_l * 6 + 2] - bvel_l_c[pos_l_km1 * 6 + 2]) + b2 * (bvel_l_c[pos_l_kp1 * 6 + 2] - bvel_l_c[pos_l_km2 * 6 + 2])) / DY;

							bivxx = (b1 * (bvel_l_c[pos_l * 6 + 3] - bvel_l_c[pos_l_im1 * 6 + 3]) + b2 * (bvel_l_c[pos_l_ip1 * 6 + 3] - bvel_l_c[pos_l_im2 * 6 + 3])) / DX;
							bivxy = (b1 * (bvel_l_c[pos_l_jp1 * 6 + 3] - bvel_l_c[pos_l * 6 + 3]) + b2 * (bvel_l_c[pos_l_jp2 * 6 + 3] - bvel_l_c[pos_l_jm1 * 6 + 3])) / DY;
							bivxz = (b1 * (bvel_l_c[pos_l_kp1 * 6 + 3] - bvel_l_c[pos_l * 6 + 3]) + b2 * (bvel_l_c[pos_l_kp2 * 6 + 3] - bvel_l_c[pos_l_km1 * 6 + 3])) / DZ;
							bivyx = (b1 * (bvel_l_c[pos_l_ip1 * 6 + 4] - bvel_l_c[pos_l * 6 + 4]) + b2 * (bvel_l_c[pos_l_ip2 * 6 + 4] - bvel_l_c[pos_l_im1 * 6 + 4])) / DX;
							bivyy = (b1 * (bvel_l_c[pos_l * 6 + 4] - bvel_l_c[pos_l_jm1 * 6 + 4]) + b2 * (bvel_l_c[pos_l_jp1 * 6 + 4] - bvel_l_c[pos_l_jm2 * 6 + 4])) / DY;
							bivyz = (b1 * (bvel_l_c[pos_l_kp1 * 6 + 4] - bvel_l_c[pos_l * 6 + 4]) + b2 * (bvel_l_c[pos_l_kp2 * 6 + 4] - bvel_l_c[pos_l_km1 * 6 + 4])) / DZ;
							bivzx = (b1 * (bvel_l_c[pos_l_ip1 * 6 + 5] - bvel_l_c[pos_l * 6 + 5]) + b2 * (bvel_l_c[pos_l_ip2 * 6 + 5] - bvel_l_c[pos_l_im1 * 6 + 5])) / DX;
							bivzy = (b1 * (bvel_l_c[pos_l_jp1 * 6 + 5] - bvel_l_c[pos_l * 6 + 5]) + b2 * (bvel_l_c[pos_l_jp2 * 6 + 5] - bvel_l_c[pos_l_jm1 * 6 + 5])) / DY;
							bivzz = (b1 * (bvel_l_c[pos_l * 6 + 5] - bvel_l_c[pos_l_km1 * 6 + 5]) + b2 * (bvel_l_c[pos_l_kp1 * 6 + 5] - bvel_l_c[pos_l_km2 * 6 + 5])) / DY;


							gradlam=0.0;
							gradlam=(fvxx+fvyy+fvzz)*(bvxx+bvyy+bvzz)+(fivxx+fivyy+fivzz)*(bivxx+bivyy+bivzz);
							gradlam=-gradlam*fdummy;

							gradmu=0.0;
							gradmu= 2*fvxx*bvxx+2*fvyy*bvyy+2*fvzz*bvzz+2*fivxx*bivxx+2*fivyy*bivyy+2*fivzz*bivzz+(fvxy+fvyx)*(bvxy+bvyx)+(fivxy+fivyx)*(bivxy+bivyx)+(fvxz+fvzx)*(bvxz+bvzx)+(fivxz+fivzx)*(bivxz+bivzx)+(fvyz+fvzy)*(bvyz+bvzy)+(fivyz+fivzy)*(bivyz+bivzy);
							gradmu=-gradmu*fdummy;

							gradrho=0.0;
							gradrho=(fvel_l_c[pos_l * 6] * bvel_l_c[pos_l * 6] + fvel_l_c[pos_l * 6 + 1] * bvel_l_c[pos_l * 6 + 1] + fvel_l_c[pos_l * 6 + 2] * bvel_l_c[pos_l * 6 + 2]) + (fvel_l_c[pos_l * 6 + 3] * bvel_l_c[pos_l * 6 + 3] + fvel_l_c[pos_l * 6 + 4] * bvel_l_c[pos_l * 6 + 4] + fvel_l_c[pos_l * 6 + 5] * bvel_l_c[pos_l * 6 + 5]);

							//gradrho=gradrho;

							//parametrisation vp, vs, rho

							grad_l[pos_grad * 3 + 0] += sqrt(rho_l[pos_pi]*pi_l[pos_pi])*2*gradlam; //gradient vp

							grad_l[pos_grad * 3 + 1] += -4*sqrt(rho_l[pos_pi]*u_l[pos_pi])*gradlam+2*sqrt(rho_l[pos_pi]*u_l[pos_pi])*gradmu; //gradient vs

							grad_l[pos_grad * 3 + 2] += gradrho+u_l[pos_pi]/rho_l[pos_pi]*gradmu+(pi_l[pos_pi]-2*u_l[pos_pi])/rho_l[pos_pi]*gradlam; //gradient rho

						}  // iz
					}  // ix

					//get next plane fv, fiv, bv, biv
	                fvel_l_tmp = fvel_l + plane_comp * ystep_l6;
	                bvel_l_tmp = bvel_l + plane_comp * ystep_l6;
					fvel_s_tmp = fvel_s1 + (iy + wy + y1) * ystep_6 - x0 * xstep_6 - z0 * 6;
					bvel_s_tmp = bvel_s1 + (iy + wy + y1) * ystep_6 - x0 * xstep_6 - z0 * 6;
					for (ix = -x0; ix < ixn + x1; ix++) {
						get_reply = 0;
						athread_get(PE_MODE, fvel_s_tmp, fvel_l_tmp, (wz + z0 + z1) * 6 * sizeof(float), &get_reply, 0, 0, 0);
						athread_get(PE_MODE, bvel_s_tmp, bvel_l_tmp, (wz + z0 + z1) * 6 * sizeof(float), &get_reply, 0, 0, 0);
						while(get_reply!=2);
						bvel_s_tmp += xstep_6;
						bvel_l_tmp += xstep_l6;
						fvel_s_tmp += xstep_6;
						fvel_l_tmp += xstep_l6;

					}

					plane_comp = addn(plane_comp, 1, wy + y0 + y1);

					//put the gradient data to main memory
					grad_l_tmp = grad_l;
					grad_s1_tmp = grad_s1 + iy * ystep_grad_3;

					for (ix = 0; ix < ixn; ix++) {

						put_reply = 0;
						athread_put(PE_MODE, grad_l_tmp, grad_s1_tmp, wz * 3 * sizeof(float), &put_reply, 0, 0);
						while(put_reply!=1);
						grad_s1_tmp += xstep_grad_3;
						grad_l_tmp += wz * 3;

					}

				} // end of loop iy

			} // end of loop iix

		} // end of loop iiz


	} // end of loop l



}















