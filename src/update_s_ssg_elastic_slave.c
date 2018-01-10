#include <slave.h>
#include <dma.h>


#define wx 4
#define wy 1
#define wz 36
#define MX 64
//#define DEBUG

typedef struct Param_str {
	int dim_x;
	int dim_y;
	int dim_z;
	int slice;
	int strip;
	int slice_up;
	int strip_up;
	int slice_pi;
	int strip_pi;
	float DT;
	float DX;
	float DY;
	float DZ;
	int FDCOEFF;
	float *vel;
	float *stress;
	float *up;
	float *pi;
	float *u;
} Param_str;


inline  addn(int x, int n, int x1) {
	return x + n < x1 ? x + n : x + n - x1;
}

inline void dma_get_set(dma_desc *p_dma_get, volatile int *p_replyget, int zlen) {
	dma_set_op(p_dma_get, DMA_GET);
	dma_set_mode(p_dma_get, PE_MODE);
	dma_set_reply(p_dma_get, p_replyget);

	dma_set_size(p_dma_get, zlen * sizeof(float));
	dma_set_bsize(p_dma_get, zlen * sizeof(float));
	dma_set_stepsize(p_dma_get, 0);
}

inline void dma_put_set(dma_desc *p_dma_put, volatile int *p_replyput, int zlen) {
	dma_set_op(p_dma_put, DMA_PUT);
	dma_set_mode(p_dma_put, PE_MODE);
	dma_set_reply(p_dma_put, p_replyput);

	dma_set_size(p_dma_put, zlen * sizeof(float));
	dma_set_bsize(p_dma_put, zlen * sizeof(float));
	dma_set_stepsize(p_dma_put, 0);
}


void update_s_kernel_fusion_slave(Param_str *param) {

	int id = athread_get_id(-1);

	int dim_x = param->dim_x;
	int dim_y = param->dim_y;
	int dim_z = param->dim_z;
	int ystep = param->slice;
	int xstep = param->strip;
	int ystep_up = param->slice_up;
	int xstep_up = param->strip_up;
	int ystep_pi = param->slice_pi;
	int xstep_pi = param->strip_pi;
	float DT = param->DT;
	float DX = param->DX;
	float DY = param->DY;
	float DZ = param->DZ;
	int FDCOEFF = param->FDCOEFF;

	float *vel_s0 = param->vel;
	float *strs_s0 = param->stress;
	float *up_s0 = param->up;
	float *pi_s0 = param->pi;
	float *u_s0 = param->u;


	//stencil operator
	const int x0 = 2;
	const int x1 = 2;
	const int y0 = 2;
	const int y1 = 2;
	const int z0 = 2;
	const int z1 = 2;

	int ix, iy, iz, iix, iiz;
	int txbeg, txend, txn;
	int izbeg, izend, izn;



	int xstep_6 = xstep * 6;
	int ystep_6 = ystep * 6;
	int xstep_3 = xstep * 3;
	int ystep_3 = ystep * 3;
	int xstep_l = wz + z0 + z1;
	int ystep_l = (wx + x0 + x1) * (wz + z0 + z1);
	int xstep_l2 = xstep_l * 2;
	int ystep_l3 = ystep_l * 3;
	int xstep_l3 = xstep_l * 3;
	int ystep_l6 = ystep_l * 6;
	int xstep_l6 = xstep_l * 6;



	int NX = (dim_x + wx * MX - 1) / (wx * MX);
	int NZ = (dim_z + wz - 1) / wz;

	float vel[(wy + y0 + y1) * (wx + x0 + x1) * (wz + z0 + z1) * 3];
	float strs[wx * wz * 6];
	float up[wx * wz * 3];
	float pi[wx * wz];
	float u[wx * wz ];

	float *vel_s1, *strs_s1, *up_s1, *pi_s1, *u_s1;
    float *vel_temp, *vel_s1_temp, *strs_temp, *strs_s1_temp,*up_temp,*up_s1_temp,*pi_temp,*pi_s1_temp,*u_temp,*u_s1_temp;
	float *vel_c = vel + x0 * xstep_l3 + z0 * 3;

	float vxx, vxy, vxz, vyx, vyy, vyz, vzx, vzy, vzz;
	float vxyyx, vyzzy, vxzzx, vxxyyzz, vyyzz, vxxzz, vxxyy;
	float g, f, fipjp, fjpkp, fipkp;
	int plane_comp = 0;
	float b1, b2;

	dma_desc  dma_get_wz, dma_get_3wz, dma_get_3wz_z0, dma_put_6wz, dma_get_6wz;
	volatile int get_wz = 0, get_3wz = 0, get_3wz_z0 = 0, put_6wz = 0, get_6wz = 0;
	dma_get_set(&dma_get_wz, &get_wz, wz);
	dma_get_set(&dma_get_3wz, &get_3wz, wz * 3);
	dma_get_set(&dma_get_6wz, &get_6wz, wz * 6);
	dma_put_set(&dma_put_6wz, &put_6wz, wz * 6);
	dma_get_set(&dma_get_3wz_z0, &get_3wz_z0, (wz + z0 + z1) * 3);


#ifdef DEBUG
	int ldm_vel = sizeof(float) * (wx + x0 + x1) * (wy + y0 + y1) * (wz + z0 + z1) * 3;
	int ldm_strs = sizeof(float) * wz * wz * 6;
	int ldm_up = sizeof(float) * wx * wz * 3;
	int ldm_pi = sizeof(float) * wx * wz * 3;
	int ldm_u = sizeof(float) * wx * wz * 3;
	int ldm_total = ldm_vel + ldm_strs + ldm_up + ldm_pi + ldm_u;
	if(id == 0) {
		//printf("LDM vel = %d Bytes\n", ldm_vel);
		//printf("LDM stress = %d Bytes\n", ldm_strs);
		//printf("LDM up = %d Bytes\n", ldm_up);
		//printf("LDM pi = %d Bytes\n", ldm_pi);
		//printf("LDM u = %d Bytes\n", ldm_u);
		printf("LDM total of update_s = %d Bytes\n", ldm_total);
  }
#endif




	b1 = 9.0 / 8.0; b2 = -1.0 / 24.0; /* Taylor coefficients*/
	if (FDCOEFF == 2) {
		b1 = 1.1382; b2 = -0.046414;
	} /* Holberg coefficients E=0.1 %*/


	for (iiz = 0 ; iiz < NZ ; iiz ++) {

		izbeg = wz * iiz;
		izend = wz * (iiz + 1);
		izend = izend < dim_z ? izend : dim_z;
		izn = izend - izbeg;
		dma_set_bsize(&dma_put_6wz, 6 *  izn * sizeof(float));

		for (iix = 0 ; iix < NX ; iix ++) {

			txbeg = wx * MX * iix + wx * id;
			txend = wx * MX * iix + wx * (id + 1);
			txend = txend < dim_x ? txend : dim_x;
			txn = txend - txbeg;


			vel_s1 = vel_s0 + (txbeg * xstep + izbeg) * 3;
			strs_s1 = strs_s0 + (txbeg * xstep + izbeg) * 6;
			up_s1 = up_s0 + (txbeg * xstep_up + izbeg) * 3;
			pi_s1 = pi_s0 + txbeg * xstep_pi + izbeg ;
			u_s1 = u_s0 + txbeg * xstep_pi + izbeg;


			vel_temp = vel;
			for (iy = -y0 ; iy < wy + y1 ; iy ++) {

				vel_s1_temp = vel_s1 + iy * ystep_3 - x0 * xstep_3 - z0 * 3;

				for (ix = -x0 ; ix < wx + x1 ; ix ++) {

					dma(dma_get_3wz_z0, (long)(vel_s1_temp), (long)(vel_temp));
					dma_wait(&get_3wz_z0, 1);
					get_3wz_z0 = 0;

					/* next x line*/
					vel_s1_temp += xstep_3;
					vel_temp += xstep_l3;
				}
			}

			plane_comp = 0;
			/*compute velocity from iy plane to iy++ plane*/
			for (iy = 0 ; iy < dim_y ; iy ++) {

				/*
					1. get velocity data from vel[iy][0] to vel[iy][wx]
					2. every geted data size is wz
				*/

				strs_temp = strs;
				strs_s1_temp = strs_s1 + iy * ystep * 6;
				for (ix = 0 ; ix < txn ; ix ++) {

					dma(dma_get_6wz, (long)(strs_s1_temp), (long)(strs_temp));
					dma_wait(&get_6wz, 1);
					get_6wz = 0;

					strs_s1_temp += xstep * 6;
					strs_temp += wz * 6;

				}

				up_temp = up;
				up_s1_temp = up_s1 + iy * ystep_up * 3;
				for (ix = 0 ; ix < txn ; ix ++) {

					dma(dma_get_3wz, (long)(up_s1_temp), (long)(up_temp));
					dma_wait(&get_3wz, 1);
					get_3wz = 0;

					up_s1_temp += xstep_up * 3;
					up_temp += wz * 3;

				}

				pi_temp = pi;
				pi_s1_temp = pi_s1 + iy * ystep_pi;
				for (ix = 0 ; ix < txn ; ix ++) {
					dma(dma_get_wz, (long)(pi_s1_temp), (long)(pi_temp));
					dma_wait(&get_wz, 1);
					get_wz = 0;
					pi_s1_temp += xstep_pi;
					pi_temp += wz;
				}


				u_temp = u;
				u_s1_temp = u_s1 + iy * ystep_pi;
				for (ix = 0 ; ix < txn ; ix ++) {

					dma(dma_get_wz, (long)(u_s1_temp), (long)(u_temp));
					dma_wait(&get_wz, 1);
					get_wz = 0;

					u_s1_temp += xstep_pi;
					u_temp += wz;

				}


				/*compute velocity at iy plane*/
				for (ix = 0 ; ix < txn ; ix ++) {
					for (iz = 0 ; iz < izn ; iz ++) {

						int ym2 = plane_comp;
						int ym1 = addn(plane_comp, 1, wy + y0 + y1);
						int y = addn(plane_comp, y0 , wy + y0 + y1);
						int yp1 = addn(plane_comp, y0 + 1, wy + y0 + y1);
						int yp2 = addn(plane_comp, y0 + 2, wy + y0 + y1);

						/*computing positon*/
						int pos = y * ystep_l + ix * xstep_l + iz;

						/* z direction */
						int pos_km2_l = pos - 2;
						int pos_km1_l = pos - 1;
						int pos_kp1_l = pos + 1;
						int pos_kp2_l = pos + 2;

						/* x direction */
						int pos_im2_l = pos - xstep_l2;
						int pos_im1_l = pos - xstep_l;
						int pos_ip1_l = pos + xstep_l;
						int pos_ip2_l = pos + xstep_l2;

						/* y direction */
						int pos_jm2_l = ym2 * ystep_l + ix * xstep_l + iz;
						int pos_jm1_l = ym1 * ystep_l + ix * xstep_l + iz;
						int pos_jp1_l = yp1 * ystep_l + ix * xstep_l + iz;
						int pos_jp2_l = yp2 * ystep_l + ix * xstep_l + iz;

						int pos_s = ix * wz + iz ;

						vxx = (b1 * (vel_c[pos * 3 + 0] - vel_c[pos_im1_l * 3 + 0]) + b2 * (vel_c[pos_ip1_l * 3 + 0] - vel_c[pos_im2_l * 3 + 0])) / DX;
						vxy = (b1 * (vel_c[pos_jp1_l * 3 + 0] - vel_c[pos * 3 + 0]) + b2 * (vel_c[pos_jp2_l * 3 + 0] - vel_c[pos_jm1_l * 3 + 0])) / DY;
						vxz = (b1 * (vel_c[pos_kp1_l * 3 + 0] - vel_c[pos * 3 + 0]) + b2 * (vel_c[pos_kp2_l * 3 + 0] - vel_c[pos_km1_l * 3 + 0])) / DZ;

						vyx = (b1 * (vel_c[pos_ip1_l * 3 + 1] - vel_c[pos * 3 + 1]) + b2 * (vel_c[pos_ip2_l * 3 + 1] - vel_c[pos_im1_l * 3 + 1])) / DX;
						vyy = (b1 * (vel_c[pos * 3 + 1] - vel_c[pos_jm1_l * 3 + 1]) + b2 * (vel_c[pos_jp1_l * 3 + 1] - vel_c[pos_jm2_l * 3 + 1])) / DY;
						vyz = (b1 * (vel_c[pos_kp1_l * 3 + 1] - vel_c[pos * 3 + 1]) + b2 * (vel_c[pos_kp2_l * 3 + 1] - vel_c[pos_km1_l * 3 + 1])) / DZ;

						vzx = (b1 * (vel_c[pos_ip1_l * 3 + 2] - vel_c[pos * 3 + 2]) + b2 * (vel_c[pos_ip2_l * 3 + 2] - vel_c[pos_im1_l * 3 + 2])) / DX;
						vzy = (b1 * (vel_c[pos_jp1_l * 3 + 2] - vel_c[pos * 3 + 2]) + b2 * (vel_c[pos_jp2_l * 3 + 2] - vel_c[pos_jm1_l * 3 + 2])) / DY;
						vzz = (b1 * (vel_c[pos * 3 + 2] - vel_c[pos_km1_l * 3 + 2]) + b2 * (vel_c[pos_kp1_l * 3 + 2] - vel_c[pos_km2_l * 3 + 2])) / DZ;



						// updating components of the stress tensor, partially
						fipjp = up[pos_s * 3 + 0] * DT;
						fjpkp = up[pos_s * 3 + 1] * DT;
						fipkp = up[pos_s * 3 + 2] * DT;
						g = pi[pos_s];
						f = 2.0 * u[pos_s];

                        vxyyx = vxy + vyx;
                        vyzzy = vyz + vzy;
                        vxzzx = vxz + vzx;
                        vxxyyzz = vxx + vyy + vzz;
                        vyyzz = vyy + vzz;
                        vxxzz = vxx + vzz;
                        vxxyy = vxx + vyy;



						strs[pos_s * 6 + 0] += DT * ((g * vxxyyzz) - (f * vyyzz));
						//strs[pos_s * 6 + 0] += vxz;
						strs[pos_s * 6 + 1] += DT * ((g * vxxyyzz) - (f * vxxzz));
						//strs[pos_s * 6 + 1] += vyx+vyy+vyz;
						strs[pos_s * 6 + 2] += DT * ((g * vxxyyzz) - (f * vxxyy));
					//	strs[pos_s * 6 + 2] += vzz;
						strs[pos_s * 6 + 3] += (fipjp * vxyyx);
						strs[pos_s * 6 + 4] += (fjpkp * vyzzy);
						strs[pos_s * 6 + 5] += (fipkp * vxzzx);


					}/*iz*/
				}/*ix*/

				/*get next plane vel data*/
			    vel_temp = vel + plane_comp * ystep_l3;
				vel_s1_temp = vel_s1 + (iy + wy + y1) * ystep_3 - x0 * xstep_3 - z0 * 3;
				for (ix = -x0 ; ix < txn + x1 ; ix ++) {

					dma(dma_get_3wz_z0, (long)(vel_s1_temp), (long)(vel_temp));
					dma_wait(&get_3wz_z0, 1);
					get_3wz_z0 = 0;

					vel_s1_temp += xstep_3;
					vel_temp += xstep_l3;

				}


				plane_comp = addn(plane_comp, 1, wy + y0 + y1);

				/*put stress data*/
				strs_temp = strs;
				strs_s1_temp = strs_s1 + iy * ystep * 6;
				for (ix = 0 ; ix < txn ; ix ++) {

					dma(dma_put_6wz, (long)(strs_s1_temp), (long)(strs_temp));
					dma_wait(&put_6wz, 1);
					put_6wz = 0;

					strs_s1_temp += xstep * 6;
					strs_temp += wz * 6;

				}

			}/*iy*/

		}/*iix*/

	}/*iiz*/

}
