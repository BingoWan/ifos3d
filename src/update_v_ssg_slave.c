#include <stdio.h>
#include <slave.h>
#include <dma.h>
#define wx 16
#define wy 1
#define wz 16
#define MX 64
//#define DEBUG

typedef struct Param_vel {
	int dim_x;
	int dim_y;
	int dim_z;
	int slice;
	int strip;
	int slice_rp;
	int strip_rp;
	float  dx;
	float  dy;
	float dz;
	int FDCOEFF;
	float *vel;
	float *stress;
	float * rp;
} Param_vel;


inline  addn(int x, int n, int x1) {
	return x + n < x1 ? x + n : x + n - x1;
}


inline int id_map(int id) {
	int i = id / 8;
	if (i % 2 == 1) {
		int beg = i * 8;
		int end = (i + 1) * 8 - 1;
		return beg + end - id;
	}
	return id;
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



void update_v_kernel_fusion_slave(Param_vel *param) {

	int id = athread_get_id(-1);
	int dim_x = param->dim_x;
	int dim_y = param->dim_y;
	int dim_z = param->dim_z;
	int xstep = param->strip;
	int ystep = param->slice;
	int xstep_rp = param->strip_rp;
	int ystep_rp = param->slice_rp;
	float dx = param->dx;
	float dy = param->dy;
	float dz = param->dz;
	int FDCOEFF = param->FDCOEFF;

	float *vel_s0 = param->vel;
	float *str_s0 = param->stress;
	float *rp_s0 = param->rp;

	const int x0 = 2;
	const int x1 = 2;
	const int y0 = 2;
	const int y1 = 2;
	const int z0 = 2;
	const int z1 = 2;

	int ix, iy, iz, iix, iiz;
	int izbeg, izend, izn;
	int ixbeg, ixend, ixn;

	int xstep_6 = xstep * 6;
	int ystep_6 = ystep * 6;

	int xstep_l = wz + z0 + z1;
	int ystep_l = (wx + x0 + x1) * (wz + z0 + z1);
	int xstep_l2 = xstep_l * 2;
	int ystep_l6 = ystep_l * 6;
	int xstep_l6 = xstep_l * 6;



	int NX = (dim_x + wx * MX - 1) / (wx * MX);
	int NZ = (dim_z + wz - 1) / wz;

	float strs[(wy + y0 + y1) * (wx + x0 + x1) * (wz + z0 + z1) * 6];
	float vel[wx * wz * 3];
	float rp[wx * wz * 3];
	float *vel_s1, * str_s1, *rp_s1;
	float *str_temp, *str_s1_temp;
	float sxx_x, sxy_y, sxz_z, syy_y, sxy_x, syz_z, szz_z, sxz_x, syz_y;
	float *strs_c = strs + x0 * xstep_l6 + z0 * 6;



	dma_desc  dma_get_3wz, dma_put_3wz, dma_get_6wz, dma_get_6wz_1;
	volatile int get_3wz = 0, put_3wz = 0, get_6wz = 0, get_6wz_1 = 0;
	dma_get_set(&dma_get_3wz, &get_3wz, wz * 3);
	dma_put_set(&dma_put_3wz, &put_3wz, wz * 3);
	dma_get_set(&dma_get_6wz, &get_6wz, (wz + z0 + z1) * 6);
	dma_get_set(&dma_get_6wz_1, &get_6wz_1, (wz + z0 + z1) * 6);
	int plane_comp = 0;
	float b1, b2;
	b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients*/
	if(FDCOEFF==2){
	b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/


#ifdef DEBUG
	int ldm_vel = sizeof(float) * wx * wz * 3;
	int ldm_strs = sizeof(float) * (wx + x0 + x1) * (wy + y0 + y1) * (wz + z0 + z1)  * 6;
	int ldm_rp = sizeof(float) * wx * wz * 3;

	int ldm_total = ldm_vel + ldm_strs + ldm_rp;
	if(id == 0) {
		//printf("LDM vel = %d Bytes\n", ldm_vel);
		//printf("LDM stress = %d Bytes\n", ldm_strs);
		//printf("LDM up = %d Bytes\n", ldm_rp);
		printf("LDM total of update_v = %d Bytes\n", ldm_total);
	}
#endif


	for (iiz = 0; iiz < NZ; iiz++) {
		izbeg = wz * iiz;
		izend = wz * (iiz + 1);
		izend = izend < dim_z ? izend : dim_z;
		izn = izend - izbeg;

		for (iix = 0; iix < NX; iix++) {

			ixbeg = wx * MX * iix + wx * id;
			ixend = wx * MX * iix + wx * (id + 1);
			ixend = ixend < dim_x ? ixend : dim_x;
			ixn = ixend - ixbeg;

			vel_s1 = vel_s0 + (ixbeg * xstep + izbeg) * 3;
			str_s1 = str_s0 + (ixbeg * xstep + izbeg) * 6;
			rp_s1 = rp_s0 + (ixbeg * xstep_rp + izbeg) * 3;

			//1. get stress data from stress[-y0][-x0][-wz0] to stress[wy+y1][wx+x1][wz+x1]
			//2. Data size obtained each time is wz * 3
			str_temp = strs;
			for (iy = -y0; iy < wy + y1; iy++) {

				str_s1_temp = str_s1 + iy * ystep_6 - x0 * xstep_6 - z0 * 6;
				for (ix = -x0; ix < wx + x1; ix++) {

					dma(dma_get_6wz, (long)(str_s1_temp), (long)(str_temp));
					dma_wait(&get_6wz, 1);
					get_6wz = 0;

					// next x line
					str_s1_temp += xstep_6;
					str_temp += xstep_l6;

				}
			}


			plane_comp = 0;
			// compute velocity from iy plane to iy++ plane
			for (iy = 0; iy < dim_y; iy++) {

				//1. get velocity data from vel[iy][0] to vel[iy][wx]
				//2. Data size obtained each time is wz * 3
				float *vel_temp = vel;
				float *vel_s1_temp = vel_s1 + iy * ystep * 3;
				for (ix = 0; ix < ixn; ix++) {

					dma(dma_get_3wz, (long)(vel_s1_temp), (long)(vel_temp));
					dma_wait(&get_3wz, 1);
					get_3wz = 0;
					vel_s1_temp += xstep * 3;
					vel_temp += wz * 3;

				}

				float *rp_temp = rp;
				float *rp_s1_temp = rp_s1 + iy * ystep_rp * 3;
				for (ix = 0 ; ix < ixn ; ix ++) {

					dma(dma_get_3wz, (long)(rp_s1_temp), (long)(rp_temp));
					dma_wait(&get_3wz, 1);
					get_3wz = 0;

					rp_s1_temp += xstep_rp * 3;
					rp_temp += wz * 3;

				}

				for (ix = 0; ix < ixn; ix++) {

					for (iz = 0; iz < izn; iz++) {
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

						int pos_v = ix * wz + iz;

						sxx_x = dx * (b1 * (strs_c[pos_ip1_l * 6 + 0] - strs_c[pos * 6 + 0]) + b2 * (strs_c[pos_ip2_l * 6 + 0] - strs_c[pos_im1_l * 6 + 0]));
						sxy_y = dy * (b1 * (strs_c[pos * 6 + 3] - strs_c[pos_jm1_l * 6 + 3]) + b2 * (strs_c[pos_jp1_l * 6 + 3] - strs_c[pos_jm2_l * 6 + 3]));
						sxz_z = dz * (b1 * (strs_c[pos * 6 + 5] - strs_c[pos_km1_l * 6 + 5]) + b2 * (strs_c[pos_kp1_l * 6 + 5] - strs_c[pos_km2_l * 6 + 5]));

						vel[pos_v * 3 + 0] += (sxx_x + sxy_y + sxz_z) / rp[pos_v * 3 + 0];

						syy_y = dy * (b1 * (strs_c[pos_jp1_l * 6 + 1] - strs_c[pos * 6 + 1]) + b2 * (strs_c[pos_jp2_l * 6 + 1] - strs_c[pos_jm1_l * 6 + 1]));
						sxy_x = dx * (b1 * (strs_c[pos * 6 + 3] - strs_c[pos_im1_l * 6 + 3]) + b2 * (strs_c[pos_ip1_l * 6 + 3] - strs_c[pos_im2_l * 6 + 3]));
						syz_z = dz * (b1 * (strs_c[pos * 6 + 4] - strs_c[pos_km1_l * 6 + 4]) + b2 * (strs_c[pos_kp1_l * 6 + 4] - strs_c[pos_km2_l * 6 + 4]));

						vel[pos_v * 3 + 1] += (syy_y + sxy_x + syz_z) / rp[pos_v * 3 + 1];

						szz_z = dz * (b1 * (strs_c[pos_kp1_l * 6 + 2] - strs_c[pos * 6 + 2]) + b2 * (strs_c[pos_kp2_l * 6 + 2] - strs_c[pos_km1_l * 6 + 2]));
						sxz_x = dx * (b1 * (strs_c[pos * 6 + 5] - strs_c[pos_im1_l * 6 + 5]) + b2 * (strs_c[pos_ip1_l * 6 + 5] - strs_c[pos_im2_l * 6 + 5]));
						syz_y = dy * (b1 * (strs_c[pos * 6 + 4] - strs_c[pos_jm1_l * 6 + 4]) + b2 * (strs_c[pos_jp1_l * 6 + 4] - strs_c[pos_jm2_l * 6 + 4]));

						vel[pos_v * 3 + 2] += (szz_z + sxz_x + syz_y) / rp[pos_v * 3 + 2];
						//printf("------>\n");
						//printf("---###----->rp[pos_v * 3 + 0]: %f, rp[pos_v * 3 + 1]: %f, rp[pos_v * 3 + 2]: %f\n",rp[pos_v * 3 + 0], rp[pos_v * 3 + 1], rp[pos_v * 3 + 2]);


					}  // iz
				}  //ix

				//get next plane stress data
				str_temp = strs + plane_comp * ystep_l6;
				str_s1_temp = str_s1 + (iy + wy + y1) * ystep_6 - x0 * xstep_6 - z0 * 6;
				for (ix = -x0 ; ix < ixn + x1 ; ix ++) {

					dma(dma_get_6wz_1, (long)(str_s1_temp), (long)(str_temp));
					dma_wait(&get_6wz_1, 1);
					get_6wz_1 = 0;

					str_s1_temp += xstep_6;
					str_temp += xstep_l6;

				}


				plane_comp = addn(plane_comp, 1, wy + y0 + y1);

				//put vel data
				vel_temp = vel;
				vel_s1_temp = vel_s1 + iy * ystep * 3;
				for (ix = 0 ; ix < ixn ; ix ++) {

					dma(dma_put_3wz, (long)(vel_s1_temp), (long)(vel_temp));
					dma_wait(&put_3wz, 1);
					put_3wz = 0;

					vel_s1_temp += xstep * 3;
					vel_temp += wz * 3;

				}
			}
		}
	}





}

