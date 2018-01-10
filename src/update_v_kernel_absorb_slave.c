/*
 * update_v_kernel_fusion_slave_2.c
 *
 *  Created on: Jan 3, 2018
 *      Author: bingo
 */
#include <slave.h>
#include <dma.h>
#define wy 1
#define wx 4
#define wz 32
#define MX 64
//#define DEBUG



typedef struct {
	float * vel;
	float * stress;
	float * absorb;
	int slice;
	int strip;
	int slice_a;
	int strip_a;
	int dimx;
	int dimy;
	int dimz;
} Param_absorb;


void update_v_kernel_absorb_salve(Param_absorb * param) {

	int id = athread_get_id(-1);
	int dimx = param->dimx;
	int dimy = param->dimy;
	int dimz = param->dimz;
	int slice = param->slice;
	int strip = param->strip;
	int slice_a = param->slice_a;
	int strip_a = param->strip_a;
	float * vel = param->vel;
	float * stress = param->stress;
	float * absorb = param->absorb;

	float vel_l[wy * wx * wz * 3];
	float stress_l[wy * wx * wz * 6];
	float absorb_l[wy * wx * wz];


	int NX = (dimx + MX * wx - 1) / (MX * wx);
	int NY = (dimy + wy - 1) / wy;
	int NZ = (dimz + wz - 1) / wz;


	int ix, iy, iz, i, j, k, ixbeg, ixend, iybeg, iyend, izbeg, izend, ixn, iyn, izn, size;

	float * vel_s, * stress_s, * absorb_s, * vel_l_s, * stress_l_s, * absorb_l_s;
	/*
	stress_l_s = stress_l;
	vel_l_s = vel_l;
	absorb_l_s = absorb_l;*/


	volatile unsigned long get_reply, put_reply;

#ifdef DEBUG
	int ldm_vel = sizeof(float) * wy  * wx * wz * 3;
	int ldm_strs = sizeof(float) * wy  * wx * wz * 6;
	int ldm_absorb = sizeof(float) * wy * wx * wz;

	int ldm_total = ldm_vel + ldm_strs + ldm_absorb;
	if(id == 0) {
		//printf("LDM vel = %d Bytes\n", ldm_vel);
		//printf("LDM stress = %d Bytes\n", ldm_strs);
		//printf("LDM absorb = %d Bytes\n", ldm_absorb);
		printf("LDM total of update_v = %d Bytes\n", ldm_total);
	}
#endif

	// To obtain the date from main memory
	for (iy = 0; iy < NY; iy ++) {
		iybeg = wy * iy;
		iyend = wy * (iy + 1);
		iyend = iyend < dimy ? iyend : dimy;
		iyn = iyend - iybeg;

		for (ix = 0; ix < NX; ix ++) {
			ixbeg = wx * MX * ix + wx * id;
			ixend = wx * MX * ix + wx * (id + 1);
			ixend = ixend < dimx ? ixend : dimx;
			ixn = ixend - ixbeg;
			if (ixend < ixbeg)
				break;

			for (iz = 0; iz < NZ; iz ++) {
				izbeg = wz * iz;
				izend = wz * (iz + 1);
				izend = izend < dimz ? izend : dimz;
				izn = izend - izbeg;

				vel_s = vel + (iybeg * slice + ixbeg * strip + izbeg) * 3;
				stress_s = stress + (iybeg * slice + ixbeg * strip+ izbeg) * 6;
				absorb_s = absorb + iybeg * slice_a + ixbeg * strip_a + izbeg;

				size = (iyn - 1) * wx * wz + (ixn - 1) * wz + izn;

				get_reply = 0;
				athread_get(PE_MODE, stress_s, stress_l, size * 6 * sizeof(float), &get_reply, 0, 0, 0);
				athread_get(PE_MODE, vel_s, vel_l, size * 3 * sizeof(float), &get_reply, 0, 0, 0);
				athread_get(PE_MODE, absorb_s, absorb_l, size * sizeof(float), &get_reply, 0, 0, 0);
				while(get_reply != 3);

				for (j = 0; j < iyn; j ++) {
					for (i = 0; i < ixn; i ++) {
						for (k = 0; k < izn; k ++) {

							int pos = j * wx * wz + i * wz + k;

							vel_l[pos * 3 + 0] *= absorb_l[pos];
							vel_l[pos * 3 + 1] *= absorb_l[pos];
							vel_l[pos * 3 + 2] *= absorb_l[pos];
							stress_l[pos * 6 + 0] *= absorb_l[pos];
							stress_l[pos * 6 + 1] *= absorb_l[pos];
							stress_l[pos * 6 + 2] *= absorb_l[pos];
							stress_l[pos * 6 + 3] *= absorb_l[pos];
							stress_l[pos * 6 + 4] *= absorb_l[pos];
							stress_l[pos * 6 + 5] *= absorb_l[pos];


						}
					}
				}

				put_reply = 0;
				athread_put(PE_MODE,stress_l, stress_s, size * 6 * sizeof(float), &put_reply, 0, 0);
				athread_put(PE_MODE, vel_l, vel_s, size * 3 * sizeof(float), &put_reply, 0, 0);
				while(put_reply != 2);

			}
		}
	}

}



