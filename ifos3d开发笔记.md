### ifos3d开发笔记
#### v1版本
修改了vx, vy, vz, sxy, syz在POS[2]>0的条件的数组长度
例如：
原来是vx  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
修改后是 vx  =  f3tensor(0-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);

#### v2版本
修改了Ffvx, Ffvy, Ffvz在POS[2]>0的条件的数组长度，
例如：
原来是Ffvx  =  f4tensor(1,NFMAX,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
修改后是Ffvx  =  f4tensor(1,NFMAX,0-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
修改数组的长度的目的主要是为了方便数组融合，若在不同的进程的数组长度不一致，数组融合会有问题，这里为了方便直接POS[2]>0的条件的数组扩充和POS[2]==0条件下的长度。
三维数组Ffvx, Ffvy, Ffvz涉及的函数有zero_invers、discfourier、exchange_Fv、gradient_F、hess_F。

修改Ffvx长度后发现进程号 2,3,6,7 segmenttation Fault的错误，进程号 0,1,4,5没有问题，既POS[2]==0条件下进程没有问题。那就判断数组扩展问题了
后来索性将所有的数组长度扩充为f4tensor(1,NFMAX,0-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);


#### v3版本

##### exchange_fv.c 的改进版，exchange_Ffv.c
考虑到整个ifos3d.c中 Ffvx,Ffvy,Ffvz,Ffivx,Ffivy,Ffivz,在函数间(zero_invers,discfourier, gradient_F, hess_F)的传递都是一起的，但是在exchange_Fv函数中Ffvx,Ffvy,Ffvz和Ffivx,Ffivy,Ffivz是分在两个函数中。
所以有意将上述的两个函数合为一体，就开发了exchange_Ffv函数。
1.将 Ffvx,Ffvy,Ffvz,Ffivx,Ffivy,Ffivz六个四维素组融合成一维数组F_Fv；
2.将原来的三个四维数组迭代，变成了六个四维数组迭代；


开发exchange_Ffv函数，增加了fbufferlef_to_rig, fbufferrig_to_lef, fbuffertop_to_bot, fbufferbot_to_top, fbufferfro_to_bac, fbufferbac_to_fro 六个变量。
exchange_Ffv函数也可以合并下面的两个函数
exchange_Fv(Fbvx,Fbvy,Fbvz,nf, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top,
												bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec,ntr_hess);
												
exchange_Fv(Fbivx,Fbivy,Fbivz,nf, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top,
												bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec,ntr_hess);


1.将 Ffvx,Ffvy,Ffvz,Ffivx,Ffivy,Ffivz六个四维素组融合成一维数组F_Fv；
2.将原来的三个四维数组迭代，变成了六个四维数组迭代；

#### v4版本			
将disc_fourier.c、zero_invers.c, gradient_F，hess_F中的Ffvx,Ffvy,Ffvz,Ffivx,Ffivy,Ffivz六个数组进行了融合	；


#### v5版本
将surface_elastic(1,u,pi,sxx,syy,szz,sxy,vx,vy,vz); 改成 surface_elastic(1,u,pi,sxx,syy,szz,sxy,syz,vx,vy,vz)，
即多传递了个sxz参数，以方便融合；

将void psource(int nt, float *** sxx, float *** syy, float *** szz,  float **  srcpos_loc, float ** signals, int nsrc)
改成
void psource(int nt, float *** sxx, float *** syy, float *** szz, float *** sxy, float *** syz, float *** sxz, float **  srcpos_loc, float ** signals, int nsrc)








sxx, syy, szz, sxy, syz, sxz需要做融合函数有
update_v_CPML  ****
update_v    ----
update_s_elastic     ---------
update_s_CPML_elastic
exchange_s     ---------
exchange_v     ---------
surface_elastic  -------
snap    -----------
zero_wavefield   --------
psource       ---------
seismo        ------------*******

vx,vy,vz函数
seismo     ----------******
zero_wavefield   --------
update_v         --------
update_s_elastic    -------
discfourier         -------
surface_elastic     -------
snap     -----------******
hess_F   -----------******


融合后的由变量而需要改变的函数
hess1,hess2,hess3: readhess, outgrad, hess_apply, hess_F;

pi: readmod, model



sxx, vx... 已经修改了
exchange_v
zero_wavefield
update_v_ssg
update_s_ssg_elastic
psource
exchange_s
surface_ssg_elastic
snap_ssg
disc_fourier
fd.h
ifos3d
共11个函数。


Ffvx 已经修改了
zero_invers
discfourier
exchange_Ffv
gradient_F
hess_F

bufferlef_to_rig 已经修改了
exchange_v
exchange_s
exchange_Ffv

uipjp.... 已经修改了
av_mat 



grad1, grad2, grad3  已经修改了
//zero_grad 
//gradient_F
//outgrad
//hess_apply
precon_grad
lbfgs
lbfgs_savegrad
//conjugate
//modelupdate

没有使用hess矩阵
所以注释了1095行的outgrad




bufferlef_to_rig
buffertop_to_bot...
修改了 exchange_v

sbufferlef_to_rig
修改了 exchange_s 


fbufferlef_to_rig
修改了 exchange_Ffv



pi
修改了
//readmod
//model
//outmod
//checkfd
//matcopy
//constant_boundary
//av_mat
//update_s_elastic
//surface_elastic
//seismo
//snap
//gradient_F
//cpmodel
//modelupdate

rip, rjp, rkp
修改了
//av_mat
//update_v


absorb_coeff 修改了
absorb
update_v


float ***  uipjp, float ***  ujpkp, float ***  uipkp 修改了
update_s_elastic
								
