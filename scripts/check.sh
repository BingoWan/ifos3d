#! /bin/bash
cp plot_seismogram.sh su_obs
cd su_obs
if [ ! -f obs_vx_it1.su.shot1.rsf ]; then
./plot_seismogram.sh obs_vx_it1.su.shot1
fi
sfattr < obs_vx_it1.su.shot1.rsf > rsf.dat
diff ./rsf.dat  /home/export/online1/wanwb/lasefis/ifos3d/par/marmousi1/su_obs/rsf.dat
cd ..
sfdisfil < ./su_obs/obs_vx_it1.su.shot1.rsf number=y col=2 format="%20.10g" | head -n 20000 |grep 19998:  |tee  1.dat
sfdisfil < ~/online1/lasefis/ifos3d/par/marmousi1/su_obs/obs_vx_it1.su.shot1.rsf number=y col=2 format="%20.10g" | head -n 20000|grep 19998: |tee 2.dat
diff 1.dat 2.dat

