#!/usr/bin/python
import os
import sys
def create_rsf(n1, n2, n3, in_file, rsf_file):
  crt_com='sfmath n1=%d n2=%d n3=%d output=0 > %s' % (n1, n2, n3, rsf_file)
  sub_com='sed -i "s/in=.*/in=%s/g" %s' % (in_file, rsf_file)
  os.popen(crt_com)
  os.popen(sub_com)

#help
if len(sys.argv) < 3 or (len(sys.argv) == 3 and (sys.argv[2][0] != 'n' or sys.argv[2][2] != '=')) or (len(sys.argv) == 4 and (sys.argv[2][0] != 'n' or sys.argv[2][2] != '=')):
  print 'Usage: ./draw.py in_file n1=%d [n2=%d]'
  print 'Example: ./draw.py toy.vp_it60 n1=68 n2=32'
  exit()

#arguments
prefix='toy'
in_file=sys.argv[1]
in_true_file='toy.vp.true'    #It needs to change!
rsf_file=prefix+'.rsf'
rsf_win_file=prefix+'_win.rsf'
rsf_true_file='%s.rsf' % in_true_file
rsf_true_win_file='%s_win.rsf' % in_true_file
vpl_file=prefix+'.vpl'
title="VP"
label1="Depth"
label2="Distance"
unit1="m"
unit2="m"
n1=184  #y
n2=160  #x
n3=160  #z

#create rsf file 
create_rsf(n1, n2, n3, in_file, rsf_file)
create_rsf(n1, n2, n3, in_true_file, rsf_true_file)

#get average bias
avg_com="sfattr want=mean < %s | awk '{print $3}'" % (rsf_file)
mean=float(os.popen(avg_com).read().strip())

#sfwindow and sfgrey
index1=int(sys.argv[2][1])
first1=int(sys.argv[2][3:])
pen=''
if len(sys.argv) == 3:
  win_com='sfwindow f%d=%d n%d=1 < %s' % (index1, first1, index1, rsf_file)
  plot_com='sfgrey color=j title=%s label1=%s label2=%s unit1=%s unit2=%s bias=%d scalebar=y barreverse=y > %s' % (title, label1, label2, unit1, unit2, mean, vpl_file)
  pen_com='sfpen %s' % (vpl_file)
  t_com = '%s | %s' % (win_com, plot_com)
  os.popen(t_com)
elif len(sys.argv) == 4:
  index2=int(sys.argv[3][1])
  first2=int(sys.argv[3][3:])
  win_com='sfwindow f%d=%d n%d=1 f%d=%d n%d=1 < %s > %s' % (index1, first1, index1, index2, first2, index2, rsf_file, rsf_win_file)
  win_com_true='sfwindow f%d=%d n%d=1 f%d=%d n%d=1 < %s > %s' % (index1, first1, index1, index2, first2, index2, rsf_true_file, rsf_true_win_file)

  plot_com='sfcat axis=2 %s %s | sfgraph title="%s" dash=0 plotfat=5  grid=y yreverse=n > %s' % (rsf_win_file, rsf_true_win_file, title, vpl_file)
  pen_com='sfpen %s' % (vpl_file)
  os.popen(win_com)
  os.popen(win_com_true)
  os.popen(plot_com)

#sfpen
os.popen(pen_com)
