#This plot was made using commands from Seismic Unix!!

psimage < wave_bicgstab_100.bin n1=201 n2=201 d1=25 d2=25 f1=-2500 f2=-2500 perc=99 wbox=6 hbox=6 label1='X(m)' label2='Z(m)' title='(a) BiCGStab-100'> wave_bicgstab_100.eps

psimage < wave_bicgstab_300.bin n1=201 n2=201 d1=25 d2=25 f1=-2500 f2=-2500 perc=99 wbox=6 hbox=6 label1='X(m)' label2='Z(m)' title='(b) BiCGStab-300'> wave_bicgstab_300.eps

psimage < wave_gmres_100.bin n1=201 n2=201 d1=25 d2=25 f1=-2500 f2=-2500 perc=99 wbox=6 hbox=6 label1='X(m)' label2='Z(m)' title='(c) GMRES-100'> wave_gmres_100.eps

psimage < wave_gmres_300.bin n1=201 n2=201 d1=25 d2=25 f1=-2500 f2=-2500 perc=99 wbox=6 hbox=6 label1='X(m)' label2='Z(m)' title='(d) GMRES-300'> wave_gmres_300.eps

psmerge in=wave_bicgstab_100.eps translate=0.,0. in=wave_bicgstab_300.eps translate=8,0 \
	in=wave_gmres_100.eps translate=0.,-7.5 in=wave_gmres_300.eps translate=8,-7.5 > final.eps


epstopdf final.eps 
rm *.eps
