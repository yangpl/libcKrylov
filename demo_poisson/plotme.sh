#This plot was made using commands from Seismic Unix!!

n1=129
n2=129
d1=0.007825
d2=0.007825
mylegend=' legend=1 lstyle=vertright lwidth=0.2 lheight=6'
cmap='wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0'

psimage < true.bin n1=$n1 n2=$n2 d1=$d1 d2=$d2 $mylegend $cmap wbox=6 hbox=6 label1='X' label2='Z' title='(a) True solution'> true.eps

psimage < cg_error.bin n1=$n1 n2=$n2 d1=$d1 d2=$d2 $mylegend $cmap wbox=6 hbox=6 label1='X' label2='Z' title='(b) Error by CG'> cg_error.eps

psimage < pcg_error.bin n1=$n1 n2=$n2 d1=$d1 d2=$d2 $mylegend $cmap wbox=6 hbox=6 label1='X' label2='Z' title='(c) Error by PCG'> pcg_error.eps

psmerge in=true.eps translate=0.,0. in=cg_error.eps translate=7.8,0 \
	in=pcg_error.eps translate=15.6,0 > final.eps


epstopdf final.eps 
rm *.eps
