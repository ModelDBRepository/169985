#!/bin/bash

# $1 (optional) is the path of the folder containing results to be plotted
# After shift:
# $1 is the index of the neuron to be plotted
# $2 is the time of the second IV curve
# $4 is the time of the third IV curve

if (($#>4)); then
	DATAPATH=$1
	shift
else
	DATAPATH=results
fi

ts=( `awk 'END {print $1" "NR-2}' $DATAPATH/funcoutfile2.dsv` )
if [ `echo ${ts[0]} | awk -v FS=. 'NF>1'` ]; then
	tmax=`echo ${ts[0]} | awk -v FS=. 'NF>1 {print $1}'`
	tmax=$((tmax+1))
else
	tmax=`echo ${ts[0]} | awk -v FS=. 'NF==1'`
fi
Ntsteps=${ts[1]}

i=1; t1=$((tmax/3)); t2=$((tmax*2/3))
if (($#>2)); then
	i=$1; t1=$2; t2=$3
elif (($#>0)); then
	i=$1
fi



dt1=$((t1*Ntsteps/tmax+3))
dt2=$((t2*Ntsteps/tmax+3))
i1=$((i*4-2)); i2=$((i*4-1)); i3=$((i*4)); i4=$((i*4+1))
trans1=( `awk 'NR==1203||NR==t1||NR==t2 {print $1" "$i1" "$i2" "$i3" "$i4}' \
	t1=$dt1 t2=$dt2 i1=$i1 i2=$i2 i3=$i3 i4=$i4 $DATAPATH/funcoutfile2.dsv` )



/usr/local/bin/gnuplot -persist <<EOF
x11=0
if (x11==1) set term x11 size 650/2,450 title "I-V curve for Neuron $i"; \
else set term postscript color; set out 'figures/IVcurve.eps'





gNMDA(V)=1/(1+0.15*exp(-0.08*V))#*(1+0.15) 
gGABAB(V)=1/(1+exp(-(-90-V-10)/10))#*(1+exp(1)) 
#gGABAB(V)=1/(1+exp(-(-80-V-10)/10))#*(1+exp(1)) 
IV(V,GABAA,AMPA,GABAB,NMDA)=\
	-(GABAA*(-70-V)+AMPA*(-V)+GABAB*(-90-V)*gGABAB(V)+NMDA*(-V)*gNMDA(V))
#	-(GABAA*(-80-V)+AMPA*(-V)+GABAB*(-80-V)*gGABAB(V)+NMDA*(-V)*gNMDA(V))


# the integral is calculated as the sum of f(x_n)*delta 
#   do this x/delta times (from x down to 0)
delta = 0.02
#
# integral_f(x) takes one variable, the upper limit.  0 is the lower limit.
# calculate the integral of function f(t) from 0 to x
integral_IV(x,Ga,A,Gb,N) = (x>0)?integral1a(x,Ga,A,Gb,N):-integral1b(x,Ga,A,Gb,N)
integral1a(x,Ga,A,Gb,N) = (x<=0)?0:(integral1a(x-delta,Ga,A,Gb,N)+delta*IV(x,Ga,A,Gb,N))
integral1b(x,Ga,A,Gb,N) = (x>=0)?0:(integral1b(x+delta,Ga,A,Gb,N)+delta*IV(x,Ga,A,Gb,N))
#
# integral2_f(x,y) takes two variables; x is the lower limit, and y the upper.
# claculate the integral of function f(t) from x to y
integral2_f(x,y) = (xy)?0:(integral2(x+delta,y)+delta*f(x))




energy=1;	# if energy==1, script plots energy function, otherwise script plots IV curve
if (energy==1) plot_IV(a,b,c,d,e)=integral_IV(a,b,c,d,e); Vmin=-100; Vmax=20;\
else plot_IV(a,b,c,d,e)=-IV(a,b,c,d,e); Vmin=-90; Vmax=0

set xlabel "Voltage (mV)"
if (energy==1) set title "Energy Landscape for Neuron $i (Local Minima are Stable Voltages)"; set ylabel "Energy (nA to return voltage to 0)"; \
else set title "I-V Curve for Neuron $i"; set ylabel "Current (nA)"
plot [V=Vmin:Vmax] plot_IV(V,${trans1[3]},${trans1[1]},${trans1[4]},${trans1[2]}) \
			title "t=${trans1[0]} ms" lt 1 lw 4, \
		plot_IV(V,${trans1[8]},${trans1[6]},${trans1[9]},${trans1[7]}) \
			title "t=${trans1[5]} ms" lt 2 lw 4, \
		plot_IV(V,${trans1[13]},${trans1[11]},${trans1[14]},${trans1[12]}) \
			title "t=${trans1[10]} ms" lt 3 lw 4, \
			0 lt 0 notitle





unset multiplot
EOF
