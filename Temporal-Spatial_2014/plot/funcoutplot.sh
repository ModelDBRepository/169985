#!/bin/sh
# $1 is the path of the folder containing results to be plotted
# $2 is the index of the neuron for which data should be plotted
if [ "-2" == $1 ] 2>/dev/null; then
	filename=funcoutfile2.dsv
	shift
else
	filename=funcoutfile1.dsv
fi

i=1; j=300
if (($#>2)); then
    DATAPATH=$1
    i=$2
    j=$3
elif (($#>1)); then
  if (($1>=1)); then    # if the first arg is a number
    DATAPATH=results
    i=$1
    j=$2
  else
    DATAPATH=$1
    i=$2; j=$i
  fi
elif (($#>0)); then
    DATAPATH=$1
else
#    echo no path given
    DATAPATH=results
fi

# echo reading from $DATAPATH/$filename


/usr/local/bin/gnuplot -persist <<EOF
x11=1
if (x11==1) set term x11 size 650,450*1.5; \
else set term postscript portrait color; set out 'figures/funcout.eps'
set multiplot

set origin 0,0.75
set size 1,0.25
set bmargin at screen 0.75
set lmargin at screen 0.1
set format x ""

#set key right center
#set ylabel "Dendritic Voltage (mV)"
#set yrange [-85:*]
#plot   '$DATAPATH/AV_d.dsv' using 1:$((i+1)) title "Neuron $i" with lines lt 1,\
#        '$DATAPATH/AV_d.dsv' using 1:$((j+1)) title "Neuron $j" with lines lt 2

set key left top
set ylabel "AMPA current"
plot 	'$DATAPATH/$filename' using 1:$((i*4-2)) title "Neuron $i" with lines lt 1, \
	'$DATAPATH/$filename' using 1:$((j*4-2)) title "Neuron $j" with lines lt 2



set origin 0,0.5
set bmargin at screen 0.5
set tmargin at screen 0.75
unset key

#set ylabel "Somatic Voltage (mV)"
#plot   '$DATAPATH/AV_s.dsv' using 1:$((i+1)) title "Neuron $i" with lines lt 1,\
#        '$DATAPATH/AV_s.dsv' using 1:$((j+1)) title "Neuron $j" with lines lt 2
#set yrange [*:*]

set ylabel "GABA-A current"
plot 	'$DATAPATH/$filename' using 1:$((i*4+0)) title "Neuron $i" with lines lt 1, \
	'$DATAPATH/$filename' using 1:$((j*4+0)) title "Neuron $j" with lines lt 2



set size 1,0.25
set origin 0,0.25
set bmargin at screen 0.25
set tmargin at screen 0.5
set ylabel "NMDA current"
plot 	'$DATAPATH/$filename' using 1:$((i*4-1)) title "Neuron $i" with lines lt 1, \
	'$DATAPATH/$filename' using 1:$((j*4-1)) title "Neuron $j" with lines lt 2


set size 1,0.25
set origin 0,0
set bmargin at screen 0
set tmargin at screen 0.25
set bmargin
set ylabel "GABA-B current"
set format x
set xlabel "Time (ms)"
plot 	'$DATAPATH/$filename' using 1:$((i*4+1)) title "Neuron $i" with lines lt 1, \
	'$DATAPATH/$filename' using 1:$((j*4+1)) title "Neuron $j" with lines lt 2
unset multiplot






EOF

exit
