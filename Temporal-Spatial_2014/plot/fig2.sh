#!/bin/sh
# $1 is the path of the folder containing results to be plotted
# $2 is the index of the neuron for which data should be plotted
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
cwd=$PWD
cd $DATAPATH

# echo reading from $DATAPATH/funcoutfile.dsv


# rearrange neurons for network plot
a=50; b=200
awk '$2<a {print $1"\t"$2-a} $2>=a {print}' a=$a Aspikecoords.dsv > Aspikecoordsa.dsv
awk '$2>b&&$2<(a+b) {print $1"\t"$2-b+a} $2<b {print} $2>a+b {print}' a=$a b=$b Aspikecoordsa.dsv > Aspikecoordsb.dsv
awk '$2<0 {print $1"\t"($2+b+a)} $2>=0 {print}' a=$a b=$b Aspikecoordsb.dsv > Aspikecoordsc.dsv

# set term x11 size 650,450*2

/usr/local/bin/gnuplot -persist <<EOF
x11=1
if (x11==1) set term x11 size 650,650; \
else set term postscript portrait color size 7,7; set out '$cwd/figures/fig2.eps'

set multiplot

set origin 0,0.65
set size 1,0.35
set border 3; set xtics nomirror; set ytics nomirror
set bmargin at screen 0.65
set lmargin at screen 0.1
unset key
#set title "Network Spiking Activity"
set ylabel "Neuron Index"
set format x ""
set yrange [0:*]
#set ytics (120,320)
plot 'Aspikecoordsc.dsv' using 1:2 with points 

set origin 0,0.35
set size 1,0.3
set tmargin at screen 0.65
set bmargin at screen 0.35
set lmargin at screen 0.1
set ylabel "Dendrite (mV)"
#set ylabel "Dendritic Voltage (mV)"
set yrange [-85:*]
set key right center
plot   'AV_d.dsv' using 1:$((i+1)) title "Neuron $i" with lines lt 1,\
        'AV_d.dsv' using 1:$((j+1)) title "Neuron $j" with lines lt 2

set origin 0,0
set size 1,0.35
set tmargin at screen 0.35
set lmargin at screen 0.1
unset bmargin
set ylabel "Soma (mV)"
#set ylabel "Somatic Voltage (mV)"
set yrange [-90:50]
set format x
set xlabel "Time (ms)"
unset key
plot   'AV_s.dsv' using 1:$((i+1)) title "Neuron $i" with lines lt 1,\
        'AV_s.dsv' using 1:$((j+1)) title "Neuron $j" with lines lt 2
EOF

cd $cwd

exit
