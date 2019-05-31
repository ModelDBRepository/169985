#!/bin/sh
	#$1 is the path of the folder containing results to be plotted
if (("$#">"0")); then
  path=$1
else
  #echo no path given
  path=results
fi
# echo reading from $PATH/Aspikecoords.dsv

pwd=$PWD
cd $path

switch=0		# =0 for normal, =1 for mixed indices

# move Neurons 1-50 to 201-250
if [ $switch==1 ]; then 
	a=50; b=200
	awk '$2<a {print $1"\t"$2-a} $2>=a {print}' a=$a Aspikecoords.dsv > Aspikecoordsa.dsv
	awk '$2>b&&$2<(a+b) {print $1"\t"$2-b} $2<b {print} $2>a+b {print}' a=$a b=$b Aspikecoordsa.dsv > Aspikecoordsb.dsv
	awk '$2<0 {print $1"\t"($2+b+a)} $2>=0 {print}' a=$a b=$b Aspikecoordsb.dsv > Aspikecoordsc.dsv
fi

:<<EOC
a=0; b=$((a+200));
awk '$2==a {print $1"\t-1"} $2!=a {print}' a=$a b=$b Aspikecoords.dsv > Aspikecoordsa.dsv
  # mark
awk '$2==b {print $1"\t"a} $2!=b {print}' a=$a b=$b Aspikecoordsa.dsv > Aspikecoordsb.dsv
awk '$2==-1 {print $1"\t"b} $2!=-1 {print}' a=$a b=$b Aspikecoordsb.dsv > Aspikecoordsc.dsv

for n in {1..50}; do
  a=$n; b=$((a+200));
  awk '$2==a {print $1"\t-1"} $2!=a {print}' a=$a b=$b Aspikecoordsc.dsv > Aspikecoordsa.dsv
  awk '$2==b {print $1"\t"a} $2!=b {print}' a=$a b=$b Aspikecoordsa.dsv > Aspikecoordsb.dsv
  awk '$2==-1 {print $1"\t"b} $2!=-1 {print}' a=$a b=$b Aspikecoordsb.dsv > Aspikecoordsc.dsv
done
EOC

/usr/local/bin/gnuplot -persist <<EOF
set mouse
x11=1
if (x11==1) set term x11 title "Network Spiking Activity"; \
else set term postscript color; set out '../figures/network.eps'
set size 1,1
set origin 0,0
unset key
#set title "Network Spiking Activity"
set ylabel "Neuron Index"
set xlabel "Time (ms)"
set yrange [0:*]
set xrange [0:*]
#set ytics (120,320)
if ($switch==1) plot 'Aspikecoordsc.dsv' using 1:2 with points ls 0; \
	else plot 'Aspikecoords.dsv' using 1:2 with points ls 0
EOF

if [ $switch==1 ]; then
	rm Aspikecoordsa.dsv Aspikecoordsb.dsv Aspikecoordsc.dsv
fi
cd $pwd
exit
