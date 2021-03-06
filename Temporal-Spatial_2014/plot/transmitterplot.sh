#!/bin/sh
# $1 is the path of the folder containing results to be plotted
# $2 is the index of the first neuron to be plotted
# $3 is the index of the second neuron to be plotted
# $4 is the index of the third neuron to be plotted

i=1; j=300; k=400
DATAPATH=results
if (($#>0)); then
	if ! [[ $1 =~ ^[0-9]+([.][0-9]+)?$ ]]; then
		DATAPATH=$1
		shift
	fi
	if (($#>2)); then
		i=$1; j=$2; k=$3
	elif (($#>1)); then
		i=$1; j=$2;
	elif (($#>0)); then
		i=$1;
	fi
fi



/usr/local/bin/gnuplot -persist <<EOF
set term x11 "Synaptic activations"
set size 1,1
set origin 0,0
set multiplot

set size 1,0.3
set origin 0,0.7
set bmargin at screen 0.7
set format x ""
set ylabel "Neuron $i"
plot '$DATAPATH/Atransmit.dsv' using 1:$((i*2)) with lines title "AMPA", \
	'$DATAPATH/Atransmit.dsv' using 1:$((i*2+1))  with lines title "NMDA" 

set size 1,0.3
set origin 0,0.4
set tmargin at screen 0.7
set bmargin at screen 0.4
set ylabel "Neuron $j"
plot '$DATAPATH/Atransmit.dsv' using 1:$((j*2)) with lines title "AMPA", \
	'$DATAPATH/Atransmit.dsv' using 1:$((j*2+1)) with lines title "NMDA"
set ylabel "Fraction of Receptors Activated"

set size 1,0.4
set origin 0,0
set tmargin at screen 0.4
unset bmargin
set ylabel "Neuron $k"
set format x
set xlabel "Time (ms)"
plot '$DATAPATH/Atransmit.dsv' using 1:$((k*2)) with lines title "GABA-A", \
	'$DATAPATH/Atransmit.dsv' using 1:$((k*2+1)) with lines title "GABA-B"
unset multiplot
EOF

exit
