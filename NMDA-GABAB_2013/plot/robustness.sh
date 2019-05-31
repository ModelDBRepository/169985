#!/bin/sh
i=0; j=128

path=results

file=sortedsweep.dsv

if [[ $1 == "-f" ]]; then
  file=$2
  path=$PWD
  shift 2
fi
if (($#>1)); then
  i=$1; j=$2;
elif (($#>0)); then
  i=$1
fi

pwd=$PWD
cd $path

bardata=$file
spacedata=$file

awk '$7>50&&$8<20&&$4==i' i=$i $spacedata > plotsweep0.dsv 
awk '$7>50&&$8<20&&$4==j' j=$j $spacedata > plotsweep1.dsv 

test=1
if (($test)); then
  awk '$7<10&&$4==i' i=$i $spacedata > plotsweep0a.dsv 
  awk '$9>0.5&&$9<5&&$4==i' i=$i $spacedata > plotsweep0b.dsv 
  awk '$7<10&&$4==j' j=$j $spacedata > plotsweep1a.dsv 
  awk '$9>0.5&&$9<5&&$4==j' j=$j $spacedata > plotsweep1b.dsv 
fi


barfile=barplot.dsv
rm -f $barfile


ks=0; base=2; xmax=10
for ((OOMk=0; OOMk<=xmax; OOMk++)); do     #gGABAB
    ks="$ks $(((base*1)**OOMk))"
#    ks="$ks $(((base*1)**OOMk))"
done
xmax=$((base**xmax))

for k in $ks; do
total=`awk '$4==k' k=$k $bardata | awk 'END {print NR}'` 
awk '$4==k&&$7>50&&$8<20' k=$k $bardata | awk 'END {print k"\t"NR/tot*100}' k=$k tot=$total >> $barfile
done

cent=100



/usr/local/bin/gnuplot -persist 2>/dev/null <<EOF
set mouse
x11=0
if (x11==1) set term x11 title "$file"; \
else set term postscript color; set out "$pwd/figures/robustness.eps"
set multiplot


set size 0.5,0.45
set origin 0,0.54
set bmargin at screen 0.54
set rmargin at screen 0.5
set format x ""
set ylabel "GABA-A"
set log xy
plot "plotsweep0a.dsv" using 1:2 title "GABA-B=$i" ls 2 lc 3 ps 2, \
	"plotsweep0b.dsv" using 1:2 notitle ls 1 lc 1 ps 2, \
	"plotsweep0.dsv" using 1:2 notitle ls 7 lc 0 ps 1

set origin 0,0
set size 0.5,0.54
set tmargin at screen 0.54
set bmargin
set format x
set xlabel "NMDA"
unset title
plot "plotsweep1a.dsv" using 1:2 title "GABA-B=$j" ls 2 lc 3 ps 2, \
	"plotsweep1b.dsv" using 1:2 notitle ls 1 lc 1 ps 2, \
	"plotsweep1.dsv" using 1:2 notitle ls 7 lc 0 ps 1


unset key
unset log
set origin 0.65,0
set size 0.35,0.5
set lmargin 0
set rmargin
set tmargin
if (x11==0) set border 1
set xrange [0.7:$((2*xmax))]
set yrange [0:*] writeback
#unset ytics
set xtics nomirror
set format y ""
unset ylabel
set log x
set xlabel "GABA-B conductance"
plot '$barfile' using 1:2 with boxes

set origin 0.5,0
set size 0.13,0.5
set rmargin 0
set lmargin
if (x11==0) set border 3
set ylabel "Proportion of Succesful Simulations"
set xlabel " "
unset logscale x
set xrange [-0.1:0.1]
set yrange restore
set format y
set ytics nomirror
#set yrange [0:0.05]
set xtics (0)
plot '$barfile' using 1:2 with boxes
unset xtics; set xtics

EOF

rm plotsweep0.dsv
rm plotsweep1.dsv
if (($test)); then
  rm plotsweep0a.dsv
  rm plotsweep0b.dsv
  rm plotsweep1a.dsv
  rm plotsweep1b.dsv
fi

cd $pwd

exit

#
