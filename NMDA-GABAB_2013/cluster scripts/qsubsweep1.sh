#!/bin/bash
DIRY=/home/honi/lismanlab/newnetwork
pwd=$PWD; cd $DIRY
SUMMARY_FILE=results/sweepresults.dsv

queue="a"	# "n" for neuro, "a" for all

rm simul.* simuls.* 2>/dev/null

if [[ $1 == "-c" ]]; then
  if (( 0<`qstat | awk 'END {print NR}'` )); then
    qdel -u honi 1>/dev/null
  fi
  g++ -g main.cpp functions.cpp -o newnetwork -O2
  {
  echo  "#Target (TEC) and Non-Target (NTEC) Excitatory Cell Firing Rates (FRs) at Different Times During Simulation for Different Parameters" 
  echo -e "#gNMDA\tN_TEC/N_EC\tInput TECFR\tPost-Input TECFR\tPost-Input FR ratio"
  echo -e "#\tgGABAie\tgGABAB\tInput NTECFR\tPost-Input NTECFR"
  } > $SUMMARY_FILE
fi


date > log

base=2
if [ $queue = "a" ]; then
  MAXSIMULS=500
elif [ $queue = "n" ]; then
  MAXSIMULS=50
fi




for ((OOMi=8; OOMi<=24; OOMi++)); do    #gNMDAee
    i=$(((4**OOMi + 3**OOMi/2)/(3**OOMi)))


while (( `qstat | awk 'END {print NR}'`>MAXSIMULS )); do
  sleep 300
done


for ((OOMj=1; OOMj<=12; OOMj++)); do    #gGABAA
    if ((OOMj==0)); then j=0.3
    elif ((OOMj==2)); then j=1.3
    else j=$(((4**OOMj + 3**OOMj/2)/(3**OOMj)))
    fi

: '
for Ns in {1..7}{,5}; do                  #N_prf
    N=0.$Ns
#    if [ $N == 0.45 ]; then
#        N=0.05
#    fi

'

: '
ks=
for ((OOMk=1; OOMk<=13; OOMk++)); do     #gGABAB
    ks="$ks $(((base*1)**$OOMk))"
done

for k in $ks; do
'


if [ $queue = "a" ]; then
  qsub -cwd -ckpt reloc newnetwork.submit1 $i $j 1>>log
elif [ $queue = "n" ]; then
  qsub -cwd -l neuro newnetwork.submit1 $i $j 1>>log
fi
#qsub -cwd -ckpt reloc newnetwork.submit $i $j $N $k 1>>log




#done	# gGABAB

#done	# N_prf

done	# gGABAA

done	# gNMDAee
		

cd $pwd
exit
