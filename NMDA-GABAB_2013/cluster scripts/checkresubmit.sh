#!/bin/bash
DIRY=/home/honi/lismanlab/newnetwork
pwd=$PWD; cd $DIRY
rm simul.* simuls.* 2>/dev/null
echo "**********MAIN.CPP**********" > results/sortedsweep.dsv
cat main.cpp >> results/sortedsweep.dsv
echo -e "\n**********FUNCTIONS.CPP**********" >> results/sortedsweep.dsv
cat functions.cpp >> results/sortedsweep.dsv
echo -e "\n**********NEWNETWORK.H**********" >> results/sortedsweep.dsv
cat newnetwork.h >> results/sortedsweep.dsv


queue="a"       # "n" for neuro, "a" for all


base=2
if [ $queue = "a" ]; then
  MAXSIMULS=600
elif [ $queue = "n" ]; then
  MAXSIMULS=50
fi


for ((OOMi=8; OOMi<=24; OOMi++)); do    #gNMDAee
    i=$(((4**OOMi + 3**OOMi/2)/(3**OOMi)))




for ((OOMj=1; OOMj<=12; OOMj++)); do    #gGABAA
    if ((OOMj==0)); then j=0.3
    elif ((OOMj==2)); then j=1.3
    else j=$(((4**OOMj + 3**OOMj/2)/(3**OOMj)))
    fi


while (( `qstat | awk 'END {print NR}'` > MAXSIMULS )); do
  sleep 100
done





for Ns in {1..7}{,5}; do                  #N_prf
    N=0.$Ns
#    if [ $N == 0.45 ]; then
#	N=0.05
#    fi





ks=0
for ((OOMk=0; OOMk<=10; OOMk++)); do     #gGABAB
    ks="$ks $((base**$OOMk))"
done

for k in $ks; do






awk 'BEGIN {test=0}
$1==a && $2==b && $3==c && $4==d && NF<14 {print; test=1; exit 0}
END {if (test==0) {
        if (q=="n") {
          command="xargs qsub -cwd -l neuro /home/honi/lismanlab/newnetwork/newnetwork.submit >/dev/null" 
        }
        else if (q=="a") {
          command="xargs qsub -cwd -ckpt reloc /home/honi/lismanlab/newnetwork/newnetwork.submit >/dev/null" 
        }
        else {
          command="xargs qsub -cwd -ckpt reloc /home/honi/lismanlab/newnetwork/newnetwork.submit >/dev/null" 
	}
#command="xargs echo "
        print a, b, c, d | command
        close(command)
       print a"\t"b"\t"c"\t"d"\tNot yet run."
    }
}
' a=$i b=$j c=$N d=$k q=$queue \
results/sweepresults.dsv >> results/sortedsweep.dsv #2>/dev/null





done    # gGABAB

done    # N_prf

done    # gGABAA

done    # gNMDAee



cd $PWD

