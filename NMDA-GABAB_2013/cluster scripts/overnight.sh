#!/bin/bash

run_qsubsweep () {	# run qsubsweep and wait til jobs are done
  ./qsubsweep1.sh -c	# compiles newnetwork and submits jobs

  # check to make sure that the jobs do not produce any errors
  Nrun=`ls simuls.e* 2>/dev/null |awk 'END {print NR}'`
  while (( Nrun == 0 )); do
    sleep 300
    Nrun=`ls simuls.e* 2>/dev/null |awk 'END {print NR}'`
  done
  errorfile=`ls simuls.e* |awk 'NR==1'`	
  if (( `cat $errorfile | awk 'END {print NR}'` > 0 )); then
    exit 200
  fi

  # sleep while there are still queued jobs
  Nsimuls=`qstat | awk 'END {print NR}'`
  while (( Nsimuls > 0 )); do
    Nrunning=`qstat | grep r | awk 'END {print NR}'`
#    sleep $((Nsimuls*300/(Nrunning+10)))
    sleep 300
    if (( `qstat |grep Eqw | awk 'END {print NR}'` > 0 )); then
      qstat |grep Eqw | awk '{print $1}' | xargs qdel
    fi
    Nsimuls=`qstat | awk 'END {print NR}'`
  done
}



run_checkresubmit () {	# run checkresubmit until it has no effect
  # if the last checkresubmit produces output files, run checkresubmit again
  ndone=1
  while (( ndone > 0 )); do
    ./checkresubmit.sh
    Nsimuls=`qstat | awk 'END {print NR}'`	# wait til jobs are done
    while (( Nsimuls > 0 )); do
      Nrunning=`qstat | grep r | awk 'END {print NR}'`
      sleep $((Nsimuls*30/Nrunning))
      if (( `qstat |grep Eqw | awk 'END {print NR}'` > 0 )); then
        qstat |grep Eqw | awk '{print $1}' | xargs qdel
      fi
      Nsimuls=`qstat | awk 'END {print NR}'`
    done
    ndone=`ls simul* 2>/dev/null | awk 'END {print NR}'`

    # check to make sure that the jobs do not produce any errors
    if (( ndone > 0 )); then
      errorfile=`ls simul.e* |awk 'NR==1'`
      if (( `cat $errorfile | awk 'END {print NR}'` > 0 )); then
        exit 200
      fi
    fi
  done
}



######## Commands start here ########

:<<EOC
# No longer necessary because of #ifdef __APPLE__ protection on writing
# ensure  that output is suppressed
mv main.cpp main1.cpp
sed 's:\tAV_sfile.o://AV_sfile.o:' < main1.cpp > main.cpp
mv main.cpp main1.cpp
sed 's:\tAV_dfile.o://AV_dfile.o:' < main1.cpp > main.cpp
mv main.cpp main1.cpp
sed 's:\tAspikesfile.o://Aspikesfile.o:' < main1.cpp > main.cpp
mv main.cpp main1.cpp
sed 's:\tAtransmitfile.o://Atransmitfile.o:' < main1.cpp > main.cpp
rm main1.cpp
EOC





run_qsubsweep
run_checkresubmit



:<<EOC
./qsubsweep1.sh -c
Nsimuls=`qstat | awk 'END {print NR}'`
while (( Nsimuls > 0 )); do
  Nrunning=`qstat | grep r | awk 'END {print NR}'`
  sleep $((Nsimuls*300/Nrunning))
done
./checkresubmit.sh
Nsimuls=`qstat | awk 'END {print NR}'`
while (( Nsimuls > 0 )); do
  Nrunning=`qstat | grep r | awk 'END {print NR}'`
  sleep $((Nsimuls*30/Nrunning))
done
./checkresubmit.sh
Nsimuls=`qstat | awk 'END {print NR}'`
while (( Nsimuls > 0 )); do
  Nrunning=`qstat | grep r | awk 'END {print NR}'`
  sleep $((Nsimuls*30/Nrunning))
done


cp functions.cpp functions1.cpp
sed 's|//w_in\[s\] = connect|w_in\[s\] = connect|' < functions1.cpp > functions.cpp

mv main.cpp main1.cpp
sed 's/eiAMPAscale=32.0/eiAMPAscale=25.0/' < main1.cpp > main.cpp

./qsubsweep1.sh -c
sleep 7200
./checkresubmit.sh
sleep 2400
./checkresubmit.sh
sleep 2400
mv results/sortedsweep.dsv results/noinputtoinhibitsortedsweep.dsv

mv main.cpp main1.cpp
sed 's/eiAMPAscale=32.0/eiAMPAscale=25.0/' < main1.cpp > main.cpp

./qsubsweep1.sh -c
sleep 7200
./checkresubmit.sh
sleep 2400
./checkresubmit.sh
sleep 2400
mv results/sortedsweep.dsv results/loweiAMPAsortedsweep.dsv
EOC

echo -e "\a"

if [ -x nohup.out ]; then
  if (( `cat nohup.out | awk 'END {print NR}'` < 2 )); then
    rm nohup.out
  fi
fi

exit 0
