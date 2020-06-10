dataFolder="dat"
  instanceFolder=${dataFolder}/Instances_345
  resultsFolder=${dataFolder}/results_Instances345_Outer
  logFile=${resultsFolder}/results.csv
  numThreads=4
  mkdir ${resultsFolder}
: '
  for ex in $(ls ${instanceFolder}/*.json); do
        instanceName=${ex/.json/}
        instanceNumber=${ex//[!0-9]/}
        echo "----------------------Running instance $ex----------------------"
        #Skip full enumeration with the first argument
        if [ ! $1 ]; then
          #echo "--------Full Enumeration"
          #./EPEC -i ${ex//.json/} -t ${numThreads} -w 2 -l ${logFile} -s ${resultsFolder}/Solution_Full_${instanceNumber} --timelimit 1800 -a 0
          #./EPEC -i ${ex//.json/} -t ${numThreads} -w 2 -l ${logFile} -s ${resultsFolder}/Solution_FullPure_${instanceNumber} --timelimit 1800 -a 0 --pure 1
          #./EPEC -i ${ex//.json/} -t ${numThreads} -w 2 -l ${logFile} -s ${resultsFolder}/Solution_${instanceNumber} --timelimit 1800 -a 1 --pure 1 --aggr 3 --add 1 --recover 1
          #./EPEC -i ${ex//.json/} -t ${numThreads} -w 2 -l ${logFile} -s ${resultsFolder}/Solution_${instanceNumber} --timelimit 1800 -a 2 --pure 1
          ./EPEC -i ${ex//.json/} -t ${numThreads} -w 2 -l ${logFile} -s ${resultsFolder}/Solution_OuterPreliminary_${instanceNumber} --timelimit 1800 -a 3
        fi
        for Aggressiveness in {1,3,5}; do
          for AddPolyMethod in {0,1,2}; do
            #echo "--------Inner Approximation"
            echo ""
            #./EPEC -i ${ex//.json/} -t ${numThreads} -w 2 -l ${logFile} -s ${resultsFolder}/Solution_Inner_${instanceNumber}_aggr${Aggressiveness}_method${AddPolyMethod} --timelimit 1800 -a 1 --aggr $Aggressiveness --add $AddPolyMethod
          done
        done
        printf "\n\n"
  done'
for (( c=1; c<=50; c++ ))
do
  echo "Running ${instanceFolder}/Instance_$c"
	./EPEC -i ${instanceFolder}/Instance_$c -t ${numThreads} -w 2 -l ${logFile} -s ${resultsFolder}/Solution_OuterPreliminary_$c --timelimit 1800 -a 3 -m 1
done