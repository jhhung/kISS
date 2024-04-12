#!/bin/bash
cd ../build && make -j
OUTPUT="/home/dave/kISS/experiment/sa_20240409.csv"
# EXP="zebrafish_protein_437"
EXP="chm13v2.0 mouse zebrafish human_protein_5640 mouse_protein_589 zebrafish_protein_437"
echo "algo,test,k,num-threads,time,space" > $OUTPUT

# pardss

for j in 1 2 3; do
    for exp in $EXP; do
        file="/home/dave/data/${exp}.fa"
        num_threads=32
        log="./sa/pardss_${exp}_${num_threads}_${j}.txt"

        if [ $num_threads -eq 128 ]; then
            cpu_list="0-127"
        else
            end=$(( $num_threads + $num_threads - 2 ))
            cpu_list="0-$end:2"
        fi

        is_general=0

        taskset -c $cpu_list ./experiment/pardss $file $num_threads $is_general > $log 2>&1
        space=`tac $log | sed -n 1p | awk '{ print $7 }'`
        t=`tac $log | sed -n 2p | awk '{ print $7 }'`
        echo "pardss,${exp},${k},${num_threads},${t},${space}" >> $OUTPUT
    done
done

# saca-k

for j in 1 2 3; do
    for exp in $EXP; do
        file="/home/dave/data/${exp}.fa"
        num_threads=1
        log="./sa/sacak_${exp}_${num_threads}_${j}.txt"

        cpu_list="0-0"

        is_general=0

        taskset -c $cpu_list ./experiment/saca-k $file $num_threads $is_general > $log 2>&1
        space=`tac $log | sed -n 1p | awk '{ print $7 }'`
        t=`tac $log | sed -n 2p | awk '{ print $7 }'`
        echo "saca-k,${exp},${k},${num_threads},${t},${space}" >> $OUTPUT
    done
done