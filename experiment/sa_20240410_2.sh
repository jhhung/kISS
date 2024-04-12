#!/bin/bash
cd ../build && make -j
OUTPUT="/home/dave/kISS/experiment/sa_20240410_2.csv"
EXP="zebrafish human_protein_5640 mouse_protein_589 zebrafish_protein_437"
echo "algo,test,k,num-threads,time,space" > $OUTPUT

# pardss

for j in 1 2 3; do
    for exp in $EXP; do
        file="/home/dave/data/${exp}.fa"
        num_threads=32
        log="./sa/pardss32_${exp}_${num_threads}_${j}.txt"

        if [ $num_threads -eq 128 ]; then
            cpu_list="0-127"
        else
            end=$(( $num_threads + $num_threads - 2 ))
            cpu_list="0-$end:2"
        fi

        # human_protein_5640 mouse_protein_589 zebrafish_protein_437
        if [ ${exp} = "human_protein_5640" ] || [ ${exp} = "mouse_protein_589" ] || [ ${exp} = "zebrafish_protein_437" ]; then
            is_general="1"
        else
            is_general="0"
        fi

        taskset -c $cpu_list ./experiment/pardss-32 $file $num_threads $is_general > $log 2>&1
        space=`tac $log | sed -n 1p | awk '{ print $7 }'`
        t=`tac $log | sed -n 2p | awk '{ print $7 }'`
        echo "pardss-32,${exp},${k},${num_threads},${t},${space}" >> $OUTPUT
    done
done

# libsais

for j in 1 2 3; do
    for exp in $EXP; do
        num_threads=32
        file="/home/dave/data/${exp}.fa"
        k=-1
        log="./sa/libsais32_${exp}_${k}_${num_threads}_${j}.txt"

        if [ $num_threads -eq 128 ]; then
            cpu_list="0-127"
        else
            end=$(( $num_threads + $num_threads - 2 ))
            cpu_list="0-$end:2"
        fi

        # human_protein_5640 mouse_protein_589 zebrafish_protein_437
        if [ ${exp} = "human_protein_5640" ] || [ ${exp} = "mouse_protein_589" ] || [ ${exp} = "zebrafish_protein_437" ]; then
            is_general="1"
        else
            is_general="0"
        fi

        taskset -c $cpu_list ./experiment/libsais-test-32 $file $num_threads $is_general > $log 2>&1
        space=`tac $log | sed -n 1p | awk '{ print $7 }'`
        t=`tac $log | sed -n 2p | awk '{ print $7 }'`
        echo "libsais-32,${exp},${k},${num_threads},${t},${space}" >> $OUTPUT
    done
done