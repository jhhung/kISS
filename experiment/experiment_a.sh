#!/bin/bash
SCRIPT_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
DATA_DIR=$(realpath $1)
LOG_DIR=$(realpath $2)
cd "$SCRIPT_DIR/../build/" && make -j
OUTPUT="$SCRIPT_DIR/experiment_32_threads.csv"

EXP="chm13v2.0 mouse zebrafish human_protein_5640 mouse_protein_589 zebrafish_protein_437"
EXP_32="zebrafish human_protein_5640 mouse_protein_589 zebrafish_protein_437"
echo "algo,test,k,num-threads,time,space" > $OUTPUT

# kiss-1
for j in 1 2 3; do
    for exp in $EXP; do
        for k in 2 4 8 16 32 64 128 256 -1; do
            file="$DATA_DIR/${exp}.fa"
            num_threads=32
            log="$LOG_DIR/kiss_1_${exp}_${k}_${num_threads}_${j}.txt"

            if [ $num_threads -eq 128 ]; then
                cpu_list="0-127"
            else
                end=$(( $num_threads + $num_threads - 2 ))
                cpu_list="0-$end:2"
            fi

            if [[ ${exp} == *"protein"* ]]; then
                is_general="1"
            else
                is_general="0"
            fi

            taskset -c $cpu_list ./experiment/kiss-1 $file $num_threads $is_general $k > $log 2>&1
            space=`tac $log | sed -n 1p | awk '{ print $7 }'`
            t=`tac $log | sed -n 2p | awk '{ print $7 }'`
            echo "kiss-1,${exp},${k},${num_threads},${t},${space}" >> $OUTPUT
        done 
    done
done

# kiss-2
for j in 1 2 3; do
    for exp in $EXP; do
        for k in 2 4 8 16 32 64 128 256 -1; do
            file="$DATA_DIR/${exp}.fa"
            num_threads=32
            log="$LOG_DIR/kiss_2_${exp}_${k}_${num_threads}_${j}.txt"

            if [ $num_threads -eq 128 ]; then
                cpu_list="0-127"
            else
                end=$(( $num_threads + $num_threads - 2 ))
                cpu_list="0-$end:2"
            fi

            if [[ ${exp} == *"protein"* ]]; then
                is_general="1"
            else
                is_general="0"
            fi

            taskset -c $cpu_list ./experiment/kiss-2 $file $num_threads $is_general $k > $log 2>&1
            space=`tac $log | sed -n 1p | awk '{ print $7 }'`
            t=`tac $log | sed -n 2p | awk '{ print $7 }'`
            echo "kiss-2,${exp},${k},${num_threads},${t},${space}" >> $OUTPUT
        done 
    done
done

# libsais-64
for j in 1 2 3; do
    for exp in $EXP; do
        num_threads=32
        file="$DATA_DIR/${exp}.fa"
        k=-1
        log="$LOG_DIR/libsais64_${exp}_${k}_${num_threads}_${j}.txt"

        if [ $num_threads -eq 128 ]; then
            cpu_list="0-127"
        else
            end=$(( $num_threads + $num_threads - 2 ))
            cpu_list="0-$end:2"
        fi

        if [[ ${exp} == *"protein"* ]]; then
            is_general="1"
        else
            is_general="0"
        fi

        taskset -c $cpu_list ./experiment/libsais-test-64 $file $num_threads $is_general > $log 2>&1
        space=`tac $log | sed -n 1p | awk '{ print $7 }'`
        t=`tac $log | sed -n 2p | awk '{ print $7 }'`
        echo "libsais-64,${exp},${k},${num_threads},${t},${space}" >> $OUTPUT
    done
done

# libsais-32
for j in 1 2 3; do
    for exp in $EXP_32; do
        num_threads=32
        file="$DATA_DIR/${exp}.fa"
        k=-1
        log="$LOG_DIR/libsais32_${exp}_${k}_${num_threads}_${j}.txt"

        if [ $num_threads -eq 128 ]; then
            cpu_list="0-127"
        else
            end=$(( $num_threads + $num_threads - 2 ))
            cpu_list="0-$end:2"
        fi

        if [[ ${exp} == *"protein"* ]]; then
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

# pardss-64
for j in 1 2 3; do
    for exp in $EXP; do
        num_threads=32
        file="$DATA_DIR/${exp}.fa"
        k=-1
        log="$LOG_DIR/pardss64_${exp}_${k}_${num_threads}_${j}.txt"

        if [ $num_threads -eq 128 ]; then
            cpu_list="0-127"
        else
            end=$(( $num_threads + $num_threads - 2 ))
            cpu_list="0-$end:2"
        fi

        is_general=0

        taskset -c $cpu_list ./experiment/pardss-64 $file $num_threads $is_general > $log 2>&1
        space=`tac $log | sed -n 1p | awk '{ print $7 }'`
        t=`tac $log | sed -n 2p | awk '{ print $7 }'`
        echo "pardss-64,${exp},${k},${num_threads},${t},${space}" >> $OUTPUT
    done
done

# pardss-32
for j in 1 2 3; do
    for exp in $EXP_32; do
        num_threads=32
        file="$DATA_DIR/${exp}.fa"
        k=-1
        log="$LOG_DIR/pardss32_${exp}_${k}_${num_threads}_${j}.txt"

        if [ $num_threads -eq 128 ]; then
            cpu_list="0-127"
        else
            end=$(( $num_threads + $num_threads - 2 ))
            cpu_list="0-$end:2"
        fi

        is_general=0

        taskset -c $cpu_list ./experiment/pardss-32 $file $num_threads $is_general > $log 2>&1
        space=`tac $log | sed -n 1p | awk '{ print $7 }'`
        t=`tac $log | sed -n 2p | awk '{ print $7 }'`
        echo "pardss-32,${exp},${k},${num_threads},${t},${space}" >> $OUTPUT
    done
done

# saca-k
for j in 1 2 3; do
    for exp in $EXP; do
        file="/home/dave/data/${exp}.fa"
        num_threads=1
        k=-1
        log="./sa/sacak_${exp}_${k}_${num_threads}_${j}.txt"

        cpu_list="0-0"

        is_general=1

        taskset -c $cpu_list ./experiment/saca-k $file $num_threads $is_general > $log 2>&1
        space=`tac $log | sed -n 1p | awk '{ print $7 }'`
        t=`tac $log | sed -n 2p | awk '{ print $7 }'`
        echo "saca-k,${exp},${k},${num_threads},${t},${space}" >> $OUTPUT
    done
done