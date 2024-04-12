#!/bin/bash
cd ../build && make -j
OUTPUT="/home/dave/kISS/experiment/sa_20240407.csv"
# EXP="example example_protein"
EXP="chm13v2.0 mouse zebrafish"
echo "algo,test,k,num-threads,time,space" > $OUTPUT

# kiss-1
for j in 1 2 3; do
    for exp in $EXP; do
        for k in 2 4 8 16 32 64 128 256 -1; do
            file="/home/dave/data/${exp}.fa"
            num_threads=32
            log="./sa/log_kiss_1_${exp}_${k}_${num_threads}_${j}.txt"

            if [ $num_threads -eq 128 ]; then
                cpu_list="0-127"
            else
                end=$(( $num_threads + $num_threads - 2 ))
                cpu_list="0-$end:2"
            fi

            if [ ${exp} = "example_protein" ] || [ ${exp} = "protein" ]; then
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
            file="/home/dave/data/${exp}.fa"
            num_threads=32
            log="./sa/log_kiss_2_${exp}_${k}_${num_threads}_${j}.txt"

            if [ $num_threads -eq 128 ]; then
                cpu_list="0-127"
            else
                end=$(( $num_threads + $num_threads - 2 ))
                cpu_list="0-$end:2"
            fi

            if [ ${exp} = "example_protein" ] || [ ${exp} = "protein" ]; then
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

# kiss-1 different threads
for j in 1 2 3; do
    for exp in $EXP; do
        for num_threads in 1 2 4 8 16 32 64 128; do
            file="/home/dave/data/${exp}.fa"
            k=256
            log="./sa/log_kiss_1_${exp}_${k}_${num_threads}_${j}.txt"

            if [ $num_threads -eq 128 ]; then
                cpu_list="0-127"
            else
                end=$(( $num_threads + $num_threads - 2 ))
                cpu_list="0-$end:2"
            fi

            if [ ${exp} = "example_protein" ] || [ ${exp} = "protein" ]; then
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

# kiss-2 different threads
for j in 1 2 3; do
    for exp in $EXP; do
        for num_threads in 1 2 4 8 16 32 64 128; do
            file="/home/dave/data/${exp}.fa"
            k=256
            log="./sa/log_kiss_2_${exp}_${k}_${num_threads}_${j}.txt"

            if [ $num_threads -eq 128 ]; then
                cpu_list="0-127"
            else
                end=$(( $num_threads + $num_threads - 2 ))
                cpu_list="0-$end:2"
            fi

            if [ ${exp} = "example_protein" ] || [ ${exp} = "protein" ]; then
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

# libsais
for j in 1 2 3; do
    for exp in $EXP; do
        num_threads=32
        file="/home/dave/data/${exp}.fa"
        k=-1
        log="./sa/libsais_${exp}_${k}_${num_threads}_${j}.txt"

        if [ $num_threads -eq 128 ]; then
            cpu_list="0-127"
        else
            end=$(( $num_threads + $num_threads - 2 ))
            cpu_list="0-$end:2"
        fi

        if [ ${exp} = "example_protein" ] || [ ${exp} = "protein" ]; then
            is_general="1"
        else
            is_general="0"
        fi

        taskset -c $cpu_list ./experiment/libsais-test $file $num_threads $is_general > $log 2>&1
        space=`tac $log | sed -n 1p | awk '{ print $7 }'`
        t=`tac $log | sed -n 2p | awk '{ print $7 }'`
        echo "libsais,${exp},${k},${num_threads},${t},${space}" >> $OUTPUT
    done
done