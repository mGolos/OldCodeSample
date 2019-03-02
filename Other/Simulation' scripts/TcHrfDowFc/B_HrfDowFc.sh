#!/bin/bash

rm oar_*/*
nodes=(  01 02 03 04 05 06 07 08 09 10 11 12 17 18 19 20 13 14 15 16 21 22 23 24)
nbcores=( 8  8  8  8  8  8  8  8  8  8  8  8  8  8  8  8 12 12 12 12 16 16 16 16)

# Start from this point
strtnode=8
node=$strtnode

ncpu=1
ncore=2
property="core<='240'"
ressource="/cpu=$ncpu/core=$ncore,walltime=15:00:00"
stdout="oar_out/%jobid%.out"
stderr="oar_err/%jobid%.err"
wait=0
dt=200

rcore=1
for p1 in $(seq 0.05 0.05 0.50); do
    for p2 in $(seq 0.05 0.05 0.45); do
        program="python B_HrfDowFc.py -p1 $p1 -p2 $p2 -w $wait"
        property="host='n${nodes[$(($node-1))]}'"
        # echo oarsub -p "$property" -l "$ressource"
        oarsub -p "$property" -l "$ressource" --stdout="$stdout" --stderr="$stderr" "$program"

        if [ $node -lt ${nodes[23]} ]; then
            node=$(( $node + 1 ))
        else
            rcore=$(( $rcore + 1 ))
            wait=$(( $wait + $dt ))
            node=$strtnode
            while [ $rcore -gt $(( ${nbcores[$(( $node-1 ))]} / $ncore )) ]; do
                node=$(( $node + 1 ))
                if [ $node -gt ${nodes[23]} ]; then
                    echo "problem"
                fi
            done
        fi

#         if [ $rcore -lt $((${nbcores[$(($node-1))]} / $ncore)) ]; then
#             rcore=$(($rcore + 1))
#         else
#             rcore=1
#             node=$(($node + 1))
#         fi
    done
done

# rm oar_*/*
# property="core<='240'"
# ressource="/cpu=1/core=2,walltime=15:00:00"
# stdout="oar_out/%jobid%.out"
# stderr="oar_err/%jobid%.err"
# for p1 in $(seq 0.05 0.05 0.50); do
#     for p2 in $(seq 0.05 0.05 0.45); do
#         program="python B_HrfDowFc.py -p1 $p1 -p2 $p2"
#         oarsub -p "$property" -l "$ressource" --stdout="$stdout" --stderr="$stderr" "$program"
#     done
# done