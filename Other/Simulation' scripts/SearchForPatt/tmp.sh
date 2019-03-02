#!/bin/bash

#rm oar_*/*
property="core>='209'"
ressource="walltime=04:00:00"
stdout="oar_out/%jobid%.out"
stderr="oar_err/%jobid%.err"

# for p1 in $(seq 0.7 0.02 1.3); do
#     for p2 in $(seq 0.02 0.03 1); do
for p1 in 1.00; do
    for p2 in 10 20 30 50; do
        program="python PreClustering.py -p1 $p1 -p2 $p2 -p3 OtherConn_103414"
        oarsub -p "$property" -l "$ressource" --stdout="$stdout" --stderr="$stderr" "$program"
    done
done
