#!/bin/bash

#rm oar_*/*
property="core>='1'"
ressource="walltime=02:00:00"
stdout="oar_out/%jobid%.out"
stderr="oar_err/%jobid%.err"

for sT in $(seq 0.3 0.15 1); do
    sTdir=D3_sT_$sT 
    mkdir $sTdir
    
    for P in $(seq 0.25 0.3 2); do
        for sx in $(seq 0.04 0.06 0.34); do
            program="python TCSimu.py -P $P -sx $sx -dir $sTdir -sT $sT"
            oarsub -p "$property" -l "$ressource" --stdout="$stdout" --stderr="$stderr" "$program"
        done
    done
done