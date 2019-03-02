#!/bin/bash

#rm oar_*/*
property="core>='1'"
ressource="walltime=03:00:00"
stdout="oar_out/%jobid%.out"
stderr="oar_err/%jobid%.err"

for G in $(seq 10 1 14); do
    for sx in 0.001 0.0032 0.01 0.0316 0.1 0.3162; do
        for sT in 0.001 0.0032 0.01 0.0316 0.1 0.3162; do    
            program="python TCSimu_SL.py -G $G -sx $sx -sT $sT -dir $1"
            oarsub -p "$property" -l "$ressource" --stdout="$stdout" --stderr="$stderr" "$program"
        done
    done
done