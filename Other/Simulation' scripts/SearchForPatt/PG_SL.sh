#!/bin/bash

#rm oar_*/*
property="core>='209'"
ressource="walltime=100:00:00"
stdout="oar_out/%jobid%.out"
stderr="oar_err/%jobid%.err"

for G in $(seq 10 10 150); do
    for revert in 0 1; do
        program="python PG_SL.py -G $G -revert $revert -dir $1 -con $2"
        oarsub -p "$property" -l "$ressource" --stdout="$stdout" --stderr="$stderr" "$program"
    done
done
