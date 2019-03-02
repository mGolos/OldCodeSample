#!/bin/bash

rm oar_*/*
property="core<='240'"
ressource="walltime=04:00:00"
stdout="oar_out/%jobid%.out"
stderr="oar_err/%jobid%.err"

for p1 in $(seq 0.05 0.05 0.50); do
    for p2 in $(seq 0.05 0.05 0.45); do
        program="python A_Tc.py -p1 $p1 -p2 $p2"
        oarsub -p "$property" -l "$ressource" --stdout="$stdout" --stderr="$stderr" "$program"
    done
done