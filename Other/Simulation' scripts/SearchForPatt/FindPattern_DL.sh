#!/bin/bash

rm oar_*/*
property="core>'42'"
ressource="walltime=04:00:00"
stdout="oar_out/%jobid%.out"
stderr="oar_err/%jobid%.err"

for p1 in $(seq 4 0.5 16); do
    for p2 in $(seq 0.02 0.03 1); do
        program="python FindPattern_DL.py -p1 $p1 -p2 $p2"
        oarsub -p "$property" -l "$ressource" --stdout="$stdout" --stderr="$stderr" "$program"
    done
done