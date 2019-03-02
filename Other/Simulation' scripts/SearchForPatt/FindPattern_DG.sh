#!/bin/bash

rm oar_*/*
property="core>='13'"
ressource="walltime=06:00:00"
stdout="oar_out/%jobid%.out"
stderr="oar_err/%jobid%.err"

# for p1 in $(seq 0.25 0.25 8); do
# for p1 in 0.60; do
#      for p2 in $(seq 0.02 0.03 1); do
for p1 in 1.00; do
   for p2 in $(seq 10 10 100); do
        program="python FindPattern_DG.py -p1 $p1 -p2 $p2 -dir $1 -con $2"
        oarsub -p "$property" -l "$ressource" --stdout="$stdout" --stderr="$stderr" "$program"
    done
done
