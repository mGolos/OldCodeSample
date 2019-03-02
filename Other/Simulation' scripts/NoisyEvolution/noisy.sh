#!/bin/bash

rm oar_*/*
property="core>='42'"
ressource="walltime=01:00:00"
stdout="oar_out/%jobid%.out"
stderr="oar_err/%jobid%.err"

for p1 in $(seq 0 1 9); do
    program="python noisy$1.py -p1 $p1"
    oarsub -p "$property" -l "$ressource" --stdout="$stdout" --stderr="$stderr" "$program"
done

#for p1 in $(seq 0 8 100); do
#    for p2 in $(seq 0 1 6); do
#        if [ $(( $p1 + $p2 )) -lt 100 ]; then
#            python noisy$1.py -p1 $(( $p1 + $p2 )) &
#        fi
#    done
#   if [ $(( $p1 + 7 )) -lt 100 ]; then
#        python noisy$1.py -p1 $(( $p1 + 7 ))
#    fi
#done
