#!/bin/bash

SEQ="./ballAlg"
OMP="./ballAlg-omp"
THREADS=(1 2 4 8)

for f in $(ls tests/*.in); do
    echo "============="
    echo "File: $f"
    ARGS=$(cat $f)
    
    echo -n "Seq: "
    seq_time=$($SEQ $ARGS 2>&1 > /dev/null)

    echo $seq_time

    for t in ${THREADS[@]}; do
        echo -n "$t threads: "
        par_time=$(OMP_NUM_THREADS=$t $OMP $ARGS 2>&1 > /dev/null)

        echo $par_time
    done
done
