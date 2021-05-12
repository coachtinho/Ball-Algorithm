#!/bin/bash

# USAGE: ./speedup.sh NUM_THREADS N_DIMS N_POINTS RANDOM_SEED

[[ $# -eq 4 ]] || exit 1;

NUM_THREADS=$1;

[[ $NUM_THREADS -gt 0 ]] || exit 1;
shift;

ARGS=$*;

SEQ_TIME=$(./ballAlg $ARGS 2>&1 > /dev/null);
echo "Seq time: $SEQ_TIME";

OMP_TIME=$(OMP_NUM_THREADS=$NUM_THREADS ./ballAlg-omp $ARGS 2>&1 > /dev/null);
echo "OpenMP time: $OMP_TIME";

OMP_SPEEDUP=$(echo "scale=6; $SEQ_TIME/$OMP_TIME" | bc);
echo "OpenMP speedup: $OMP_SPEEDUP";

MPI_TIME=$(mpirun --use-hwthread-cpus -n $NUM_THREADS ballAlg-mpi $ARGS 2>&1 > /dev/null | grep -E '^[0-9]+\.[0-9]+$');
echo "MPI time: $MPI_TIME";

MPI_SPEEDUP=$(echo "scale=6; $SEQ_TIME/$MPI_TIME" | bc);
echo "MPI speedup: $MPI_SPEEDUP";
