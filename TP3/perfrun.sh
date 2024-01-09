#!/bin/bash
#SBATCH --ntasks=40
#SBATCH --time=00:10:00
#SBATCH --partition=cpar
#SBATCH --exclusive

threads=(1 4 8 12 14 18 22 26 28 32 36 40)

for nthreads in "${threads[@]}"
do
    export OMP_NUM_THREADS=${nthreads}
    echo "Running with ${OMP_NUM_THREADS} threads"
    perf record -o perf.data -e cpu-clock ./MDpar.exe <inputdata.txt >lixo
    perf report -i perf.data --stdio
done