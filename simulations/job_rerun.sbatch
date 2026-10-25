#!/bin/bash
# One array task = one simulation. Submitted by submit_rerun.sh, which sets --array, -t and --mem.
#SBATCH -J rerun
#SBATCH -p mit_normal
#SBATCH -c 1
#SBATCH -o logs/rerun-%A-%a.out

MODE=${1:?give MODE}
source "$SLURM_SUBMIT_DIR/julia_env.sh"

# one thread per task; 2D loops are single-threaded, so extra threads would only oversubscribe
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OUT_DIR=$SLURM_SUBMIT_DIR/results

cd "$SLURM_SUBMIT_DIR"
echo "task $SLURM_ARRAY_TASK_ID  mode=$MODE  node=$(hostname)  julia=$(command -v julia)"
julia --startup-file=no rerun.jl "$MODE" "$SLURM_ARRAY_TASK_ID"
