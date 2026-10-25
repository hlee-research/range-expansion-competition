#!/bin/bash
# Submit all tasks of one MODE as a Slurm job array. Run on the login node, from this directory.
#   ./submit_rerun.sh MODE [extra sbatch options]
#   MODE: fig2a | fig2b | fig2c | fig3 | conv | ic | strip
#   e.g. ./submit_rerun.sh fig2b
#        ./submit_rerun.sh conv --mail-type=END --mail-user=you@mit.edu
set -euo pipefail
MODE=${1:?give MODE (fig2a|fig2b|fig2c|fig3|conv|ic|strip)}; shift
source ./julia_env.sh

N=$(julia --startup-file=no rerun.jl count "$MODE")
case $MODE in
    conv) TIME=06:00:00; MEM=8G ;;   # includes dx = 0.5 and T = 1000 runs
    *)    TIME=02:00:00; MEM=4G ;;
esac
mkdir -p logs results
echo "submitting $N tasks for $MODE (time $TIME, mem $MEM)"
sbatch --array=0-$((N-1)) -t "$TIME" --mem="$MEM" -J "rerun-$MODE" "$@" job_rerun.sbatch "$MODE"
