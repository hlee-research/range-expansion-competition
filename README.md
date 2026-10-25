# Velocity differences have distinct effects on exclusion, coexistence, and bistability during range expansions

Code and data for

> H. Lee, D. W. Swartz, and K. S. Korolev, *Velocity differences have distinct effects on exclusion, coexistence, and bistability during range expansions* (2026). <!-- TODO: journal, DOI -->

The repository contains the reaction-diffusion solver for two interacting types in a planar range expansion (Eq 4 of the paper), the scripts used to run the parameter sweeps on a Slurm cluster, the summarized simulation results, and the analysis code for the figures.

## Layout

```
simulations/
  two_dimensional_module_interaction.jl   reference solver (Eq 4): explicit Euler, operator splitting, 5-point Laplacian
  rerun.jl                                 driver used for Figs 2, 3C and S1 Fig: fused kernel, on-the-fly front tracking,
                                           velocity fits; `julia rerun.jl selftest` checks it against the reference solver
  submit_rerun.sh, job_rerun.sbatch,       Slurm job-array submission (one task = one simulation)
  julia_env.sh
  collect.py                               merges per-task results into results/MODE/summary.csv
analysis/
  theory.py                                Eqs (7), (8), (9) and (11): theory curves for Figs 2 and 3
  plot_S1.py                               S1 Fig (numerical convergence)
  fig5_community.ipynb                     Fig 5 (community-level competition under the trade-off)
data/
  MODE/summary.csv                         one row per simulation (see data/README.md for the columns)
figures/                                   output directory for the plotting scripts
```

## Which data and code produce which figure

| Figure | Data | Code |
|---|---|---|
| Fig 1 | phase boundaries are analytical (Eq 8); insets from earlier simulations with D = 0.1 (S1 Text, Section S4.1) | `analysis/F1.ipynb` |
| Fig 2A, 2B, 2C | `data/fig2a`, `data/fig2b`, `data/fig2c` | `analysis/theory.py`, `analysis/F2.ipynb` |
| Fig 3B | earlier simulations with D = 0.1 | `data/fig3b` |
| Fig 3C | `data/fig3` (labels `traj1`, `traj2`, `traj3`) | `analysis/theory.py`,`analysis/F3.ipynb` |
| Fig 4 | sampled from the trade-off, no simulation data | `analysis/F4.ipynb` |
| Fig 5 | computed in the notebook | `analysis/F5.ipynb` |
| S1 Fig | `data/conv` | `analysis/plot_S1.py` |
| S1 Text, Section S4.5 | `data/ic`, `data/strip` | — |

## Requirements

- Julia ≥ 1.6. `rerun.jl` uses only the standard library. The reference module (and therefore `julia rerun.jl selftest`) needs `Distributions.jl`; install it with `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
- Python ≥ 3.9 with the packages in `requirements.txt` (`pip install -r requirements.txt`).

## Reproducing the simulations

From `simulations/`:

```bash
julia --project=.. rerun.jl selftest          # kernel vs. reference solver; should print PASS
julia rerun.jl list fig2b                      # parameter table of a sweep
julia rerun.jl fig2b 60                        # run one task locally (a few minutes)
./submit_rerun.sh fig2b                        # submit the whole sweep as a Slurm job array
python collect.py fig2b                        # merge results into results/fig2b/summary.csv
```

Sweeps: `fig2a`, `fig2b`, `fig2c`, `fig3`, `conv` (S1 Fig), `ic` and `strip` (S1 Text, Section S4.5). The Slurm scripts were written for the MIT Engaging cluster (partition `mit_normal`); change the partition in `job_rerun.sbatch` and the module name in `julia_env.sh` for other clusters.

All simulations in Figs 2 and 3C use r_f = 1, r_s = 0.8, D = 1, Δx = 1, Δt = 0.025, a 1600 × 950 lattice and T = 400, in units in which lengths are measured in √(D/r_f) and times in 1/r_f. Velocities are fitted over 200 ≤ t ≤ 400.

## Reproducing the figures

From the repository root:

```bash
python analysis/plot_S1.py                     # writes figures/S1.pdf
jupyter notebook analysis/fig5_community.ipynb # Fig 5
```

## License

MIT; see `LICENSE`.

## Citation

See `CITATION.cff`. <!-- TODO: add the Zenodo DOI badge after the first release -->
