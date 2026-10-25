# Simulation summaries

Each `MODE/summary.csv` has one row per simulation. Velocities are in simulation units (lengths in √(D/r_f), times in 1/r_f) and are **not** normalized.

| Column | Meaning |
|---|---|
| `mode`, `task`, `label` | sweep name, task index, curve or trajectory label |
| `regime` | predicted regime: `dent`, `circ` (circular bulge), `comp` (composite bulge), `coex` |
| `rf`, `rs`, `afs`, `asf`, `D` | model parameters (Δr = rf − rs) |
| `dx`, `dt`, `T`, `Lx`, `Ly` | lattice spacing, time step, run length, domain size |
| `init`, `h0`, `offset`, `width` | initial condition (`half`: fast type left, slow type right; `strip`: fast strip of width `width`) |
| `u_theory`, `u2_theory` | Eq (8) for exclusion; for coexistence `u` = u_ms (mixture into slow) and `u2` = u_mf (mixture into fast) |
| `vf_theory`, `vs_theory` | 2√(D r) for the pure types |
| `t_valid_end` | last time at which the front and boundaries were away from the domain edges |
| `u_late`, `u2_late` | linear fit of the boundary position over the second half of the valid window (**used in the figures**) |
| `u_mid`, `u2_mid` | linear fit over the second quarter of the valid window |
| `u_sqrt`, `u2_sqrt` | fit x = x0 + u t + c√t over the last 80% of the valid window (unreliable near transitions) |
| `vf_meas`, `vs_meas` | measured expansion velocities of the pure types |
| `u_late_thr50`, `u2_late_thr50` | as `u_late`, with the front defined by n > 0.5 instead of n > 0.99 |
| `final_width` | width of the region occupied by the fast type at the end (`strip` sweep) |
| `runtime_s` | wall-clock time of the run |

Sign convention: u > 0 when the fast type advances into the slow type; u < 0 when the slow type advances into the fast type. For `strip`, `u_late` is the speed of each strip boundary.

`data/conv/task_0005_traj.csv`, `task_0011_traj.csv` and `task_0017_traj.csv` (boundary trajectories of the T = 1000 runs) are needed for panel B of S1 Fig.
