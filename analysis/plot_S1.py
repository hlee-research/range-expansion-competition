"""S1 Fig: numerical convergence of invasion velocities.
usage (from the repository root):  python analysis/plot_S1.py
reads data/conv/summary.csv (panel A) and data/conv/task_00{05,11,17}_traj.csv (panel B)
writes figures/S1.pdf
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

names = {"circ": "Circular bulge", "comp": "Composite bulge", "dent": "Dent"}
colors = {"circ": "C0", "comp": "C1", "dent": "C2"}
long_runs = {"dent": 5, "circ": 11, "comp": 17}          # T = 1000 runs in the conv task table

c = pd.read_csv("data/conv/summary.csv")
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7, 2.8))

for lab, g in c[c.label.isin(names)].sort_values("dt").groupby("label"):
    axA.plot(g.dt, g.u_late / g.u_theory, "o-", color=colors[lab], label=names[lab])
axA.axhline(1, ls="--", lw=0.8, color="k")
axA.set_xlabel(r"Time step, $\Delta t$")
axA.set_ylabel(r"$u_{\mathrm{sim}}/u_{\mathrm{theory}}$")
axA.legend(frameon=False, fontsize=7)

w = 20                                                    # half-window, in time units (records every 1.0)
for lab, task in long_runs.items():
    path = f"data/conv/task_{task:04d}_traj.csv"
    if not os.path.exists(path):
        print(f"missing {path}: panel B will lack the {names[lab]} curve")
        continue
    uth = abs(c.loc[c.task == task, "u_theory"].item())
    tr = pd.read_csv(path).query("valid_99 == 1")
    x, t = tr.xL_99.to_numpy(), tr.t.to_numpy()
    u_t = np.abs(x[2 * w:] - x[:-2 * w]) / (t[2 * w:] - t[:-2 * w])
    axB.plot(t[w:-w], u_t / uth, color=colors[lab], lw=1)
axB.axhline(1, ls="--", lw=0.8, color="k")
axB.set_xlabel(r"Time, $t$")
axB.set_ylabel(r"$|u(t)|/u_{\mathrm{theory}}$")

for ax, l in zip((axA, axB), "AB"):
    ax.text(-0.18, 1.02, l, transform=ax.transAxes, fontweight="bold")
fig.tight_layout()
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/S1.pdf")
print("saved figures/S1.pdf")
