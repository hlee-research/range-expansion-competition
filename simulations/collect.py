"""Merge per-task summaries and compare measured invasion velocities with theory.
usage: python collect.py MODE [MODE ...]      (reads results/MODE/task_*.csv)
writes results/MODE/summary.csv and results/MODE/summary.png
"""
import glob, os, sys
import pandas as pd
import matplotlib.pyplot as plt

for mode in sys.argv[1:]:
    files = sorted(glob.glob(f"results/{mode}/task_[0-9][0-9][0-9][0-9].csv"))
    if not files:
        print(f"{mode}: no results yet"); continue
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True).sort_values("task")
    df["label"] = df["label"].fillna("")
    df.to_csv(f"results/{mode}/summary.csv", index=False)
    df["rel_err"] = df["u_late"] / df["u_theory"] - 1
    print(f"\n{mode}: {len(df)} tasks")
    cols = ["task", "label", "regime", "afs", "asf", "dx", "dt", "T", "u_theory", "u_late",
            "u_sqrt", "u2_theory", "u2_late", "vf_meas", "vf_theory", "t_valid_end", "rel_err"]
    with pd.option_context("display.width", 200, "display.max_rows", 500):
        print(df[cols].to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    if mode == "conv":                      # u / u_theory vs dt, one line per regime
        fig, ax = plt.subplots(figsize=(4.5, 3.5))
        for lab, g in df[df.label.isin(["dent", "circ", "comp"])].groupby("label"):
            ax.plot(g.dt, g.u_late / g.u_theory, "o-", label=lab)
        ax.axhline(1, color="k", lw=0.8); ax.set_xscale("log")
        ax.set_xlabel("dt"); ax.set_ylabel("u_measured / u_theory"); ax.legend()
    elif mode in ("fig2a", "fig2b", "fig2c", "fig3"):
        x = "afs" if mode in ("fig2a", "fig2c") else "asf"
        fig, ax = plt.subplots(figsize=(5, 3.5))
        for lab, g in df.groupby("label"):
            gx = g.sort_values("afs" if lab == "traj3" else x)
            xx = gx["afs"] if lab == "traj3" else gx[x]
            ax.plot(xx, gx.u_theory, "-", lw=1)
            ax.plot(xx, gx.u_late, "o", ms=3, label=lab or "u")
            if gx.u2_theory.notna().any():
                ax.plot(xx, gx.u2_theory, "--", lw=1)
                ax.plot(xx, gx.u2_late, "s", ms=3, mfc="none")
        ax.axhline(0, color="k", lw=0.5); ax.set_xlabel(x); ax.set_ylabel("u"); ax.legend(fontsize=7)
    else:
        continue
    fig.tight_layout(); fig.savefig(f"results/{mode}/summary.png", dpi=150)
    print(f"saved results/{mode}/summary.png")
