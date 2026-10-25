#!/usr/bin/env python3
"""
N5_aggregate_nested_runs.py

Purpose:
    Aggregate held-out AUC/BACC across multiple independent nested-validation
    runs (different random train/held-out splits) to check whether the
    leakage-free classification result is stable, or was a lucky/unlucky
    artefact of a single ~20-sample held-out split.

Outputs:
    results/nested_validation/nested_validation_summary.csv
    results/nested_validation/nested_validation_summary.png
"""

import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_DIR = Path.cwd()
NESTED_DIR = PROJECT_DIR / "results/nested_validation"

run_dirs = sorted([d for d in NESTED_DIR.iterdir() if d.is_dir() and d.name.startswith("run")])
print(f"Found {len(run_dirs)} runs: {[d.name for d in run_dirs]}")

rows = []
for run_dir in run_dirs:
    f = run_dir / "tables" / "nested_binary_holdout_metrics.csv"
    if not f.exists():
        print(f"  [SKIP] {run_dir.name}: no binary holdout metrics found")
        continue
    df = pd.read_csv(f)
    df["run"] = run_dir.name
    rows.append(df)

if not rows:
    raise SystemExit("No completed runs found — nothing to aggregate yet.")

all_binary = pd.concat(rows, ignore_index=True)
all_binary.to_csv(NESTED_DIR / "nested_validation_all_runs_binary.csv", index=False)

print("\n" + "=" * 60)
print("Held-out AUC by run and model")
print("=" * 60)
pivot = all_binary.pivot(index="run", columns="model", values="holdout_AUC")
print(pivot.round(3))

summary = all_binary.groupby("model")["holdout_AUC"].agg(["mean", "std", "min", "max", "count"]).round(3)
summary_bacc = all_binary.groupby("model")["holdout_BACC"].agg(["mean", "std", "min", "max"]).round(3)
print("\n" + "=" * 60)
print("Summary across runs — held-out AUC")
print("=" * 60)
print(summary)
print("\n" + "=" * 60)
print("Summary across runs — held-out Balanced Accuracy")
print("=" * 60)
print(summary_bacc)

summary_out = summary.reset_index().rename(columns={
    "mean": "AUC_mean_across_runs", "std": "AUC_std_across_runs",
    "min": "AUC_min", "max": "AUC_max", "count": "n_runs"
})
summary_out.to_csv(NESTED_DIR / "nested_validation_summary.csv", index=False)

# Multiclass
multi_rows = []
for run_dir in run_dirs:
    f = run_dir / "tables" / "nested_multiclass_holdout_metrics.csv"
    if f.exists():
        df = pd.read_csv(f)
        df["run"] = run_dir.name
        multi_rows.append(df)
if multi_rows:
    all_multi = pd.concat(multi_rows, ignore_index=True)
    all_multi.to_csv(NESTED_DIR / "nested_validation_all_runs_multiclass.csv", index=False)
    print("\n" + "=" * 60)
    print("Held-out 3-class BACC across runs")
    print("=" * 60)
    print(all_multi[["run", "holdout_BACC"]].round(3))
    print(f"\nMean: {all_multi['holdout_BACC'].mean():.3f}  SD: {all_multi['holdout_BACC'].std():.3f}")

# Plot: AUC by run and model
fig, ax = plt.subplots(figsize=(8, 5))
colors = {"RandomForest": "#8B0000", "LogisticRegression": "#012345", "SVM_RBF": "#2D6A4F"}
for model in all_binary["model"].unique():
    sub = all_binary[all_binary["model"] == model].sort_values("run")
    ax.plot(sub["run"], sub["holdout_AUC"], marker="o", label=model.replace("_", " "),
             color=colors.get(model, "grey"), linewidth=1.5)
ax.axhline(0.5, color="grey", linestyle="--", linewidth=0.8, label="Chance")
ax.set_ylabel("Held-out AUC (leakage-free)")
ax.set_xlabel("Run (independent random train/held-out split)")
ax.set_title("Stability of leakage-free held-out AUC across independent splits\n(T2D vs ND, hub genes selected from TRAIN samples only)")
ax.set_ylim(0, 1.05)
ax.legend()
fig.tight_layout()
fig.savefig(NESTED_DIR / "nested_validation_summary.png", dpi=300)
plt.close(fig)

print("\n" + "=" * 60)
print("N5 COMPLETE")
print("=" * 60)
print(f"Outputs: {NESTED_DIR / 'nested_validation_summary.csv'}")
print(f"         {NESTED_DIR / 'nested_validation_summary.png'}")
