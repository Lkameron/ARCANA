#!/usr/bin/env python3
"""
N4_nested_ml_validation.py

Purpose:
    Leakage-free re-evaluation of the hub-gene classification claim from
    script 08. The original pipeline selected hub-gene features using
    module-trait correlation and gene significance computed across ALL 98
    samples (including the labels of every 5-fold CV "held-out" fold), which
    inflates the reported AUC. This script:
      1. Uses hub genes selected ONLY from the 78 TRAIN samples (N1 split,
         N0 train-only DE, N2 train-only WGCNA, N3 train-only hub scoring —
         none of which ever saw the 20 held-out samples' labels).
      2. Fits classifiers on the 78 TRAIN samples only.
      3. Evaluates ONCE on the 20 held-out TEST samples, which had zero
         influence on feature selection, network construction, or model
         fitting at any stage.

    Expression VALUES (VST-normalised, blind=TRUE) are reused from the full
    98-sample matrix — VST-blind is unsupervised (does not use status labels)
    so this does not reintroduce leakage; only gene IDENTITY selection needed
    to be train-only.

Outputs:
    results/nested_validation/tables/nested_binary_holdout_metrics.csv
    results/nested_validation/tables/nested_multiclass_holdout_metrics.csv
    results/nested_validation/plots/nested_holdout_ROC_curves.png
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score, roc_curve, balanced_accuracy_score

warnings.filterwarnings("ignore")

PROJECT_DIR = Path.cwd()
RUN_ID = sys.argv[1] if len(sys.argv) > 1 else "run1"

# Train-only hub genes (leakage-free feature selection)
HUB_PATH_TRAIN = PROJECT_DIR / f"results/nested_validation/{RUN_ID}/wgcna_train/hubgenes/tables/wgcna_top20_hub_genes_per_focus_module.csv"

# Full-dataset VST expression (unsupervised transform, safe to reuse) + labels
EXPR_PATH   = PROJECT_DIR / "results/wgcna/tables/wgcna_expression_input.tsv"
TRAITS_PATH = PROJECT_DIR / "results/wgcna/tables/wgcna_traits_input.csv"

SPLIT_PATH  = PROJECT_DIR / f"results/nested_validation/{RUN_ID}/tables/train_test_split.csv"

OUT_DIR    = PROJECT_DIR / f"results/nested_validation/{RUN_ID}"
OUT_TABLES = OUT_DIR / "tables"
OUT_PLOTS  = OUT_DIR / "plots"
OUT_TABLES.mkdir(parents=True, exist_ok=True)
OUT_PLOTS.mkdir(parents=True, exist_ok=True)

STATUS_COLORS = {"ND": "#012345", "IGT": "#345678", "T2D": "#8B0000"}

print("=" * 60)
print(f"N4 - Nested (leakage-free) held-out validation — {RUN_ID}")
print("=" * 60)

expr_df   = pd.read_csv(EXPR_PATH, sep="\t", index_col="sampleid")
traits_df = pd.read_csv(TRAITS_PATH, index_col="sampleid")
hub_df    = pd.read_csv(HUB_PATH_TRAIN)
split_df  = pd.read_csv(SPLIT_PATH).set_index("sampleid")

common = expr_df.index.intersection(traits_df.index)
expr_df, traits_df = expr_df.loc[common], traits_df.loc[common]

train_ids = split_df.index[split_df["split"] == "train"]
test_ids  = split_df.index[split_df["split"] == "test"]
print(f"Train samples: {len(train_ids)}  |  Held-out test samples: {len(test_ids)}")

# Feature set: top 10 hub genes per module, selected using TRAIN-only ranking
hub_genes = (
    hub_df.sort_values("hub_rank_all_genes")
    .groupby("module_color")
    .head(10)["gene_id"]
    .unique()
    .tolist()
)
hub_genes = [g for g in hub_genes if g in expr_df.columns]
print(f"Hub gene features (train-derived): {len(hub_genes)}")

label_map = {}
symbol_col = next((c for c in ["gene_symbol", "gene_symbol.x", "gene_symbol.y"] if c in hub_df.columns), None)
if symbol_col:
    sym = hub_df.dropna(subset=[symbol_col]).drop_duplicates("gene_id")
    label_map = dict(zip(sym["gene_id"], sym[symbol_col]))

def readable(g):
    s = label_map.get(g, g)
    return s if s and s != g else g

X_full = expr_df[hub_genes].copy().astype(float)
X_full.columns = [readable(g) for g in hub_genes]
feature_names = X_full.columns.tolist()

# ──────────────────────────────────────────────
# Binary: T2D vs ND, single held-out evaluation
# ──────────────────────────────────────────────
print("\n" + "-" * 40)
print("Binary held-out evaluation: T2D vs ND")
print("-" * 40)

mask_bin = traits_df["status"].isin(["T2D", "ND"])
bin_ids = traits_df.index[mask_bin]

train_bin = [s for s in bin_ids if s in train_ids]
test_bin  = [s for s in bin_ids if s in test_ids]
print(f"  Train (T2D/ND): {len(train_bin)}  |  Held-out test (T2D/ND): {len(test_bin)}")
print(f"  Held-out composition: {traits_df.loc[test_bin, 'status'].value_counts().to_dict()}")

X_train = X_full.loc[train_bin].values
y_train = (traits_df.loc[train_bin, "status"] == "T2D").astype(int).values
X_test  = X_full.loc[test_bin].values
y_test  = (traits_df.loc[test_bin, "status"] == "T2D").astype(int).values

models = {
    "RandomForest": Pipeline([
        ("scaler", StandardScaler()),
        ("clf", RandomForestClassifier(n_estimators=300, max_features="sqrt",
                                        class_weight="balanced", random_state=42, n_jobs=-1)),
    ]),
    "LogisticRegression": Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(C=1.0, max_iter=2000, class_weight="balanced",
                                    solver="lbfgs", random_state=42)),
    ]),
    "SVM_RBF": Pipeline([
        ("scaler", StandardScaler()),
        ("clf", SVC(kernel="rbf", C=1.0, probability=True, class_weight="balanced", random_state=42)),
    ]),
}

results = {}
roc_curves = {}
for name, model in models.items():
    model.fit(X_train, y_train)
    y_prob = model.predict_proba(X_test)[:, 1]
    y_pred = model.predict(X_test)
    auc = roc_auc_score(y_test, y_prob) if len(set(y_test)) > 1 else float("nan")
    bacc = balanced_accuracy_score(y_test, y_pred)
    results[name] = {"holdout_AUC": auc, "holdout_BACC": bacc, "n_test": len(test_bin)}
    roc_curves[name] = roc_curve(y_test, y_prob) if len(set(y_test)) > 1 else None
    print(f"  {name:22s}  held-out AUC={auc:.3f}  BACC={bacc:.3f}  (n_test={len(test_bin)})")

metrics_df = pd.DataFrame(results).T.reset_index().rename(columns={"index": "model"})
metrics_df.to_csv(OUT_TABLES / "nested_binary_holdout_metrics.csv", index=False)

print(
    f"\n  NOTE: with only {len(test_bin)} held-out binary samples "
    f"({(y_test==1).sum()} T2D / {(y_test==0).sum()} ND), AUC resolution is "
    f"coarse (steps of ~1/{(y_test==1).sum()*(y_test==0).sum()}) and confidence "
    f"intervals are wide. This single-split estimate is the leakage-free "
    f"answer, not a low-variance one — report it with that caveat."
)

# ROC plot
fig, ax = plt.subplots(figsize=(7, 6))
ax.plot([0, 1], [0, 1], "k--", lw=0.8, label="Random (AUC = 0.50)")
colors = {"RandomForest": "#8B0000", "LogisticRegression": "#012345", "SVM_RBF": "#2D6A4F"}
for name, rc in roc_curves.items():
    if rc is None:
        continue
    fpr, tpr, _ = rc
    ax.plot(fpr, tpr, lw=2, color=colors.get(name, "grey"),
             label=f"{name.replace('_', ' ')} (held-out AUC = {results[name]['holdout_AUC']:.3f})")
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title(f"Nested held-out ROC — T2D vs ND (n_test={len(test_bin)})\nLeakage-free: features + model both derived from TRAIN only")
ax.legend(loc="lower right", framealpha=0.9)
fig.tight_layout()
fig.savefig(OUT_PLOTS / "nested_holdout_ROC_curves.png", dpi=300)
plt.close(fig)

# ──────────────────────────────────────────────
# 3-class: ND / IGT / T2D, single held-out evaluation
# ──────────────────────────────────────────────
print("\n" + "-" * 40)
print("3-class held-out evaluation: ND / IGT / T2D")
print("-" * 40)

le = LabelEncoder()
le.fit(traits_df["status"])
y_train_multi = le.transform(traits_df.loc[train_ids, "status"])
y_test_multi  = le.transform(traits_df.loc[test_ids, "status"])
X_train_multi = X_full.loc[train_ids].values
X_test_multi  = X_full.loc[test_ids].values

multi_model = Pipeline([
    ("scaler", StandardScaler()),
    ("clf", RandomForestClassifier(n_estimators=300, max_features="sqrt",
                                    class_weight="balanced", random_state=42, n_jobs=-1)),
])
multi_model.fit(X_train_multi, y_train_multi)
y_pred_multi = multi_model.predict(X_test_multi)
bacc_multi = balanced_accuracy_score(y_test_multi, y_pred_multi)
print(f"  Held-out 3-class BACC = {bacc_multi:.3f}  (n_test={len(test_ids)})")

pd.DataFrame({
    "model": ["RandomForest_3class"],
    "holdout_BACC": [bacc_multi],
    "n_test": [len(test_ids)],
}).to_csv(OUT_TABLES / "nested_multiclass_holdout_metrics.csv", index=False)

print("\n" + "=" * 60)
print("N4 COMPLETE")
print("=" * 60)
print("Outputs written to results/nested_validation/")
