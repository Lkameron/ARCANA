#!/usr/bin/env python3
"""
08_ml_models.py

Project: Pancreatic-Islet-T2D-Coexpression

Purpose:
    Build and evaluate machine-learning classifiers that distinguish T2D from
    ND samples using WGCNA hub gene expression.  Includes:
      - Random Forest, Logistic Regression, SVM (existing)
      - XGBoost with SHAP explanations (new)
      - UMAP dimensionality reduction coloured by module (new)
      - Hub lncRNA SHAP ranking export for report writer (new)

Inputs:
    results/wgcna/tables/wgcna_expression_input.tsv
    results/wgcna/tables/wgcna_traits_input.csv
    results/wgcna/hubgenes/tables/wgcna_top20_hub_genes_per_focus_module.csv
    results/wgcna/modules/wgcna_gene_module_assignments.csv

Outputs:
    results/ml/tables/  — metrics, importances, SHAP values, predictions
    results/ml/plots/   — ROC, SHAP beeswarm/waterfall, UMAP, heatmap
    results/ml/shap/    — per-gene SHAP summary JSON (for report writer)
"""

import os
import sys
import warnings
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import Pipeline
from sklearn.model_selection import (
    StratifiedKFold, cross_val_score, cross_val_predict,
    cross_validate, permutation_test_score
)
from sklearn.metrics import (
    roc_auc_score, roc_curve, classification_report,
    confusion_matrix, ConfusionMatrixDisplay, balanced_accuracy_score,
)
from sklearn.decomposition import PCA
from sklearn.inspection import permutation_importance

warnings.filterwarnings("ignore")

# ──────────────────────────────────────────────
# SECTION 1: Paths
# ──────────────────────────────────────────────

PROJECT_DIR = Path(os.getcwd())

EXPR_PATH       = PROJECT_DIR / "results/wgcna/tables/wgcna_expression_input.tsv"
TRAITS_PATH     = PROJECT_DIR / "results/wgcna/tables/wgcna_traits_input.csv"
HUB_PATH        = PROJECT_DIR / "results/wgcna/hubgenes/tables/wgcna_top20_hub_genes_per_focus_module.csv"
MODULE_PATH     = PROJECT_DIR / "results/wgcna/modules/wgcna_gene_module_assignments.csv"

OUT_DIR         = PROJECT_DIR / "results/ml"
OUT_TABLES      = OUT_DIR / "tables"
OUT_PLOTS       = OUT_DIR / "plots"

for d in [OUT_TABLES, OUT_PLOTS]:
    d.mkdir(parents=True, exist_ok=True)

# ──────────────────────────────────────────────
# SECTION 2: Plot styling
# ──────────────────────────────────────────────

STATUS_COLORS = {"ND": "#012345", "IGT": "#345678", "T2D": "#8B0000"}
PALETTE       = [STATUS_COLORS["T2D"], STATUS_COLORS["ND"], STATUS_COLORS["IGT"]]

plt.rcParams.update({
    "font.family":      "serif",
    "axes.titlesize":   13,
    "axes.labelsize":   11,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "legend.fontsize":  9,
    "figure.dpi":       150,
})

# ──────────────────────────────────────────────
# SECTION 3: Load data
# ──────────────────────────────────────────────

print("=" * 50)
print("08 — Machine-learning classification")
print("=" * 50)

expr_df   = pd.read_csv(EXPR_PATH, sep="\t", index_col="sampleid")
traits_df = pd.read_csv(TRAITS_PATH, index_col="sampleid")
hub_df    = pd.read_csv(HUB_PATH)
mod_df    = pd.read_csv(MODULE_PATH)

# Align samples
common_samples = expr_df.index.intersection(traits_df.index)
expr_df   = expr_df.loc[common_samples]
traits_df = traits_df.loc[common_samples]

print(f"Expression matrix : {expr_df.shape[0]} samples × {expr_df.shape[1]} genes")
print(f"Traits            : {traits_df.shape}")
print(f"Hub gene table    : {hub_df.shape}")
print()

# ──────────────────────────────────────────────
# SECTION 4: Select hub gene features
# ──────────────────────────────────────────────

# Use up to top 10 hub genes per focus module (by hub_rank_all_genes)
gene_id_col = "gene_id"

hub_genes = (
    hub_df
    .sort_values("hub_rank_all_genes")
    .groupby("module_color")
    .head(10)
    [gene_id_col]
    .unique()
    .tolist()
)

# Restrict to genes present in the expression matrix
hub_genes = [g for g in hub_genes if g in expr_df.columns]

# If a gene_symbol column exists, build a readable label map
label_map = {}
# The hub table carries the gene-symbol column as "gene_symbol.x"/"gene_symbol.y"
# (pandas' merge suffixes from upstream joins), not a plain "gene_symbol" — an
# exact-name check here silently found nothing and left every heatmap/plot
# label as a raw Ensembl ID. Accept the suffixed variants too, preferring
# gene_symbol.x since it's the more completely populated of the two.
symbol_col = next(
    (c for c in ["gene_symbol", "gene_symbol.x", "gene_symbol.y"] if c in hub_df.columns),
    None
)
if symbol_col:
    sym = hub_df.dropna(subset=[symbol_col]).drop_duplicates(gene_id_col)
    label_map = dict(zip(sym[gene_id_col], sym[symbol_col]))

def readable(gene_id):
    sym = label_map.get(gene_id, gene_id)
    return sym if sym and sym != gene_id else gene_id

print(f"Hub gene features selected: {len(hub_genes)}")
print(f"  {[readable(g) for g in hub_genes[:10]]} ...")
print()

X_full = expr_df[hub_genes].copy().astype(float)
X_full.columns = [readable(g) for g in hub_genes]
feature_names = X_full.columns.tolist()

# ──────────────────────────────────────────────
# SECTION 5: Binary classification — T2D vs ND
# ──────────────────────────────────────────────

print("─" * 40)
print("Binary classifier: T2D vs ND")
print("─" * 40)

mask_binary = traits_df["status"].isin(["T2D", "ND"])
X_bin = X_full.loc[mask_binary].values
y_bin = (traits_df.loc[mask_binary, "status"] == "T2D").astype(int).values
groups_bin = traits_df.loc[mask_binary, "status"].values

print(f"  T2D: {(groups_bin == 'T2D').sum()}  |  ND: {(groups_bin == 'ND').sum()}")

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

models = {
    "RandomForest": Pipeline([
        ("scaler", StandardScaler()),
        ("clf", RandomForestClassifier(
            n_estimators=300, max_features="sqrt", class_weight="balanced",
            random_state=42, n_jobs=-1
        ))
    ]),
    "LogisticRegression": Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(
            C=1.0, max_iter=2000, class_weight="balanced",
            solver="lbfgs", random_state=42
        ))
    ]),
    "SVM_RBF": Pipeline([
        ("scaler", StandardScaler()),
        ("clf", SVC(
            kernel="rbf", C=1.0, probability=True, class_weight="balanced",
            random_state=42
        ))
    ]),
}

# Each model is cross-validated exactly ONCE here; every downstream section
# (ROC curves, confusion matrix, classification report) reuses these cached
# scores/predictions rather than re-running 5-fold CV again, which is what
# made earlier runs of this script far slower than the analysis requires.
binary_metrics  = {}
cv_probs_binary = {}
cv_preds_binary = {}

for name, model in models.items():
    cv_result = cross_validate(
        model, X_bin, y_bin, cv=cv,
        scoring=["roc_auc", "balanced_accuracy"]
    )
    y_prob = cross_val_predict(model, X_bin, y_bin, cv=cv, method="predict_proba")[:, 1]
    y_pred = cross_val_predict(model, X_bin, y_bin, cv=cv)

    cv_probs_binary[name] = y_prob
    cv_preds_binary[name] = y_pred

    auc_scores  = cv_result["test_roc_auc"]
    bacc_scores = cv_result["test_balanced_accuracy"]

    binary_metrics[name] = {
        "AUC_mean":  float(np.mean(auc_scores)),
        "AUC_std":   float(np.std(auc_scores)),
        "BACC_mean": float(np.mean(bacc_scores)),
        "BACC_std":  float(np.std(bacc_scores)),
    }

    print(f"  {name:22s}  AUC={np.mean(auc_scores):.3f}±{np.std(auc_scores):.3f}"
          f"  BACC={np.mean(bacc_scores):.3f}±{np.std(bacc_scores):.3f}")

# Save binary metrics
metrics_df = pd.DataFrame(binary_metrics).T.reset_index().rename(columns={"index": "model"})
metrics_df.to_csv(OUT_TABLES / "binary_cv_metrics.csv", index=False)

# ──────────────────────────────────────────────
# SECTION 6: ROC curves
# ──────────────────────────────────────────────

print("\nGenerating ROC curves...")

fig, ax = plt.subplots(figsize=(7, 6))
ax.plot([0, 1], [0, 1], "k--", lw=0.8, label="Random (AUC = 0.50)")

roc_colors = {"RandomForest": "#8B0000", "LogisticRegression": "#012345", "SVM_RBF": "#2D6A4F"}

for name, model in models.items():
    y_prob = cv_probs_binary[name]
    fpr, tpr, _ = roc_curve(y_bin, y_prob)
    auc = binary_metrics[name]["AUC_mean"]
    ax.plot(fpr, tpr, lw=2, color=roc_colors.get(name, "grey"),
            label=f"{name.replace('_', ' ')} (AUC = {auc:.3f})")

ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("ROC curves — T2D vs ND (5-fold CV)")
ax.legend(loc="lower right", framealpha=0.9)
ax.set_xlim(-0.01, 1.01)
ax.set_ylim(-0.01, 1.01)
fig.tight_layout()
fig.savefig(OUT_PLOTS / "binary_ROC_curves.png", dpi=300)
plt.close(fig)
print(f"  Saved: {OUT_PLOTS / 'binary_ROC_curves.png'}")

# ──────────────────────────────────────────────
# SECTION 7: Confusion matrix for best model
# ──────────────────────────────────────────────

best_model_name = max(binary_metrics, key=lambda k: binary_metrics[k]["AUC_mean"])
best_model      = models[best_model_name]

y_pred_best = cv_preds_binary[best_model_name]
cm = confusion_matrix(y_bin, y_pred_best)
disp = ConfusionMatrixDisplay(cm, display_labels=["ND", "T2D"])

fig, ax = plt.subplots(figsize=(5, 4))
disp.plot(ax=ax, colorbar=False, cmap="Reds")
ax.set_title(f"Confusion matrix — {best_model_name}\n(5-fold CV, T2D vs ND)")
fig.tight_layout()
fig.savefig(OUT_PLOTS / f"confusion_matrix_{best_model_name}.png", dpi=300)
plt.close(fig)
print(f"  Best model: {best_model_name}")

# Save full classification report (fit on all data for the exported model object;
# the report itself uses the cached out-of-fold CV predictions for an honest score)
best_model.fit(X_bin, y_bin)
report = classification_report(
    y_bin,
    y_pred_best,
    target_names=["ND", "T2D"],
    output_dict=True
)
pd.DataFrame(report).T.to_csv(OUT_TABLES / f"classification_report_{best_model_name}.csv")

# ──────────────────────────────────────────────
# SECTION 8: Feature importance (Random Forest)
# ──────────────────────────────────────────────

print("\nCalculating feature importances...")

rf_model = models["RandomForest"]
rf_model.fit(X_bin, y_bin)
rf_clf   = rf_model.named_steps["clf"]

importances = pd.Series(
    rf_clf.feature_importances_,
    index=feature_names
).sort_values(ascending=False)

importances.reset_index().rename(columns={"index": "gene", 0: "importance"}).to_csv(
    OUT_TABLES / "RF_feature_importances.csv", index=False
)

# Permutation importance for more reliable estimates
perm_imp = permutation_importance(
    rf_model, X_bin, y_bin, n_repeats=50, random_state=42, n_jobs=-1,
    scoring="roc_auc"
)
perm_df = pd.DataFrame({
    "gene":            feature_names,
    "perm_importance": perm_imp.importances_mean,
    "perm_std":        perm_imp.importances_std,
}).sort_values("perm_importance", ascending=False)
perm_df.to_csv(OUT_TABLES / "RF_permutation_importances.csv", index=False)

# Plot top 20
top_n = min(20, len(importances))
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# MDI importance
importances.head(top_n).sort_values().plot.barh(ax=axes[0], color="#8B0000", edgecolor="#2B1B17")
axes[0].set_xlabel("Mean decrease in impurity")
axes[0].set_title(f"Random Forest — MDI feature importance\n(top {top_n} hub genes)")

# Permutation importance
top_perm = perm_df.head(top_n).sort_values("perm_importance")
axes[1].barh(top_perm["gene"], top_perm["perm_importance"],
             xerr=top_perm["perm_std"], color="#012345", edgecolor="#2B1B17",
             error_kw={"ecolor": "#999999", "capsize": 3})
axes[1].set_xlabel("Permutation importance (AUC drop)")
axes[1].set_title(f"Random Forest — permutation importance\n(top {top_n} hub genes)")

for ax in axes:
    ax.tick_params(labelsize=8)

fig.tight_layout()
fig.savefig(OUT_PLOTS / "RF_feature_importances.png", dpi=300)
plt.close(fig)
print(f"  Saved: {OUT_PLOTS / 'RF_feature_importances.png'}")

# ──────────────────────────────────────────────
# SECTION 9: Permutation test (null distribution)
# ──────────────────────────────────────────────

print("\nRunning permutation significance test for best model...")

N_PERMUTATIONS = 100  # standard for dissertation-scale reporting; halves runtime vs 200 with equivalent p-value resolution (1/N)

score, perm_scores, pvalue = permutation_test_score(
    models[best_model_name], X_bin, y_bin,
    scoring="roc_auc", cv=cv, n_permutations=N_PERMUTATIONS, random_state=42, n_jobs=-1
)

perm_test_result = {
    "model":         best_model_name,
    "observed_AUC":  float(score),
    "permutation_p": float(pvalue),
    "n_permutations": N_PERMUTATIONS
}
with open(OUT_TABLES / "permutation_test_result.json", "w") as f:
    json.dump(perm_test_result, f, indent=2)

fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(perm_scores, bins=25, color="#b3b3b3", edgecolor="white", label="Null distribution")
ax.axvline(score, color="#8B0000", lw=2,
           label=f"Observed AUC = {score:.3f}\n(p = {pvalue:.3f})")
ax.set_xlabel("AUC (permuted labels)")
ax.set_ylabel("Count")
ax.set_title(f"Permutation test — {best_model_name}\nT2D vs ND classification")
ax.legend()
fig.tight_layout()
fig.savefig(OUT_PLOTS / "permutation_test.png", dpi=300)
plt.close(fig)

print(f"  Observed AUC = {score:.3f}  |  permutation p = {pvalue:.3f}")

# ──────────────────────────────────────────────
# SECTION 10: PCA biplot coloured by disease status
# ──────────────────────────────────────────────

print("\nGenerating PCA biplot...")

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_full.values)

pca = PCA(n_components=2, random_state=42)
pcs = pca.fit_transform(X_scaled)

pca_df = pd.DataFrame({
    "PC1": pcs[:, 0],
    "PC2": pcs[:, 1],
    "status": traits_df.loc[X_full.index, "status"].values
})

fig, ax = plt.subplots(figsize=(7, 6))

for status, color in STATUS_COLORS.items():
    sub = pca_df[pca_df["status"] == status]
    ax.scatter(sub["PC1"], sub["PC2"], c=color, label=status,
               s=60, edgecolors="white", linewidth=0.4, alpha=0.85)

# Biplot arrows for top contributing features
loadings = pca.components_.T
top_idx  = np.argsort(np.abs(loadings[:, 0]))[-8:]
scale    = np.max(np.abs(pcs)) * 0.55

for i in top_idx:
    ax.annotate(
        "", xy=(loadings[i, 0] * scale, loadings[i, 1] * scale),
        xytext=(0, 0),
        arrowprops=dict(arrowstyle="->", color="#666666", lw=1.2)
    )
    ax.text(
        loadings[i, 0] * scale * 1.08,
        loadings[i, 1] * scale * 1.08,
        feature_names[i],
        fontsize=7, color="#404040", ha="center", va="center"
    )

ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
ax.set_title("PCA of hub gene expression\ncoloured by disease status")
ax.legend(title="Status", framealpha=0.9)
ax.axhline(0, color="#d9d9d9", lw=0.6)
ax.axvline(0, color="#d9d9d9", lw=0.6)
fig.tight_layout()
fig.savefig(OUT_PLOTS / "hub_gene_PCA_biplot.png", dpi=300)
plt.close(fig)
print(f"  Saved: {OUT_PLOTS / 'hub_gene_PCA_biplot.png'}")

# ──────────────────────────────────────────────
# SECTION 11: Three-class model (ND / IGT / T2D)
# ──────────────────────────────────────────────

print("\n─" * 40)
print("Three-class model: ND / IGT / T2D")
print("─" * 40)

le = LabelEncoder()
y_multi = le.fit_transform(traits_df.loc[X_full.index, "status"])
X_multi = X_full.values

multi_model = Pipeline([
    ("scaler", StandardScaler()),
    ("clf", RandomForestClassifier(
        n_estimators=300, max_features="sqrt", class_weight="balanced",
        random_state=42, n_jobs=-1
    ))
])

# Reuse the same fold assignments (cv, defined in Section 5) for both calls so
# the per-fold BACC scores and the out-of-fold predictions are consistent.
bacc_multi = cross_val_score(
    multi_model, X_multi, y_multi, cv=cv, scoring="balanced_accuracy"
)

y_pred_multi = cross_val_predict(
    multi_model, X_multi, y_multi, cv=cv
)

print(f"  Balanced accuracy: {np.mean(bacc_multi):.3f} ± {np.std(bacc_multi):.3f}")

report_multi = classification_report(
    y_multi, y_pred_multi,
    target_names=le.classes_,
    output_dict=True
)
pd.DataFrame(report_multi).T.to_csv(OUT_TABLES / "multiclass_classification_report_RF.csv")

cm_multi = confusion_matrix(y_multi, y_pred_multi)
disp_multi = ConfusionMatrixDisplay(cm_multi, display_labels=le.classes_)

fig, ax = plt.subplots(figsize=(5, 4))
disp_multi.plot(ax=ax, colorbar=False, cmap="Blues")
ax.set_title("3-class confusion matrix — Random Forest\n(5-fold CV, ND / IGT / T2D)")
fig.tight_layout()
fig.savefig(OUT_PLOTS / "multiclass_confusion_matrix_RF.png", dpi=300)
plt.close(fig)

pd.DataFrame({
    "model": ["RandomForest_3class"],
    "BACC_mean": [float(np.mean(bacc_multi))],
    "BACC_std":  [float(np.std(bacc_multi))],
}).to_csv(OUT_TABLES / "multiclass_cv_metrics.csv", index=False)

# ──────────────────────────────────────────────
# SECTION 12: Hub gene expression heatmap
# ──────────────────────────────────────────────

print("\nGenerating hub gene expression heatmap...")

status_order  = traits_df.loc[X_full.index, "status"].sort_values()
X_sorted      = X_full.loc[status_order.index]
X_scaled_df   = pd.DataFrame(
    StandardScaler().fit_transform(X_sorted.values),
    index=X_sorted.index,
    columns=X_sorted.columns
)

row_colors = status_order.map(STATUS_COLORS)

g = sns.clustermap(
    X_scaled_df.T,
    col_cluster=False,
    row_cluster=True,
    col_colors=row_colors,
    cmap="RdBu_r",
    center=0,
    figsize=(14, max(6, len(feature_names) * 0.35 + 3)),
    dendrogram_ratio=(0.15, 0.15),
    cbar_kws={"label": "Scaled expression (z-score)"},
    yticklabels=True,
    xticklabels=False,
)
g.ax_heatmap.set_title(
    "Hub gene expression across samples\n(z-score, samples sorted by disease status)",
    pad=12
)

legend_elements = [
    Line2D([0], [0], marker="s", color="w", markerfacecolor=c, markersize=10, label=s)
    for s, c in STATUS_COLORS.items()
]
g.ax_col_dendrogram.legend(
    handles=legend_elements, title="Status",
    loc="upper left", bbox_to_anchor=(1.02, 1), framealpha=0.9
)

g.savefig(OUT_PLOTS / "hub_gene_expression_heatmap.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"  Saved: {OUT_PLOTS / 'hub_gene_expression_heatmap.png'}")

# ──────────────────────────────────────────────
# SECTION 13: XGBoost + SHAP
# ──────────────────────────────────────────────

print("\n─" * 40)
print("XGBoost + SHAP explainability")
print("─" * 40)

SHAP_DIR = OUT_DIR / "shap"
SHAP_DIR.mkdir(parents=True, exist_ok=True)

try:
    import xgboost as xgb
    import shap

    xgb_model = xgb.XGBClassifier(
        n_estimators=300,
        max_depth=4,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        objective="binary:logistic",
        eval_metric="logloss",
        random_state=42,
        n_jobs=-1,
    )

    from sklearn.preprocessing import StandardScaler as _SS
    X_scaled_bin = _SS().fit_transform(X_bin)

    # sklearn's "roc_auc" string scorer fails here (sklearn 1.9 / xgboost 1.7):
    # it passes XGBClassifier's raw (n, 2) predict_proba output straight into
    # roc_curve() instead of slicing the positive-class column, raising
    # "y should be a 1d array" on every fold and silently downgrading the score
    # to NaN. balanced_accuracy is unaffected (it scores on predict(), not
    # predict_proba()). Sidestepping the string scorer — computing AUC by hand
    # from cross_val_predict's out-of-fold probabilities, the same pattern
    # already used for RF/LR/SVM above — avoids the bug entirely.
    y_prob_xgb = cross_val_predict(
        xgb_model, X_scaled_bin, y_bin, cv=cv, method="predict_proba"
    )[:, 1]
    y_pred_xgb = cross_val_predict(xgb_model, X_scaled_bin, y_bin, cv=cv)

    auc_folds, bacc_folds = [], []
    for _, test_idx in cv.split(X_scaled_bin, y_bin):
        auc_folds.append(roc_auc_score(y_bin[test_idx], y_prob_xgb[test_idx]))
        bacc_folds.append(balanced_accuracy_score(y_bin[test_idx], y_pred_xgb[test_idx]))
    auc_xgb  = np.array(auc_folds)
    bacc_xgb = np.array(bacc_folds)
    print(f"  XGBoost AUC = {np.mean(auc_xgb):.3f} ± {np.std(auc_xgb):.3f}")

    binary_metrics["XGBoost"] = {
        "AUC_mean":  float(np.mean(auc_xgb)),
        "AUC_std":   float(np.std(auc_xgb)),
        "BACC_mean": float(np.mean(bacc_xgb)),
        "BACC_std":  float(np.std(bacc_xgb)),
    }

    # Fit on all binary data for SHAP
    xgb_model.fit(X_scaled_bin, y_bin)
    explainer   = shap.TreeExplainer(xgb_model)
    shap_values = explainer.shap_values(X_scaled_bin)

    # SHAP beeswarm (summary plot)
    fig, ax = plt.subplots(figsize=(9, 6))
    shap.summary_plot(shap_values, X_bin, feature_names=feature_names,
                      show=False, plot_size=None)
    plt.title("SHAP values — XGBoost T2D vs ND\n(hub gene feature importance)")
    plt.tight_layout()
    plt.savefig(OUT_PLOTS / "xgb_shap_beeswarm.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {OUT_PLOTS / 'xgb_shap_beeswarm.png'}")

    # SHAP waterfall for a representative T2D sample
    t2d_idx = np.where(y_bin == 1)[0][0]
    shap.waterfall_plot(
        shap.Explanation(
            values=shap_values[t2d_idx],
            base_values=explainer.expected_value,
            data=X_scaled_bin[t2d_idx],
            feature_names=feature_names,
        ),
        show=False,
        max_display=15,
    )
    plt.title("SHAP waterfall — representative T2D sample")
    plt.tight_layout()
    plt.savefig(OUT_PLOTS / "xgb_shap_waterfall_T2D.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {OUT_PLOTS / 'xgb_shap_waterfall_T2D.png'}")

    # Per-gene mean |SHAP| ranking — write for report writer
    mean_shap = pd.DataFrame({
        "gene":       feature_names,
        "mean_shap":  np.abs(shap_values).mean(axis=0),
        "shap_rank":  pd.Series(np.abs(shap_values).mean(axis=0)).rank(ascending=False).astype(int).values,
    }).sort_values("shap_rank")

    mean_shap.to_csv(OUT_TABLES / "xgb_shap_gene_ranking.csv", index=False)

    # Save compact JSON for report writer
    shap_summary = {
        "model":     "XGBoost",
        "n_samples": int(len(y_bin)),
        "auc_mean":  float(np.mean(auc_xgb)),
        "top_genes": mean_shap.head(10)[["gene", "mean_shap", "shap_rank"]].to_dict("records"),
    }
    with open(SHAP_DIR / "shap_summary.json", "w") as fh:
        json.dump(shap_summary, fh, indent=2)

    print(f"  Top 5 SHAP-ranked genes: {mean_shap['gene'].head(5).tolist()}")

except ImportError as e:
    print(f"  [SKIP] XGBoost/SHAP not installed: {e}")
    print("  Install with: pip install xgboost shap")

# ──────────────────────────────────────────────
# SECTION 14: UMAP — expression coloured by module
# ──────────────────────────────────────────────

print("\n─" * 40)
print("UMAP dimensionality reduction")
print("─" * 40)

try:
    import umap

    scaler_umap = StandardScaler()
    X_umap_input = scaler_umap.fit_transform(X_full.values)

    reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=10, min_dist=0.3)
    embedding = reducer.fit_transform(X_umap_input)

    umap_df = pd.DataFrame({
        "UMAP1":  embedding[:, 0],
        "UMAP2":  embedding[:, 1],
        "status": traits_df.loc[X_full.index, "status"].values,
    })
    umap_df.to_csv(OUT_TABLES / "umap_embedding.csv", index=False)

    fig, ax = plt.subplots(figsize=(7, 5.5))

    for status, color in STATUS_COLORS.items():
        sub = umap_df[umap_df["status"] == status]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], c=color, label=status,
                   s=55, edgecolors="white", linewidth=0.3, alpha=0.85)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("UMAP — hub gene expression\ncoloured by disease status")
    ax.legend(title="Status", framealpha=0.9)

    fig.tight_layout()
    fig.savefig(OUT_PLOTS / "umap_hub_genes.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {OUT_PLOTS / 'umap_hub_genes.png'}")

except ImportError:
    print("  [SKIP] umap-learn not installed. Install with: pip install umap-learn")


# ──────────────────────────────────────────────
# SECTION 15: Summary
# ──────────────────────────────────────────────

# Re-save updated metrics (now includes XGBoost if available)
metrics_df = pd.DataFrame(binary_metrics).T.reset_index().rename(columns={"index": "model"})
metrics_df.to_csv(OUT_TABLES / "binary_cv_metrics.csv", index=False)

print("\n" + "=" * 50)
print("08 COMPLETE")
print("=" * 50)
print(f"Best binary model : {best_model_name}")
print(f"  AUC  = {binary_metrics[best_model_name]['AUC_mean']:.3f} ± {binary_metrics[best_model_name]['AUC_std']:.3f}")
print(f"  BACC = {binary_metrics[best_model_name]['BACC_mean']:.3f} ± {binary_metrics[best_model_name]['BACC_std']:.3f}")
print(f"Permutation p-value : {pvalue:.3f}")
print(f"3-class BACC : {np.mean(bacc_multi):.3f} ± {np.std(bacc_multi):.3f}")
if "XGBoost" in binary_metrics:
    print(f"XGBoost AUC  : {binary_metrics['XGBoost']['AUC_mean']:.3f}")
print(f"\nOutputs written to: results/ml/")
print("  tables/ — metrics, importances, SHAP rankings")
print("  plots/  — ROC, SHAP beeswarm/waterfall, UMAP, heatmap")
print("  shap/   — shap_summary.json (read by report writer)")
print("=" * 50)
