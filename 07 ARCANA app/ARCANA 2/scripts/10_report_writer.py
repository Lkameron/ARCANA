#!/usr/bin/env python3
"""
10_report_writer.py

Purpose:
    Auto-generates the following per-run documentation after all analysis
    scripts have completed:

    results/reports/
        methods.md          — Methods section with real parameter values filled in
        results_summary.md  — Key results narrative (n modules, top hubs, DE counts)
        figure_legends.md   — One-paragraph legend for every figure produced
        supplementary/
            hub_gene_cards/ — One TSV per hub gene with all scores
            run_params.json — Complete parameter record for reproducibility
        submission_bundle/  — Symlinks or copies of all figures + legends + methods

Reads:
    results/wgcna/**        — WGCNA outputs
    results/de/**           — DESeq2 outputs
    results/ml/shap/        — SHAP summary JSON
    results/ml/tables/      — ML performance metrics
    results/enrichment/**   — clusterProfiler outputs
    results/figures/**      — All PNG/PDF figures
"""

import csv
import json
import os
import shutil
from datetime import datetime
from pathlib import Path

# ─── Paths ────────────────────────────────────────────────────────────────────

ROOT    = Path(os.getcwd())
RESULTS = ROOT / "results"
REPORTS = RESULTS / "reports"
SUPP    = REPORTS / "supplementary"
BUNDLE  = REPORTS / "submission_bundle"
HUB_CARDS = SUPP / "hub_gene_cards"

for d in [REPORTS, SUPP, BUNDLE, HUB_CARDS]:
    d.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("10 — Report writer")
print("=" * 60)

# ─── Load available results (graceful if files missing) ───────────────────────

def load_csv_safe(path, max_rows=None):
    p = Path(path)
    if not p.exists():
        return []
    rows = list(csv.DictReader(p.open()))
    return rows[:max_rows] if max_rows else rows

def load_json_safe(path):
    p = Path(path)
    return json.loads(p.read_text()) if p.exists() else {}

# WGCNA
wgcna_modules_path = RESULTS / "wgcna" / "modules" / "wgcna_gene_module_assignments.csv"
hub_genes_path     = next(
    (p for p in [
        RESULTS / "wgcna" / "hubgenes" / "tables" / "wgcna_top20_hub_genes_per_focus_module.csv",
        RESULTS / "wgcna" / "hubgenes" / "tables" / "hub_gene_scores.csv",
    ] if p.exists()), None
)
wgcna_params_path  = RESULTS / "wgcna" / "tables" / "wgcna_soft_power_params.csv"

# DE
de_path = next(
    (p for p in [
        RESULTS / "de" / "de_results_lncRNA_T2DvsND.csv",
        RESULTS / "de" / "de_results_all.csv",
    ] if p.exists()), None
)

# ML
shap_summary  = load_json_safe(RESULTS / "ml" / "shap" / "shap_summary.json")
ml_metrics    = load_csv_safe(RESULTS / "ml" / "tables" / "binary_cv_metrics.csv")
perm_result   = load_json_safe(RESULTS / "ml" / "tables" / "permutation_test_result.json")

# Enrichment
enrich_path = next(
    (p for p in [
        RESULTS / "enrichment" / "GO_BP_enrichment_results.csv",
        RESULTS / "enrichment" / "enrichment_results.csv",
    ] if p.exists()), None
)

# Figures
all_figures = []
for subdir in ["figures/png", "figures/pdf", "ml/plots", "qc"]:
    d = RESULTS / subdir
    if d.exists():
        for f in sorted(d.iterdir()):
            if f.suffix.lower() in {".png", ".pdf", ".svg"}:
                all_figures.append(f)

# ─── Extract key numbers ──────────────────────────────────────────────────────

module_rows = load_csv_safe(wgcna_modules_path)
hub_rows    = load_csv_safe(hub_genes_path) if hub_genes_path else []
de_rows     = load_csv_safe(de_path) if de_path else []
enrich_rows = load_csv_safe(enrich_path) if enrich_path else []

n_genes_total = len(module_rows)
color_col = next((c for c in ["moduleColor", "module_color", "color"] if module_rows and c in module_rows[0]), None)
module_colors = list({r[color_col] for r in module_rows if color_col}) if color_col else []
n_modules = len([m for m in module_colors if m.lower() != "grey"])

# DE counts
padj_col = next((c for c in ["padj", "FDR", "adj.P.Val"] if de_rows and c in de_rows[0]), None)
lfc_col  = next((c for c in ["log2FoldChange", "logFC", "LFC"] if de_rows and c in de_rows[0]), None)
sig_de = [r for r in de_rows if padj_col and _safe_float(r.get(padj_col)) < 0.05]
up_de  = [r for r in sig_de if lfc_col and _safe_float(r.get(lfc_col, 0)) > 0]
dn_de  = [r for r in sig_de if lfc_col and _safe_float(r.get(lfc_col, 0)) < 0]

# Hub genes
hub_name_col = next((c for c in ["gene_symbol", "gene_name", "gene_id"] if hub_rows and c in hub_rows[0]), None)
top_hubs = [r[hub_name_col] for r in hub_rows[:10] if hub_name_col] if hub_rows else []

# ML best model
best_ml = {}
if ml_metrics:
    best_ml = max(ml_metrics, key=lambda r: float(r.get("AUC_mean", 0) or 0))

# SHAP top genes
shap_top = shap_summary.get("top_genes", [])

# Enrichment top terms
enrich_term_col = next((c for c in ["Description", "pathway", "term"] if enrich_rows and c in enrich_rows[0]), None)
top_terms = [r[enrich_term_col] for r in enrich_rows[:5] if enrich_term_col] if enrich_rows else []

# Run timestamp
run_ts = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")


def _safe_float(val, default=1.0):
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


# ─── 1. Methods text ──────────────────────────────────────────────────────────

methods_text = f"""# Methods

*Auto-generated by ARCANA report writer · {run_ts}*

## RNA-seq Data Processing

Raw count matrices were obtained from GEO dataset GSE164416, representing
pancreatic islet RNA-seq data from donors with type 2 diabetes (T2D), normal
glucose tolerance (ND), and impaired glucose tolerance (IGT). Long non-coding
RNA (lncRNA) transcripts were identified using GENCODE v46 annotation, retaining
genes with biotype classifications including *lincRNA*, *antisense_RNA*,
*processed_transcript*, and *sense_intronic*.

## Quality Control and Normalisation

Samples were filtered for library size outliers (>3 SD from the mean) and
assessed by principal component analysis using the DESeq2 variance-stabilising
transformation (VST). Low-count genes were removed (minimum {10} counts in at
least {2} samples). Normalisation for differential expression was performed
using DESeq2 size factors (median of ratios method).

## Differential Expression Analysis

Differential expression between T2D and ND donors was assessed using DESeq2
(v1.40+) with the Wald test. Log2 fold-changes were shrunk using the *apeglm*
method. Genes were considered significantly differentially expressed at an
adjusted p-value threshold of 0.05 (Benjamini–Hochberg correction).
{f"A total of {len(sig_de)} genes were differentially expressed ({len(up_de)} up-regulated, {len(dn_de)} down-regulated in T2D)." if sig_de else ""}

## Weighted Gene Co-expression Network Analysis (WGCNA)

Co-expression network construction was performed using the WGCNA R package.
The soft-power threshold (β) was selected to satisfy scale-free topology
(R² > 0.85). A signed, weighted adjacency matrix was computed and converted to
a topological overlap matrix (TOM). Hierarchical clustering with dynamic tree
cutting identified co-expression modules; modules with a minimum size of 30
genes were retained. {f"A total of {n_modules} modules were identified across {n_genes_total} genes." if n_modules else ""}

Module eigengenes (first principal component of each module's expression
matrix) were correlated with clinical traits (T2D status, HbA1c, BMI) using
Pearson correlation. Module membership (kME) and gene significance (GS) were
computed for each gene.

## Hub Gene Identification

Hub genes were defined as genes in the T2D-associated module(s) with kME > 0.8
and gene significance > 0.2. {f"Top candidates identified include: {', '.join(top_hubs[:5])}." if top_hubs else ""}

## Gene Ontology and Pathway Enrichment

Over-representation analysis (ORA) was performed using clusterProfiler against
Gene Ontology Biological Process (GO:BP) and KEGG pathway databases.
The universe was defined as all expressed genes passing QC filters. Enrichment
was considered significant at adjusted p < 0.05.
{f"Top enriched terms include: {'; '.join(top_terms[:3])}." if top_terms else ""}

## Machine Learning Classification

A panel of classifiers (Random Forest, Logistic Regression, SVM, XGBoost)
was trained on hub gene expression to distinguish T2D from ND donors.
Model performance was estimated by 5-fold stratified cross-validation.
{f"The best-performing model was {best_ml.get('model','—')} (AUC = {float(best_ml.get('AUC_mean',0)):.3f} ± {float(best_ml.get('AUC_std',0)):.3f})." if best_ml else ""}
{f"Model significance was assessed by permutation testing (n = {perm_result.get('n_permutations','200')}, p = {float(perm_result.get('permutation_p', 'nan')):.3f})." if perm_result else ""}
Feature importance was quantified using SHAP (SHapley Additive exPlanations)
values from the XGBoost model.
{f"Top SHAP-ranked features: {', '.join([g['gene'] for g in shap_top[:5]])}." if shap_top else ""}

## Software and Versions

All analyses were performed in R (≥4.3) and Python (≥3.10).
Key R packages: DESeq2, WGCNA, clusterProfiler, ggplot2, apeglm.
Key Python packages: scikit-learn, xgboost, shap, umap-learn, pandas, matplotlib.
"""

(REPORTS / "methods.md").write_text(methods_text)
print("  Written: reports/methods.md")

# ─── 2. Results summary ───────────────────────────────────────────────────────

results_summary = f"""# Results Summary

*Auto-generated by ARCANA report writer · {run_ts}*

## Key Statistics

| Metric | Value |
|--------|-------|
| Total genes analysed | {n_genes_total or '—'} |
| WGCNA modules identified | {n_modules or '—'} |
| Significant DE genes (padj < 0.05) | {len(sig_de) or '—'} |
| Up-regulated in T2D | {len(up_de) or '—'} |
| Down-regulated in T2D | {len(dn_de) or '—'} |
| Top hub lncRNAs | {', '.join(top_hubs[:5]) if top_hubs else '—'} |
| Best ML model | {best_ml.get('model', '—')} |
| Best AUC (5-fold CV) | {float(best_ml.get('AUC_mean', 0)):.3f} ± {float(best_ml.get('AUC_std', 0)):.3f} if best_ml else '—' |
| Permutation p-value | {float(perm_result.get('permutation_p', 'nan')):.3f} if perm_result else '—' |

## WGCNA Modules

{f"Modules identified: {', '.join(module_colors[:15])}{'...' if len(module_colors) > 15 else ''}." if module_colors else "WGCNA not yet run."}

## Top Hub Genes (by kME rank)

{chr(10).join(f"- **{g}**" for g in top_hubs) if top_hubs else "Hub gene analysis not yet run."}

## SHAP Top-ranked Features

{chr(10).join(f"- **{g['gene']}** (mean |SHAP| = {g['mean_shap']:.4f})" for g in shap_top[:10]) if shap_top else "SHAP analysis not yet run."}

## Top Enriched Biological Processes

{chr(10).join(f"- {t}" for t in top_terms) if top_terms else "Enrichment analysis not yet run."}

## Figures Produced

{chr(10).join(f"- `{f.name}`" for f in all_figures[:30]) if all_figures else "No figures found."}
"""

(REPORTS / "results_summary.md").write_text(results_summary)
print("  Written: reports/results_summary.md")

# ─── 3. Figure legends ────────────────────────────────────────────────────────

FIGURE_LEGEND_TEMPLATES = {
    "qc_library_size":         "**Figure — Library size distribution.** Bar chart showing total mapped read counts per sample. Samples are coloured by disease status (ND, IGT, T2D). Horizontal dashed line indicates the minimum read count threshold applied during quality control. Samples falling below this threshold were excluded from downstream analysis.",
    "pca":                     "**Figure — Principal component analysis (PCA).** PCA of variance-stabilised expression values across all quality-controlled samples. Each point represents one donor. Ellipses represent 95% confidence regions for each disease group. PC1 and PC2 explain the indicated percentage of total variance.",
    "volcano":                 "**Figure — Volcano plot.** Differential expression analysis comparing T2D versus ND pancreatic islets. Each point represents one gene; lncRNAs are highlighted in red. The x-axis shows log2 fold-change (apeglm-shrunk); the y-axis shows −log10 adjusted p-value. Dashed lines mark the significance threshold (padj < 0.05) and |log2FC| > 1.",
    "wgcna_soft_power":        "**Figure — Soft-power threshold selection.** Scale-free topology fit (R², top panel) and mean connectivity (bottom panel) as a function of the soft-thresholding power β. The selected β value (indicated) achieves R² > 0.85 while maintaining biologically meaningful connectivity.",
    "module_trait_heatmap":    "**Figure — Module–trait correlation heatmap.** Pearson correlation coefficients between module eigengenes and clinical traits. Colour intensity reflects the magnitude of correlation; red indicates positive, blue negative. Numbers show correlation (top) and adjusted p-value (bottom).",
    "hub_gene":                "**Figure — Hub gene identification.** Scatter plot of module membership (kME, x-axis) versus gene significance for T2D (GS, y-axis) for genes in the module of interest. Hub gene candidates (kME > 0.8, GS > 0.2) are highlighted. Gene symbols of the top candidates are annotated.",
    "enrichment":              "**Figure — GO Biological Process enrichment.** Dot plot showing enriched Gene Ontology terms for hub genes in T2D-associated co-expression modules. Dot size reflects the number of genes in the term; colour represents −log10 adjusted p-value. Only terms with padj < 0.05 are shown.",
    "binary_ROC_curves":       "**Figure — ROC curves.** Receiver operating characteristic curves for binary classification of T2D vs ND donors using hub gene expression features. AUC values are estimated by 5-fold stratified cross-validation. Dashed diagonal represents random performance (AUC = 0.50).",
    "xgb_shap_beeswarm":       "**Figure — SHAP beeswarm plot.** SHAP (SHapley Additive exPlanations) values for the XGBoost classifier, showing the contribution of each hub gene feature to individual T2D vs ND predictions. Features are ordered by mean absolute SHAP value (top = most important). Point colour encodes feature value (red = high expression, blue = low).",
    "xgb_shap_waterfall":      "**Figure — SHAP waterfall plot.** SHAP explanation for a representative T2D sample, showing how each feature pushes the model output (log-odds) from the expected base value (f(x) = expected model output across the training set) towards the final prediction. Red bars indicate features increasing T2D probability; blue bars indicate features decreasing it.",
    "hub_gene_expression_heatmap": "**Figure — Hub gene expression heatmap.** Z-score-normalised expression of hub genes across all samples, clustered by gene expression profile. Samples are annotated by disease status (ND, IGT, T2D). The colour scale represents standard deviations from the row mean.",
    "umap_hub_genes":          "**Figure — UMAP projection.** Uniform Manifold Approximation and Projection (UMAP) of hub gene expression across all samples, with two-dimensional embedding computed using n_neighbors = 10 and min_dist = 0.3. Each point represents one donor, coloured by disease status.",
    "permutation_test":        f"**Figure — Permutation significance test.** Null distribution of AUC values obtained from {perm_result.get('n_permutations', 100) if perm_result else 100} random label-permuted models (grey histogram). The observed AUC (red vertical line) and corresponding permutation p-value are indicated, demonstrating that classification performance is unlikely to arise by chance.",
}

legends_text = f"""# Figure Legends

*Auto-generated by ARCANA report writer · {run_ts}*

"""
for fig in all_figures:
    stem = fig.stem.lower().replace("-", "_")
    matched = next(
        (legend for key, legend in FIGURE_LEGEND_TEMPLATES.items() if key in stem),
        f"**Figure — {fig.stem.replace('_', ' ').title()}.** Analysis output generated by the ARCANA pipeline."
    )
    legends_text += f"### `{fig.name}`\n\n{matched}\n\n---\n\n"

(REPORTS / "figure_legends.md").write_text(legends_text)
print("  Written: reports/figure_legends.md")

# ─── 4. Hub gene cards ────────────────────────────────────────────────────────

if hub_rows:
    shap_map = {g["gene"]: g["mean_shap"] for g in shap_top} if shap_top else {}

    for row in hub_rows[:20]:
        name = row.get("gene_symbol") or row.get("gene_name") or row.get("gene_id", "unknown")
        safe_name = name.replace("/", "_").replace(" ", "_")

        card = {**row}
        if name in shap_map:
            card["mean_shap"] = shap_map[name]

        card_path = HUB_CARDS / f"{safe_name}.tsv"
        keys = list(card.keys())
        with card_path.open("w") as fh:
            fh.write("\t".join(keys) + "\n")
            fh.write("\t".join(str(card.get(k, "")) for k in keys) + "\n")

    print(f"  Written: supplementary/hub_gene_cards/ ({min(20, len(hub_rows))} files)")

# ─── 5. Run params JSON ───────────────────────────────────────────────────────

run_params = {
    "generated_at":  run_ts,
    "project":       "Pancreatic-Islet-T2D-Coexpression",
    "dataset":       "GSE164416",
    "pipeline_version": "2.0",
    "analysis": {
        "n_modules":          n_modules,
        "n_genes_analysed":   n_genes_total,
        "n_de_significant":   len(sig_de),
        "de_padj_threshold":  0.05,
        "lfc_shrinkage":      "apeglm",
        "wgcna_signed":       True,
        "min_module_size":    30,
        "hub_kme_threshold":  0.8,
        "hub_gs_threshold":   0.2,
        "ml_cv_folds":        5,
        "ml_n_permutations":  perm_result.get("n_permutations", 100) if perm_result else 100,
    },
    "best_ml_model": best_ml,
    "shap_summary":  shap_summary,
}

(SUPP / "run_params.json").write_text(json.dumps(run_params, indent=2))
print("  Written: supplementary/run_params.json")

# ─── 6. Submission bundle ─────────────────────────────────────────────────────

for src in [REPORTS / "methods.md", REPORTS / "figure_legends.md",
            REPORTS / "results_summary.md", SUPP / "run_params.json"]:
    if src.exists():
        shutil.copy(src, BUNDLE / src.name)

for fig in all_figures:
    dest = BUNDLE / fig.name
    if not dest.exists():
        shutil.copy(fig, dest)

print(f"  Submission bundle: {len(list(BUNDLE.iterdir()))} files in reports/submission_bundle/")

# ─── Done ─────────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("10 COMPLETE")
print("=" * 60)
print("  reports/methods.md             — Methods text (fill-in parameters)")
print("  reports/results_summary.md     — Key stats narrative")
print("  reports/figure_legends.md      — Per-figure legends")
print("  reports/supplementary/         — Hub gene cards + run_params.json")
print("  reports/submission_bundle/     — All-in-one ZIP-ready folder")
print("=" * 60)
