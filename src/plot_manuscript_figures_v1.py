"""Render first-draft manuscript figures from public aggregate outputs only.

This script does not access individual response values, expression matrices,
sample-level predictions, or model checkpoints. It creates editable SVGs and
600 dpi PNG previews under figures/manuscript_v1/.
"""
from __future__ import annotations

import json
import hashlib
import platform
from pathlib import Path
from datetime import datetime, timezone

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch
import numpy as np
import pandas as pd
import matplotlib


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "manuscript_v1"
DRUGS = [
    (1931, "MIRA-1"),
    (1036, "PLX-4720"),
    (1191, "Bortezomib"),
    (1930, "Telomerase Inhibitor IX"),
    (1199, "Tamoxifen"),
    (1091, "BMS-536924"),
    (1010, "Gefitinib"),
    (1922, "Cediranib"),
]
NAVY = "#17324D"
BLUE = "#2878B5"
TEAL = "#1B998B"
ORANGE = "#E58B2A"
GRAY = "#73808C"
LIGHT = "#F1F5F8"
RED = "#B84A4A"
SOURCE_FILES = [
    "results/cohort_v2/public_cohort_snapshot.json",
    "results/ordinal_repeated/paired_effect_intervals.csv",
    "results/p6_preflight/final_model_manifest.jsonl",
    "results/p6_external/evaluation_summary.json",
    "docs/ORDINAL_REPEATED_V1_PROTOCOL.md",
    "docs/EXTERNAL_VALIDATION_PROTOCOL_V1_1.md",
    "results/tissue_loto/contrasts_by_drug_tissue.csv",
    "config/tissue_sensitivity_loto_v1.json",
    "docs/TISSUE_SENSITIVITY_LOTO_PROTOCOL_V1.md",
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 12,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#CBD3DA",
            "axes.linewidth": 0.8,
            "savefig.facecolor": "white",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def _save(fig: plt.Figure, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    svg_path = OUT / f"{name}.svg"
    fig.savefig(svg_path, bbox_inches="tight")
    # Matplotlib wraps some SVG path data with trailing spaces; strip them so
    # generated text assets pass repository whitespace checks without changing
    # SVG geometry or rendering.
    svg_lines = svg_path.read_text(encoding="utf-8").splitlines()
    svg_path.write_text("\n".join(line.rstrip() for line in svg_lines) + "\n",
                        encoding="utf-8", newline="\n")
    fig.savefig(OUT / f"{name}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def _box(ax, xy, width, height, title, body, color, title_color="white"):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.02,rounding_size=0.025",
        linewidth=0.8, edgecolor=color, facecolor=LIGHT,
    )
    ax.add_patch(patch)
    ax.add_patch(FancyBboxPatch(
        (x, y + height - 0.105), width, 0.105,
        boxstyle="round,pad=0.02,rounding_size=0.025",
        linewidth=0, facecolor=color,
    ))
    ax.text(x + width / 2, y + height - 0.053, title, ha="center", va="center",
            color=title_color, weight="bold", fontsize=9)
    ax.text(x + width / 2, y + (height - 0.105) / 2, body, ha="center", va="center",
            color=NAVY, fontsize=7.4, linespacing=1.16)


def _arrow(ax, start, end, color=GRAY):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12,
                                 linewidth=1.2, color=color, connectionstyle="arc3"))


def fig1_design() -> None:
    cohort = json.loads((ROOT / SOURCE_FILES[0]).read_text(encoding="utf-8"))["cohort_counts"]
    final_rows = [json.loads(line) for line in (ROOT / SOURCE_FILES[2]).read_text(encoding="utf-8").splitlines() if line.strip()]
    ext_summary = json.loads((ROOT / SOURCE_FILES[3]).read_text(encoding="utf-8"))
    development_groups = cohort["development_unique_groups"]
    development_observations = cohort["development_drug_model_observations"]
    external_groups = cohort["external_unique_groups"]
    external_responses = ext_summary["response_audit"]["response_valid_drug_model_rows"]
    final_models = len(final_rows)
    final_fits = sum(int(row["fit_count"]) for row in final_rows)
    prediction_rows = sum(int(row["predicted"]) for row in ext_summary["prediction_coverage"])
    fig, ax = plt.subplots(figsize=(12, 5.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.5, 0.955, "Study design and evaluation boundaries", ha="center",
            va="center", fontsize=15, weight="bold", color=NAVY)

    _box(ax, (0.025, 0.60), 0.19, 0.27, "Development data",
         f"GDSC2 continuous AUC\n{development_groups:,} independent groups\n{development_observations:,} drug–model records\n{len(DRUGS)} prespecified compounds", BLUE)
    _box(ax, (0.025, 0.17), 0.19, 0.27, "Expression input",
         "CCLE gene TPM\nShared by development and\nexternal prediction\nNo expression-platform test", TEAL)
    _box(ax, (0.275, 0.60), 0.23, 0.27, "Nested development",
         "5 outer folds × 3 repeats\n3 inner folds\nGrouped by cell-line identity\nAll selection inside training", NAVY)
    _box(ax, (0.275, 0.17), 0.23, 0.27, "Eight model arms",
         "Tissue reference + 2,000-gene\nexpression baseline\n2 candidate rules × 3 forms:\nexpression / rank / gene pairs", TEAL)
    _box(ax, (0.565, 0.60), 0.18, 0.27, "Development readout",
         "Out-of-fold predictions\nDrug-wise Spearman and RMSE\nPaired contrasts and intervals\nZero solutions retained", BLUE)
    _box(ax, (0.79, 0.60), 0.18, 0.27, "Final model audit",
         f"{final_models} models / {final_fits:,} fits\nIndependent reload audit\n{prediction_rows:,} predictions\nmade before response access", NAVY)
    _box(ax, (0.79, 0.17), 0.18, 0.27, "Locked CTRPv2",
         f"{external_groups:,} non-overlapping groups\n{external_responses:,} drug–model responses\nOne-time cross-screen ranking\nSame CCLE expression source", ORANGE)

    _arrow(ax, (0.215, 0.735), (0.275, 0.735), BLUE)
    _arrow(ax, (0.215, 0.305), (0.275, 0.685), TEAL)
    _arrow(ax, (0.39, 0.60), (0.39, 0.44), NAVY)
    _arrow(ax, (0.505, 0.305), (0.565, 0.685), TEAL)
    _arrow(ax, (0.745, 0.735), (0.79, 0.735), BLUE)
    _arrow(ax, (0.88, 0.60), (0.88, 0.44), ORANGE)
    ax.text(0.88, 0.08, "Retrospective screen transfer; not clinical or cross-platform validation",
            ha="center", va="center", color=RED, fontsize=8, style="italic")
    _save(fig, "Figure1_study_design")


def fig2_development_primary() -> None:
    path = ROOT / "results" / "ordinal_repeated" / "paired_effect_intervals.csv"
    df = pd.read_csv(path)
    df["gdsc_drug_id"] = df["gdsc_drug_id"].astype(str)
    rows = df[(df["contrast"] == "primary:association_pair-strong_expression") &
              (df["metric"] == "spearman") &
              (df["family_name"] == "primary_fixed9")].copy()
    order = [str(x[0]) for x in DRUGS] + ["macro_equal_fixed_eight_drugs"]
    rows = rows.set_index("gdsc_drug_id").loc[order].reset_index()
    names = [x[1] for x in DRUGS] + ["Equal-drug macro"]
    y = np.arange(len(rows))[::-1]
    point = rows["point"].to_numpy(float)
    low = rows["family_adjusted_conditional_lower"].to_numpy(float)
    high = rows["family_adjusted_conditional_upper"].to_numpy(float)

    fig, ax = plt.subplots(figsize=(8.8, 5.4))
    ax.axvline(0, color=GRAY, linewidth=1, linestyle="--", zorder=0)
    ax.errorbar(point[:-1], y[:-1], xerr=[point[:-1] - low[:-1], high[:-1] - point[:-1]],
                fmt="o", color=BLUE, ecolor=BLUE, capsize=2.8, markersize=5, linewidth=1.5)
    ax.errorbar(point[-1], y[-1], xerr=[[point[-1] - low[-1]], [high[-1] - point[-1]]],
                fmt="D", color=ORANGE, ecolor=ORANGE, capsize=3.5, markersize=6.3, linewidth=1.8)
    ax.set_yticks(y, names)
    ax.set_xlabel("Δ Spearman ρ (REO pair − strong expression)")
    ax.set_title("Repeated nested development evaluation", loc="left", color=NAVY, weight="bold")
    ax.text(0.01, -0.16,
            "Family-9 adjusted conditional intervals; bootstrap is conditional on recorded OOF predictions.",
            transform=ax.transAxes, color=GRAY, fontsize=8)
    ax.grid(axis="x", color="#E7ECF0", linewidth=0.7)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.09, 0.09)
    _save(fig, "Figure2_development_primary")


def fig3_representation_and_zero() -> None:
    paired = pd.read_csv(ROOT / "results" / "ordinal_repeated" / "paired_effect_intervals.csv")
    paired["gdsc_drug_id"] = paired["gdsc_drug_id"].astype(str)
    rep = paired[(paired["metric"] == "spearman") &
                 (paired["family_name"] == "representation_fixed18") &
                 (paired["contrast"].isin(["mad:pair_full-expression", "response_association:pair_full-expression"]))].copy()
    order = [str(x[0]) for x in DRUGS] + ["macro_equal_fixed_eight_drugs"]
    rep["drug_order"] = rep["gdsc_drug_id"].map({d: i for i, d in enumerate(order)})
    rep["rule_order"] = rep["contrast"].map({"mad:pair_full-expression": 0, "response_association:pair_full-expression": 1})
    rep = rep.sort_values(["drug_order", "rule_order"])
    labels = []
    ys = []
    for i, drug_id in enumerate(order):
        for rule, yoff in [("mad:pair_full-expression", 0.15), ("response_association:pair_full-expression", -0.15)]:
            row = rep[(rep["gdsc_drug_id"] == drug_id) & (rep["contrast"] == rule)].iloc[0]
            labels.append((row, rule))
            ys.append(len(order) * 2 - 1 - 2 * i - (0 if yoff > 0 else 1))
    yvals = np.array(ys, dtype=float)

    manifest = [json.loads(line) for line in (ROOT / SOURCE_FILES[2]).read_text(encoding="utf-8").splitlines() if line.strip()]
    matrix = []
    model_names = ["expression_strong", "mad_expression", "mad_rank", "mad_pair_full",
                   "response_association_expression", "response_association_rank", "response_association_pair_full"]
    for drug_id, _ in DRUGS:
        lookup = {x["task"]["model"]: x["n_nonzero_molecular_coeffs"] for x in manifest if x["task"]["drug_id"] == drug_id}
        matrix.append([1 if lookup[m] > 0 else 0 for m in model_names])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.4, 6.1), gridspec_kw={"width_ratios": [1.35, 1]})
    ax1.axvline(0, color=GRAY, linestyle="--", linewidth=1)
    for idx, (row, rule) in enumerate(labels):
        color = BLUE if rule.startswith("mad") else TEAL
        pt = float(row["point"])
        lo = float(row["family_adjusted_conditional_lower"])
        hi = float(row["family_adjusted_conditional_upper"])
        ax1.errorbar(pt, yvals[idx], xerr=[[pt - lo], [hi - pt]], fmt="o", color=color,
                     ecolor=color, capsize=2, markersize=4.4, linewidth=1.15)
    ylabels = []
    for _, name in DRUGS:
        ylabels.extend([name + " · MAD", name + " · response"])
    ylabels.extend(["Equal-drug macro · MAD", "Equal-drug macro · response"])
    ax1.set_yticks(np.arange(len(ylabels))[::-1], ylabels)
    ax1.set_xlabel("Δ Spearman ρ (pair − expression, same 200 genes)")
    ax1.set_title("A  Matched candidate-set representation", loc="left", color=NAVY, weight="bold")
    ax1.grid(axis="x", color="#E7ECF0", linewidth=0.7)
    ax1.spines[["top", "right", "left"]].set_visible(False)
    ax1.tick_params(axis="y", length=0, labelsize=7)
    ax1.set_xlim(-0.12, 0.12)

    cmap = plt.matplotlib.colors.ListedColormap(["#DDE3E8", TEAL])
    ax2.imshow(np.asarray(matrix), aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
    short = ["Expr\n2,000", "MAD\nExpr", "MAD\nRank", "MAD\nPair", "Assoc.\nExpr", "Assoc.\nRank", "Assoc.\nPair"]
    ax2.set_xticks(range(len(short)), short)
    ax2.set_yticks(range(len(DRUGS)), [n for _, n in DRUGS])
    ax2.set_title("B  Full-development model coefficients", loc="left", color=NAVY, weight="bold")
    ax2.set_xticks(np.arange(-.5, len(short), 1), minor=True)
    ax2.set_yticks(np.arange(-.5, len(DRUGS), 1), minor=True)
    ax2.grid(which="minor", color="white", linewidth=2)
    ax2.tick_params(which="minor", bottom=False, left=False)
    ax2.tick_params(axis="y", length=0, labelsize=7)
    ax2.tick_params(axis="x", length=0, labelsize=7)
    ax2.text(0.02, -0.17, "Gray: zero molecular coefficients   Teal: at least one non-zero coefficient",
             transform=ax2.transAxes, color=GRAY, fontsize=7.5)
    fig.suptitle("Representation contrasts and molecular model solutions", x=0.04, ha="left",
                 fontsize=14, weight="bold", color=NAVY)
    fig.subplots_adjust(top=0.88, bottom=0.18, wspace=0.42)
    _save(fig, "Figure3_representation_and_zero_solutions")


def fig4_external_primary() -> None:
    summary = json.loads((ROOT / SOURCE_FILES[3]).read_text(encoding="utf-8"))
    contrast = summary["evaluation"]["contrasts"]["primary"]["inference"]
    per_drug = contrast["per_drug"]
    point = []
    low = []
    high = []
    ns = []
    for drug_id, _ in DRUGS:
        row = per_drug[str(drug_id)]
        interval = row["delta_interval"]
        point.append(float(row["delta_rho"]))
        low.append(float(interval["family_adjusted_lower"]))
        high.append(float(interval["family_adjusted_upper"]))
        ns.append(int(row["n_common_groups"]))
    macro = contrast["strict_macro"]["delta_interval"]
    point.append(float(macro["point"]))
    low.append(float(macro["family_adjusted_lower"]))
    high.append(float(macro["family_adjusted_upper"]))
    ns.append("8 drugs")

    names = [n for _, n in DRUGS] + ["Equal-drug macro"]
    y = np.arange(len(names))[::-1]
    p, lo, hi = np.asarray(point), np.asarray(low), np.asarray(high)
    loto_config = json.loads((ROOT / SOURCE_FILES[7]).read_text(encoding="utf-8"))
    display_tissues = loto_config["expected_primary_display_tissues"]
    loto = pd.read_csv(ROOT / SOURCE_FILES[6])
    loto["drug_id"] = loto["drug_id"].astype(str)
    loto = loto[loto["heldout_tissue"].isin(display_tissues)]
    defined = loto.pivot(index="drug_id", columns="heldout_tissue",
                         values="pair_minus_expression_strong_defined").reindex(
                             index=[str(d) for d, _ in DRUGS], columns=display_tissues)
    defined_matrix = defined.fillna(False).astype(bool).to_numpy(dtype=int)
    defined_counts = defined_matrix.sum(axis=1)
    tissue_labels = ["Esoph./\nStomach" if x == "Esophagus/Stomach" else
                     "Ovary/\nFallopian" if x == "Ovary/Fallopian Tube" else
                     "CNS/\nBrain" if x == "CNS/Brain" else x
                     for x in display_tissues]

    fig, (ax_loto, ax_ext) = plt.subplots(
        1, 2, figsize=(13.6, 6.0), gridspec_kw={"width_ratios": [0.95, 1.45]})
    loto_cmap = plt.matplotlib.colors.ListedColormap(["#DDE3E8", TEAL])
    ax_loto.imshow(defined_matrix, aspect="auto", interpolation="nearest",
                   cmap=loto_cmap, vmin=0, vmax=1)
    ax_loto.set_xticks(range(len(display_tissues)), tissue_labels, rotation=45, ha="right")
    ax_loto.set_yticks(range(len(DRUGS)),
                       [f"{name} ({count}/10)" for (_, name), count in zip(DRUGS, defined_counts)])
    ax_loto.set_title("A  LOTO defined contrasts", loc="left", color=NAVY, weight="bold")
    ax_loto.set_xticks(np.arange(-.5, len(display_tissues), 1), minor=True)
    ax_loto.set_yticks(np.arange(-.5, len(DRUGS), 1), minor=True)
    ax_loto.grid(which="minor", color="white", linewidth=1.5)
    ax_loto.tick_params(which="minor", bottom=False, left=False)
    ax_loto.tick_params(axis="both", length=0, labelsize=7)
    ax_loto.set_xlabel("Prespecified held-out tissue")
    ax_loto.legend(handles=[
        Patch(facecolor=TEAL, label="Defined Δρ"),
        Patch(facecolor="#DDE3E8", label="Undefined Δρ"),
    ], loc="lower left", bbox_to_anchor=(0.0, -0.31), frameon=False, ncol=2, fontsize=7)

    ax_ext.axvline(0, color=GRAY, linewidth=1, linestyle="--", zorder=0)
    ax_ext.errorbar(p[:-1], y[:-1], xerr=[p[:-1] - lo[:-1], hi[:-1] - p[:-1]],
                    fmt="o", color=BLUE, ecolor=BLUE, capsize=2.7, markersize=5, linewidth=1.35)
    ax_ext.errorbar(p[-1], y[-1], xerr=[[p[-1] - lo[-1]], [hi[-1] - p[-1]]],
                    fmt="D", color=ORANGE, ecolor=ORANGE, capsize=3.5, markersize=6.3, linewidth=1.8)
    ax_ext.set_yticks(y, names)
    ax_ext.set_xlabel("Δ Spearman ρ (REO pair − strong expression)")
    ax_ext.set_title("B  Locked CTRPv2 cross-screen evaluation", loc="left", color=NAVY, weight="bold")
    ax_ext.grid(axis="x", color="#E7ECF0", linewidth=0.7)
    ax_ext.spines[["top", "right", "left"]].set_visible(False)
    ax_ext.tick_params(axis="y", length=0, labelsize=8)
    ax_ext.set_xlim(-0.34, 0.34)
    for yi, n in zip(y, ns):
        label = n if isinstance(n, str) else f"N={n}"
        ax_ext.text(0.355, yi, label, va="center", ha="left", fontsize=7, color=GRAY, clip_on=False)
    fig.text(0.01, 0.025,
             "LOTO strict eight-drug macro is NA (34/80 cells defined). CTRPv2 intervals are family-9 adjusted and conditional on locked models; three drug contrasts are identically zero because both models are molecular zero solutions.",
             color=GRAY, fontsize=7.2, ha="left", va="bottom")
    fig.subplots_adjust(left=0.14, right=0.93, top=0.86, bottom=0.25, wspace=0.34)
    _save(fig, "Figure4_loto_and_external")


def main() -> None:
    _style()
    fig1_design()
    fig2_development_primary()
    fig3_representation_and_zero()
    fig4_external_primary()
    output_files = sorted(p for p in OUT.iterdir() if p.suffix in {".png", ".svg"})
    manifest = {
        "schema_version": "1.0",
        "status": "first_draft_rendered_from_aggregate_outputs",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script": "src/plot_manuscript_figures_v1.py",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "individual_response_values_read": False,
        "script_sha256": _sha256(ROOT / "src" / "plot_manuscript_figures_v1.py"),
        "source_sha256": {name: _sha256(ROOT / name) for name in SOURCE_FILES},
        "output_sha256": {p.name: _sha256(p) for p in output_files},
    }
    (OUT / "figure_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "rendered", "output": str(OUT), "figures": 4, "manifest": "figure_manifest.json"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
