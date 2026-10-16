"""M08 figure layout revision using the original frozen aggregate inputs only.

No training, scoring, model loading, individual response or prediction access.
Old assets are preserved. Numerical plotted values are compared against artists
created by the original figure script in memory, without saving its figures.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.container import ErrorbarContainer
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
from matplotlib.path import Path as PlotPath
from matplotlib.text import Text

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "manuscript_m08_v1"
WIDTH_MM = 170.0
WORD_WIDTH_MM = 590 / 96 * 25.4
MIN_FONT_PT = 9.5
NAVY, BLUE, TEAL = "#17324D", "#2878B5", "#1B998B"
ORANGE, GRAY, LIGHT, RED = "#E58B2A", "#73808C", "#F1F5F8", "#B84A4A"
DRUGS = [(1931, "MIRA-1"), (1036, "PLX-4720"), (1191, "Bortezomib"),
         (1930, "Telomerase Inhibitor IX"), (1199, "Tamoxifen"),
         (1091, "BMS-536924"), (1010, "Gefitinib"), (1922, "Cediranib")]
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
TITLES = {
    "Figure1_study_design": "Study design and evaluation boundaries",
    "Figure2_development_primary": "Primary REO-pair versus strong-expression contrast in repeated nested development evaluation",
    "Figure3_representation_and_zero_solutions": "Candidate-set representation contrasts and full-development molecular solutions",
    "Figure4_loto_and_external": "Defined leave-one-tissue-out contrasts and locked CTRPv2 cross-screen evaluation",
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11, "axes.labelsize": 10, "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5, "legend.fontsize": 9.5,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "axes.edgecolor": "#CBD3DA", "axes.linewidth": .8,
        "savefig.facecolor": "white", "svg.fonttype": "none", "pdf.fonttype": 42})


def model_manifest():
    return [json.loads(line) for line in (ROOT / SOURCE_FILES[2]).read_text(encoding="utf-8").splitlines() if line.strip()]


def numerical_artists(fig):
    """Only x point/interval values and heatmaps; layout y positions may change."""
    payload = []
    for ax in fig.axes:
        points = []
        for item in ax.containers:
            if isinstance(item, ErrorbarContainer):
                xs = np.asarray(item.lines[0].get_xdata(), dtype=float).ravel()
                segments = item.lines[2][0].get_segments()
                points.extend([[float(x), float(seg[0][0]), float(seg[1][0])]
                               for x, seg in zip(xs, segments)])
        if points:
            payload.append({"point_lower_upper": points})
        for im in ax.images:
            payload.append({"matrix": np.asarray(im.get_array()).astype(int).tolist()})
    return payload


def original_numerical_artists():
    path = ROOT / "src" / "plot_manuscript_figures_v1.py"
    spec = importlib.util.spec_from_file_location("original_figures_m08_audit", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    captured = {}
    def capture(fig, name):
        captured[name] = numerical_artists(fig)
        plt.close(fig)
    module._save = capture
    module._style()
    module.fig2_development_primary()
    module.fig3_representation_and_zero()
    module.fig4_external_primary()
    return captured


def box(ax, x, y, w, h, title, body, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.015",
                               linewidth=.8, edgecolor=color, facecolor=LIGHT))
    ax.add_patch(FancyBboxPatch((x, y+h-.046), w, .046, boxstyle="round,pad=0.008,rounding_size=0.015",
                               linewidth=0, facecolor=color))
    ax.text(x+w/2, y+h-.023, title, color="white", ha="center", va="center", fontsize=10, weight="bold")
    ax.text(x+w/2, y+(h-.046)/2, body, ha="center", va="center", color=NAVY,
            fontsize=9.5, linespacing=1.12)


def arrow(ax, start, end, color=GRAY):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=11,
                               linewidth=1.2, color=color))


def routed_arrow(ax, points, color=GRAY):
    path = PlotPath(points, [PlotPath.MOVETO] + [PlotPath.LINETO] * (len(points)-1))
    ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=11,
                                linewidth=1.2, color=color))


def figure1():
    cohort = json.loads((ROOT / SOURCE_FILES[0]).read_text(encoding="utf-8"))["cohort_counts"]
    final = model_manifest()
    external = json.loads((ROOT / SOURCE_FILES[3]).read_text(encoding="utf-8"))
    fig = plt.figure(figsize=(WIDTH_MM/25.4, 6.5))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis("off")
    left, right, w, h = .035, .535, .43, .193
    ys = [.765, .535, .305, .075]
    box(ax, left, ys[0], w, h, "Development data",
        f"GDSC2 continuous AUC\n{cohort['development_unique_groups']:,} independent groups\n{cohort['development_drug_model_observations']:,} drug–model records\n8 prespecified compounds", BLUE)
    box(ax, right, ys[0], w, h, "Expression input",
        "CCLE gene TPM\nShared development/external input\nNo expression-platform test", TEAL)
    box(ax, left, ys[1], w, h, "Nested development",
        "5 outer folds × 3 repeats\n3 inner folds; cell-line groups\nAll selection within training", NAVY)
    box(ax, right, ys[1], w, h, "Eight model arms",
        "Tissue + expression 2,000\n2 candidate rules × 3 forms\nExpression / rank / gene pairs", TEAL)
    box(ax, left, ys[2], w, h, "Development readout",
        "Out-of-fold predictions\nSpearman ρ and RMSE\nPaired intervals; zeros retained", BLUE)
    box(ax, right, ys[2], w, h, "Final model audit",
        f"{len(final)} models / {sum(int(x['fit_count']) for x in final):,} fits\nIndependent reload audit\n{sum(int(x['predicted']) for x in external['prediction_coverage']):,} predictions before responses", NAVY)
    box(ax, left, ys[3], .93, h, "Locked CTRPv2 cross-screen evaluation",
        f"{cohort['external_unique_groups']:,} non-overlapping groups; {external['response_audit']['response_valid_drug_model_rows']:,} drug–model responses\nOne-time retrospective ranking; same CCLE expression source\nNot clinical or cross-expression-platform validation", ORANGE)
    arrow(ax, (.25, ys[0]), (.25, ys[1]+h), BLUE)
    routed_arrow(ax, [(.75, ys[0]), (.75, .747), (.405, .747), (.405, ys[1]+h+.003)], TEAL)
    arrow(ax, (.465, ys[1]+h/2), (.535, ys[1]+h/2), NAVY)
    routed_arrow(ax, [(.75, ys[1]), (.75, .517), (.25, .517), (.25, ys[2]+h+.003)], TEAL)
    arrow(ax, (.465, ys[2]+h/2), (.535, ys[2]+h/2), BLUE)
    arrow(ax, (.75, ys[2]), (.75, ys[3]+h), ORANGE)
    return fig


def forest_style(ax, xlim, xlabel):
    ax.axvline(0, color=GRAY, linewidth=1, linestyle="--", zorder=0)
    ax.set_xlabel(xlabel)
    ax.set_xlim(*xlim)
    ax.grid(axis="x", color="#E7ECF0", linewidth=.7)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)


def primary_rows():
    frame = pd.read_csv(ROOT / SOURCE_FILES[1])
    frame["gdsc_drug_id"] = frame["gdsc_drug_id"].astype(str)
    selected = frame[(frame["contrast"] == "primary:association_pair-strong_expression") &
                     (frame["metric"] == "spearman") & (frame["family_name"] == "primary_fixed9")]
    return selected.set_index("gdsc_drug_id").loc[[str(d) for d, _ in DRUGS] + ["macro_equal_fixed_eight_drugs"]]


def draw_primary(ax, points, lower, upper, ys):
    p, lo, hi = np.asarray(points), np.asarray(lower), np.asarray(upper)
    ax.errorbar(p[:-1], ys[:-1], xerr=[p[:-1]-lo[:-1], hi[:-1]-p[:-1]], fmt="o",
                color=BLUE, ecolor=BLUE, capsize=2.8, markersize=5, linewidth=1.4)
    ax.errorbar(p[-1], ys[-1], xerr=[[p[-1]-lo[-1]], [hi[-1]-p[-1]]], fmt="D",
                color=ORANGE, ecolor=ORANGE, capsize=3.5, markersize=6.2, linewidth=1.8)


def figure2():
    rows = primary_rows(); y = np.arange(9)[::-1]
    fig = plt.figure(figsize=(WIDTH_MM/25.4, 4.5))
    ax = fig.add_axes([.35, .19, .61, .76])
    draw_primary(ax, rows["point"].to_numpy(float), rows["family_adjusted_conditional_lower"].to_numpy(float),
                 rows["family_adjusted_conditional_upper"].to_numpy(float), y)
    ax.set_yticks(y, [n for _, n in DRUGS] + ["Equal-drug macro"])
    forest_style(ax, (-.09, .09), "Δ Spearman ρ (REO pair − strong expression)")
    ax.set_xticks([-.08, -.04, 0, .04, .08])
    fig.text(.35, .045, "Family-9 adjusted conditional intervals", color=GRAY, fontsize=9.5)
    return fig


def heatmap_style(ax, matrix):
    cmap = matplotlib.colors.ListedColormap(["#DDE3E8", TEAL])
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(np.arange(-.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def figure3():
    df = pd.read_csv(ROOT / SOURCE_FILES[1]); df["gdsc_drug_id"] = df["gdsc_drug_id"].astype(str)
    rep = df[(df["metric"] == "spearman") & (df["family_name"] == "representation_fixed18")]
    order = [str(d) for d, _ in DRUGS] + ["macro_equal_fixed_eight_drugs"]
    fig = plt.figure(figsize=(WIDTH_MM/25.4, 6.8))
    ax1 = fig.add_axes([.35, .595, .61, .305])
    for i, drug in enumerate(order):
        for rule, offset, color in [("mad:pair_full-expression", .14, BLUE),
                                    ("response_association:pair_full-expression", -.14, TEAL)]:
            row = rep[(rep["gdsc_drug_id"] == drug) & (rep["contrast"] == rule)].iloc[0]
            p, lo, hi = (float(row[x]) for x in ["point", "family_adjusted_conditional_lower", "family_adjusted_conditional_upper"])
            ax1.errorbar(p, 8-i+offset, xerr=[[p-lo], [hi-p]], fmt="o", color=color,
                         ecolor=color, capsize=2, markersize=4.3, linewidth=1.1)
    ax1.set_yticks(np.arange(9)[::-1], [n for _, n in DRUGS] + ["Equal-drug macro"])
    ax1.set_ylim(-.6, 8.6)
    forest_style(ax1, (-.12, .12), "Δ Spearman ρ (pair − expression; same 200 genes)")
    ax1.set_xticks([-.10, -.05, 0, .05, .10])
    ax1.set_title("A  Matched candidate-set contrasts", loc="left", color=NAVY, weight="bold", pad=32)
    ax1.legend(handles=[plt.Line2D([], [], color=BLUE, marker="o", linestyle="", label="MAD"),
                        plt.Line2D([], [], color=TEAL, marker="o", linestyle="", label="Response-associated")],
               loc="lower left", bbox_to_anchor=(0, 1.005), frameon=False, ncol=2, columnspacing=1,
               handletextpad=.4, borderaxespad=0)
    models = ["expression_strong", "mad_expression", "mad_rank", "mad_pair_full",
              "response_association_expression", "response_association_rank", "response_association_pair_full"]
    final = model_manifest(); matrix = []
    for drug, _ in DRUGS:
        rows = {x["task"]["model"]: x["n_nonzero_molecular_coeffs"] for x in final if x["task"]["drug_id"] == drug}
        matrix.append([int(rows[model] > 0) for model in models])
    ax2 = fig.add_axes([.35, .13, .61, .28])
    heatmap_style(ax2, np.asarray(matrix))
    ax2.set_yticks(range(8), [n for _, n in DRUGS])
    ax2.set_xticks(range(7), ["Expr\n2,000", "MAD\nExpr", "MAD\nRank", "MAD\nPair",
                              "Assoc.\nExpr", "Assoc.\nRank", "Assoc.\nPair"])
    ax2.set_title("B  Full-development molecular coefficients", loc="left", color=NAVY, weight="bold", pad=13)
    fig.legend(handles=[Patch(facecolor="#DDE3E8", label="Zero molecular coefficients"),
                        Patch(facecolor=TEAL, label="At least one nonzero coefficient")],
               loc="lower left", bbox_to_anchor=(.345, .005), frameon=False, ncol=1,
               handlelength=1.3, handletextpad=.5, labelspacing=.5, borderaxespad=0)
    return fig


def figure4():
    summary = json.loads((ROOT / SOURCE_FILES[3]).read_text(encoding="utf-8"))
    contrast = summary["evaluation"]["contrasts"]["primary"]["inference"]
    points, lower, upper, ns = [], [], [], []
    for drug, _ in DRUGS:
        row = contrast["per_drug"][str(drug)]; interval = row["delta_interval"]
        points.append(float(row["delta_rho"])); lower.append(float(interval["family_adjusted_lower"]))
        upper.append(float(interval["family_adjusted_upper"])); ns.append(int(row["n_common_groups"]))
    macro = contrast["strict_macro"]["delta_interval"]
    points.append(float(macro["point"])); lower.append(float(macro["family_adjusted_lower"]))
    upper.append(float(macro["family_adjusted_upper"])); ns.append("8 drugs")
    tissues = json.loads((ROOT / SOURCE_FILES[7]).read_text(encoding="utf-8"))["expected_primary_display_tissues"]
    loto = pd.read_csv(ROOT / SOURCE_FILES[6]); loto["drug_id"] = loto["drug_id"].astype(str)
    matrix = loto[loto["heldout_tissue"].isin(tissues)].pivot(index="drug_id", columns="heldout_tissue",
        values="pair_minus_expression_strong_defined").reindex(index=[str(d) for d, _ in DRUGS], columns=tissues)
    defined = matrix.fillna(False).astype(bool).to_numpy(dtype=int)
    labels = ["Esoph./\nStomach" if t == "Esophagus/Stomach" else "Ovary/\nFallopian" if t == "Ovary/Fallopian Tube"
              else "CNS/\nBrain" if t == "CNS/Brain" else t for t in tissues]
    fig = plt.figure(figsize=(WIDTH_MM/25.4, 7.0))
    ax1 = fig.add_axes([.395, .635, .57, .265])
    heatmap_style(ax1, defined)
    ax1.set_yticks(range(8), [f"{n} ({c}/10)" for (_, n), c in zip(DRUGS, defined.sum(axis=1))])
    ax1.set_xticks(range(10), labels, rotation=55, ha="right", rotation_mode="anchor")
    ax1.set_title("A  LOTO contrast definedness", loc="left", color=NAVY, weight="bold", pad=12)
    fig.legend(handles=[Patch(facecolor=TEAL, label="Defined Δρ"), Patch(facecolor="#DDE3E8", label="Undefined Δρ")],
               loc="lower left", bbox_to_anchor=(.39, .475), frameon=False, ncol=2,
               columnspacing=1, handletextpad=.5, handlelength=1.1, borderaxespad=0)
    fig.text(.395, .46, "34/80 defined; strict equal-drug macro: NA", color=GRAY, fontsize=9.5)
    ax2 = fig.add_axes([.35, .125, .50, .255]); y = np.arange(9)[::-1]
    draw_primary(ax2, points, lower, upper, y)
    ax2.set_yticks(y, [n for _, n in DRUGS] + ["Equal-drug macro"])
    forest_style(ax2, (-.34, .34), "Δ Spearman ρ (REO pair − strong expression)")
    ax2.set_xticks([-.3, -.15, 0, .15, .3])
    ax2.set_title("B  Locked CTRPv2 cross-screen contrasts", loc="left", color=NAVY, weight="bold", pad=13)
    for yi, n in zip(y, ns):
        ax2.text(1.045, yi, n if isinstance(n, str) else f"N={n}", transform=ax2.get_yaxis_transform(),
                 va="center", fontsize=9.5, color=GRAY, clip_on=False)
    fig.text(.35, .025, "Family-9 adjusted conditional intervals", color=GRAY, fontsize=9.5)
    return fig


def save(fig, name, baseline):
    fig.canvas.draw()
    text_objects = [x for x in fig.findobj(Text) if x.get_visible() and x.get_text().strip()]
    min_font = min(x.get_fontsize() for x in text_objects)
    if min_font < MIN_FONT_PT:
        raise ValueError(f"{name}: {min_font} pt below {MIN_FONT_PT} pt design minimum")
    renderer = fig.canvas.get_renderer()
    bounds = fig.bbox
    clipped = []
    for item in text_objects:
        extent = item.get_window_extent(renderer)
        if (extent.x0 < bounds.x0-1 or extent.y0 < bounds.y0-1 or
                extent.x1 > bounds.x1+1 or extent.y1 > bounds.y1+1):
            clipped.append(item.get_text())
    if clipped:
        raise ValueError(f"{name}: visible text outside figure bounds: {clipped}")
    numeric = numerical_artists(fig)
    equal = name not in baseline or numeric == baseline[name]
    if not equal:
        raise AssertionError(f"{name}: plotted numerical values differ from original figures")
    svg = OUT / f"{name}.svg"
    fig.savefig(svg, metadata={"Date": None})
    svg.write_text("\n".join(x.rstrip() for x in svg.read_text(encoding="utf-8").splitlines())+"\n",
                   encoding="utf-8", newline="\n")
    fig.savefig(OUT / f"{name}.png", dpi=600)
    # A display proof at approximately 150 dpi. Final source PNG remains 600 dpi.
    fig.savefig(OUT / f"{name}_proof.png", dpi=150)
    result = {"title_in_caption_only": TITLES[name], "width_mm": WIDTH_MM,
        "height_mm": float(fig.get_figheight()*25.4), "minimum_font_pt_at_170mm": min_font,
        "minimum_font_pt_at_word_590px": min_font*WORD_WIDTH_MM/WIDTH_MM,
        "word_590px_width_mm": WORD_WIDTH_MM, "png_dpi": 600,
        "same_numerical_artists_as_original": equal,
        "numerical_artists": numeric, "visible_text_count": len(text_objects),
        "visible_text_outside_canvas": clipped}
    if name == "Figure3_representation_and_zero_solutions":
        matrix = np.asarray(numeric[1]["matrix"])
        result["full_development_molecular_zero_models"] = int((matrix == 0).sum())
        result["full_development_molecular_models"] = int(matrix.size)
    if name == "Figure4_loto_and_external":
        matrix = np.asarray(numeric[0]["matrix"])
        result["loto_defined_contrasts"] = int(matrix.sum())
        result["loto_display_contrasts"] = int(matrix.size)
        result["strict_loto_macro"] = "NA"
    plt.close(fig)
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    hashes_before = {p: sha256(ROOT / p) for p in SOURCE_FILES}
    baseline = original_numerical_artists(); style()
    figures = {}
    for function, name in [(figure1, "Figure1_study_design"), (figure2, "Figure2_development_primary"),
                           (figure3, "Figure3_representation_and_zero_solutions"), (figure4, "Figure4_loto_and_external")]:
        figures[name] = save(function(), name, baseline)
    hashes_after = {p: sha256(ROOT / p) for p in SOURCE_FILES}
    if hashes_before != hashes_after:
        raise AssertionError("Frozen aggregate inputs changed while rendering")
    manifest = {"schema_version": "2.0", "status": "M08_layout_revision_numerical_artists_exactly_preserved",
        "generated_utc": datetime.now(timezone.utc).isoformat(), "script": "src/plot_manuscript_figures_m08_v1.py",
        "script_sha256": sha256(Path(__file__).resolve()), "source_sha256": hashes_before,
        "original_script_sha256": sha256(ROOT / "src/plot_manuscript_figures_v1.py"),
        "original_assets_manifest_sha256": sha256(ROOT / "figures/manuscript_v1/figure_manifest.json"),
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__, "individual_response_values_read": False,
        "sample_level_predictions_read": False, "models_loaded": False, "training_or_scoring_performed": False,
        "font_minimum_is_project_legibility_target_not_journal_requirement": True,
        "layout_policy": "No overall titles inside figures; panel headings only. Figures 3 and 4 use vertical panels.",
        "na_policy": "LOTO undefined cells remain gray; no missing result converted to zero. Matrix encodes definedness only.",
        "figures": figures,
        "output_sha256": {p.name: sha256(p) for p in sorted(OUT.iterdir()) if p.suffix in {".png", ".svg"}},
    }
    (OUT / "figure_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"figures": len(figures), "all_plotted_numbers_exactly_preserved": True,
                      "minimum_font_at_word_width_pt": min(x["minimum_font_pt_at_word_590px"] for x in figures.values()),
                      "output": str(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
