"""Display existing aggregate LOTO effects and locked external contrasts.

This is a presentation-only revision. No model, individual response, sample-level
prediction, training, scoring, or resampling is accessed or performed. The input
tables and earlier figures remain unchanged. Numeric audit records are emitted
beside the new figure so that NA cells cannot be mistaken for zero effects.
"""
from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.container import ErrorbarContainer
from matplotlib.patches import Patch
from matplotlib.text import Text


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "manuscript_scope_diagnostics_v1"
NAME = "Figure4_loto_and_external"
WIDTH_MM, HEIGHT_MM = 170.0, 180.0
MIN_FONT_PT = 9.5
NAVY, BLUE, ORANGE, GRAY = "#17324D", "#2878B5", "#E58B2A", "#647382"
NA_COLOR = "#E0E4E8"
DRUGS = [(1931, "MIRA-1"), (1036, "PLX-4720"), (1191, "Bortezomib"),
         (1930, "Telomerase Inhibitor IX"), (1199, "Tamoxifen"),
         (1091, "BMS-536924"), (1010, "Gefitinib"), (1922, "Cediranib")]
SOURCE_FILES = [
    "results/tissue_loto/contrasts_by_drug_tissue.csv",
    "config/tissue_sensitivity_loto_v1.json",
    "results/p6_external/evaluation_summary.json",
    "figures/manuscript_m08_v1/Figure4_loto_and_external.png",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def text_number(value: float) -> str:
    """Three decimal places; suppress leading zero only in narrow heatmap cells."""
    rounded = f"{abs(value):.3f}"
    if rounded == "0.000":
        return "0" if value == 0 else ("−<.001" if value < 0 else "<.001")
    return ("−" if value < 0 else "") + rounded.removeprefix("0")


def load_aggregates():
    config = json.loads((ROOT / SOURCE_FILES[1]).read_text(encoding="utf-8"))
    tissues = config["expected_primary_display_tissues"]
    source = pd.read_csv(ROOT / SOURCE_FILES[0], float_precision="round_trip")
    selected = source[source["heldout_tissue"].isin(tissues)].copy()
    assert len(selected) == 80
    assert not selected.duplicated(["drug_id", "heldout_tissue"]).any()
    index = [d for d, _ in DRUGS]
    matrix = selected.pivot(index="drug_id", columns="heldout_tissue",
                            values="pair_minus_expression_strong_delta_rho").reindex(index=index, columns=tissues)
    defined = selected.pivot(index="drug_id", columns="heldout_tissue",
                             values="pair_minus_expression_strong_defined").reindex(index=index, columns=tissues)
    ns = selected.pivot(index="drug_id", columns="heldout_tissue", values="n_samples").reindex(index=index, columns=tissues)
    values = matrix.to_numpy(dtype=float)
    mask = defined.to_numpy(dtype=bool)
    assert np.array_equal(np.isfinite(values), mask)
    assert int(mask.sum()) == 34 and int((~mask).sum()) == 46
    assert int(ns.min().min()) >= 20

    summary = json.loads((ROOT / SOURCE_FILES[2]).read_text(encoding="utf-8"))
    contrast = summary["evaluation"]["contrasts"]["primary"]["inference"]
    external = []
    for drug, name in DRUGS:
        row = contrast["per_drug"][str(drug)]
        interval = row["delta_interval"]
        external.append({"drug_id": drug, "compound": name,
                         "point": float(row["delta_rho"]),
                         "lower": float(interval["family_adjusted_lower"]),
                         "upper": float(interval["family_adjusted_upper"]),
                         "n_common_groups": int(row["n_common_groups"])})
    macro = contrast["strict_macro"]["delta_interval"]
    external.append({"drug_id": "macro", "compound": "Equal-drug macro",
                     "point": float(macro["point"]),
                     "lower": float(macro["family_adjusted_lower"]),
                     "upper": float(macro["family_adjusted_upper"]),
                     "n_common_groups": None})
    exact_zeros = [row["drug_id"] for row in external[:-1]
                   if row["point"] == row["lower"] == row["upper"] == 0.0]
    assert exact_zeros == [1931, 1036, 1199]
    return tissues, values, mask, ns.to_numpy(dtype=int), external


def draw_figure(tissues, values, mask, counts, external):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": MIN_FONT_PT,
        "axes.titlesize": 11.5, "axes.labelsize": MIN_FONT_PT,
        "xtick.labelsize": MIN_FONT_PT, "ytick.labelsize": MIN_FONT_PT,
        "legend.fontsize": MIN_FONT_PT, "figure.facecolor": "white",
        "axes.facecolor": "white", "axes.edgecolor": "#CBD3DA",
        "axes.linewidth": .8, "savefig.facecolor": "white",
        "svg.fonttype": "none", "pdf.fonttype": 42})
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    fig.text(.035, .970, "A  Within-tissue contrasts after tissue holdout",
             color=NAVY, weight="bold", fontsize=11)
    fig.text(.035, .942, "Cell values: Δρ; row counts: defined contrasts / 10", color=GRAY)
    ax1 = fig.add_axes([.285, .590, .690, .330])
    cmap = LinearSegmentedColormap.from_list("signed_effect", ["#C56B30", "#FCFCFC", "#347FA3"])
    cmap.set_bad(NA_COLOR)
    norm = TwoSlopeNorm(vmin=-.5, vcenter=0, vmax=.5)
    im = ax1.imshow(np.ma.masked_invalid(values), cmap=cmap, norm=norm,
                    interpolation="nearest", aspect="auto")
    ax1.set_xticks(np.arange(-.5, 10, 1), minor=True)
    ax1.set_yticks(np.arange(-.5, 8, 1), minor=True)
    ax1.grid(which="minor", color="white", linewidth=1.5)
    ax1.tick_params(which="both", length=0)
    for spine in ax1.spines.values():
        spine.set_visible(False)
    row_labels = []
    for (_, name), n in zip(DRUGS, mask.sum(axis=1)):
        label = "Telomerase\nInhibitor IX" if name == "Telomerase Inhibitor IX" else name
        row_labels.append(f"{label} ({n}/10)")
    ax1.set_yticks(range(8), row_labels)
    labels = []
    for j, tissue in enumerate(tissues):
        label = {"Esophagus/Stomach": "Esoph./\nStomach",
                 "Ovary/Fallopian Tube": "Ovary/\nFallopian", "CNS/Brain": "CNS/Brain"}.get(tissue, tissue)
        labels.append(label)
    ax1.set_xticks(range(10), labels, rotation=55, ha="right", rotation_mode="anchor")
    for i in range(8):
        for j in range(10):
            if mask[i, j]:
                ax1.text(j, i, text_number(values[i, j]), ha="center", va="center",
                         color="white" if abs(values[i, j]) > .36 else NAVY, fontsize=MIN_FONT_PT)
            else:
                ax1.text(j, i, "×", ha="center", va="center", color="#7C8791", fontsize=12)

    cax = fig.add_axes([.285, .455, .415, .014])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[-.5, -.25, 0, .25, .5])
    cb.ax.set_xticklabels(["−0.50", "−0.25", "0", "0.25", "0.50"])
    cb.set_label("Δ Spearman ρ (pair − expression 2,000)", labelpad=3)
    fig.legend(handles=[Patch(facecolor=NA_COLOR, label="×  Undefined")], loc="center left",
               bbox_to_anchor=(.720, .463), frameon=False, borderaxespad=0,
               handlelength=1.0, handletextpad=.4)
    fig.text(.285, .378, "34/80 defined; strict eight-drug macro: NA", color=GRAY)

    fig.text(.035, .337, "B  Locked CTRP v2 cross-screen contrasts",
             color=NAVY, weight="bold", fontsize=11)
    ax2 = fig.add_axes([.335, .110, .515, .205])
    for i, row in enumerate(external):
        p, lo, hi = row["point"], row["lower"], row["upper"]
        is_macro = i == 8
        ax2.errorbar(p, 8 - i, xerr=[[p - lo], [hi - p]], fmt="D" if is_macro else "o",
                     color=ORANGE if is_macro else BLUE, capsize=3.5 if is_macro else 2.8,
                     markersize=6.2 if is_macro else 5, linewidth=1.8 if is_macro else 1.4)
    ax2.axvline(0, color=GRAY, linewidth=1, linestyle="--", zorder=0)
    ax2.set(xlim=(-.34, .34), ylim=(-.55, 8.55))
    ax2.grid(axis="x", color="#E7ECF0", linewidth=.7)
    ax2.spines[["top", "right", "left"]].set_visible(False)
    ax2.tick_params(axis="y", length=0)
    ax2.set_yticks(np.arange(9)[::-1], [n for _, n in DRUGS] + ["Equal-drug macro"])
    ax2.set_xticks([-.3, -.15, 0, .15, .3], ["−0.30", "−0.15", "0", "0.15", "0.30"])
    ax2.set_xlabel("Δ Spearman ρ (pair − expression 2,000)", labelpad=4)
    for i, row in enumerate(external):
        ntext = "8 drugs" if i == 8 else f"N={row['n_common_groups']}"
        ax2.text(1.055, 8 - i, ntext, transform=ax2.get_yaxis_transform(), va="center",
                 color=GRAY, clip_on=False)
    fig.text(.335, .018, "Family-9 adjusted conditional intervals", color=GRAY)
    return fig, ax1, ax2


def audit_and_save(fig, ax1, ax2, tissues, values, mask, counts, external, hashes):
    renderer = fig.canvas.get_renderer()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [item for item in fig.findobj(Text) if item.get_visible() and item.get_text().strip()]
    minimum_font = min(item.get_fontsize() for item in texts)
    assert minimum_font >= MIN_FONT_PT
    clipped = []
    for item in texts:
        bounds = item.get_window_extent(renderer)
        if bounds.x0 < -1 or bounds.y0 < -1 or bounds.x1 > fig.bbox.x1 + 1 or bounds.y1 > fig.bbox.y1 + 1:
            clipped.append(item.get_text())
    if clipped:
        raise AssertionError(f"Visible text outside canvas: {clipped}")
    plotted = np.ma.asarray(ax1.images[0].get_array())
    assert np.array_equal(np.ma.getmaskarray(plotted), ~mask)
    assert np.array_equal(plotted.data[mask], values[mask])
    artist_values = []
    for item in ax2.containers:
        if isinstance(item, ErrorbarContainer):
            point = float(np.asarray(item.lines[0].get_xdata()).ravel()[0])
            segment = item.lines[2][0].get_segments()[0]
            artist_values.append([point, float(segment[0, 0]), float(segment[1, 0])])
    expected = [[row["point"], row["lower"], row["upper"]] for row in external]
    assert np.allclose(artist_values, expected, rtol=0, atol=1e-16)
    assert {name: sha256(ROOT / name) for name in SOURCE_FILES} == hashes
    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf", "svg"):
        kwargs = {"dpi": 600} if extension == "png" else {"metadata": {"Date": None}} if extension == "svg" else {}
        fig.savefig(OUT / f"{NAME}.{extension}", **kwargs)
    fig.savefig(OUT / f"{NAME}_proof.png", dpi=180)
    rows = []
    for i, (drug, compound) in enumerate(DRUGS):
        for j, tissue in enumerate(tissues):
            rows.append({"drug_id": drug, "compound": compound, "heldout_tissue": tissue,
                         "n_groups": int(counts[i, j]), "delta_rho_defined": bool(mask[i, j]),
                         "delta_rho": float(values[i, j]) if mask[i, j] else None,
                         "display_text": text_number(values[i, j]) if mask[i, j] else "×"})
    pd.DataFrame(rows).to_csv(OUT / "figure4_loto_source_value_audit.csv", index=False, na_rep="NA")
    pd.DataFrame(external).to_csv(OUT / "figure4_external_source_value_audit.csv", index=False, na_rep="NA")
    manifest = {"schema_version": "1.0", "status": "presentation_only_existing_aggregate_values",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script": "src/plot_scope_diagnostics_figure4_v1.py", "script_sha256": sha256(Path(__file__)),
        "source_sha256": hashes, "width_mm": WIDTH_MM, "height_mm": HEIGHT_MM,
        "figure_only_height_at_170_mm": HEIGHT_MM,
        "height_remaining_for_caption_under_225_mm_limit": 225.0 - HEIGHT_MM,
        "tissue_sample_counts_location": "source-value audit and manuscript Table S1; omitted from rotated labels to preserve readability",
        "minimum_font_pt_at_170_mm": minimum_font, "minimum_font_pt_at_156_mm": minimum_font * 156 / 170,
        "project_legibility_target_not_journal_requirement": True,
        "loto_display_cells": 80, "loto_defined_cells": int(mask.sum()), "loto_undefined_cells": int((~mask).sum()),
        "loto_defined_per_drug": {name: int(n) for (_, name), n in zip(DRUGS, mask.sum(axis=1))},
        "strict_loto_macro": "NA", "undefined_cells_not_replaced_by_zero": True,
        "heatmap_artist_values_match_source_exactly": True,
        "external_points_and_intervals_preserved": True,
        "external_artist_max_absolute_error": float(np.max(np.abs(np.asarray(artist_values) - np.asarray(expected)))),
        "external_exact_zero_compound_ids": [1931, 1036, 1199],
        "visible_text_outside_canvas": clipped,
        "individual_response_values_read": False, "sample_level_predictions_read": False,
        "models_loaded": False, "training_or_scoring_performed": False, "resampling_performed": False,
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "output_sha256": {path.name: sha256(path) for path in sorted(OUT.iterdir())
                          if path.suffix in {".png", ".pdf", ".svg", ".csv"}},
    }
    (OUT / "figure4_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUT), "defined": 34, "NA": 46,
                      "minimum_font_at_156_mm_pt": minimum_font * 156 / 170,
                      "external_max_artist_error": manifest["external_artist_max_absolute_error"]}, ensure_ascii=False))
    plt.close(fig)


def main():
    hashes = {name: sha256(ROOT / name) for name in SOURCE_FILES}
    tissues, values, mask, counts, external = load_aggregates()
    fig, ax1, ax2 = draw_figure(tissues, values, mask, counts, external)
    audit_and_save(fig, ax1, ax2, tissues, values, mask, counts, external, hashes)


if __name__ == "__main__":
    main()
