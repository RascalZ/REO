"""Demonstrate actual frozen zero/null policy on generated, label-free input."""
from __future__ import annotations
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from repro_inference_v1 import load_bundle, load_json_response_free, predict_bundle


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_new(path, obj):
    with Path(path).open("x", encoding="utf-8", newline="\n") as f:
        json.dump(obj, f, ensure_ascii=False, allow_nan=False, indent=2)
        f.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "reproducibility_v1/local_assets/inference_examples/zero_state")
    args = parser.parse_args()
    base = ROOT / "reproducibility_v1/local_assets"
    out = args.out.resolve()
    if not out.is_relative_to(base.resolve()):
        raise ValueError("Actual fitted-model demonstration remains local pending release review")
    if out.exists():
        raise FileExistsError("Refusing to overwrite an existing demonstration")
    paths = {"zero_molecular": base / "models/d1931_response_association_pair_full.json.gz", "tissue_reference": base / "models/d1931_tissue_null.json.gz"}
    expected = json.loads((base / "model_export_manifest.json").read_text(encoding="utf-8"))
    expected_sha = {r["path"]: r["portable_sha256"] for r in expected["models"]}
    for p in paths.values():
        if sha(p) != expected_sha[p.relative_to(base).as_posix()]:
            raise ValueError("Portable source hash mismatch")
    bundles = {key: load_bundle(path) for key, path in paths.items()}
    if not bundles["zero_molecular"]["declared_molecular_arm"] or bundles["zero_molecular"]["n_nonzero_molecular_coeffs"] != 0:
        raise ValueError("Expected existing zero molecular solution")
    sample = load_json_response_free(base / "synthetic_inputs/full_gene_order.synthetic.json")
    normal = sample["raw_tpm"][0]
    selected = set(bundles["zero_molecular"]["candidate_gene_indices"])
    unselected = next(i for i in range(len(normal)) if i not in selected)
    rows = [list(normal) for _ in range(3)]
    rows[1][unselected] = None
    sample = {"schema_version": "reo-inference-input-1.0", "sample_ids": ["SIM_VALID", "SIM_MISSING_UNSELECTED_GENE", "SIM_MISSING_TISSUE"], "biological_group_ids": ["SIM_GROUP_VALID", "SIM_GROUP_MISSING_GENE", "SIM_GROUP_MISSING_TISSUE"], "tissues": ["Lung", "Lung", None], "gene_ids": sample["gene_ids"], "raw_tpm": rows}
    results = {key: predict_bundle(bundle, sample) for key, bundle in bundles.items()}
    if results["zero_molecular"]["statuses"] != ["predicted", "abstained", "abstained"] or results["tissue_reference"]["statuses"] != ["predicted", "predicted", "abstained"]:
        raise ValueError("Zero/null abstention contract failed")
    if results["zero_molecular"]["predictions"][0] != results["tissue_reference"]["predictions"][0]:
        raise ValueError("Valid zero molecular solution differs from same-policy tissue reference")
    out.mkdir(parents=True)
    save_new(out / "input_generated_only.json", sample)
    for key, result in results.items(): save_new(out / (key + "_output.json"), result)
    manifest = {"scope": "existing frozen models, generated numerical input only; not a predictive performance result", "new_training_or_scoring": False, "source_CTRP_labels_indices_or_individual_predictions_read": False, "invalid_unselected_gene_index": unselected, "source_model_sha256": {key: sha(path) for key, path in paths.items()}, "files": {p.name: sha(p) for p in out.iterdir() if p.is_file()}, "statuses": {key: value["statuses"] for key, value in results.items()}}
    save_new(out / "demo_manifest.json", manifest)
    print(json.dumps({"status": "passed", "case_count": 3, "zero_molecular_statuses": results["zero_molecular"]["statuses"], "tissue_reference_statuses": results["tissue_reference"]["statuses"], "labels_used": False}))


if __name__ == "__main__":
    main()
