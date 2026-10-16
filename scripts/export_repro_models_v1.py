"""Export trusted, frozen final models and check response-free numerical parity.

No fitting, real expression row, response label, original sentinel or scoring
function is used. Complete derived assets remain local pending release terms.
"""
from __future__ import annotations

import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[_key] = "1"
import sys
sys.dont_write_bytecode = True
import argparse
import csv
import gc
import gzip
import hashlib
import importlib.metadata
import io
import json
from pathlib import Path
import platform
import time

ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import joblib
import numpy as np
from threadpoolctl import threadpool_limits
import ordinal_diagnostic_models as diagnostic
import modeling.nested_cv as nested
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from repro_inference_v1 import predict_bundle

RUN = ROOT / "runs/gdsc2-cancer-v2-fixed8-full-development-final-v1_1"
FIRST = ROOT / "config/external_validation_v1_1.lock.json"
SECOND = ROOT / "config/external_execution_v1_1.lock.json"
MANIFEST = ROOT / "results/p6_preflight/final_model_manifest.jsonl"
CONTEXT = RUN / "worker_context.json"
METADATA = RUN / "final_model_metadata.json"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _canonical(obj):
    return (json.dumps(obj, ensure_ascii=False, allow_nan=False, separators=(",", ":")) + "\n").encode("utf-8")


def _write_new(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as f:
        f.write(data)


def _forbidden_fit(*args, **kwargs):
    raise RuntimeError("Training/scoring is forbidden during artifact export")


def _guard(allowed, out, summary):
    allowed = {Path(p).absolute() for p in allowed}
    opened, blocked = set(), []
    def hook(event, args):
        if event in {"socket.connect", "socket.getaddrinfo", "subprocess.Popen"}:
            blocked.append(event)
            raise PermissionError("Network/subprocess disabled in response-free export")
        if event != "open" or not args or not isinstance(args[0], (str, bytes, os.PathLike)):
            return
        p = Path(os.fsdecode(args[0])).absolute()
        mode = args[1] if len(args) > 1 else None
        flags = args[2] if len(args) > 2 and isinstance(args[2], int) else 0
        writing = (isinstance(mode, str) and any(x in mode for x in "wax+")) or bool(flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND))
        if writing:
            ok = p.is_relative_to(out) or p == summary
        else:
            # Libraries were imported before installation. Remaining project
            # imports may read code, but cannot open an unlisted data asset.
            ok = p in allowed or p.is_relative_to(out) or p.suffix.lower() in {".py", ".pyc", ".pyd", ".dll"}
        if not ok:
            blocked.append(str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else "outside_allowed_assets")
            raise PermissionError("Unlisted file access disabled during response-free export")
        if p in allowed:
            opened.add(p.relative_to(ROOT).as_posix())
    sys.addaudithook(hook)
    return opened, blocked


def _bundle(artifact, manifest_row):
    m = artifact["model"]
    if type(m) is not diagnostic.DiagnosticRegressor:
        raise ValueError("Unexpected native model class")
    stable = list(artifact["stable_gene_ids"])
    if len(stable) != 19836 or len(set(stable)) != len(stable) or m.n_features_in_ != 19838:
        raise ValueError("Unexpected frozen full-expression schema")
    b = m.builder_
    genes = [] if b is None else list(b.gene_ids_)
    features = [] if b is None else b.feature_names_.tolist()
    p = len(features)
    beta = np.asarray(m.regressor_.coef_, dtype=np.float64) if hasattr(m, "regressor_") else np.zeros(p)
    if p != len(beta):
        raise ValueError("Feature/coefficient dimension mismatch")
    active = int(np.count_nonzero(beta))
    if active != manifest_row["n_nonzero_molecular_coeffs"]:
        raise ValueError("Native coefficient state differs from frozen manifest")
    scaled = hasattr(m, "scaler_")
    if scaled != (m.fitted_kind_ in {"expression", "rank", "pair_full"}):
        raise ValueError("Unexpected preprocessor state")
    payload = {
        "schema_version": "reo-portable-model-1.0", "task": artifact["task"],
        "source": {"scientific_source_commit": "ca0f837ac10665ab8d4d6ab97691e58e5bed4577", "native_artifact_path": "runs/gdsc2-cancer-v2-fixed8-full-development-final-v1_1/" + manifest_row["artifact_path"], "native_artifact_sha256": manifest_row["artifact_sha256"], "first_lock_sha256": artifact["freeze_sha256"], "expression_genes_sha256": artifact["expression_genes_sha256"], "expression_tpm_sha256": artifact["expression_tpm_sha256"], "scope": "existing final model, derived parameters only; no training/test labels or original predictions"},
        "input_gene_ids": stable, "tissue_code": artifact["tissue_code"],
        "model_kind": m.model_kind, "fitted_kind": m.fitted_kind_, "selection_rule": m.selection_rule,
        "candidate_gene_indices": [] if b is None else b.gene_indices_.tolist(), "candidate_gene_ids": genes,
        "feature_names": features, "pair_indices": [] if b is None else b.pair_indices_.tolist(),
        "feature_mean": m.scaler_.mean_.tolist() if scaled else [], "feature_scale": m.scaler_.scale_.tolist() if scaled else [],
        "feature_tissue_means": m.feature_tissue_means_.tolist() if scaled else [[] for _ in m.training_tissue_codes_],
        "target_tissue_means": m.target_tissue_means_.tolist(), "training_tissue_codes": m.training_tissue_codes_,
        "target_mean": m.target_mean_, "target_sd": m.target_sd_, "coefficients": beta.tolist(),
        "n_nonzero_molecular_coeffs": active, "declared_molecular_arm": artifact["task"]["model"] != "tissue_null",
        "fallback_reason": m.fallback_reason_, "solver_reason": m.solver_reason_, "chosen_params": manifest_row["chosen_params"],
    }
    return payload


def _synthetic_input(genes):
    rng = np.random.Generator(np.random.PCG64(2026101601))
    tpm = np.exp2(rng.uniform(0, 7, (6, len(genes)))) - 1.0
    tpm[1] = np.floor(tpm[1] / 2.0) * 2.0
    tpm[2] = 0.0
    tpm[3] = tpm[0]
    tpm[4] = np.arange(len(genes), dtype=np.float64) % 13
    tpm[5] = np.arange(len(genes), dtype=np.float64) * 1e-12
    return {"schema_version": "reo-inference-input-1.0", "sample_ids": [f"SYNTHETIC_CASE_{i}" for i in range(6)], "biological_group_ids": [f"SYNTHETIC_GROUP_{i}" for i in range(6)], "tissues": ["Lung", "Skin", "Lung", "SYNTHETIC_UNSEEN_TISSUE", "Skin", "Lung"], "gene_ids": genes, "raw_tpm": tpm.tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "reproducibility_v1/local_assets")
    parser.add_argument("--summary", type=Path, default=ROOT / "docs/reproducibility_assets_v1/model_export_public_summary.json")
    args = parser.parse_args()
    out, summary = args.out.absolute(), args.summary.absolute()
    if not out.is_relative_to(ROOT / "reproducibility_v1/local_assets"):
        raise ValueError("Full assets must remain in the ignored local-assets namespace")
    if (out / "model_export_manifest.json").exists() or (out / "models").exists() or summary.exists():
        raise FileExistsError("Refusing to overwrite an export; use a new local output and summary path")
    # Setup reads only two published locks and one published model manifest.
    first, second = _read_json(FIRST), _read_json(SECOND)
    rows = [json.loads(line) for line in MANIFEST.read_text(encoding="utf-8").splitlines() if line]
    if len(rows) != 64 or second["model_count"] != 64 or second["first_lock_sha256"] != sha(FIRST) or second["final_model_manifest_sha256"] != sha(MANIFEST):
        raise ValueError("Frozen final-model manifest/lock binding mismatch")
    environment = {"python": platform.python_version(), "packages": {k: importlib.metadata.version(k) for k in first["environment"]["packages"]}}
    if environment != first["environment"]:
        raise ValueError("Use the exact recorded final-model environment for native artifact export")
    artifacts = [RUN / row["artifact_path"] for row in rows]
    source_files = [ROOT / name for name in first["execution_files"]]
    opened, blocked = _guard([FIRST, SECOND, MANIFEST, CONTEXT, METADATA, *artifacts, *source_files], out, summary)
    for p in source_files:
        if sha(p) != first["files"][p.relative_to(ROOT).as_posix()]:
            raise ValueError("Frozen scientific source changed: " + p.name)
    context = _read_json(CONTEXT)
    metadata = _read_json(METADATA)
    if sha(METADATA) != second["final_model_metadata_sha256"] or sha(CONTEXT) != metadata["prepared_input_hashes"]["worker_context.json"] or context["freeze_sha256"] != sha(FIRST):
        raise ValueError("Final worker context lock mismatch")
    # Explicitly prevent fitting/scoring even if a future export edit calls one.
    for cls in (diagnostic.DiagnosticRegressor, diagnostic.DiagnosticFeatureBuilder, nested.TrainOnlyRegressor, ElasticNet, StandardScaler, OneHotEncoder):
        cls.fit = _forbidden_fit
    diagnostic.DiagnosticRegressor.score = _forbidden_fit
    start = time.perf_counter()
    output_rows, checks = [], []
    tables = out / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    candidate_path, coeff_path, tissue_path = tables / "candidate_genes.csv", tables / "all_coefficients.csv.gz", tables / "tissue_effects.csv"
    for p in [candidate_path, coeff_path, tissue_path]:
        if p.exists():
            raise FileExistsError(p.name)
    candidate_count = coefficient_count = tissue_count = 0
    candidate_fields = ["drug_id", "model", "candidate_position", "original_gene_index", "stable_gene_id", "selection_score", "signed_selection_association"]
    coefficient_fields = ["drug_id", "model", "feature_index", "feature_name", "coefficient_standardized", "coefficient_original_response_per_encoded_feature", "scaler_mean", "scaler_scale", "zero_coefficient"]
    tissue_fields = ["drug_id", "model", "tissue_code", "tissue_label", "standardized_target_tissue_mean", "tissue_reference_original_scale"]
    with candidate_path.open("x", encoding="utf-8", newline="") as cg, coeff_path.open("xb") as cb, tissue_path.open("x", encoding="utf-8", newline="") as th:
        with gzip.GzipFile(filename="", mode="wb", fileobj=cb, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as cf:
                gw, cw, tw = (csv.DictWriter(f, fieldnames=fields, lineterminator="\n") for f, fields in [(cg, candidate_fields), (cf, coefficient_fields), (th, tissue_fields)])
                for w in [gw, cw, tw]: w.writeheader()
                sample, full_ids = None, None
                for i, (row, path) in enumerate(zip(rows, artifacts)):
                    relative = path.relative_to(ROOT).as_posix()
                    expected = row["artifact_sha256"]
                    if sha(path) != expected or second["model_artifact_hashes"].get(relative) != expected:
                        raise ValueError("Trusted model hash mismatch")
                    artifact = joblib.load(path)
                    if artifact["task"] != row["task"] or artifact["freeze_sha256"] != sha(FIRST) or artifact["tissue_code"] != context["tissue_code"]:
                        raise ValueError("Native artifact context mismatch")
                    payload = _bundle(artifact, row)
                    if full_ids is None:
                        full_ids = payload["input_gene_ids"]
                        sample = _synthetic_input(full_ids)
                        _write_new(out / "synthetic_inputs/full_gene_order.synthetic.json", _canonical(sample))
                        _write_new(out / "schemas/input_gene_ids.json", _canonical(full_ids))
                        _write_new(out / "schemas/tissue_code.json", _canonical(payload["tissue_code"]))
                    if payload["input_gene_ids"] != full_ids:
                        raise ValueError("Input gene order changed between models")
                    # Round-trip all float64 parameters through canonical JSON,
                    # not the live model object, before numerical checking.
                    raw = _canonical(payload)
                    portable = json.loads(raw)
                    with threadpool_limits(limits=1):
                        result = predict_bundle(portable, sample)
                        tpm = np.asarray(sample["raw_tpm"], dtype=np.float64)
                        codes = np.array([artifact["tissue_code"].get(t, max(artifact["tissue_code"].values()) + 1) for t in sample["tissues"]], dtype=np.float64)
                        native = artifact["model"].predict(np.column_stack((tpm, codes, np.arange(len(tpm), dtype=np.float64) + 100000)))
                    predicted = np.asarray(result["predictions"], dtype=np.float64)
                    error = float(np.max(np.abs(predicted - native)))
                    if result["statuses"] != ["predicted"] * len(tpm) or not np.allclose(predicted, native, atol=1e-12, rtol=1e-12):
                        raise ValueError("Portable/native inference disagreement")
                    b, m, task = artifact["model"].builder_, artifact["model"], row["task"]
                    base = {"drug_id": task["drug_id"], "model": task["model"]}
                    if b is not None:
                        for j, (index, gid) in enumerate(zip(payload["candidate_gene_indices"], payload["candidate_gene_ids"])):
                            gw.writerow({**base, "candidate_position": j, "original_gene_index": index, "stable_gene_id": gid, "selection_score": float(b.selected_scores_[j]), "signed_selection_association": float(b.selected_signed_associations_[j])}); candidate_count += 1
                    for j, feature in enumerate(payload["feature_names"]):
                        beta = payload["coefficients"][j]
                        cw.writerow({**base, "feature_index": j, "feature_name": feature, "coefficient_standardized": beta, "coefficient_original_response_per_encoded_feature": payload["target_sd"] * beta / payload["feature_scale"][j], "scaler_mean": payload["feature_mean"][j], "scaler_scale": payload["feature_scale"][j], "zero_coefficient": beta == 0.0}); coefficient_count += 1
                    labels = {v: k for k, v in payload["tissue_code"].items()}
                    for code, a in zip(payload["training_tissue_codes"], payload["target_tissue_means"]):
                        tw.writerow({**base, "tissue_code": code, "tissue_label": labels[code], "standardized_target_tissue_mean": a, "tissue_reference_original_scale": payload["target_mean"] + payload["target_sd"] * a}); tissue_count += 1
                    compressed_bytes = gzip.compress(raw, mtime=0)
                    target = out / "models" / (f"d{task['drug_id']}_{task['model']}.json.gz")
                    _write_new(target, compressed_bytes)
                    record = {**base, "source_artifact_sha256": expected, "path": target.relative_to(out).as_posix(), "portable_sha256": hashlib.sha256(compressed_bytes).hexdigest(), "canonical_json_sha256": hashlib.sha256(raw).hexdigest(), "portable_bytes": len(compressed_bytes), "candidate_genes": len(payload["candidate_gene_ids"]), "features": len(payload["feature_names"]), "nonzero_coefficients": payload["n_nonzero_molecular_coeffs"], "declared_molecular_arm": payload["declared_molecular_arm"], "zero_molecular_solution": payload["n_nonzero_molecular_coeffs"] == 0 if payload["declared_molecular_arm"] else None, "fallback_reason": payload["fallback_reason"], "synthetic_case_count": len(tpm), "max_absolute_native_portable_difference": error}
                    output_rows.append(record)
                    checks.append({**base, "max_absolute_difference": error, "within_frozen_numeric_tolerance": True, "unknown_tissue_cases": sum(result["unknown_tissue"]), "synthetic_input_only": True})
                    if sha(path) != expected:
                        raise ValueError("Original model artifact changed during read-only export")
                    del artifact, payload, portable, raw, compressed_bytes, m, b
                    gc.collect()
                    print(json.dumps({"completed_models": i + 1, "total_models": 64, "drug": task["drug_id"], "arm": task["model"], "max_abs_difference": error}), flush=True)
    if blocked:
        raise ValueError("Export attempted a prohibited operation")
    if sum(r["zero_molecular_solution"] is True for r in output_rows) != 29:
        raise ValueError("Final zero-solution count changed")
    table_rows = []
    for p, n in [(candidate_path, candidate_count), (coeff_path, coefficient_count), (tissue_path, tissue_count)]:
        table_rows.append({"path": p.relative_to(out).as_posix(), "sha256": sha(p), "rows": n})
    common = {"schema_version": "reo-local-assets-1.0", "date": "2026-10-16", "status": "local_technical_assets_prepared_public_access_and_rights_pending", "first_lock_sha256": sha(FIRST), "second_lock_sha256": sha(SECOND), "source_manifest_sha256": sha(MANIFEST), "environment": environment, "threads": 1, "model_count": len(output_rows), "molecular_models": 56, "zero_molecular_models": 29, "nonzero_molecular_models": 27, "tissue_reference_models": 8, "synthetic_numeric_checks": 384, "max_absolute_native_portable_difference": max(r["max_absolute_native_portable_difference"] for r in output_rows), "absolute_tolerance": 1e-12, "relative_tolerance": 1e-12, "metric_or_model_selection_scoring": False, "real_expression_rows_read": False, "new_training": False, "CTRP_labels_indices_or_individual_predictions_read": False, "original_sentinels_used": False, "fit_methods_disabled": True, "wall_seconds": time.perf_counter() - start, "models": output_rows, "tables": table_rows, "read_allowlist_opened": sorted(opened), "blocked_accesses": blocked, "interpretation": "Equality of numerical implementations on generated inputs; no predictive performance or scientific replication claim."}
    _write_new(out / "model_export_manifest.json", _canonical(common))
    # Public summary publishes only counts, hashes and synthetic numerical
    # agreement. No real mapping, fitted parameter vector or tissue mean.
    summary.parent.mkdir(parents=True, exist_ok=True)
    _write_new(summary, (json.dumps(common, ensure_ascii=False, allow_nan=False, indent=2) + "\n").encode("utf-8"))
    print(json.dumps({"completed": 64, "candidate_rows": candidate_count, "coefficient_rows": coefficient_count, "zero_molecular_models": 29, "public_summary": summary.relative_to(ROOT).as_posix()}), flush=True)


if __name__ == "__main__":
    main()
