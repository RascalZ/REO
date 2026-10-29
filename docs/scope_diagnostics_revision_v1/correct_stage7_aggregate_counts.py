"""One-time Stage 7 display correction using saved aggregate counts only.

Does not invoke the original analyzer or access observations, OOF predictions,
model checkpoints, response labels, or CTRP inputs.
"""

import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
BACKUP = OUT / "local_private" / "stage7_before"
sys.path.insert(0, str(ROOT / "src"))
from analyze_development_diagnostics import aggregate_zero_counts


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def rel(path):
    return path.relative_to(ROOT).as_posix()


def main():
    receipt_path = OUT / "stage7_correction_audit.json"
    if receipt_path.exists():
        raise RuntimeError("Correction receipt already exists; do not overwrite provenance.")
    source = ROOT / "results/ordinal_repeated/molecular_and_failure_summary.csv"
    s6 = ROOT / "tables/manuscript_v1/TableS6_zero_solution_audit.csv"
    profile = ROOT / "results/development_diagnostics/drug_heterogeneity.csv"
    report = ROOT / "reports/阶段7_零解与药物异质性诊断.md"
    manifest_path = ROOT / "results/development_diagnostics/diagnostics_manifest.json"
    analyzer = ROOT / "src/analyze_development_diagnostics.py"
    old_manifest = json.loads((BACKUP / rel(manifest_path)).read_text(encoding="utf-8"))
    for path in (profile, report, manifest_path):
        if digest(path) != digest(BACKUP / rel(path)):
            raise RuntimeError(f"Unexpected post-backup change: {rel(path)}")
    counts = aggregate_zero_counts(read_csv(source))
    s6_counts = {
        (int(row["GDSC2_compound_ID"]), row["Model_arm"]): int(row["zero_molecular_fits"])
        for row in read_csv(s6)
    }
    for key, count in counts.items():
        if s6_counts[key] != count:
            raise RuntimeError(f"S6 count mismatch: {key}")
    old_rows = read_csv(profile)
    rows = [dict(row) for row in old_rows]
    count_columns = [name for name in rows[0] if name.endswith("_zero_fits")]
    if len(rows) != 8 or len(count_columns) != 7:
        raise RuntimeError("Unexpected profile dimensions")
    changed = []
    for row, previous in zip(rows, old_rows):
        for column in count_columns:
            model = column.removesuffix("_zero_fits")
            row[column] = str(counts[(int(row["gdsc_drug_id"]), model)])
            if row[column] != previous[column]:
                changed.append({"drug": row["drug_name"], "column": column,
                                "before": int(previous[column]), "after": int(row[column])})
        if any(row[key] != previous[key] for key in row if key not in count_columns):
            raise RuntimeError("A non-count profile field changed")
    if sum(int(row[column]) for row in rows for column in count_columns) != 430:
        raise RuntimeError("Corrected molecular zero counts do not sum to 430")
    csv_buffer = io.StringIO(newline="")
    writer = csv.DictWriter(csv_buffer, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    profile.write_bytes(csv_buffer.getvalue().encode("utf-8-sig"))

    by_name = {row["drug_name"]: row for row in rows}
    report_lines = report.read_text(encoding="utf-8").splitlines(keepends=True)
    in_table = False
    corrected_rows = 0
    for index, line in enumerate(report_lines):
        if line.startswith("|药物|N|组织数|"):
            in_table = True
            continue
        if in_table and not line.startswith("|"):
            in_table = False
        fields = line.rstrip("\r\n").split("|")
        if in_table and len(fields) == 11 and fields[1] in by_name:
            fields[-2] = by_name[fields[1]]["response_association_pair_full_zero_fits"] + "/15"
            report_lines[index] = "|".join(fields) + "\n"
            corrected_rows += 1
    if corrected_rows != 8:
        raise RuntimeError("Expected eight report rows")
    note = (
        "\n勘误说明（2026-10-29）：第二张药物异质性表的零解列及对应CSV七个零解列，"
        "已由误取repeat 0修正为三个repeat合计（每臂15次外层拟合）。"
        "本次仅使用既有聚合计数并与Table S6核对；未重跑分析、训练或评分。"
        "下方运行时间、执行提交及原分析程序SHA仍为2026-10-14原执行记录。"
        "原文件保留，详见[勘误记录](../docs/scope_diagnostics_revision_v1/STAGE7_CORRECTION.md)。\n"
    )
    report_lines.insert(3, note)
    report.write_text("".join(report_lines), encoding="utf-8", newline="\n")

    created = datetime.now(timezone.utc).isoformat()
    manifest = dict(old_manifest)
    manifest["original_execution_output_sha256"] = dict(old_manifest["output_sha256"])
    manifest["output_sha256"] = dict(old_manifest["output_sha256"])
    manifest["output_sha256"][profile.name] = digest(profile)
    manifest["artifact_correction"] = {
        "created_utc": created,
        "receipt": rel(receipt_path),
        "description": "Aggregate-only correction: sum all three repeats in profile zero-fit columns.",
        "original_manifest_sha256": digest(BACKUP / rel(manifest_path)),
        "original_manifest_backup": rel(BACKUP / rel(manifest_path)),
        "original_execution_metadata_unchanged": True,
        "historical_fields_scope": "analysis_script_sha256, execution_commit, created_utc, runtime and input hashes describe the original 2026-10-14 run, not this correction",
        "output_hashes_scope": "output_sha256 binds current files; original_execution_output_sha256 preserves the original output bindings",
        "corrected_analyzer_sha256": digest(analyzer),
        "new_model_fits": 0, "new_scoring": False,
        "development_label_accessed": False, "external_response_accessed": False,
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    for name, sha in manifest["output_sha256"].items():
        if digest(profile.parent / name) != sha:
            raise RuntimeError(f"Current manifest output checksum mismatch: {name}")

    receipt = {
        "status": "completed", "created_utc": created,
        "scope": "Stage7 zero-count display correction only; historical model results and numerical performance unchanged",
        "procedure": rel(Path(__file__)), "procedure_sha256": digest(Path(__file__)),
        "inputs": {rel(path): digest(path) for path in (source, s6)},
        "before_sha256": {rel(path): digest(BACKUP / rel(path)) for path in (analyzer, profile, report, manifest_path)},
        "after_sha256": {rel(path): digest(path) for path in (analyzer, profile, report, manifest_path)},
        "backup_root": rel(BACKUP),
        "changed_count_cells": changed, "count_cells_checked": 56,
        "corrected_report_rows_checked": corrected_rows,
        "non_count_profile_fields_unchanged": True,
        "all_counts_match_table_s6": True,
        "molecular_outer_fits": 840, "zero_molecular_outer_fits": 430,
        "new_model_fits": 0, "new_scoring": False,
        "development_label_accessed": False, "external_response_accessed": False,
        "regression_test": {
            "command": ".venv/Scripts/python.exe -m unittest discover -s tests -p test_development_zero_count_aggregation.py -v",
            "passed": 4,
            "initial_environment_attempt": "Bundled Python lacked psutil; no tests executed there; repository .venv completed all four tests",
        },
    }
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "completed", "changed_cells": len(changed), "checked_cells": 56,
                      "total_zero": 430, "receipt": rel(receipt_path)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
