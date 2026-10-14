from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/synthetic_mechanisms"
INPUT = RESULTS / "performance_by_outer_task.csv"
AUDIT_CSV = RESULTS / "representation_gene_set_audit.csv"
REPORT = ROOT / "reports/阶段9b_候选基因集一致性审计.md"
EXPECTED_MODELS = {
    "tissue_null",
    "expression_strong",
    "mad_expression",
    "mad_rank",
    "mad_pair_full",
    "response_association_expression",
    "response_association_rank",
    "response_association_pair_full",
}
RULES = {
    "MAD": ("mad_expression", "mad_rank", "mad_pair_full"),
    "response_association": (
        "response_association_expression",
        "response_association_rank",
        "response_association_pair_full",
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def main() -> None:
    with INPUT.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 960:
        raise ValueError(f"Expected 960 outer-task rows, observed {len(rows)}")

    grouped: dict[tuple[str, int, int], dict[str, dict]] = defaultdict(dict)
    for row in rows:
        key = (row["scenario"], int(row["replicate"]), int(row["outer_fold"]))
        model = row["model"]
        if model in grouped[key]:
            raise ValueError(f"Duplicate model row for {key}: {model}")
        grouped[key][model] = row

    if len(grouped) != 120:
        raise ValueError(f"Expected 120 scenario × replicate × outer-fold contexts, observed {len(grouped)}")

    audit_rows: list[dict] = []
    failures: list[str] = []
    context_n_values = set()
    for key in sorted(grouped):
        scenario, replicate, outer_fold = key
        models = grouped[key]
        if set(models) != EXPECTED_MODELS:
            failures.append(f"{key}: model set mismatch {sorted(models)}")
            continue
        sample_ns = {int(row["n"]) for row in models.values()}
        context_n_values.update(sample_ns)
        if len(sample_ns) != 1:
            failures.append(f"{key}: outer test sample sizes differ: {sorted(sample_ns)}")

        if int(models["expression_strong"]["selected_gene_count"]) != 220:
            failures.append(f"{key}: expression_strong did not use the prespecified 220 genes")
        if int(models["expression_strong"]["pair_feature_count"]) != 220:
            failures.append(f"{key}: expression_strong feature count differs from its 220-gene budget")

        for rule_name, rule_models in RULES.items():
            selected_counts = {int(models[m]["selected_gene_count"]) for m in rule_models}
            selected_hashes = {models[m]["selected_gene_ids_sha256"] for m in rule_models}
            feature_counts = {m: int(models[m]["pair_feature_count"]) for m in rule_models}
            expected_feature_counts = {
                rule_models[0]: 200,
                rule_models[1]: 200,
                rule_models[2]: 19900,
            }
            if feature_counts != expected_feature_counts:
                failures.append(f"{key}/{rule_name}: representation feature budget mismatch: {feature_counts}")
            matched = len(selected_counts) == 1 and len(selected_hashes) == 1 and "" not in selected_hashes
            if not matched:
                failures.append(
                    f"{key}/{rule_name}: candidate sets differ: counts={sorted(selected_counts)}, hashes={sorted(selected_hashes)}"
                )
            if selected_counts != {200}:
                failures.append(f"{key}/{rule_name}: expected a fixed 200-gene candidate set, observed {sorted(selected_counts)}")
            audit_rows.append(
                {
                    "scenario": scenario,
                    "replicate": replicate,
                    "outer_fold": outer_fold,
                    "selection_rule": rule_name,
                    "models_compared": ";".join(rule_models),
                    "selected_gene_count": next(iter(selected_counts)) if len(selected_counts) == 1 else "",
                    "unique_selected_gene_hashes": len(selected_hashes),
                    "selected_gene_ids_sha256": next(iter(selected_hashes)) if len(selected_hashes) == 1 else "",
                    "feature_counts_expression_rank_pair": json.dumps(feature_counts, sort_keys=True),
                    "outer_test_n_shared": len(sample_ns) == 1,
                    "outer_test_n": next(iter(sample_ns)) if len(sample_ns) == 1 else "",
                    "candidate_gene_sets_identical": matched,
                }
            )

    if context_n_values != {80}:
        failures.append(f"Expected 80 samples in every outer test fold, observed {sorted(context_n_values)}")

    write_csv(
        AUDIT_CSV,
        audit_rows,
        [
            "scenario",
            "replicate",
            "outer_fold",
            "selection_rule",
            "models_compared",
            "selected_gene_count",
            "unique_selected_gene_hashes",
            "selected_gene_ids_sha256",
            "feature_counts_expression_rank_pair",
            "outer_test_n_shared",
            "outer_test_n",
            "candidate_gene_sets_identical",
        ],
    )

    matched_rows = sum(bool(row["candidate_gene_sets_identical"]) for row in audit_rows)
    lines = [
        "# 阶段9b：表示臂候选基因集一致性独立审计",
        "",
        "该审计在合成机制模型拟合全部结束后运行，不读取真实药物开发集或CTRP，也不新增模型拟合。",
        "",
        f"- 输入外层任务行：{len(rows)}（SHA256 `{sha256(INPUT)}`）。",
        f"- 场景×重复×外层折上下文：{len(grouped)}；每个测试折样本数：{sorted(context_n_values)}。",
        f"- MAD 与响应关联两种筛选规则下的表达/秩/pair 三表示候选集核查：{matched_rows}/{len(audit_rows)} 完全一致。",
        "- 每个表示臂使用同一筛选规则产生的同一组200个基因；表达和秩各编码200个特征，完整pair展开为19,900个基因对特征；单独强表达基线实际编码220个特征。",
        "- 单独的强表达基线按冻结方案使用220个基因，不属于200基因三表示配对比较。",
        "- 此审计证明的是候选基因集合及外层测试样本量一致；它不替代对底层划分、标签或模型选择过程的主协议校验，也不提供真实药物证据。",
        f"- 审计源文件SHA256：`{sha256(Path(__file__))}`。",
        f"- 审计表SHA256：`{sha256(AUDIT_CSV)}`。",
    ]
    if failures:
        lines.extend(["", "## 未通过项", ""])
        lines.extend(f"- {item}" for item in failures)
    else:
        lines.extend(["", "## 结论", "", "全部预设场景、重复、外层折和模型臂均满足检查条件。"])
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(
        json.dumps(
            {
                "contexts": len(grouped),
                "rule_context_checks": len(audit_rows),
                "matched_candidate_sets": matched_rows,
                "outer_test_n_values": sorted(context_n_values),
                "failures": failures,
                "audit_csv_sha256": sha256(AUDIT_CSV),
                "report_sha256": sha256(REPORT),
            },
            ensure_ascii=False,
        )
    )
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
