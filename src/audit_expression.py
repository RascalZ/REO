"""Audit and prepare historical CCLE gene TPM without changing original data.

Inputs are CCLE 2018-09-29 RSEM gene TPM and matching gene annotation.
Output expression_log2tpm.npy is samples x genes in the exact TSV index order.
No outcome-dependent operations, sample exclusion, or variable-gene selection.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path
import re
import urllib.request

import numpy as np
import pandas as pd
from project_paths import REPOSITORY, ORIGINAL_ROOT, original_input


ROOT = ORIGINAL_ROOT
BUILD = REPOSITORY
LOCAL = ROOT / "1.数据/2.CCLE-耐药敏感表达数据/0.CCLE原始数据"
TPM = original_input('ccle_tpm')
COUNTS = original_input('ccle_counts')
ANNOTATION = original_input('ccle_annotations')
BROAD = "https://data.broadinstitute.org/ccle/"
FIGSHARE = "https://api.figshare.com/v2/articles/27993248/versions/1"
GTF_NAME = "gencode.v19.genes.v7_model.patched_contigs.gtf.gz"


def digests(path):
    h = hashlib.sha256()
    md = hashlib.md5()
    lines = 0
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            h.update(block)
            md.update(block)
            lines += block.count(b"\n")
    return {"path": str(path), "bytes": path.stat().st_size,
            "sha256": h.hexdigest(), "md5": md.hexdigest(), "newline_count": lines}


def download(url, path):
    """Download only named small reference assets, never a full release."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        req = urllib.request.Request(url, headers={"User-Agent": "CCLE-expression-audit/1.0"})
        with urllib.request.urlopen(req, timeout=60) as response, path.open("wb") as out:
            while block := response.read(1024 * 1024):
                out.write(block)
    return {**digests(path), "source_url": url}


def references(allow_download):
    meta_dir = BUILD / "data_raw/cell_metadata"
    expression_dir = BUILD / "data_raw/expression"
    meta_dir.mkdir(parents=True, exist_ok=True)
    expression_dir.mkdir(parents=True, exist_ok=True)
    meta_path = meta_dir / "DepMap_24Q4_Public_v1_figshare_metadata.json"
    if allow_download:
        download(FIGSHARE, meta_path)
    metadata = json.loads(meta_path.read_text(encoding="utf-8"))
    provenance = []
    for asset in metadata["files"]:
        if asset["name"] not in {"Model.csv", "README.txt"}:
            continue
        path = meta_dir / ("DepMap_24Q4_v1_" + asset["name"])
        item = download(asset["download_url"], path) if allow_download else digests(path)
        item.update(source_url=asset["download_url"],
                    expected_md5=asset["computed_md5"],
                    expected_bytes=asset["size"],
                    source_doi=metadata["doi"], release="DepMap 24Q4 Public v1",
                    license=metadata["license"])
        if item["md5"] != asset["computed_md5"] or item["bytes"] != asset["size"]:
            raise ValueError(f"Official checksum/size mismatch: {path}")
        item["official_md5_verified"] = True
        provenance.append(item)
    gtf_path = expression_dir / GTF_NAME
    if allow_download:
        provenance.append(download(BROAD + GTF_NAME, gtf_path))
    elif gtf_path.exists():
        provenance.append({**digests(gtf_path), "source_url": BROAD + GTF_NAME})
    return meta_dir / "DepMap_24Q4_v1_Model.csv", gtf_path, provenance


def gtf_annotations(path):
    rows = {}
    if not path.exists():
        return rows
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)";', fields[8]))
            gene = attrs["gene_id"]
            if gene in rows:
                raise ValueError(f"Duplicated full gene ID in GTF: {gene}")
            rows[gene] = {"gene_symbol": attrs.get("gene_name", ""),
                          "gene_type": attrs.get("gene_type", ""),
                          "chromosome": fields[0]}
    return rows


def header(path):
    with path.open(encoding="utf-8-sig") as stream:
        return stream.readline().rstrip("\r\n").split("\t")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--download-references", action="store_true",
                        help="Fetch only fixed Model/README and matched 7 MB GTF if absent")
    args = parser.parse_args()
    processed = BUILD / "data_processed"
    reports = BUILD / "reports"
    processed.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    model_path, gtf_path, provenance = references(args.download_references)
    print("Auditing historical TPM and counts headers", flush=True)
    tpm_hash = digests(TPM)
    counts_hash = digests(COUNTS)
    annotation_hash = digests(ANNOTATION)
    tpm_header = header(TPM)
    count_header = header(COUNTS)
    if tpm_header[:2] != ["gene_id", "transcript_ids"]:
        raise ValueError("Unexpected TPM header")
    sample_names = tpm_header[2:]
    if len(set(sample_names)) != len(sample_names):
        raise ValueError("Duplicated expression sample header")
    if sample_names != count_header[2:]:
        raise ValueError("TPM/count sample names or ordering differ")

    # Keep original five-decimal TPM precision for pairwise ordering and ties.
    dtype = {name: np.float64 for name in sample_names}
    dtype.update(gene_id=str, transcript_ids=str)
    data = pd.read_csv(TPM, sep="\t", dtype=dtype, float_precision="round_trip")
    gene_ids = data.pop("gene_id").astype(str)
    transcript_ids = data.pop("transcript_ids").astype(str)
    matrix = data.to_numpy(dtype=np.float64, copy=False).T.copy()
    del data
    if not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError("TPM contains negative or nonfinite entries")
    if gene_ids.duplicated().any():
        raise ValueError("Duplicated full Ensembl gene IDs; explicit resolution required")
    ensembl_ids = gene_ids.str.replace(r"\.\d+$", "", regex=True)
    tpm_sums = matrix.sum(axis=1, dtype=np.float64)
    zeros = (matrix == 0).sum(axis=1)
    min_value, max_value = float(matrix.min()), float(matrix.max())
    # The original float64 TPM is retained so REO need not depend on a transform.
    np.save(processed / "expression_tpm.npy", matrix)
    raw_tpm = np.load(processed / "expression_tpm.npy", mmap_mode="r")
    matrix += 1.0
    np.log2(matrix, out=matrix)
    # Checking adjacent original values after sorting proves all pairwise weak
    # and strict order relations for each sample, without quadratic pair arrays.
    order_checks = {"method": "all sorted adjacent TPM values per sample; suffices for all pairwise order relations",
                    "samples_checked": len(sample_names),
                    "adjacent_comparisons": 0, "original_equal_adjacent_pairs": 0,
                    "original_strict_adjacent_pairs": 0,
                    "tie_relation_violations": 0, "strict_order_violations": 0}
    for sample_index in range(len(sample_names)):
        order = np.argsort(raw_tpm[sample_index], kind="stable")
        raw_delta = np.diff(raw_tpm[sample_index, order])
        log_delta = np.diff(matrix[sample_index, order])
        original_equal = raw_delta == 0
        order_checks["adjacent_comparisons"] += len(raw_delta)
        order_checks["original_equal_adjacent_pairs"] += int(original_equal.sum())
        order_checks["original_strict_adjacent_pairs"] += int((~original_equal).sum())
        order_checks["tie_relation_violations"] += int((original_equal & (log_delta != 0)).sum())
        order_checks["strict_order_violations"] += int((~original_equal & (log_delta <= 0)).sum())
    if order_checks["tie_relation_violations"] or order_checks["strict_order_violations"]:
        raise ValueError(f"TPM-to-log transform failed rank preservation: {order_checks}")
    del raw_tpm
    np.save(processed / "expression_log2tpm.npy", matrix)

    hist = pd.read_csv(ANNOTATION, sep="\t", dtype=str, keep_default_na=False)
    modern = pd.read_csv(model_path, dtype=str, keep_default_na=False)
    if hist["CCLE_ID"].duplicated().any() or modern["ModelID"].duplicated().any():
        raise ValueError("Nonunique cell annotation primary keys")
    hist = hist.set_index("CCLE_ID")
    modern_by_id = modern.set_index("ModelID")
    # No fuzzy aliases: direct CCLE exact names first, historical ACH fallback.
    ccle_groups = modern.loc[modern["CCLEName"] != ""].groupby("CCLEName")
    direct = {name: group.iloc[0] for name, group in ccle_groups if len(group) == 1}
    sample_rows = []
    conflicts = []
    for i, ccle in enumerate(sample_names):
        old_ach = hist.loc[ccle, "depMapID"] if ccle in hist.index else ""
        row = direct.get(ccle)
        strategy = "exact_modern_CCLEName"
        if row is None and old_ach in modern_by_id.index:
            row = modern_by_id.loc[old_ach]
            strategy = "historical_CCLE_ID_to_ACH"
        if row is None:
            row = pd.Series(dtype=str)
            strategy = "historical_ACH_only" if old_ach else "unmapped"
        ach = row.get("ModelID", old_ach)
        if ach and old_ach and ach != old_ach:
            conflicts.append({"ccle_id": ccle, "historical_depmap_id": old_ach,
                              "modern_depmap_id": ach})
            strategy = "CONFLICT_REQUIRES_REVIEW"
        sample_rows.append({"expression_row_index": i, "ccle_id": ccle,
                            "depmap_id": ach, "historical_depmap_id": old_ach,
                            "mapping_strategy": strategy,
                            "patient_id": row.get("PatientID", ""),
                            "cell_line_name": row.get("CellLineName", ""),
                            "stripped_cell_line_name": row.get("StrippedCellLineName", ""),
                            "oncotree_lineage": row.get("OncotreeLineage", ""),
                            "oncotree_primary_disease": row.get("OncotreePrimaryDisease", ""),
                            "tcga_code_historical": hist.loc[ccle, "tcga_code"] if ccle in hist.index else "",
                            "sanger_model_id": row.get("SangerModelID", ""),
                            "cosmic_id": row.get("COSMICID", ""),
                            "rrid": row.get("RRID", ""),
                            "tpm_sum": float(tpm_sums[i]),
                            "zero_gene_fraction": float(zeros[i] / len(gene_ids))})
    sample_df = pd.DataFrame(sample_rows)
    sample_df.to_csv(processed / "expression_samples.tsv", sep="\t", index=False)
    gene_annotation = gtf_annotations(gtf_path)
    genes = pd.DataFrame({"expression_column_index": np.arange(len(gene_ids)),
                          "gene_id_full": gene_ids, "ensembl_id": ensembl_ids,
                          "transcript_ids": transcript_ids})
    for key in ["gene_symbol", "gene_type", "chromosome"]:
        genes[key] = genes["gene_id_full"].map(lambda g: gene_annotation.get(g, {}).get(key, ""))
    genes["annotation_exact_version_match"] = genes["gene_id_full"].isin(gene_annotation)
    genes.to_csv(processed / "expression_genes.tsv", sep="\t", index=False)
    summary = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "expression_source_version": "CCLE 2019 release, RNA-seq RSEM gene TPM 20180929",
        "original_files_modified": False,
        "large_expression_downloaded": False,
        "local_tpm_byte_identity_to_official_remote_verified": False,
        "matrix_orientation": "samples x genes",
        "matrix_dtype": str(matrix.dtype),
        "matrix_shape": list(matrix.shape),
        "transformation": "log2(TPM + 1); float64 input, calculation and storage",
        "raw_tpm_retained_as_float64_npy": True,
        "log_transform_pairwise_order_verification": order_checks,
        "feature_selection_applied": False,
        "samples_excluded": 0,
        "gene_id_version_collapsing_applied": False,
        "tpm_count_sample_headers_identical": sample_names == count_header[2:],
        "count_format": "Tab-delimited table; #1.2 and dimension lines absent; not intact standard GCT",
        "tpm_min": min_value, "tpm_max": max_value,
        "tpm_sum_min": float(tpm_sums.min()), "tpm_sum_max": float(tpm_sums.max()),
        "tpm_sum_median": float(np.median(tpm_sums)),
        "nonfinite_values": 0, "negative_values": 0,
        "unique_full_gene_ids": int(gene_ids.nunique()),
        "duplicate_unversioned_gene_ids": int(ensembl_ids.duplicated().sum()),
        "genes_exact_gtf_version_matches": int(genes.annotation_exact_version_match.sum()),
        "protein_coding_genes": int((genes.gene_type == "protein_coding").sum()),
        "samples_with_depmap_id": int((sample_df.depmap_id != "").sum()),
        "samples_with_patient_id": int((sample_df.patient_id != "").sum()),
        "duplicate_depmap_id_rows": int(sample_df.loc[sample_df.depmap_id != "", "depmap_id"].duplicated().sum()),
        "mapping_conflicts": conflicts,
        "modern_metadata_rows": len(modern),
        "modern_metadata_has_parent_model_column": False,
        "patient_id_definition": "Broad 24Q4 README: unique identifier for models derived from the same tissue sample; not a complete patient genealogy",
        "source_inputs": [dict(tpm_hash, source_url=BROAD + "CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz"),
                          dict(counts_hash, source_url=BROAD + "CCLE_RNAseq_genes_counts_20180929.gct.gz"),
                          dict(annotation_hash, source_url=BROAD + "Cell_lines_annotations_20181226.txt")],
        "downloaded_references": provenance,
    }
    summary["derived_files"] = [digests(processed / name) for name in
                                ["expression_samples.tsv", "expression_genes.tsv", "expression_tpm.npy", "expression_log2tpm.npy"]]
    (processed / "expression_audit.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    (BUILD / "data_raw/expression/expression_provenance.json").write_text(
        json.dumps({"retrieved_utc": summary["created_utc"], "local_inputs": summary["source_inputs"],
                    "downloaded_references": provenance}, ensure_ascii=False, indent=2), encoding="utf-8")
    with (reports / "expression_audit.md").open("w", encoding="utf-8") as out:
        out.write("# CCLE表达数据审计\n\n")
        out.write("本地完整历史RSEM gene TPM可用，无需下载新表达矩阵；原始文件只读。\n\n")
        out.write(f"- TPM：{matrix.shape[0]}个样本 × {matrix.shape[1]}个基因；样本列与counts完全一致。\n")
        out.write(f"- 全基因TPM列和：最小{tpm_sums.min():.3f}，中位{np.median(tpm_sums):.3f}，最大{tpm_sums.max():.3f}；未发现负值或非有限值。\n")
        out.write(f"- 所有{gene_ids.nunique()}个完整ENSG.version唯一；去版本重复数{ensembl_ids.duplicated().sum()}。未合并基因、未筛选变异基因。\n")
        out.write(f"- matched GENCODE19 GTF精确版本匹配{genes.annotation_exact_version_match.sum()}；protein_coding {sum(genes.gene_type == 'protein_coding')}。所有基因保留。\n")
        out.write(f"- ACH映射{summary['samples_with_depmap_id']}，有Broad PatientID {summary['samples_with_patient_id']}，ACH重复行{summary['duplicate_depmap_id_rows']}，新旧映射冲突{len(conflicts)}。全部未映射样本也保留。\n")
        out.write("- 本地counts已去除标准GCT的版本和维度前导行，是加工过的tab表，不能仅按扩展名当作完整原始GCT。它只用于header审计，不作为本次REO输入。\n")
        out.write("- 原TPM和log2(TPM+1)均以float64保存，NPY行为样本、列为基因，对应TSV中0基索引；REO可直接读取原TPM以免依赖转换。输入用round-trip浮点解析，未先转成float32。\n")
        out.write(f"- 对数变换秩序检查覆盖全部{len(sample_names)}个样本、{order_checks['adjacent_comparisons']}个排序后相邻比较；原相等{order_checks['original_equal_adjacent_pairs']}对、严格不等{order_checks['original_strict_adjacent_pairs']}对，ties改变{order_checks['tie_relation_violations']}、严格次序改变{order_checks['strict_order_violations']}。排序后相邻检查足以验证所有基因对的相等和大小关系，无需二次方枚举。此项仅保证本次双精度数据及转换，后续不要再转float32构造REO。\n")
        out.write("- 仅核验本地TPM的结构、数值完整性、版本和官方来源路径，未重新下载官方TPM来逐字节比对，因此不声称本地与远端字节相同。\n")
        out.write("- Model.csv固定DepMap24Q4 v1，来自Broad Figshare并验证官方MD5及大小；含PatientID但无ParentModelID。PatientID在官方README定义为源于相同组织样本的模型标识，并非完整患者谱系；需联合Sanger parent/patient信息后冻结分组。\n\n")
        out.write("## 输出\n\n")
        for name in ["expression_tpm.npy", "expression_log2tpm.npy", "expression_samples.tsv", "expression_genes.tsv", "expression_audit.json"]:
            out.write(f"- `rebuild_v2/data_processed/{name}`\n")
        out.write("\n下载来源、MD5/SHA256和输入哈希：`rebuild_v2/data_raw/expression/expression_provenance.json`。\n\n")
        out.write("## 官方来源\n\n")
        out.write("- [Broad CCLE历史文件目录](https://data.broadinstitute.org/ccle/)：同一20180929 gene TPM、counts及配套GENCODE19 GTF。\n")
        out.write("- [Ghandi et al., Nature 2019](https://doi.org/10.1038/s41586-019-1186-3)：CCLE资源研究。\n")
        out.write("- [DepMap24Q4 v1固定记录](https://doi.org/10.25452/figshare.plus.27993248.v1)：Model.csv与README，CC BY4.0。\n")
        out.write("- 2018表达与2024身份注释属于不同用途的固定版本；未混入2024表达重新量化结果。\n")
    print(json.dumps({k: v for k, v in summary.items() if k not in ["source_inputs", "downloaded_references", "derived_files"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
