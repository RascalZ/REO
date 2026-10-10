"""Repository-relative outputs and configurable, read-only original input paths."""
import json
import os
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
ORIGINAL_ROOT = Path(os.environ.get('REO_ORIGINAL_ROOT', str(REPOSITORY.parent))).resolve()
CONFIG_PATH = Path(os.environ.get('REO_INPUT_CONFIG', str(REPOSITORY / 'config/input_paths.local.json')))
OVERRIDES = json.loads(CONFIG_PATH.read_text(encoding='utf-8')) if CONFIG_PATH.exists() else {}
DEFAULTS = {
    'gdsc2': '3.Project/Pan-Drugs-REOs/data/GDSC2_fitted_dose_response_27Oct23.xlsx',
    'legacy_cell_annotations': '3.Project/Pan-Drugs-REOs/data/1_Cell_lines_annotations_20181226.txt',
    'legacy_ctrp_bridge': '3.Project/Pan-Drugs-REOs/data/1.CTRP-CCLE-TCGACode-名字对照表.txt',
    'legacy_ctrp_response': '3.Project/Pan-Drugs-REOs/data/CTRP--Drug--Cellline--AUC.txt',
    'legacy_project_counts': '3.Project/Pan-Drugs-REOs/data/1_CCLE_RNAseq_genes_counts_20180929.gct',
    'ccle_tpm': '1.数据/2.CCLE-耐药敏感表达数据/0.CCLE原始数据/1.CCLE_RNAseq_rsem_genes_tpm_20180929.txt',
    'ccle_counts': '1.数据/2.CCLE-耐药敏感表达数据/0.CCLE原始数据/1.CCLE_RNAseq_genes_counts_20180929.gct',
    'ccle_annotations': '1.数据/2.CCLE-耐药敏感表达数据/0.CCLE原始数据/1.Cell_lines_annotations_20181226.txt',
}

def original_input(name):
    """Resolve a named original file; overrides may be absolute or ORIGINAL_ROOT-relative."""
    value = Path(OVERRIDES.get(name, DEFAULTS[name]))
    return value.resolve() if value.is_absolute() else (ORIGINAL_ROOT / value).resolve()

def original_label(path):
    """Use a relative label where possible, otherwise retain the local absolute provenance."""
    path = Path(path).resolve()
    return path.relative_to(ORIGINAL_ROOT).as_posix() if path.is_relative_to(ORIGINAL_ROOT) else str(path)

def prepare_output_directories():
    for name in ['config', 'data_raw', 'data_processed', 'data_manifest', 'results', 'reports', 'logs']:
        (REPOSITORY / name).mkdir(parents=True, exist_ok=True)
