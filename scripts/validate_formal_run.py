"""Stream independent integrity audit of frozen resumable formal CV artifacts.

Reads prepared/run artifacts and identity manifests, never original GDSC/CTRP
responses, CTRP indices or expression arrays. Use --allow-partial for progress
checks; a partial check never certifies formal completion or inferential results.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not bool(condition):
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def names_sha(values):
    return hashlib.sha256("\0".join(values).encode()).hexdigest()


def key(task):
    return (int(task["drug_id"]), int(task["repeat"]), int(task["outer_fold"]), task["model"])


def filename(task):
    return f"d{task['drug_id']}_r{task['repeat']}_f{task['outer_fold']}_{task['model']}.json"


def grid(kind, cfg):
    if kind == "mean":
        return [{}]
    if kind == "tissue":
        return [{"alpha": value} for value in cfg["ridge_alphas"]]
    if kind.endswith("_rf"):
        return [{"max_features": features, "min_samples_leaf": leaf}
                for features, leaf in itertools.product(cfg["random_forest"]["max_features"], cfg["random_forest"]["min_samples_leaf"])]
    return [{"alpha_ratio": value, "l1_ratio": l1} for value, l1 in
            itertools.product(cfg["elastic_net"]["alpha_ratios"], cfg["elastic_net"]["l1_ratios"])]


def feature_audit(payload, event, cfg, model):
    if event["n_selected_genes"] is None:
        require(payload is None, "Unexpected feature payload on unfitted-feature fallback")
        return
    require(payload is not None, "Missing outer feature schema")
    genes, indices, names = payload["gene_ids"], payload["gene_indices"], payload["feature_names"]
    require(len(genes) == len(set(genes)) == len(indices) == event["n_selected_genes"] and len(indices) == len(set(indices)), "Outer gene schema duplicated/misaligned")
    require(payload["training_samples"] == len(event["fit_sample_ids"]), "Outer schema training count mismatch")
    require(names_sha(genes) == event["selected_gene_ids_sha256"] and names_sha(names) == event["built_feature_names_sha256"], "Outer feature schema disagrees with actual-fit hashes")
    pairs = payload["pair_indices_original_input_columns"]
    identity = dict(zip(indices, genes))
    if payload["kind"] == "reo":
        require(genes == sorted(genes) and len(pairs) == len(names), "REO orientation/schema alignment invalid")
        for pair, name in zip(pairs, names):
            require(len(pair) == 2 and pair[0] in identity and pair[1] in identity, "REO pair outside selected genes")
            a, b = identity[pair[0]], identity[pair[1]]
            require(a < b and name == f"REO:{a}>{b}", "REO name cannot be traced to canonical pair indices")
    else:
        require(names == ["log2TPM:" + gene for gene in genes], "Expression feature names disagree with stable genes")
    if model.endswith("_rf"):
        forest = payload["random_forest"]
        require(forest["n_estimators"] == cfg["random_forest"]["n_estimators"] and forest["n_jobs"] == 1, "RF fixed tree/thread budget changed")
        require(forest["max_features"] == event["params"]["max_features"] and forest["min_samples_leaf"] == event["params"]["min_samples_leaf"], "RF chosen parameter mismatch")
        require(forest["molecular_features_used"] == event["n_features"], "RF molecular feature-count mismatch")
        if model == "reo_rf":
            selection = forest["pair_selection"]
            chosen = selection["selected_feature_indices"]
            selected_names = selection["selected_feature_names"]
            require(len(chosen) == len(set(chosen)) <= cfg["random_forest"]["reo_max_pairs"] and all(0 <= i < len(names) for i in chosen), "RF supervised pair budget/index mismatch")
            require(selected_names == [names[i] for i in chosen] == forest["molecular_feature_names"], "RF selected stable names mismatch")
            require(selection["selected_pairs_original_input_gene_indices"] == [pairs[i] for i in chosen], "RF pair indices lost canonical identity")
            require(names_sha(selected_names) == event["rf_selected_pair_names_sha256"], "Actual RF pair hash differs from selected schema")
    else:
        require(event["n_features"] == len(names), "EN candidate feature-count mismatch")


def validate(run, allow_partial=False, statistics_dir=None):
    run = Path(run).resolve()
    meta = load(run / "formal_metadata.json")
    cfg = meta["configuration"]
    lock_path = ROOT / "config/formal_analysis_v1.lock.json"
    lock = load(lock_path)
    require(sha(lock_path) == meta["freeze_sha256"], "Formal pre-run lock changed")
    require(meta["config_sha256"] == lock["files"][lock["model_config"]] and meta["statistics_sha256"] == lock["files"][lock["statistics_config"]] and meta["protocol_sha256"] == lock["files"][lock["protocol"]], "Formal model/statistical protocol linkage mismatch")
    for name, value in lock["files"].items():
        require(sha(ROOT / name) == value, "Frozen scientific dependency changed: " + name)
    cohort_lock = load(ROOT / "config/cohort_freeze_v2.lock.json")
    require(meta["cohort_sha256"] == lock["cohort_freeze_sha256"] == cohort_lock["freeze_sha256"], "Cancer cohort linkage mismatch")
    require(meta["external_response_accessed"] is False, "Formal metadata reports external response access")
    if meta.get("completed"):
        require(meta["data_guard"]["blocked_attempts"] == 0, "Formal controller attempted forbidden data access")
    revision = meta["execution_commit"]
    for name, value in meta["source_hashes"].items():
        blob = subprocess.check_output(["git", "show", f"{revision}:{name.replace(chr(92), '/') }"], cwd=ROOT)
        require(hashlib.sha256(blob).hexdigest() == value, "Execution source differs from immutable Git snapshot: " + name)
    for name, value in meta["prepared_input_hashes"].items():
        require(sha(run / name) == value, "Prepared formal input changed: " + name)

    manifest = pd.read_csv(ROOT / "results/cohort_v2/development_sample_manifest.csv", dtype=str, keep_default_na=False)
    response = pd.read_csv(run / "development_responses.csv", dtype={"gdsc_drug_id": str}, keep_default_na=False)
    global_ids = manifest.drop_duplicates("canonical_model_id").set_index("canonical_model_id")
    require(len(global_ids) == global_ids.biological_group_id.nunique(), "Biological groups no longer have the frozen unique-model structure")
    folds = pd.read_csv(run / "global_fold_manifest.csv", keep_default_na=False)
    memberships = {}
    require(not folds.duplicated(["repeat", "outer_fold", "inner_fold", "role", "canonical_model_id"]).any(), "Duplicate global fold roles")
    for (repeat, outer, inner, role), block in folds.groupby(["repeat", "outer_fold", "inner_fold", "role"]):
        memberships[(int(repeat), int(outer), int(inner), role)] = set(block.canonical_model_id)
    expected_keys = {(r, f, -1, role) for r in range(cfg["outer_repeats"]) for f in range(cfg["outer_folds"]) for role in ("train", "test")}
    expected_keys |= {(r, f, i, role) for r in range(cfg["outer_repeats"]) for f in range(cfg["outer_folds"]) for i in range(cfg["inner_folds"]) for role in ("train", "validation")}
    require(set(memberships) == expected_keys, "Global fold-role schema incomplete")
    group = global_ids.biological_group_id.to_dict()
    for repeat in range(cfg["outer_repeats"]):
        coverage = Counter()
        for outer in range(cfg["outer_folds"]):
            train, test = (memberships[(repeat, outer, -1, role)] for role in ("train", "test"))
            require(not train & test and train | test == set(global_ids.index), "Incomplete/overlapping outer fold")
            require(not {group[i] for i in train} & {group[i] for i in test}, "Biological group crosses outer fold")
            coverage.update(test)
            inner_coverage = Counter()
            for inner in range(cfg["inner_folds"]):
                it, iv = (memberships[(repeat, outer, inner, role)] for role in ("train", "validation"))
                require(not it & iv and it | iv == train and not {group[i] for i in it} & {group[i] for i in iv}, "Inner split incomplete or crosses groups")
                inner_coverage.update(iv)
            require(set(inner_coverage) == train and set(inner_coverage.values()) == {1}, "Inner validation coverage is not once")
        require(set(coverage) == set(global_ids.index) and set(coverage.values()) == {1}, "Outer test coverage is not once per repeat")

    tasks = [{"drug_id": int(d), "repeat": r, "outer_fold": f, "model": m} for d in cfg["drug_ids"]
             for r in range(cfg["outer_repeats"]) for f in range(cfg["outer_folds"]) for m in cfg["models"]]
    tasks_by_name = {filename(task): task for task in tasks}
    require(len(tasks_by_name) == meta["expected_tasks"] == 840, "Formal task product/budget changed")
    completed = fit_count = prediction_count = warned_fits = fallback_fits = wholly_constant = 0
    completed_rows = []
    registered = meta["checkpoint_hashes"]
    for path in sorted((run / "checkpoints").glob("*.json")):
        require(path.name in tasks_by_name, "Unexpected or duplicate-alias checkpoint filename")
        if path.name in registered:
            require(sha(path) == registered[path.name], "Registered checkpoint was altered: " + path.name)
        elif not allow_partial:
            raise ValueError("Completed run contains unregistered checkpoint")
        obj = load(path)
        digest = obj.pop("payload_sha256", None)
        require(digest == hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(",", ":")).encode()).hexdigest(), "Checkpoint payload integrity mismatch")
        task = tasks_by_name[path.name]
        require(obj["task"] == task and obj["freeze_sha256"] == meta["freeze_sha256"], "Checkpoint scientific/task context changed")
        if obj["status"] != "completed":
            require(allow_partial, "Failed checkpoint cannot certify full analysis")
            continue
        require(obj["error"] is None and obj["data_guard"]["external_response_accessed"] is False and obj["data_guard"]["blocked_attempts"] == 0, "Completed checkpoint contains failure/access violation")
        drug, repeat, outer, model = key(task)
        data = response[response.gdsc_drug_id.eq(str(drug))].set_index("canonical_model_id")
        ids = set(data.index)
        require(ids == set(manifest[manifest.gdsc_drug_id.eq(str(drug))].canonical_model_id), "Prepared drug observation set differs from frozen manifest")
        def expected(inner, role):
            return memberships[(repeat, outer, inner, role)] & ids
        def expr(identity_set):
            return {int(data.loc[identity, "expression_row_index"]) for identity in identity_set}
        train, test = expected(-1, "train"), expected(-1, "test")
        require(obj["expected_test_ids"] == sorted(test) and set(obj["expected_training_expression_ids"]) == expr(train) and set(obj["outer_test_expression_ids"]) == expr(test), "Checkpoint declarations differ from actual frozen outer split")
        require([set(values) for values in obj["inner_training_expression_ids"]] == [expr(expected(i, "train")) for i in range(cfg["inner_folds"])], "Checkpoint declarations differ from actual frozen inner split")
        rows = obj["oof_rows"]
        require(len(rows) == len(test) and {row["canonical_model_id"] for row in rows} == test, "OOF duplicate/missing identities in checkpoint")
        for row in rows:
            identity = row["canonical_model_id"]
            require(key(row) == key(task) and row["biological_group_id"] == data.loc[identity, "biological_group_id"] and row["tissue"] == data.loc[identity, "tissue"], "OOF identity/task metadata mismatch")
            require(np.isfinite([row["observed_auc"], row["predicted_auc"]]).all() and row["observed_auc"] == float(data.loc[identity, "response_auc"]), "OOF observed outcome differs from bound prepared development outcome")
        completed_rows.extend(rows)
        y = np.array([row["observed_auc"] for row in rows]); p = np.array([row["predicted_auc"] for row in rows])
        mse = np.mean((y-p)**2)
        for name, value in (("rmse", np.sqrt(mse)), ("mae", np.mean(np.abs(y-p)))):
            require(np.isclose(obj["result"]["metrics"][name], value, atol=1e-12, rtol=0), "Checkpoint metric differs from individual OOF records")
        choices = obj["result"]
        if model != "mean":
            scores = choices["inner_scores"]
            require(len(scores) == len(grid(model, cfg)), "Inner candidate budget incomplete")
            best = min(scores, key=lambda item: item["pooled_inner_rmse"])
            require(choices["chosen_params"] == best["params"] and choices["best_inner_rmse"] == best["pooled_inner_rmse"], "HPO did not select minimum logged pooled inner RMSE")
        seen = Counter(); hashes = defaultdict(lambda: defaultdict(set)); thresholds = defaultdict(list)
        lookup = data.reset_index().set_index("expression_row_index")
        outer_event = None
        for event in obj["fit_events"]:
            stage, candidate, inner = event["stage"], event["candidate"], event["inner_fold"]
            seen.update([(stage, candidate, inner)])
            actual = event["fit_sample_ids"]
            require(len(actual) == len(set(actual)), "Duplicated actual-fit identity")
            if stage == "outer_training_refit":
                require(candidate == inner == -1 and event["params"] == choices["chosen_params"], "Outer refit metadata/parameters invalid")
                desired = expr(train); outer_event = event
            else:
                require(stage == "inner_training" and model != "mean" and 0 <= inner < cfg["inner_folds"] and 0 <= candidate < len(grid(model,cfg)), "Invalid actual-inner fit context")
                require(event["params"] == grid(model,cfg)[candidate], "Actual candidate differs from frozen grid")
                desired = expr(expected(inner, "train"))
            require(set(actual) == desired and not set(actual) & expr(test), "Actual fit contains validation/test or omits intended training identities")
            local = lookup.loc[actual]
            target = local.response_auc.to_numpy(float)
            require(np.isclose(event["target_mean"], target.mean(), atol=1e-12, rtol=0) and np.isclose(event["target_sd"], target.std() if target.std() else 1., atol=1e-12, rtol=0), "Actual target normalization uses wrong training observations")
            require(event["prediction_calls"] == 1 and event["predict_seconds"] >= 0, "Missing actual prediction timing event")
            constant = all(np.ptp(values) == 0 for _, values in local.groupby("tissue").response_auc)
            wholly_constant += int(constant)
            if model in {"expression","expression_matched","reo"} and event["alpha_max"] is not None:
                require(np.isclose(event["actual_alpha"], event["alpha_max"] * event["params"]["alpha_ratio"], atol=1e-14, rtol=1e-12), "Actual EN alpha not current-fit-relative")
                thresholds[(stage,inner)].append(event["alpha_max"] * event["params"]["l1_ratio"])
                if event["params"]["alpha_ratio"] == 1 or constant:
                    require(event["n_nonzero_molecular_coeffs"] == 0, "Theoretical/exact-tissue null fit has molecular coefficients")
            if model == "reo_rf" and event["n_features"] is not None:
                require(0 <= event["n_features"] <= cfg["random_forest"]["reo_max_pairs"] and (not constant or event["n_features"] == 0), "RF pair budget/zero-association rule violated")
            for field in ("selected_gene_ids_sha256","built_feature_names_sha256","rf_selected_pair_names_sha256"):
                if event[field] is not None:
                    hashes[(stage,inner)][field].add(event[field])
            warned_fits += bool(event["warnings"]); fallback_fits += bool(event["fallback_reason"])
        expected_events = {("outer_training_refit",-1,-1)}
        if model != "mean":
            expected_events |= {("inner_training", c, i) for c in range(len(grid(model,cfg))) for i in range(cfg["inner_folds"])}
        require(set(seen) == expected_events and set(seen.values()) == {1}, "Candidate/partition missing or duplicated actual fit")
        require(all(len(values) == 1 for fields in hashes.values() for values in fields.values()), "Training features/pairs change across HPO candidates in same partition")
        require(all(np.allclose(values,values[0],atol=1e-12,rtol=1e-10) for values in thresholds.values()), "alpha_max reused/miscomputed across HPO candidates")
        feature_audit(choices["feature_audit"], outer_event, cfg, model)
        completed += 1; fit_count += len(obj["fit_events"]); prediction_count += len(rows)

    complete = meta.get("completed") is True
    if not allow_partial:
        require(complete and completed == meta["completed_tasks"] == meta["expected_tasks"] and fit_count == meta["expected_fit_events"] == 18120 and prediction_count == meta["expected_oof_predictions"] == 104307, "Formal analysis is incomplete or has incorrect totals")
        for name, value in meta["output_hashes"].items():
            require(sha(run / name) == value, "Formal aggregate output changed: " + name)
        saved = pd.read_csv(run / "oof_predictions.csv")
        expected_frame = pd.DataFrame(completed_rows)
        columns = list(expected_frame.columns)
        sort = ["drug_id","model","repeat","canonical_model_id"]
        pd.testing.assert_frame_equal(saved[columns].sort_values(sort).reset_index(drop=True), expected_frame.sort_values(sort).reset_index(drop=True), check_dtype=False, check_exact=False, rtol=0, atol=1e-14)
    frame = pd.DataFrame(completed_rows)
    if not frame.empty:
        require(not frame.duplicated(["drug_id","model","repeat","canonical_model_id"]).any(), "An identity has duplicate OOF predictions within a repeat")
        if not allow_partial:
            for drug,model,repeat in itertools.product(cfg["drug_ids"],cfg["models"],range(cfg["outer_repeats"])):
                block=frame[frame.drug_id.eq(drug)&frame.model.eq(model)&frame["repeat"].eq(repeat)]
                expected_ids=set(manifest[manifest.gdsc_drug_id.eq(str(drug))].canonical_model_id)
                require(set(block.canonical_model_id)==expected_ids and len(block)==len(expected_ids), "Complete repeat lacks exactly one OOF per frozen identity")
    stats_check = None
    if statistics_dir:
        directory=Path(statistics_dir).resolve(); provenance=load(directory/'analysis_provenance.json')
        statistics=load(ROOT/lock['statistics_config'])
        groups=sorted(set(manifest.biological_group_id))
        saved_groups=pd.read_csv(directory/'bootstrap_group_manifest.csv').biological_group_id.tolist()
        require(saved_groups==groups, "Bootstrap group order differs from global frozen biological units")
        with np.load(directory/'bootstrap_group_multiplicities.npz') as archive:
            counts=archive['multiplicities']
        generated=np.random.Generator(np.random.PCG64(statistics['bootstrap']['seed'])).multinomial(len(groups),np.full(len(groups),1/len(groups)),size=statistics['bootstrap']['draws'])
        require(np.array_equal(counts,generated) and np.all(counts.sum(axis=1)==len(groups)), "Bootstrap does not reuse the frozen global subject draws")
        require(hashlib.sha256(counts.tobytes(order='C')).hexdigest()==provenance['bootstrap_count_matrix_sha256'], "Bootstrap count-matrix fingerprint mismatch")
        require(provenance['formal_analysis_lock_sha256']==meta['freeze_sha256'] and provenance['OOF_source_sha256']==meta['output_hashes']['oof_predictions.csv'] and provenance['bootstrap_shared_across_models_repeats_drugs'] is True, "CI source/pairing provenance mismatch")
        intervals=pd.read_csv(directory/'paired_effect_intervals.csv')
        require(intervals[intervals.comparison.eq('primary_EN')&intervals.metric.eq('spearman')].family_size.eq(9).all(), "Primary family-adjusted interval budget altered")
        require(intervals.groupby('gdsc_drug_id').n_unique_biological_groups.max().max()<=len(groups), "Repeat or fold count inflated biological N")
        stats_check={'status':'verified','draws':len(counts),'global_biological_groups':len(groups),'paired_global_subject_draws':True,'repeat_not_added_to_N':True}
    return {"status":"passed" if complete and not allow_partial else "partial_verified", "validated_utc":datetime.now(timezone.utc).isoformat(),
            "execution_commit":revision,"formal_lock_sha256":meta['freeze_sha256'], "cohort_freeze_sha256":meta['cohort_sha256'],
            "original_GDSC_CTRP_response_index_or_expression_read":False,
            "counts":{"completed_tasks_verified":completed,"expected_tasks":840,"actual_fits_verified":fit_count,"expected_fits":18120,"oof_rows_verified":prediction_count,"expected_oof_rows":104307,"unique_biological_groups":len(global_ids)},
            "diagnostics":{"warning_fits_retained":warned_fits,"documented_fallback_fits_retained":fallback_fits,"exact_tissue_constant_target_fits":wholly_constant},
            "statistics_validation":stats_check,
            "scope":"Recorded actual fit memberships and independent prepared-metadata boundaries; schema/path/provenance checks; final aggregate reconstruction. Synthetic tests independently inspect actual feature inputs and paired bootstrap estimands.",
            "limitations":["Recorded inner candidate scores are checked for selection consistency; individual inner OOF predictions are not exported for independent score reconstruction.","Conditional bootstrap does not refit training/tuning or make repeated folds independent biological samples."]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True)
    parser.add_argument('--allow-partial',action='store_true')
    parser.add_argument('--statistics-dir',type=Path)
    args=parser.parse_args(); run=args.run_dir.resolve()
    try:
        report=validate(run,args.allow_partial,args.statistics_dir)
    except Exception as exc:
        report={'status':'failed','validated_utc':datetime.now(timezone.utc).isoformat(),'error_type':type(exc).__name__,'error':str(exc).replace(str(ROOT),'<repository>')}
    target=run/('formal_validation_partial.json' if args.allow_partial else 'formal_validation.json')
    target.write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'status':report['status'],'counts':report.get('counts'),'error':report.get('error')},ensure_ascii=False))
    return 1 if report['status']=='failed' else 0


if __name__=='__main__':
    sys.exit(main())
