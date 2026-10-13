"""Additive adapter for the eight-drug repeated ordinal comparison.

All numerical feature construction and EN/FWL operations are inherited from
the frozen diagnostic implementation. The strong expression reference changes
only the training-selected gene budget to 2000; its column convention, empty
feature handling, tissue adjustment and inner tuning remain identical.
"""
from copy import deepcopy

from ordinal_diagnostic_models import evaluate_diagnostic_outer_fold


def model_spec(model, config):
    if model == 'tissue_null':
        return dict(selection_rule='none', representation='tissue_null', n_genes=0)
    if model == 'expression_strong':
        return dict(selection_rule='mad', representation='expression',
                    n_genes=int(config['feature'].get('strong_expression_genes', 2000)))
    for rule in config['selection_rules']:
        if model.startswith(rule+'_'):
            representation = model[len(rule)+1:]
            if representation not in config['representations']:
                raise ValueError('Unknown repeated representation: '+representation)
            return dict(selection_rule=rule, representation=representation,
                        n_genes=int(config['feature']['n_genes']))
    raise ValueError('Unknown repeated model: '+str(model))


def evaluate_repeated_outer_fold(X, y, groups, train_idx, test_idx, inner_splits,
                                 model_kind, config, stable_gene_ids, audit_sink=None):
    spec = model_spec(model_kind, config)
    adapted = deepcopy(config)
    if model_kind != 'tissue_null':
        adapted['feature']['n_genes'] = spec['n_genes']
        # The frozen diagnostic accepts either key; explicitly make both agree.
        adapted['feature']['diagnostic_genes'] = spec['n_genes']
    def annotate(event):
        event.update(declared_model=model_kind, candidate_gene_budget=spec['n_genes'])
        if audit_sink is not None:
            audit_sink(event)
    result = evaluate_diagnostic_outer_fold(
        X, y, groups, train_idx, test_idx, inner_splits,
        spec['representation'], spec['selection_rule'], adapted, stable_gene_ids,
        audit_sink=annotate)
    result['outer_fit_audit'].update(declared_model=model_kind, candidate_gene_budget=spec['n_genes'])
    if result['feature_audit'] is not None:
        result['feature_audit'].update(declared_model=model_kind, candidate_gene_budget=spec['n_genes'])
    return result
