"""Development-only grouped nested CV with every fitted state inside training.

Input: float64 TPM gene columns, then tissue-code and audit-only sample-ID columns.
Neither audit ID nor tissue code enters gene filtering or REO construction.
No file access is performed by this module.
"""
from __future__ import annotations
import itertools,time,warnings,hashlib
import numpy as np
from scipy.stats import spearmanr
from sklearn.base import BaseEstimator,RegressorMixin
from sklearn.ensemble import RandomForestRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNet,Ridge
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.metrics import mean_absolute_error,mean_squared_error,r2_score
from .features import FoldFeatureBuilder
from .reo_rf import TrainOnlyReoPairSelector

def make_grouped_splits(groups,tissue,n_splits,seed):
    groups=np.asarray(groups);tissue=np.asarray(tissue)
    split=list(StratifiedGroupKFold(n_splits=n_splits,shuffle=True,random_state=seed).split(np.zeros((len(groups),1)),tissue,groups))
    coverage=np.zeros(len(groups),dtype=int)
    for train,test in split:
        assert not(set(groups[train])&set(groups[test]))
        coverage[test]+=1
    assert (coverage==1).all()
    return split

def metrics(y,pred):
    y=np.asarray(y);pred=np.asarray(pred)
    constant_y=bool(np.ptp(y)==0);constant_pred=bool(np.ptp(pred)==0)
    corr=None if constant_y or constant_pred else float(spearmanr(y,pred).statistic)
    return {'n':len(y),'spearman':corr,'spearman_unavailable_reason':'constant_response' if constant_y else 'constant_prediction' if constant_pred else '',
            'mae':float(mean_absolute_error(y,pred)),'rmse':float(np.sqrt(mean_squared_error(y,pred))),
            'r2':None if constant_y else float(r2_score(y,pred)),
            'predictions_outside_0_1':int(((pred<0)|(pred>1)).sum())}

class _ZeroMolecularRegressor:
    """Explicit zero molecular solution when training alpha_max equals zero."""
    def __init__(self,n_features):
        self.coef_=np.zeros(n_features,dtype=np.float64)
        self.n_iter_=0;self.dual_gap_=0.0


class TrainOnlyRegressor(RegressorMixin,BaseEstimator):
    """Learn every state from the rows supplied to the current training fit.

    EN uses fold-fitted standardization and unpenalized tissue effects through
    FWL projection of already-created molecular features. With config mode
    ``alpha_max_relative``, the current training residuals alone set
    ``alpha_max=max(abs(Zr.T @ yr))/(n*l1_ratio)``; ``alpha_ratio`` selects a
    fraction of that value. Legacy configs retain their absolute-alpha path.
    At ratio exactly one, the analytically valid zero molecular solution is
    used explicitly: solver roundoff must not break ties within a tissue and
    create a spurious molecular ranking. We do not round weaker-path solutions.

    RF instead receives original log-expression or original three-value REO
    columns plus the same training-fitted tissue one-hot design. REO RF selects
    at most the configured number of pairs by within-tissue training Pearson
    association. The RF input itself is never residualized. Both RF variants
    use the identical fixed tree count, seed and grid/tissue policy. An unknown
    tissue has all-zero one-hot columns; molecular predictors remain active.

    The final input columns are tissue code and audit-only unique sample ID;
    neither enters gene selection. This class performs no file access.
    """
    def __init__(self,model_kind='reo',config=None,stable_gene_ids=None,alpha=0.1,l1_ratio=0.5,
                 alpha_ratio=None,max_features='sqrt',min_samples_leaf=1):
        self.model_kind=model_kind;self.config=config;self.stable_gene_ids=stable_gene_ids
        self.alpha=alpha;self.l1_ratio=l1_ratio;self.alpha_ratio=alpha_ratio
        self.max_features=max_features;self.min_samples_leaf=min_samples_leaf

    def fit(self,X,y):
        start=time.perf_counter();X=np.asarray(X);y=np.asarray(y,dtype=np.float64)
        assert X.ndim==2 and X.shape[0]==len(y) and X.shape[1]==len(self.stable_gene_ids)+2
        assert len(y)>0 and np.isfinite(y).all() and len(np.unique(X[:,-1]))==len(y)
        if self.model_kind not in {'mean','tissue','reo','expression','expression_matched','reo_rf','expression_rf'}:
            raise ValueError('Unsupported model_kind: '+str(self.model_kind))
        self.fit_sample_ids_=X[:,-1].astype(int).tolist()
        self.target_mean_=float(y.mean());self.target_sd_=float(y.std())
        self.fallback_reason_='';self.warning_messages_=[];self.builder_=None;self.n_iter_=0
        self.pair_selector_=None;self.predict_seconds_=0.0;self.last_predict_seconds_=0.0;self.prediction_calls_=0
        self.actual_alpha_=None;self.alpha_max_=None;self.path_ratio_=None;self.dual_gap_=None
        self.null_endpoint_=False;self.solver_reason_='not_applicable'
        self.n_nonzero_molecular_coeffs_=None;self.n_molecular_features_used_=None
        self.rf_molecular_feature_names_=();self.rf_tissue_policy_=None
        if self.target_sd_==0:self.target_sd_=1.0;self.fallback_reason_='constant_training_response'
        ys=(y-self.target_mean_)/self.target_sd_
        self.encoder_=OneHotEncoder(handle_unknown='ignore',sparse_output=False,dtype=np.float64)
        d=self.encoder_.fit_transform(X[:,-2].reshape(-1,1))
        self.training_tissue_codes_=self.encoder_.categories_[0].tolist()
        self.n_features_in_=X.shape[1]
        if self.model_kind=='mean' or self.fallback_reason_:self.fitted_kind_='mean'
        elif self.model_kind=='tissue':
            self.fitted_kind_='tissue';self.actual_alpha_=float(self.alpha)
            self.regressor_=Ridge(alpha=self.alpha,fit_intercept=False,solver='svd').fit(d,ys)
        else:
            feature=self.config['feature'];kind='reo' if self.model_kind in {'reo','reo_rf'} else 'expression'
            if self.model_kind=='expression_matched':
                n_genes=feature.get('matched_expression_genes',feature['reo_genes'])
            else:n_genes=feature['expression_genes'] if kind=='expression' else feature['reo_genes']
            self.builder_=FoldFeatureBuilder(kind=kind,n_genes=n_genes,min_tpm=feature['min_tpm'],min_detection_fraction=feature['min_detection_fraction'],max_pair_tie_fraction=feature['max_pair_tie_fraction'],stable_gene_ids=self.stable_gene_ids)
            try:
                z=self.builder_.fit_transform(X[:,:-2])
            except ValueError as exc:
                # Only a documented empty/constant feature set can fall back.
                if not any(text in str(exc) for text in ['No genes','No REO pairs','REO requires at least two']):raise
                self.fitted_kind_='mean';self.fallback_reason_='empty_features:'+str(exc)
            else:
                self.fitted_kind_=self.model_kind
                if self.model_kind.endswith('_rf'):
                    self._fit_random_forest(z,d,ys)
                else:
                    self._fit_elastic_net(z,d,ys)
        self.fit_seconds_=time.perf_counter()-start
        return self

    def _fit_elastic_net(self,z,d,ys):
        self.solver_reason_='elastic_net_coordinate_descent'
        self.scaler_=StandardScaler().fit(z);z=self.scaler_.transform(z)
        counts=d.sum(axis=0)
        self.feature_tissue_means_=(d.T@z)/counts[:,None]
        self.target_tissue_means_=(d.T@ys)/counts
        zr=z-d@self.feature_tissue_means_;yr=ys-d@self.target_tissue_means_
        e=self.config['elastic_net'];mode=e.get('regularization_mode','absolute')
        if mode=='alpha_max_relative':
            if not 0<float(self.l1_ratio)<=1:
                raise ValueError('alpha_max_relative requires l1_ratio in (0,1].')
            ratio=1.0 if self.alpha_ratio is None else float(self.alpha_ratio)
            if not np.isfinite(ratio) or not 0<ratio<=1:
                raise ValueError('alpha_ratio must be in (0,1].')
            self.path_ratio_=ratio
            self.alpha_max_=float(np.max(np.abs(zr.T@yr))/(len(yr)*self.l1_ratio))
            self.actual_alpha_=ratio*self.alpha_max_
            if self.alpha_max_==0.0:
                self.regressor_=_ZeroMolecularRegressor(zr.shape[1])
                self.fallback_reason_='zero_training_alpha_max_molecular_coefficients_zero'
                self.solver_reason_='zero_training_alpha_max'
                self.null_endpoint_=(ratio==1.0)
            elif ratio==1.0:
                self.regressor_=_ZeroMolecularRegressor(zr.shape[1])
                self.null_endpoint_=True
                self.solver_reason_='theoretical_alpha_max_null_endpoint'
        elif mode=='absolute':self.actual_alpha_=float(self.alpha)
        else:raise ValueError('Unknown elastic_net regularization_mode: '+str(mode))
        if self.alpha_max_!=0.0 and not self.null_endpoint_:
            self.regressor_=ElasticNet(alpha=self.actual_alpha_,l1_ratio=self.l1_ratio,fit_intercept=False,max_iter=e['max_iter'],tol=e['tol'],selection=e.get('selection','cyclic'))
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always',ConvergenceWarning);self.regressor_.fit(zr,yr)
            self.warning_messages_=[str(w.message) for w in caught]
        self.n_iter_=int(self.regressor_.n_iter_)
        self.dual_gap_=float(self.regressor_.dual_gap_)
        self.n_nonzero_molecular_coeffs_=int(np.count_nonzero(self.regressor_.coef_))
        self.n_molecular_features_used_=int(zr.shape[1])

    def _fit_random_forest(self,z,d,ys):
        rf=self.config['random_forest']
        names=tuple(self.builder_.get_feature_names_out())
        if self.model_kind=='reo_rf':
            self.pair_selector_=TrainOnlyReoPairSelector(max_pairs=rf.get('reo_max_pairs',500),association=rf.get('association','within_tissue_pearson'))
            self.pair_selector_.fit(z,ys,tissue_design=d,feature_names=names)
            z=self.pair_selector_.transform(z)
            names=self.pair_selector_.selected_feature_names_
            if self.pair_selector_.fallback_reason_:
                self.fallback_reason_='rf_tissue_only:'+self.pair_selector_.fallback_reason_
        self.rf_molecular_feature_names_=names
        self.n_molecular_features_used_=int(z.shape[1])
        self.rf_tissue_policy_='original molecular features plus training-fitted full tissue one-hot; unseen tissue columns all zero'
        seed=int(rf.get('seed',self.config.get('seed',0)))
        if rf.get('n_jobs',1)!=1:raise ValueError('Pilot RF must use n_jobs=1 for resource comparability.')
        self.regressor_=RandomForestRegressor(n_estimators=rf.get('n_estimators',300),max_features=self.max_features,min_samples_leaf=self.min_samples_leaf,n_jobs=1,random_state=seed)
        self.regressor_.fit(np.column_stack((z,d)),ys)

    def predict(self,X):
        start=time.perf_counter();X=np.asarray(X)
        assert X.shape[1]==self.n_features_in_
        if self.fitted_kind_=='mean':pred=np.full(len(X),self.target_mean_)
        else:
            d=self.encoder_.transform(X[:,-2].reshape(-1,1))
            if self.fitted_kind_=='tissue':yp=self.regressor_.predict(d)
            elif self.fitted_kind_.endswith('_rf'):
                z=self.builder_.transform(X[:,:-2])
                if self.pair_selector_ is not None:z=self.pair_selector_.transform(z)
                yp=self.regressor_.predict(np.column_stack((z,d)))
            else:
                z=self.scaler_.transform(self.builder_.transform(X[:,:-2]))
                yp=(z-d@self.feature_tissue_means_)@self.regressor_.coef_+d@self.target_tissue_means_
            pred=self.target_mean_+self.target_sd_*yp
        self.last_predict_seconds_=time.perf_counter()-start
        self.predict_seconds_+=self.last_predict_seconds_;self.prediction_calls_+=1
        return pred

    def feature_audit(self):
        if self.builder_ is None or not hasattr(self.builder_,'audit_metadata_'):return None
        result=self.builder_.get_audit_metadata(include_feature_names=True)
        if self.model_kind.endswith('_rf'):
            rf=self.config['random_forest']
            result['random_forest']={'molecular_feature_names':list(self.rf_molecular_feature_names_),
                'molecular_features_used':self.n_molecular_features_used_,'tissue_policy':self.rf_tissue_policy_,
                'n_estimators':int(rf.get('n_estimators',300)),'n_jobs':1,
                'seed':int(rf.get('seed',self.config.get('seed',0))),
                'max_features':self.max_features,'min_samples_leaf':int(self.min_samples_leaf)}
            if self.pair_selector_ is not None:
                pair_audit=self.pair_selector_.get_audit_metadata()
                pair_audit['selected_pairs_original_input_gene_indices']=self.builder_.pair_indices_[self.pair_selector_.selected_indices_].tolist()
                result['random_forest']['pair_selection']=pair_audit
        return result

    def audit(self):
        built_features=None if self.builder_ is None or not hasattr(self.builder_,'n_output_features_') else int(self.builder_.n_output_features_)
        n_features=self.n_molecular_features_used_ if self.model_kind.endswith('_rf') else built_features
        return {'fit_sample_ids':self.fit_sample_ids_,'model_kind':self.model_kind,'fitted_kind':self.fitted_kind_,'fallback_reason':self.fallback_reason_,
                'target_mean':self.target_mean_,'target_sd':self.target_sd_,'training_tissue_codes':self.training_tissue_codes_,
                'fit_seconds':self.fit_seconds_,'predict_seconds':self.predict_seconds_,'last_predict_seconds':self.last_predict_seconds_,'prediction_calls':self.prediction_calls_,
                'n_iter':self.n_iter_,'warnings':self.warning_messages_,
                'actual_alpha':self.actual_alpha_,'alpha_max':self.alpha_max_,'path_ratio':self.path_ratio_,
                'null_endpoint':self.null_endpoint_,'solver_reason':self.solver_reason_,
                'n_nonzero_molecular_coeffs':self.n_nonzero_molecular_coeffs_,'dual_gap':self.dual_gap_,
                'n_selected_genes':None if self.builder_ is None or not hasattr(self.builder_,'gene_indices_') else len(self.builder_.gene_indices_),
                'n_features':n_features,'n_candidate_molecular_features':built_features,
                'n_molecular_features_used':self.n_molecular_features_used_,'n_tissue_features':len(self.training_tissue_codes_),
                'rf_tissue_policy':self.rf_tissue_policy_,
                'rf_pair_selection_reason':None if self.pair_selector_ is None else self.pair_selector_.fallback_reason_,
                'selected_gene_ids_sha256':None if self.builder_ is None or not hasattr(self.builder_,'gene_ids_') else hashlib.sha256('\0'.join(self.builder_.gene_ids_).encode()).hexdigest(),
                'built_feature_names_sha256':None if self.builder_ is None or not hasattr(self.builder_,'feature_names_') else hashlib.sha256('\0'.join(self.builder_.feature_names_).encode()).hexdigest(),
                'rf_selected_pair_names_sha256':None if self.pair_selector_ is None else hashlib.sha256('\0'.join(self.pair_selector_.selected_feature_names_).encode()).hexdigest()}

def candidate_grid(model_kind,config):
    if model_kind=='mean':return [{}]
    if model_kind=='tissue':return [{'alpha':a} for a in config['ridge_alphas']]
    if model_kind.endswith('_rf'):
        rf=config['random_forest']
        return [{'max_features':f,'min_samples_leaf':s} for f,s in itertools.product(rf['max_features'],rf['min_samples_leaf'])]
    if config['elastic_net'].get('regularization_mode','absolute')=='alpha_max_relative':
        return [{'alpha_ratio':a,'l1_ratio':r} for a,r in itertools.product(config['elastic_net']['alpha_ratios'],config['elastic_net']['l1_ratios'])]
    return [{'alpha':a,'l1_ratio':r} for a,r in itertools.product(config['elastic_net']['alphas'],config['elastic_net']['l1_ratios'])]


def _grid_boundary_details(best_params,grid):
    """Report boundaries for each varied numerical grid dimension.

    RF categorical max_features choices are reported as categorical choices,
    not endpoints. A one-value grid cannot diagnose a search boundary.
    """
    details={}
    for key,value in best_params.items():
        values=[candidate[key] for candidate in grid]
        unique=list(dict.fromkeys(values))
        numeric=all(isinstance(x,(int,float,np.integer,np.floating)) and not isinstance(x,bool) for x in unique)
        if numeric:
            details[key]={'chosen':value,'minimum':min(unique),'maximum':max(unique),'at_boundary':bool(len(unique)>1 and value in (min(unique),max(unique))),
                          'search_extent_diagnostic':key in {'alpha','alpha_ratio','min_samples_leaf'}}
        else:details[key]={'chosen':value,'choices':unique,'at_boundary':False,'search_extent_diagnostic':False}
    return details

def evaluate_outer_fold(X,y,groups,train_idx,test_idx,inner_splits,model_kind,config,stable_gene_ids,audit_sink=None):
    """Choose parameters using inner outcomes only, then evaluate outer once.

    inner_splits index the outer training rows. Every fit has its own feature
    builder, encoder, scalers and target/tissue state. Training/test ID checks
    are enforced independently of any response. ``audit_sink`` is called once
    per fit after that model's corresponding prediction, so its event includes
    prediction time. Relative-alpha candidates recompute alpha_max independently
    in each inner-training fit and in the final outer-training refit; no value
    is shared across partitions. Parameter selection uses pooled inner RMSE.
    """
    X=np.asarray(X);y=np.asarray(y);groups=np.asarray(groups)
    train_idx=np.asarray(train_idx);test_idx=np.asarray(test_idx)
    assert not(set(groups[train_idx])&set(groups[test_idx]))
    assert not(set(X[train_idx,-1])&set(X[test_idx,-1]))
    xt=X[train_idx];yt=y[train_idx];gt=groups[train_idx]
    events=[];scores=[];best_score=float('inf');best_params=None
    def fit_one(indices,params,stage,candidate,inner_fold):
        model=TrainOnlyRegressor(model_kind,config,stable_gene_ids,**params)
        model.fit(xt[indices],yt[indices])
        event={**model.audit(),'stage':stage,'candidate':candidate,'inner_fold':inner_fold,'params':params}
        assert set(event['fit_sample_ids'])<=set(xt[:,-1].astype(int))
        assert not(set(event['fit_sample_ids'])&set(X[test_idx,-1].astype(int)))
        events.append(event)
        return model
    def finish_prediction(model):
        events[-1].update(model.audit())
        if audit_sink is not None:audit_sink(events[-1])
    grid=candidate_grid(model_kind,config)
    if model_kind=='mean':best_params={};best_score=None
    else:
        for ci,params in enumerate(grid):
            oof=np.full(len(yt),np.nan);covered=np.zeros(len(yt),dtype=int);invalid=[]
            for ii,(it,iv) in enumerate(inner_splits):
                it=np.asarray(it);iv=np.asarray(iv)
                assert not(set(gt[it])&set(gt[iv]))
                assert not(set(it)&set(iv))
                model=fit_one(it,params,'inner_training',ci,ii)
                assert set(model.fit_sample_ids_)==set(xt[it,-1].astype(int))
                oof[iv]=model.predict(xt[iv]);covered[iv]+=1
                finish_prediction(model)
                if model.fallback_reason_:invalid.append(model.fallback_reason_)
            assert (covered==1).all() and np.isfinite(oof).all()
            score=float(np.sqrt(mean_squared_error(yt,oof)))
            scores.append({'candidate':ci,'params':params,'pooled_inner_rmse':score,'fallbacks':invalid})
            if score<best_score:best_score=score;best_params=params.copy()
    model=fit_one(np.arange(len(yt)),best_params,'outer_training_refit',-1,-1)
    pred=model.predict(X[test_idx])
    finish_prediction(model)
    feature_audit=model.feature_audit()
    boundary_details=_grid_boundary_details(best_params,grid)
    # This is the first use of outer-test outcomes, after all choices and predictions.
    result={'chosen_params':best_params,'best_inner_rmse':best_score,'predictions':pred,'metrics':metrics(y[test_idx],pred),
            'inner_scores':scores,'feature_audit':feature_audit,'fit_events':events,
            'outer_unknown_tissue_samples':int((~np.isin(X[test_idx,-2],model.training_tissue_codes_)).sum()),
            'chosen_grid_boundary':bool(any(d['at_boundary'] and d['search_extent_diagnostic'] for d in boundary_details.values())),
            'chosen_grid_boundary_details':boundary_details,
            'outer_fit_audit':model.audit()}
    return result
