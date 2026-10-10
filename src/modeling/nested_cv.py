"""Development-only grouped nested CV with every fitted state inside training.

Input: float64 TPM gene columns, then tissue-code and audit-only sample-ID columns.
Neither audit ID nor tissue code enters gene filtering or REO construction.
No file access is performed by this module.
"""
from __future__ import annotations
import itertools,time,warnings
import numpy as np
from scipy.stats import spearmanr
from sklearn.base import BaseEstimator,RegressorMixin
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNet,Ridge
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.metrics import mean_absolute_error,mean_squared_error,r2_score
from .features import FoldFeatureBuilder

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

class TrainOnlyRegressor(RegressorMixin,BaseEstimator):
    """Training-only features, scaling and tissue-adjusted regression.

    Tissue terms are unpenalized for both EN representations. Full one-hot on
    centered target gives unseen tissue an explicit zero effect (global mean).
    FWL projects already-created features, never alters TPM before REO encoding.
    """
    def __init__(self,model_kind='reo',config=None,stable_gene_ids=None,alpha=0.1,l1_ratio=0.5):
        self.model_kind=model_kind;self.config=config;self.stable_gene_ids=stable_gene_ids
        self.alpha=alpha;self.l1_ratio=l1_ratio

    def fit(self,X,y):
        start=time.perf_counter();X=np.asarray(X);y=np.asarray(y,dtype=np.float64)
        assert X.ndim==2 and X.shape[0]==len(y) and X.shape[1]==len(self.stable_gene_ids)+2
        assert np.isfinite(y).all() and len(np.unique(X[:,-1]))==len(y)
        self.fit_sample_ids_=X[:,-1].astype(int).tolist()
        self.target_mean_=float(y.mean());self.target_sd_=float(y.std())
        self.fallback_reason_='';self.warning_messages_=[];self.builder_=None;self.n_iter_=0
        if self.target_sd_==0:self.target_sd_=1.0;self.fallback_reason_='constant_training_response'
        ys=(y-self.target_mean_)/self.target_sd_
        self.encoder_=OneHotEncoder(handle_unknown='ignore',sparse_output=False,dtype=np.float64)
        d=self.encoder_.fit_transform(X[:,-2].reshape(-1,1))
        self.training_tissue_codes_=self.encoder_.categories_[0].tolist()
        self.n_features_in_=X.shape[1]
        if self.model_kind=='mean' or self.fallback_reason_:self.fitted_kind_='mean'
        elif self.model_kind=='tissue':
            self.fitted_kind_='tissue';self.regressor_=Ridge(alpha=self.alpha,fit_intercept=False,solver='svd').fit(d,ys)
        else:
            feature=self.config['feature'];kind='reo' if self.model_kind=='reo' else 'expression'
            n_genes=feature['expression_genes'] if self.model_kind=='expression' else feature['reo_genes']
            self.builder_=FoldFeatureBuilder(kind=kind,n_genes=n_genes,min_tpm=feature['min_tpm'],min_detection_fraction=feature['min_detection_fraction'],max_pair_tie_fraction=feature['max_pair_tie_fraction'],stable_gene_ids=self.stable_gene_ids)
            try:
                z=self.builder_.fit_transform(X[:,:-2])
            except ValueError as exc:
                # Only a documented empty/constant feature set can fall back.
                if not any(text in str(exc) for text in ['No genes','No REO pairs','REO requires at least two']):raise
                self.fitted_kind_='mean';self.fallback_reason_='empty_features:'+str(exc)
            else:
                self.fitted_kind_=self.model_kind
                self.scaler_=StandardScaler().fit(z);z=self.scaler_.transform(z)
                counts=d.sum(axis=0)
                self.feature_tissue_means_=(d.T@z)/counts[:,None]
                self.target_tissue_means_=(d.T@ys)/counts
                zr=z-d@self.feature_tissue_means_;yr=ys-d@self.target_tissue_means_
                e=self.config['elastic_net']
                self.regressor_=ElasticNet(alpha=self.alpha,l1_ratio=self.l1_ratio,fit_intercept=False,max_iter=e['max_iter'],tol=e['tol'],selection=e.get('selection','cyclic'))
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always',ConvergenceWarning);self.regressor_.fit(zr,yr)
                self.warning_messages_=[str(w.message) for w in caught]
                self.n_iter_=int(self.regressor_.n_iter_)
        self.fit_seconds_=time.perf_counter()-start
        return self

    def predict(self,X):
        X=np.asarray(X)
        assert X.shape[1]==self.n_features_in_
        if self.fitted_kind_=='mean':return np.full(len(X),self.target_mean_)
        d=self.encoder_.transform(X[:,-2].reshape(-1,1))
        if self.fitted_kind_=='tissue':yp=self.regressor_.predict(d)
        else:
            z=self.scaler_.transform(self.builder_.transform(X[:,:-2]))
            yp=(z-d@self.feature_tissue_means_)@self.regressor_.coef_+d@self.target_tissue_means_
        return self.target_mean_+self.target_sd_*yp

    def audit(self):
        return {'fit_sample_ids':self.fit_sample_ids_,'model_kind':self.model_kind,'fitted_kind':self.fitted_kind_,'fallback_reason':self.fallback_reason_,
                'target_mean':self.target_mean_,'target_sd':self.target_sd_,'training_tissue_codes':self.training_tissue_codes_,
                'fit_seconds':self.fit_seconds_,'n_iter':self.n_iter_,'warnings':self.warning_messages_,
                'n_selected_genes':None if self.builder_ is None or not hasattr(self.builder_,'gene_indices_') else len(self.builder_.gene_indices_),
                'n_features':None if self.builder_ is None or not hasattr(self.builder_,'n_output_features_') else self.builder_.n_output_features_}

def candidate_grid(model_kind,config):
    if model_kind=='mean':return [{}]
    if model_kind=='tissue':return [{'alpha':a} for a in config['ridge_alphas']]
    return [{'alpha':a,'l1_ratio':r} for a,r in itertools.product(config['elastic_net']['alphas'],config['elastic_net']['l1_ratios'])]

def evaluate_outer_fold(X,y,groups,train_idx,test_idx,inner_splits,model_kind,config,stable_gene_ids,audit_sink=None):
    """Choose parameters using inner outcomes only, then evaluate outer once.

    inner_splits index the outer training rows. Every fit has its own feature
    builder, encoder, scalers and target/tissue state. Training/test ID checks
    are enforced independently of any response.
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
        if audit_sink is not None:audit_sink(event)
        return model
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
                if model.fallback_reason_:invalid.append(model.fallback_reason_)
            assert (covered==1).all() and np.isfinite(oof).all()
            score=float(np.sqrt(mean_squared_error(yt,oof)))
            scores.append({'candidate':ci,'params':params,'pooled_inner_rmse':score,'fallbacks':invalid})
            if score<best_score:best_score=score;best_params=params.copy()
    model=fit_one(np.arange(len(yt)),best_params,'outer_training_refit',-1,-1)
    pred=model.predict(X[test_idx])
    feature_audit=None if model.builder_ is None or not hasattr(model.builder_,'audit_metadata_') else model.builder_.get_audit_metadata(include_feature_names=True)
    # This is the first use of outer-test outcomes, after all choices and predictions.
    result={'chosen_params':best_params,'best_inner_rmse':best_score,'predictions':pred,'metrics':metrics(y[test_idx],pred),
            'inner_scores':scores,'feature_audit':feature_audit,'fit_events':events,
            'outer_unknown_tissue_samples':int((~np.isin(X[test_idx,-2],model.training_tissue_codes_)).sum()),
            'chosen_grid_boundary':bool(model_kind!='mean' and (best_params['alpha']==max(x['alpha'] for x in grid) or best_params['alpha']==min(x['alpha'] for x in grid))),
            'outer_fit_audit':model.audit()}
    return result
