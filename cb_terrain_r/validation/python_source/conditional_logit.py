"""Conditional logit for one observed endpoint per fixed-size choice stratum.

Inputs X have shape (n_steps, n_choices, n_terms); choice zero is observed.
No intercept or stratum-constant main effect is identifiable.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import norm,t

def objective(beta,X):
 eta=np.einsum('nkp,p->nk',X,beta)
 lse=logsumexp(eta,axis=1);prob=np.exp(eta-lse[:,None])
 mean=np.einsum('nk,nkp->np',prob,X)
 score=X[:,0,:]-mean
 return float(np.mean(lse-eta[:,0])),-score.mean(axis=0)

def information(beta,X):
 eta=np.einsum('nkp,p->nk',X,beta);prob=np.exp(eta-logsumexp(eta,axis=1)[:,None])
 mean=np.einsum('nk,nkp->np',prob,X)
 second=np.einsum('nk,nkp,nkq->pq',prob,X,X)
 h=second-mean.T@mean
 return (h+h.T)/2,X[:,0,:]-mean

def fit(X,clusters,start=None):
 X=np.asarray(X,dtype=float);clusters=np.asarray(clusters)
 if not np.isfinite(X).all():raise ValueError('Nonfinite model input')
 # Remove constants within strata for numerical stability; likelihood is unchanged.
 X=X-X.mean(axis=1,keepdims=True)
 opt=minimize(objective,np.zeros(X.shape[2]) if start is None else start,args=(X,),jac=True,
              method='BFGS',options={'gtol':1e-9,'maxiter':400})
 grad=float(np.max(abs(objective(opt.x,X)[1])))
 if not opt.success and grad>1e-7:raise RuntimeError(f'Fit failed: {opt.message}; mean score {grad}')
 H,scores=information(opt.x,X);ev=np.linalg.eigvalsh(H)
 if ev.min()<=0:raise RuntimeError('Nonpositive information matrix; rank deficiency or separation')
 cov=np.linalg.inv(H);groups=np.unique(clusters);G=len(groups);n=len(X);p=X.shape[-1]
 grouped=np.stack([scores[clusters==g].sum(axis=0) for g in groups])
 if G>1:
  cr1=G/(G-1)*(n-1)/(n-p);cc=cov@(grouped.T@grouped)@cov*cr1
 else:cc=np.full_like(cov,np.nan)
 return {'beta':opt.x,'cov_model':cov,'cov_cluster':cc,'n_steps':n,'n_clusters':G,
         'loglik':-objective(opt.x,X)[0]*n,'converged':bool(opt.success or grad<=1e-7),
         'mean_score_max_abs':grad,'hessian_condition':float(ev.max()/ev.min()),
         'optimizer_message':str(opt.message)}

def coefficient_rows(result,terms,model_id):
 rows=[]
 for j,term in enumerate(terms):
  b=result['beta'][j];se=np.sqrt(result['cov_model'][j,j]);sc=np.sqrt(max(0,result['cov_cluster'][j,j]));df=result['n_clusters']-1
  q=t.ppf(.975,df) if df>0 else np.nan
  rows.append({'model':model_id,'term':term,'beta':b,'se':se,'p_normal':2*norm.sf(abs(b/se)),
   'ci_lo':b-norm.ppf(.975)*se,'ci_hi':b+norm.ppf(.975)*se,'cluster_se':sc,
   'cluster_p_t':2*t.sf(abs(b/sc),df) if sc>0 and df>0 else np.nan,
   'cluster_ci_lo':b-q*sc,'cluster_ci_hi':b+q*sc,'n_steps':result['n_steps'],'n_clusters':result['n_clusters']})
 return rows
