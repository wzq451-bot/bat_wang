import numpy as np
from scipy.stats import t
def design(d,cols,group='ind'):
 ids=sorted(d[group].unique());terms=['Intercept']+cols+[str(i) for i in ids[1:]]
 return np.column_stack([np.ones(len(d))]+[d[c].to_numpy(float) for c in cols]+[d[group].eq(i).to_numpy(float) for i in ids[1:]]),terms
def ols(X,y,g,weights=None):
 n,k=X.shape;assert np.linalg.matrix_rank(X)==k;w=np.ones(n) if weights is None else np.asarray(weights,float)
 bread=np.linalg.inv(X.T@(w[:,None]*X));b=np.linalg.lstsq(X*np.sqrt(w[:,None]),y*np.sqrt(w),rcond=None)[0];e=y-X@b;ids=np.unique(g);G=len(ids)
 scores=np.array([X[g==j].T@(w[g==j]*e[g==j]) for j in ids]);cc=bread@scores.T@scores@bread*G/(G-1)*(n-1)/(n-k)
 return dict(beta=b,cov_cluster=cc,cov_model=bread*(np.sum(w*e*e)/(n-k)),fitted=X@b,residuals=e,n=n,k=k,G=G)
def stats(r,j=1):
 b=r['beta'][j];se=np.sqrt(r['cov_cluster'][j,j]);q=t.ppf(.975,r['G']-1)
 return dict(beta=float(b),cluster_se=float(se),cluster_lo=float(b-q*se),cluster_hi=float(b+q*se),cluster_p=float(2*t.sf(abs(b/se),r['G']-1)),n=r['n'],n_bats=r['G'])
