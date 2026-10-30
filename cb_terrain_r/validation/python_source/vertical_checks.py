"""Reconstructed vertical analysis; explicitly replaces, rather than reproduces, undocumented old null."""
from pathlib import Path
import json
import numpy as np,pandas as pd
from model_utils import ols,stats
R=Path(__file__).resolve().parents[1];I=R/'input';O=R/'results'
m=pd.read_csv(I/'steps.csv');m['t']=pd.to_datetime(m.t)
end=m[['ind','t','end_elev']].rename(columns={'end_elev':'ground'});start=m[['ind','t','dt_min','start_elev']].copy();start.t=(start.t-pd.to_timedelta(start.dt_min,unit='min')).dt.round('s');start=start.rename(columns={'start_elev':'ground'})[['ind','t','ground']]
g=pd.concat([end,start]);assert g.groupby(['ind','t']).ground.agg(lambda x:x.max()-x.min()).max()<1e-6;g=g.drop_duplicates(['ind','t'])
f=pd.read_csv(I/'fixes.csv');f.t=pd.to_datetime(f.t);f=f.merge(g,on=['ind','t'],how='left',validate='one_to_one');f.to_csv(O/'ground_recovery.csv.gz',index=False)
segs=pd.read_csv(I/'segments.csv');d=f[f.seg.isin(segs.seg)&f.altitude_valid&f.ground.notna()].copy().sort_values(['ind','seg','t']).reset_index(drop=True)
d['x']=d.ground-d.groupby('seg').ground.transform('mean');d['y']=d.gps_alt-d.groupby('seg').gps_alt.transform('mean');r=ols(d.x.to_numpy()[:,None],d.y.to_numpy(),d.ind.to_numpy());rows=[dict(model='within',**stats(r,0))]
means=d.groupby(['ind','seg']).agg(ground=('ground','mean'),alt=('gps_alt','mean')).reset_index();r=ols(np.column_stack([np.ones(len(means)),means.ground]),means.alt.to_numpy(),means.ind.to_numpy());rows.append(dict(model='between',**stats(r,1)))
pd.DataFrame(rows).to_csv(O/'vertical_refit.csv',index=False)
ss=[]
for (ind,seg),a in d.groupby(['ind','seg']):ss.append(dict(ind=ind,seg=seg,n=len(a),beta=np.dot(a.x,a.y)/np.dot(a.x,a.x),ground_sd=a.ground.std(ddof=1)))
pd.DataFrame(ss).to_csv(O/'vertical_segments.csv',index=False);d.to_csv(O/'vertical_analysis.csv.gz',index=False)
rng=np.random.default_rng(20261027);alt=d.gps_alt.to_numpy();ground=d.ground.to_numpy();indices=[a.index.to_numpy() for _,a in d.groupby('seg')];bats=[a.index.to_numpy() for _,a in d.groupby('ind')]
def score(a):return np.mean([np.std((a-ground)[ix],ddof=1) for ix in bats])
obs=score(alt);null=[]
for rep in range(2000):
 shuffled=alt.copy()
 for ix in indices:
  q=int(np.ceil(.1*len(ix)));lag=int(rng.integers(q,len(ix)-q+1));shuffled[ix]=np.roll(alt[ix],lag)
 null.append(dict(replicate=rep+1,mean_SD_AGL=score(shuffled),within_beta=np.dot(d.x,shuffled)/np.dot(d.x,d.x)))
pd.DataFrame(null).to_csv(O/'vertical_null.csv',index=False);v=np.array([a['mean_SD_AGL'] for a in null]);summary=dict(n_fixes=len(d),n_bats=d.ind.nunique(),n_segments=d.seg.nunique(),seed=20261027,n_randomizations=2000,observed_mean_SD_AGL=obs,null_median=np.median(v),null_lo=np.quantile(v,.025),null_hi=np.quantile(v,.975),p_lower=(1+sum(v<=obs))/2001,all_segment_slopes_positive=all(a['beta']>0 for a in ss),ground_recovered=int(f.ground.notna().sum()),ground_missing=int(f.ground.isna().sum()))
(O/'vertical_summary.json').write_text(json.dumps(summary,indent=2));print(pd.DataFrame(rows).to_string(index=False));print(summary)
