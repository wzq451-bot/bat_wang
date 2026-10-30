"""Documented additional diagnostics; these are new specifications, not old gates."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from conditional_logit import fit,coefficient_rows
R=Path(__file__).resolve().parents[1];O=R/'results/refit'
m=pd.read_csv(R/'results/master_validated.csv').sort_values('step_id').reset_index(drop=True)
q=pd.read_csv(O/'matched_choices.csv.gz');scale=json.loads((O/'scaling.json').read_text())['predictors']
x=np.stack([((q[k]-scale[k]['mean'])/scale[k]['sd']).to_numpy().reshape(-1,41) for k in ['tpi','elev','tri']],axis=2)
z=(m.start_agl.to_numpy()-scale['start_agl']['mean'])/scale['start_agl']['sd'];X=np.concatenate([x,(x[:,:,0]*z[:,None])[:,:,None]],axis=2)
bats=m.ind.to_numpy();terms=['tpi_z','elev_z','tri_z','agl_z:tpi_z'];rows=[];diags={}
def one(name,X,terms,mask=None):
 mask=np.ones(len(m),bool) if mask is None else np.asarray(mask)
 r=fit(X[mask],bats[mask]);rows.extend(coefficient_rows(r,terms,name))
 np.savez_compressed(O/(name+'_fit.npz'),beta=r['beta'],cov_model=r['cov_model'],cov_cluster=r['cov_cluster'],terms=np.asarray(terms),included_step_ids=m.loc[mask,'step_id'].to_numpy())
 diags[name]={k:v for k,v in r.items() if k not in ['beta','cov_model','cov_cluster']};print(name,np.round(r['beta'],5),flush=True)
 return r
# Clock time recorded in the supplied table; do not assert astronomical night phase.
t=pd.to_datetime(m.t);hours=t.dt.hour+t.dt.minute/60+t.dt.second/3600
nh=(hours-18)%24;hz=((nh-nh.mean())/nh.std(ddof=1)).to_numpy();h2=hz**2;h2-=h2.mean()
Xt=np.concatenate([X,x*hz[:,None,None],x*h2[:,None,None]],axis=2)
one('clock_time_terrain_interactions',Xt,terms+[f'clock_z:{v}' for v in terms[:3]]+[f'clock_z2:{v}' for v in terms[:3]])
one('before_midnight',X,terms,hours>=18);one('after_midnight',X,terms,hours<12)
# Different baseline terrain responses for every bat; common within-bat height interaction.
wc=m.start_agl-m.groupby('ind').start_agl.transform('mean');wz=wc.to_numpy()/scale['start_agl']['sd']
individuals=sorted(m.ind.unique());XX=np.concatenate([x*(bats==b)[:,None,None] for b in individuals]+[(x[:,:,0]*wz[:,None])[:,:,None]],axis=2)
tt=[f'HQ{b}:{v}' for b in individuals for v in terms[:3]]+['within_agl_z:tpi_z']
one('individual_terrain_slopes',XX,tt)
for b in individuals:one(f'individual_HQ{b}',X,terms,bats==b)
pd.DataFrame(rows).to_csv(O/'additional_coefficients.csv',index=False)
(O/'additional_diagnostics.json').write_text(json.dumps(diags,indent=2))
# Preserve choice-set duplicates: repeated sampled alternatives retain sampling weights.
c=pd.read_csv(R/'input/full_candidates.csv');v=c.groupby('step_id').apply(lambda g:len(g[['elev','tpi','tri','slope']].drop_duplicates()),include_groups=False)
v.rename('distinct_terrain_vectors').to_csv(O/'candidate_terrain_diversity.csv')
(O/'candidate_audit.json').write_text(json.dumps({'rows':len(c),'steps':c.step_id.nunique(),'every_step_has_40':bool(c.groupby('step_id').size().eq(40).all()),'missing_values':int(c.isna().sum().sum()),'duplicate_terrain_rows_within_step':int(c.duplicated().sum()),'distinct_vectors_median':float(v.median()),'constant_terrain_sets':int(v.eq(1).sum()),'coordinates_provided':False,'candidate_generation_audited':False,'duplicate_policy':'Retained exactly as supplied; without coordinates cannot distinguish identical raster attributes from duplicate sampled locations.'},indent=2))
