"""Independent, explicit Fig 3 specification; no attempt to target legacy coefficients."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import t
R=Path(__file__).resolve().parents[1];O=R/'results';O.mkdir(exist_ok=True)
def distance(la,lo,lb,lp):
 a,b,c,d=np.deg2rad([la,lo,lb,lp]);v=np.sin((a-c)/2)**2+np.cos(a)*np.cos(c)*np.sin((b-d)/2)**2
 return 2*6371000*np.arcsin(np.sqrt(np.clip(v,0,1)))
def ols(X,y,groups):
 X=np.asarray(X,float);y=np.asarray(y,float);n,k=X.shape
 if np.linalg.matrix_rank(X)!=k:raise ValueError('Rank-deficient design')
 b=np.linalg.lstsq(X,y,rcond=None)[0];e=y-X@b;bread=np.linalg.inv(X.T@X);cm=bread*(e@e)/(n-k);ids=np.unique(groups);G=len(ids)
 if G>1:
  scores=np.stack([X[groups==g].T@e[groups==g] for g in ids]);cc=bread@scores.T@scores@bread*(G/(G-1))*((n-1)/(n-k))
 else:cc=np.full_like(cm,np.nan)
 return {'beta':b,'cov_model':cm,'cov_cluster':cc,'n':n,'k':k,'G':G,'r2':1-e@e/np.sum((y-y.mean())**2),'residuals':e,'fitted':X@b}
def coef(r,j):
 b=r['beta'][j];se=np.sqrt(r['cov_model'][j,j]);q=t.ppf(.975,r['n']-r['k']);sc=np.sqrt(r['cov_cluster'][j,j]) if r['G']>1 else np.nan;qc=t.ppf(.975,r['G']-1) if r['G']>1 else np.nan
 return dict(beta=b,se=se,ci_lo=b-q*se,ci_hi=b+q*se,p=2*t.sf(abs(b/se),r['n']-r['k']),cluster_se=sc,cluster_ci_lo=b-qc*sc,cluster_ci_hi=b+qc*sc,cluster_p=2*t.sf(abs(b/sc),r['G']-1) if r['G']>1 else np.nan,n=r['n'],n_bats=r['G'])
def run():
 m=pd.read_csv(R/'input/master_steps.csv');m['time']=pd.to_datetime(m.t)
 raw=pd.read_excel(R/'input/raw_tracking.xlsx').rename(columns={'individual':'ind','时间':'time'});raw.time=pd.to_datetime(raw.time)
 rows=[]
 for _,g in raw.sort_values('time').groupby(['ind','seg']):
  g=g.reset_index(drop=True);la=g['纬度'].to_numpy();lo=g['经度'].to_numpy();ts=g.time.astype('int64').to_numpy()/1e9
  step=np.r_[0,distance(la[1:],lo[1:],la[:-1],lo[:-1])];cum=np.cumsum(step)
  for i in range(len(g)):
   a=np.searchsorted(ts,ts[i]-300);b=np.searchsorted(ts,ts[i]+300,side='right')-1;length=cum[b]-cum[a];span=(ts[b]-ts[a])/60;gap=np.max(np.diff(ts[a:b+1]))/60 if b>a else np.nan
   straight=distance(la[a],lo[a],la[b],lo[b])/length if length>0 else np.nan
   rows.append(dict(ind=g.ind.iloc[i],time=g.time.iloc[i],straight_recomputed=straight,window_span_min=span,window_max_gap_min=gap,window_n_steps=b-a,window_path_m=length))
 d=m.merge(pd.DataFrame(rows),on=['ind','time'],validate='one_to_one');d['window_valid']=(d.window_n_steps>=3)&(d.window_span_min>=5)&(d.window_max_gap_min<=2.5)&(d.window_path_m>0)&d.straight_recomputed.between(0,1+1e-10)
 d['straight_recomputed']=d.straight_recomputed.clip(upper=1) # Only floating-point excess <=1e-10 can enter.
 d['analysis_valid']=d.window_valid&d.dt_min.between(0,2.5,inclusive='right')&d.end_tri.notna()&d['运动量'].gt(0)
 d.to_csv(O/'all_steps_with_window_audit.csv',index=False);base=d[d.analysis_valid].copy()
 base['log_activity']=np.log(base['运动量']);base['log_dt']=np.log(base.dt_min);base['log_span']=np.log(base.window_span_min)
 scales={}
 for col in ['end_tri','straight_recomputed','log_dt','log_span']:
  mu=base[col].mean();sd=base[col].std(ddof=1);scales[col]={'mean':float(mu),'sd':float(sd)};base[col+'_z']=(base[col]-mu)/sd
 bats=sorted(base.ind.unique());specs={'straightness':['end_tri_z','log_dt_z','log_span_z'],'activity':['end_tri_z','straight_recomputed_z','log_dt_z','log_span_z']};ys={'straightness':'straight_recomputed','activity':'log_activity'};allrows=[];results={}
 def fitmodel(name,dat,outcome,cols,fe=True):
  bs=sorted(dat.ind.unique());terms=['Intercept']+cols+([f'HQ{b}' for b in bs[1:]] if fe else [])
  X=np.column_stack([np.ones(len(dat))]+[dat[c].to_numpy() for c in cols]+([dat.ind.eq(b).to_numpy(float) for b in bs[1:]] if fe else []));r=ols(X,dat[outcome].to_numpy(),dat.ind.to_numpy());results[name]=r
  for j,v in enumerate(terms):allrows.append({'model':name,'term':v,**coef(r,j)})
  np.savez_compressed(O/(name+'.npz'),beta=r['beta'],cov_model=r['cov_model'],cov_cluster=r['cov_cluster'],terms=terms,step_id=dat.step_id.to_numpy(),fitted=r['fitted'],residuals=r['residuals'])
  return r,terms
 for name,cols in specs.items():
  r,terms=fitmodel(name,base,ys[name],cols);print(name,coef(r,1),flush=True)
  for b in bats:fitmodel(f'{name}_individual_{b}',base[base.ind==b],ys[name],cols,False)
  for b in bats:fitmodel(f'{name}_omit_{b}',base[base.ind!=b],ys[name],cols)
  # Same complete cohort, time adjustment is a stated secondary specification.
  hour=base.time.dt.hour+base.time.dt.minute/60+base.time.dt.second/3600;base['clock_sin']=np.sin(2*np.pi*hour/24);base['clock_cos']=np.cos(2*np.pi*hour/24)
  fitmodel(name+'_clock',base,ys[name],cols+['clock_sin','clock_cos'])
  fitmodel(name+'_full_windows',base[base.window_span_min>=8],ys[name],cols)
  st=base.copy();st['end_tri_z']=(st.start_tri-scales['end_tri']['mean'])/scales['end_tri']['sd'];fitmodel(name+'_start_TRI',st,ys[name],cols)
 # Original straight5 is tested without silently clipping its 42 impossible values.
 old=base[base.straight5.between(0,1)].copy();fitmodel('straightness_supplied_metric',old,'straight5',specs['straightness'])
 old['straight_recomputed_z']=(old.straight5-scales['straight_recomputed']['mean'])/scales['straight_recomputed']['sd'];fitmodel('activity_supplied_metric',old,'log_activity',specs['activity'])
 fitmodel('activity_without_straightness',base,'log_activity',['end_tri_z','log_dt_z','log_span_z'])
 ledger=pd.DataFrame(allrows);ledger.to_csv(O/'all_coefficients.csv',index=False);ledger[(ledger.term=='end_tri_z')].to_csv(O/'TRI_effects.csv',index=False)
 # Prediction contrasts: average individual intercept with equal weight; nuisance terms at cohort means.
 grid=np.linspace(base.end_tri.quantile(.025),base.end_tri.quantile(.975),120)
 for name,cols in specs.items():
  r=results[name];L=np.zeros((len(grid),len(r['beta'])));L[:,0]=1;L[:,1]=(grid-scales['end_tri']['mean'])/scales['end_tri']['sd'];L[:,1+len(cols):]=1/len(bats)
  pred=L@r['beta'];se=np.sqrt(np.einsum('np,pq,nq->n',L,r['cov_model'],L));sc=np.sqrt(np.maximum(0,np.einsum('np,pq,nq->n',L,r['cov_cluster'],L)));qm=t.ppf(.975,r['n']-r['k']);qc=t.ppf(.975,r['G']-1)
  z=pd.DataFrame({'TRI':grid,'fit_link':pred,'model_lo_link':pred-qm*se,'model_hi_link':pred+qm*se,'cluster_lo_link':pred-qc*sc,'cluster_hi_link':pred+qc*sc})
  for col in list(z.columns[1:]):z[col.replace('_link','')]=np.exp(z[col]) if name=='activity' else z[col]
  z.to_csv(O/(name+'_predictions.csv'),index=False)
 base.to_csv(O/'analysis_rows.csv',index=False)
 audit={'raw_fixes':len(raw),'master_steps':len(m),'analysis_steps':len(base),'bats':len(bats),'excluded_steps':len(d)-len(base),'supplied_straight_missing':int(m.straight5.isna().sum()),'supplied_straight_above_one':int((m.straight5>1).sum()),'supplied_straight_max':float(m.straight5.max()),'recomputed_vs_supplied_correlation_on_analysis_rows':float(base.straight_recomputed.corr(base.straight5)),'n_by_bat':{str(k):int(v) for k,v in base.groupby('ind').size().items()},'TRI_plot_range':list(map(float,[grid.min(),grid.max()])),'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (R/'input').iterdir()}}
 (O/'audit.json').write_text(json.dumps(audit,indent=2));(O/'scaling.json').write_text(json.dumps(scales,indent=2));(O/'diagnostics.json').write_text(json.dumps({k:{v:r[v] for v in ['n','k','G','r2']} for k,r in results.items()},indent=2));print(audit,flush=True)
if __name__=='__main__':run()
