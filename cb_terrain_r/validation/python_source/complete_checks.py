"""Post hoc completion checks. Run from any directory; supplied analytical terrain is retained.
No candidate coordinates or environmental values are simulated.
"""
from pathlib import Path
import json,itertools
import numpy as np,pandas as pd
from scipy.stats import t
from model_utils import design,ols,stats
from conditional_logit import fit,coefficient_rows
R=Path(__file__).resolve().parents[1]; I=R/'input'; O=R/'results';O.mkdir(exist_ok=True)
d=pd.read_csv(I/'activity_rows.csv.gz');m=pd.read_csv(I/'steps.csv').sort_values('step_id').reset_index(drop=True)
for c in ['time','window_start','window_end']:d[c]=pd.to_datetime(d[c])
d['log_rate']=np.log(d['运动量']/d.raw_dt_min)
rows=[];individual=[];shapes=[];bins=[]
def nonoverlap(b):
 ix=[]
 for _,a in b.sort_values('time').groupby(['ind','seg']):
  end=pd.Timestamp.min
  for j,r in a.iterrows():
   if r.window_start>end:ix.append(j);end=r.window_end
 return b.loc[ix]
def one(kind,name,b,cols=None,equal=False,segment=False,wild=False):
 y='straight_recomputed' if kind=='straightness' else 'log_rate'
 cs=cols or ['end_tri_z']+(['straight_recomputed_z'] if kind=='activity' else [])+['log_dt_z','log_span_z']
 X,names=design(b,cs,group='seg' if segment else 'ind'); yy=b[y].to_numpy();g=b.ind.to_numpy()
 r=ols(X,yy,g,1/b.groupby('ind').ind.transform('size').to_numpy() if equal else None)
 rr=dict(outcome=kind,model=name,**stats(r))
 if wild:
  X0=np.delete(X,1,axis=1);r0=ols(X0,yy,g);ids=np.unique(g);ts=[];obs=abs(r['beta'][1]/np.sqrt(r['cov_cluster'][1,1]))
  for ss in itertools.product([-1.,1.],repeat=len(ids)):
   rb=ols(X,r0['fitted']+r0['residuals']*np.asarray(ss)[np.searchsorted(ids,g)],g)
   ts.append(abs(rb['beta'][1]/np.sqrt(rb['cov_cluster'][1,1])))
  rr['wild_p']=float(np.mean(np.asarray(ts)>=obs-1e-12))
 rows.append(rr)
 if name=='primary':
  np.savez_compressed(O/(kind+'_primary.npz'),**r,terms=names)
  for ind,a in b.groupby('ind'):
   xx,_=design(a,cs);individual.append(dict(outcome=kind,ind=ind,n=len(a),beta=np.linalg.lstsq(xx,a[y],rcond=None)[0][1]))
 return r,names
for kind in ['straightness','activity']:
 one(kind,'primary',d,wild=True);one(kind,'equal_bat_weight',d,equal=True);one(kind,'segment_intercepts',d,segment=True)
 one(kind,'nonoverlapping',nonoverlap(d));one(kind,'raw_duration_45_75_seconds',d[d.raw_dt_min.between(.75,1.25)])
 cs=['end_tri_z']+(['straight_recomputed_z'] if kind=='activity' else [])+['log_dt_z','log_span_z']
 one(kind,'clock_time',d,cs+['clock_sin','clock_cos'])
 d['log_hdop']=np.log(d.HDOP);one(kind,'GPS_quality_covariates',d,cs+['log_hdop','卫星'])
 one(kind,'strict_window_quality',d[(d.window_max_hdop<=5)&(d.window_min_sat>=5)])
 for cut in [6,10,20,50]:one(kind,'count_ge_'+str(cut),d[d['运动量']>=cut])
 for cut in [5,10,20]:one(kind,'rate_ge_'+str(cut),d[d['运动量']/d.raw_dt_min>=cut])
 for ind in sorted(d.ind.unique()):one(kind,'without_HQ'+str(ind),d[d.ind!=ind])
 # Fixed, exploratory natural cubic curves; covariates held at full-sample means.
 knots=np.array(json.loads((I/'shape_knots.json').read_text())['knots_z']);sc=json.loads((I/'movement_scaling.json').read_text())
 def basis(x):
  def dd(j):return (np.maximum(x-knots[j],0)**3-np.maximum(x-knots[-1],0)**3)/(knots[-1]-knots[j])
  return dd(0)-dd(2),dd(1)-dd(2)
 d['spline1'],d['spline2']=basis(d.end_tri_z.to_numpy());y='straight_recomputed' if kind=='straightness' else 'log_rate'
 cs=['end_tri_z','spline1','spline2']+(['straight_recomputed_z'] if kind=='activity' else [])+['log_dt_z','log_span_z']
 grid=np.linspace(d.end_tri.quantile(.025),d.end_tri.quantile(.975),180);z=(grid-sc['end_tri']['mean'])/sc['end_tri']['sd'];s1,s2=basis(z)
 for name,b in [('primary',d),('equal_bat_weight',d)]+[('without_HQ'+str(ind),d[d.ind!=ind]) for ind in sorted(d.ind.unique())]:
  X,names=design(b,cs);r=ols(X,b[y].to_numpy(),b.ind.to_numpy(),1/b.groupby('ind').ind.transform('size').to_numpy() if name=='equal_bat_weight' else None)
  L=np.zeros((len(grid),X.shape[1]));L[:,0]=1;L[:,1]=z;L[:,2]=s1;L[:,3]=s2;L[:,len(cs)+1:]=1/b.ind.nunique();p=L@r['beta'];se=np.sqrt(np.maximum(0,np.einsum('ij,jk,ik->i',L,r['cov_cluster'],L)));q=t.ppf(.975,b.ind.nunique()-1)
  shapes.append(pd.DataFrame(dict(outcome=kind,model=name,TRI=grid,fit=p,lo=p-q*se,hi=p+q*se,response_mean=d[y].mean(),response_sd=d[y].std(ddof=1))))
  if name=='primary':
   curve=r['beta'][0]+r['beta'][len(cs)+1:].sum()/b.ind.nunique()+X[:,1:4]@r['beta'][1:4]
   bb=pd.DataFrame(dict(TRI=b.end_tri,adjusted=b[y].to_numpy()-r['fitted']+curve));bb['bin']=pd.qcut(bb.TRI,12,labels=False,duplicates='drop');bb=bb.groupby('bin').agg(TRI=('TRI','mean'),adjusted=('adjusted','mean'),n=('TRI','size')).reset_index();bb['outcome']=kind;bins.append(bb)
# Additive and interaction models for activity, with corrected duration.
rr=[]
for name,cs in [('additive',['end_tri_z','straight_recomputed_z','log_dt_z','log_span_z']),('interaction',['end_tri_z','straight_recomputed_z','log_dt_z','log_span_z','TRI_x_straightness'])]:
 X,names=design(d,cs);r=ols(X,d.log_rate.to_numpy(),d.ind.to_numpy())
 for term in cs:rr.append(dict(model=name,term=term,**stats(r,names.index(term))))
pd.DataFrame(rr).to_csv(O/'activity_interactions.csv',index=False)
# Whole-window count rates include every original record, irrespective of GPS validity/count screen.
raw=pd.read_csv(I/'raw_records.csv.gz');raw.time=pd.to_datetime(raw.time);m['time']=pd.to_datetime(m.t)
ep=m[['ind','time','end_tri']].rename(columns={'end_tri':'terrain'});sp=m[['ind','time','dt_min','start_tri']].copy();sp.time=(sp.time-pd.to_timedelta(sp.dt_min,unit='min')).dt.round('s');sp=sp.rename(columns={'start_tri':'terrain'})[['ind','time','terrain']]
terrain=pd.concat([ep,sp]).drop_duplicates(['ind','time']);sync=[]
for ind,b in d.groupby('ind'):
 a=raw[raw.ind==ind];a=a.sort_values('time');ter=terrain[terrain.ind==ind]
 for ix,r in b.iterrows():
  z=a[(a.time>r.window_start)&(a.time<=r.window_end)];tt=ter[(ter.time>=r.window_start)&(ter.time<=r.window_end)]
  duration=z.raw_dt_min.sum();coverage=len(tt)/(r.window_n_steps_x+1)
  valid=len(z)>0 and z.raw_dt_min.gt(0).all() and z.raw_dt_min.le(2.5).all() and abs(duration-r.window_span_min)<1e-6 and coverage>=.8
  sync.append(dict(row=ix,valid=valid,duration=duration,count=z['运动量'].sum(),mean_TRI=tt.terrain.mean(),coverage=coverage))
sy=pd.DataFrame(sync);sy.to_csv(O/'synchronized_window_audit.csv',index=False);ok=sy[sy.valid].set_index('row');b=d.loc[ok.index].copy();b['log_rate']=np.log(ok['count']/ok.duration);b['end_tri_z']=(ok.mean_TRI-sc['end_tri']['mean'])/sc['end_tri']['sd'];b.to_csv(O/'synchronized_analysis.csv.gz',index=False)
one('activity','whole_window',b);one('activity','whole_window_nonoverlapping',nonoverlap(b))
pd.DataFrame(rows).to_csv(O/'movement_checks.csv',index=False);pd.DataFrame(individual).to_csv(O/'individual_slopes.csv',index=False);pd.concat(shapes).to_csv(O/'shape_predictions.csv',index=False);pd.concat(bins).to_csv(O/'shape_bins.csv',index=False)
# Choice-set composition interaction, in addition to original height interaction.
c=pd.read_csv(I/'choices.csv.gz').sort_values(['step_id','choice_index']);sc=json.loads((I/'selection_scaling.json').read_text())['predictors']
for k in ['tpi','elev','tri','start_agl']:c[k+'_z']=(c[k]-sc[k]['mean'])/sc[k]['sd']
X=c[['tpi_z','elev_z','tri_z']].to_numpy().reshape(len(m),41,3);X=np.concatenate([X,(c.start_agl_z.to_numpy().reshape(len(m),41)*X[:,:,0])[:,:,None]],axis=2)
terms=['tpi_z','elev_z','tri_z','agl_z:tpi_z'];sel=[]
support=c[c.used==0].groupby('step_id').agg(tpi_sd=('tpi','std'),elev_min=('elev','min'),elev_max=('elev','max'),tpi_min=('tpi','min'),tpi_max=('tpi','max'));support['elev_range']=support.elev_max-support.elev_min;support.to_csv(O/'candidate_support.csv')
for name,xx,tt in [('primary',X,terms),('candidate_terrain_variation',np.concatenate([X]+[(X[:,:,0]*((support[k]-support[k].mean())/support[k].std(ddof=1)).to_numpy()[:,None])[:,:,None] for k in ['tpi_sd','elev_range']],axis=2),terms+['tpi_z:available_tpi_sd','tpi_z:available_elev_range'])]:
 r=fit(xx,m.ind.to_numpy());sel.extend(coefficient_rows(r,tt,name));np.savez_compressed(O/('selection_'+name+'.npz'),beta=r['beta'],cov_cluster=r['cov_cluster'],cov_model=r['cov_model'],terms=tt)
pd.DataFrame(sel).to_csv(O/'selection_new_checks.csv',index=False)
print(pd.DataFrame(rows).query("model in ['primary','equal_bat_weight','count_ge_50','whole_window','whole_window_nonoverlapping']").to_string(index=False),flush=True)
print(pd.DataFrame(sel).query("term == 'agl_z:tpi_z'").to_string(index=False),flush=True)
