import hashlib,json,numpy as np,mpmath as m
from pathlib import Path
from scipy.spatial import HalfspaceIntersection
from scripts.fresh_h9_split_leaf_caps_20261013 import seed,PHI
m.mp.dps=40;p=(1+m.sqrt(5))/2;k=lambda x:m.mpf(x['a'])+m.mpf(x['b'])*p
j=json.load(open('data/enumerations/fresh_h9_split_endpoint_selection_20261013_scout.json'));best=min((r for r in j['history']if'selected'in r),key=lambda r:r['sampled_max_square']);d={}
for path in j['sources']:
 assert hashlib.sha256(Path(path['path']).read_bytes()).hexdigest()==path['sha256']
 d[path['path'].split('/')[-1]]=json.load(open(path['path']))['caps']
R=[];refs=[]
for ii in best['selected']:
 fn,i=j['references'][ii];r=d[fn][i];sd=m.sqrt(max(0,k(r['radicand'])));R.append([float(k(x[0])+k(x[1])*sd)for x in r['normal']]);refs.append([fn,i])
a=PHI+.5;b=-a/(12*a+1);G=seed(a,b);ev,Q=np.linalg.eigh(G);ix=list(np.where(ev< -1e-7)[0])+list(np.where(ev>1e-7)[0]);N=Q[:,ix]*np.sqrt(abs(ev[ix]));N[:,0]*=np.sign(sum(N[:,0]));W=np.vstack([N,np.array(R)@N]);H=HalfspaceIntersection(np.c_[W[:,1:],-W[:,0]],np.zeros(9));q=np.sum(H.intersections**2,axis=1);finite=np.flatnonzero(q<1-1e-7);ideal=np.flatnonzero(abs(q-1)<1e-7);outside=np.flatnonzero(q>1+1e-7)
print({'iteration':best['iteration'],'facets':len(W),'finite':len(finite),'ideal':len(ideal),'outside':len(outside),'max_square':max(q)},flush=True)
for ii in finite:print('finite',q[ii],H.dual_facets[ii],flush=True)
json.dump({'source_iteration':best['iteration'],'refs':refs,'root_coefficients':R,'physical_normals':W.tolist(),'N':N.tolist(),'finite_vertices':[{'point':H.intersections[ii].tolist(),'active':H.dual_facets[ii]}for ii in finite],'leading_escape':{'point':H.intersections[int(np.argmax(q))].tolist(),'active':H.dual_facets[int(np.argmax(q))]}},open('data/enumerations/fresh_h9_split_endpoint_selected_geometry_20261013_scout.json','w'),indent=2)
