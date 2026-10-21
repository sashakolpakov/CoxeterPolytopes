exec(open('scripts/fresh_h20_coupled_completion_20261013_next_polar_scout.py').read().split('out=[];best=[]')[0])
from scipy.optimize import linprog
import networkx as nx
D=json.load(open('data/enumerations/fresh_h20_coupled_completion_20261013_dotted_scout.json'));caps=np.array([r['normal']for r in D['good']]);pp=-caps@Q@caps.T/2;gp=(pp>=1-1e-7)|(np.min(abs(pp[...,None]-allowed[None,None,:]/2),axis=2)<1e-7);np.fill_diagonal(gp,False);g=nx.from_numpy_array(gp);print('graph',len(g),g.number_of_edges(),flush=True)
hn=h/np.sqrt(-h@Q@h)
def escape(ids,seed=202610135):
 U=np.vstack([R,caps[ids]])if ids else R;qn=np.einsum('ij,jk,ik->i',U,Q,U);A=U@Q/np.sqrt(qn[:,None]);rng=np.random.default_rng(seed)
 for k in range(80):
  ans=linprog(rng.normal(size=21),A_ub=A,b_ub=np.zeros(len(A)),A_eq=(hn@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
  if ans.success and ans.x@Q@ans.x>1e-7:return {'vector':ans.x.tolist(),'norm':float(ans.x@Q@ans.x),'active_indices':np.flatnonzero(abs(A@ans.x)<1e-7).tolist(),'attempt':k}
 return None
full=escape(list(range(len(caps))));print('full bank escape',None if full is None else(full['norm'],full['active_indices']),flush=True)
cliques=list(nx.find_cliques(g));cliques.sort(key=lambda x:(-len(x),x));print('cliques',len(cliques),'max',len(cliques[0]),flush=True)
rows=[]
for ids in cliques:
 e=escape(ids);rows.append({'cap_ids':ids,'escape':e});print('clique',ids,'escape',None if e is None else e['norm'],flush=True)
json.dump({'status':'NUMERICAL_ONLY','source':'data/enumerations/fresh_h20_coupled_completion_20261013_dotted_scout.json','full_bank_escape':full,'maximal_compatible_cliques':rows},open('data/enumerations/fresh_h20_coupled_completion_20261013_dotted_volume_scout.json','w'),indent=2)
