import sys,json,numpy as np
sys.path[:0]=['scripts','src']
from fresh_h20_gap_next_20261013 import setup,read,OUT,ROOT
from fresh_h20_edge_family_completion_20261013_bank import gate
from scipy.optimize import linprog
import networkx as nx
Q,R,h,w=setup();d=read(OUT);C=np.array([x['normal']for x in d['caps']]);M=C@Q@C.T;A=gate(-M/2);np.fill_diagonal(A,False)
g=nx.from_numpy_array(A);print('edges',g.number_of_edges(),'nodes',len(C))
cl=[]
for i in sorted(g,key=lambda i:-g.degree(i)):
 if all(g.has_edge(i,j)for j in cl):cl.append(i)
print('greedy',len(cl),flush=True)
U=np.vstack([R,C]);norms=np.einsum('ij,jk,ik->i',U,Q,U);walls=U@Q/np.sqrt(norms[:,None]);hn=h/np.sqrt(-h@Q@h);rng=np.random.default_rng(202720)
for k in range(160):
 p=linprog(rng.normal(size=21),A_ub=walls,b_ub=np.zeros(len(walls)),A_eq=(hn@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
 if p.success and p.x@Q@p.x>1e-7:
  print('allbankescape',k,p.x@Q@p.x,'active',np.flatnonzero(abs(walls@p.x)<1e-7).tolist())
  json.dump({'status':'NUMERICAL_WHOLE_FREE_FAMILY_SAMPLE_BANK_ESCAPE','vector':p.x.tolist(),'norm':float(p.x@Q@p.x),'active_indices':np.flatnonzero(abs(walls@p.x)<1e-7).tolist(),'source':OUT},open('data/enumerations/fresh_h20_gap_next_20261013_sample_bank.json','w'),indent=2);break
else:print('NOESCAPE IN160')
print('minimumcenter', sorted(enumerate(d['caps']),key=lambda z:abs(z[1]['center_pairing']))[:1])
