#!/usr/bin/env python3
import json,numpy as np,networkx as nx
from scipy.optimize import linprog
from fresh_h20_coupled_completion_20261013_dotted_exact import setup,OUT as FIRST
from fresh_h20_coupled_completion_20261013 import num,read,ROOT
from fresh_dimension_bridge_20261013_verify import vec,parse
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension
Qe,he,hc,roots=setup();c=read(FIRST)['new_cap'];f=Extension(B,parse(c['field']['radicand']));roots.append((f,tuple(map(f.parse,c['normal']))));Q=np.array([[num(B,x)for x in r]for r in Qe]);h=np.array([num(B,x)for x in he]);R=np.array([[num(g,x)for x in r]for g,r in roots]);orb=read('data/enumerations/fresh_h20_coupled_completion_20261013_orbit_scout.json');R=np.vstack([R,[x['normal']for x in orb['accepted']]]);second=read('data/enumerations/fresh_h20_coupled_completion_20261013_second_scout.json');chosen=next(c for c in second['caps']if c['type']=='finite_leaves'and c['labels']==[2,4,2]);r=np.array(chosen['normal']);ns=np.einsum('ij,jk,ik->i',R,Q,R);angles=np.array([0]+[np.cos(np.pi/m)for m in(3,4,5,6,8,10,12)])
def gate(v):return(v>=1-1e-7)|(np.min(abs(v[...,None]-angles),axis=-1)<1e-7)
base=read('data/enumerations/fresh_dimension_bridge_20261013.json');W=np.array([[num(B,x)for x in vec(row)]for row in base['old_normals']]);eps=dict(base['height_signs']);perms=read('data/enumerations/fresh_h20_coupled_completion_20261013_automorphisms.json')['automorphisms'];caps=[];rows=[];seen=set()
for k,p in enumerate(perms):
 sg=eps[p[24]]//eps[24];M=np.zeros((21,21));M[:20]=W[np.array(p[:20])];M[20,20]=sg;z=r@M;key=tuple(np.round(z,8))
 if key in seen:continue
 seen.add(key);rat=-R@Q@z/np.sqrt(2*ns);bad=np.flatnonzero(~gate(rat)).tolist();rows.append({'permutation_index':k,'normal':z.tolist(),'bad':bad})
 if not bad:caps.append({'permutation_index':k,'normal':z.tolist()})
C=np.array([x['normal']for x in caps]);good=gate(-C@Q@C.T/2);np.fill_diagonal(good,False);g=nx.from_numpy_array(good);cliques=list(nx.find_cliques(g));cliques.sort(key=lambda c:-len(c));print('orbit',len(rows),'accepted',len(caps),'graph',len(g),g.number_of_edges(),'maxclique',len(cliques[0]),'cliques',len(cliques),flush=True)
sel=cliques[0];U=np.vstack([R,C[sel]]);n=np.einsum('ij,jk,ik->i',U,Q,U);A=U@Q/np.sqrt(n[:,None]);hn=h/np.sqrt(-h@Q@h);rng=np.random.default_rng(2026101311);esc=None;trials=[]
for k in range(300):
 x=linprog(rng.normal(size=21),A_ub=A,b_ub=np.zeros(len(A)),A_eq=(hn@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
 if not x.success:continue
 norm=float(x.x@Q@x.x);trials.append(norm)
 if norm>1e-7:esc={'vector':x.x.tolist(),'norm':norm,'active_indices':np.flatnonzero(abs(A@x.x)<1e-7).tolist(),'attempt':k};break
print('escape',None if esc is None else(esc['norm'],esc['active_indices']),flush=True)
(ROOT/'data/enumerations/fresh_h20_coupled_completion_20261013_second_orbit_scout.json').write_text(json.dumps({'status':'NUMERICAL_SECOND_ORBIT','chosen':chosen,'all_orbit':rows,'accepted':caps,'maximal_cliques':cliques,'selected':sel,'escape':esc,'trials':trials},indent=2)+'\n')
