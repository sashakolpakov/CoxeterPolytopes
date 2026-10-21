#!/usr/bin/env python3
"""Numerical72-symmetry orbit of the exact new free-dotted-anchor cap."""
import json,numpy as np,networkx as nx
from scipy.optimize import linprog
from fresh_h20_coupled_completion_20261013_dotted_exact import setup,OUT as EXACT,SOURCE
from fresh_h20_coupled_completion_20261013 import num,read,ROOT
from fresh_dimension_bridge_20261013_verify import vec,parse
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension
A_SOURCE='data/enumerations/fresh_h20_coupled_completion_20261013_automorphisms.json'
OUT='data/enumerations/fresh_h20_coupled_completion_20261013_orbit_scout.json'
Qe,he,hc,roots=setup();Q=np.array([[num(B,x)for x in row]for row in Qe]);h=np.array([num(B,x)for x in he]);R=np.array([[num(f,x)for x in r]for f,r in roots]);base=read('data/enumerations/fresh_dimension_bridge_20261013.json');W=np.array([[num(B,x)for x in vec(row)]for row in base['old_normals']]);eps=dict(base['height_signs']);d=read(EXACT);cap=d['new_cap'];f=Extension(B,parse(cap['field']['radicand']));r=np.array([num(f,x)for x in map(f.parse,cap['normal'])]);R=np.vstack([R,r]);ns=np.einsum('ij,jk,ik->i',R,Q,R)
angles=np.array([0]+[np.cos(np.pi/m)for m in(3,4,5,6,8,10,12)])
def gate(v):return(v>=1-1e-7)|(np.min(abs(v[...,None]-angles),axis=-1)<1e-7)
perms=read(A_SOURCE)['automorphisms'];seen=set();caps=[];allrows=[]
for k,p in enumerate(perms):
 sg=eps[p[24]]//eps[24];assert all(eps[p[i]]==sg*eps[i]for i in eps)
 M=np.zeros((21,21));M[:20]=W[np.array(p[:20])];M[20,20]=sg
 assert np.max(abs(M@Q@M.T-Q))<1e-8
 z=r@M;key=tuple(np.round(z,8))
 if key in seen:continue
 seen.add(key);rat=-R@Q@z/np.sqrt(2*ns);bad=np.flatnonzero(~gate(rat)).tolist();allrows.append({'permutation_index':k,'normal':z.tolist(),'bad':bad})
 if not bad:caps.append({'permutation_index':k,'normal':z.tolist()})
print('orbit',len(allrows),'accepted',len(caps),'badcounts',sorted((len(x['bad']),x['bad'])for x in allrows)[:20],flush=True)
C=np.array([x['normal']for x in caps]);graph=gate(-C@Q@C.T/2);np.fill_diagonal(graph,False);g=nx.from_numpy_array(graph);cliques=list(nx.find_cliques(g));print('graph',len(g),g.number_of_edges(),'cliques',len(cliques),'max',max(map(len,cliques),default=0),flush=True)
rng=np.random.default_rng(202610139);U=np.vstack([R,C]);n=np.einsum('ij,jk,ik->i',U,Q,U);A=U@Q/np.sqrt(n[:,None]);hn=h/np.sqrt(-h@Q@h);escape=None
for k in range(100):
 x=linprog(rng.normal(size=21),A_ub=A,b_ub=np.zeros(len(A)),A_eq=(hn@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
 if x.success and x.x@Q@x.x>1e-7:escape={'vector':x.x.tolist(),'norm':float(x.x@Q@x.x),'active_indices':np.flatnonzero(abs(A@x.x)<1e-7).tolist(),'attempt':k};break
print('all accepted escape',None if escape is None else(escape['norm'],escape['active_indices']),flush=True)
(ROOT/OUT).write_text(json.dumps({'status':'NUMERICAL_ORBIT_SCOUT','all_orbit':allrows,'accepted':caps,'maximal_cliques':cliques,'all_accepted_escape':escape},indent=2)+'\n')
