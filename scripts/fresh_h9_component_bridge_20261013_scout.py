#!/usr/bin/env python3
"""Fresh whole-chamber trivalent endpoint profiles; numerical proposals only."""
import sys,json,itertools,time
from pathlib import Path
from fractions import Fraction
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import networkx as nx
from scripts import fresh_h9_prism_completion_20261013_scout as S
V=S.V
PREFIX='data/enumerations/fresh_h9_component_bridge_20261013'

def setup():
 d=S.setup();G=np.array(d['gram']);old=np.array(d['normals']);perm=[0,2,1,3,4,6,5,7,8,9]
 R=np.vstack([old,old[90:,perm]]);names=d['names']+['mirror'+n[0].upper()+n[1:]for n in d['names'][90:]]
 src=S.PREFIX+'_caps.json';r=V.read(src);y=np.array([float(Fraction(x['a']))+float(Fraction(x['b']))*((1+5**.5)/2)for x in r['escape']])
 orb=S.PREFIX+'_orbit_obstruction_point.json';ob=V.read(orb)
 print('ORBITKEYS',list(ob),flush=True)
 return d,G,R,names,y

def run():
 d,G,R,names,y=setup();PP=R@G@R.T
 graph=nx.Graph();graph.add_nodes_from(range(len(R)));graph.add_edges_from((i,j)for i in range(len(R))for j in range(i)if PP[i,j]>-2+1e-7)
 cliques=list(nx.find_cliques(graph));cand={tuple(sorted(c))for C in cliques if len(C)>=8 for c in itertools.combinations(C,8)}
 faces=[]
 for c in sorted(cand):
  M=PP[np.ix_(c,c)]
  if np.linalg.eigvalsh(M)[0]>1e-7:faces.append(c)
 out={'source':S.PREFIX+'_caps.json','source_sha256':V.sha(S.PREFIX+'_caps.json'),'normals':R.tolist(),'names':names,'gram':G.tolist(),'h':d['h'],'center':d['center'],'target':y.tolist(),
      'graph_maximal_cliques':len(cliques),'max_clique_size':max(map(len,cliques)),'clique8_candidates':len(cand),'spherical8_faces':faces}
 (V.ROOT/(PREFIX+'_faces_scout.json')).write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items()if k in('graph_maximal_cliques','max_clique_size','clique8_candidates')},'PD8',len(faces),flush=True)
if __name__=='__main__':run()
