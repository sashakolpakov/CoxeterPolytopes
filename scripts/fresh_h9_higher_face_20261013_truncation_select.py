#!/usr/bin/env python3
"""Joint geometry selection of the fresh roof polar/two-plane truncations."""
import json,sys
from itertools import combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,networkx as nx
from scripts import fresh_h9_higher_face_20261013_linked_scout as L
S=L.S

def run():
 _,_,_,_,_,names,_,_,R,G,c=S.load();d=S.V.read(S.PREFIX+'_selected_scout.json');rr=np.array([x['normal']for x in d['recipes']]);R=np.vstack([R,rr]);names=names+[f'roof{i}'for i in range(len(rr))];fixed=[names.index(n)for n in S.COMMON if n!='seed1'];pd=S.V.read(S.PREFIX+'_polar_scout.json');td=S.V.read(S.PREFIX+'_truncation_scout.json');bank=[{'kind':'polar','right':x['incident'],'normal':x['normal'],'target':x['point']}for x in pd['rows']if not x['bad']]+[dict(x,kind='two_plane')for x in td['candidates']];A=np.array([x['normal']for x in bank]);P=A@G@A.T;graph=nx.Graph();graph.add_nodes_from(range(len(bank)));graph.add_edges_from((i,j)for i,j in combinations(range(len(bank)),2)if L.gate(P[i,j]));cases=[]
 for ids in nx.find_cliques(graph):
  geo=L.geom(np.vstack([R,A[ids]]),G,c,fixed);cases.append({'indices':sorted(ids),'geometry':geo})
 cases.sort(key=lambda x:x['geometry'].get('max_norm',1e99));out={'status':'NUMERICAL_JOINT_ROOF_TRUNCATION_GEOMETRY','bank':bank,'compatible_edges':[list(x)for x in graph.edges],'cases':cases};(S.V.ROOT/(S.PREFIX+'_truncation_select.json')).write_text(json.dumps(out,indent=2)+'\n');print('bank',len(bank),'maximal',len(cases));print(cases[:5])
if __name__=='__main__':run()
