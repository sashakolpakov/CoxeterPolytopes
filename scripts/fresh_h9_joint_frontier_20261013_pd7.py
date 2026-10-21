#!/usr/bin/env python3
"""Exact PD7 face-support and existing-edge incidence handoff."""
from pathlib import Path
from itertools import combinations
from collections import Counter,defaultdict
import argparse,json
import fresh_h9_joint_frontier_20261013 as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_joint_frontier_20261013_pd7'
EDGES=ROOT/'results/fresh_h9_joint_frontier_20261013_all_edges.json';OUT=ROOT/'results'/f'{TAG}.json'
def run():
 src=json.loads(EDGES.read_text())
 for p,h in src['sources'].items():assert A.sha(ROOT/p)==h
 s,Q,W,C,c,G,adj=A.setup();cl=A.clique_list(adj,7)
 pd=[ids for ids in cl if A.pd([[G[i][j]for j in ids]for i in ids])]
 incidence=defaultdict(list)
 for edge in src['all_edge_supports']:
  for face in combinations(edge['support'],7):incidence[face].append({'edge_id':edge.get('finite_frontier_edge_id',edge.get('new_edge_id')),'kind':edge['kind'],'support':edge['support']})
 assert set(pd)==set(incidence) and len(cl)==2234 and len(pd)==1339
 rows=[{'id':k,'support':list(ids),'existing_edges':incidence[ids]}for k,ids in enumerate(pd)]
 counts=Counter(len(r['existing_edges'])for r in rows)
 return {'status':'PASS_EXACT_PD7_SUPPORT_AND_EXISTING_EDGE_INCIDENCE','sources':{str(p.relative_to(ROOT)):A.sha(p)for p in[EDGES,A.SOURCE,Path(A.__file__),Path(__file__)]},'finite_pair_seven_cliques':len(cl),'PD7_faces':len(pd),'edge_count_distribution':{str(k):v for k,v in sorted(counts.items())},'faces':rows,'scope':'Complete positive-definite rank-seven supports in the fixed 27-wall cone and their already certified geodesic-edge incidences. These are actual hyperbolic two-faces by strict-center projection; no new cap or closure is asserted.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run()
 if a.verify:assert json.loads(OUT.read_text())==out
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['PD7_faces'],out['edge_count_distribution'],flush=True)
