#!/usr/bin/env python3
"""Replayable finite-bank selection and unchanged-face diagnostics, not a new gate."""
import argparse,json
from hashlib import sha256
from pathlib import Path
import networkx as nx
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_nonright_s6_20261016_strategy'
SEL='results/fresh_nonarithmetic_nonright_s6_20261016_selection.json'
C53='results/fresh_nonarithmetic_escape52_20261016_cone.json'
C55='results/fresh_nonarithmetic_nonright_s6_20261016_cone.json'
S={0,1,2,3,5,6}
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for x,h in d.get('sources',{}).items():assert sha(x)==h,(p,x)
 return d
def run():
 s=read(SEL);a=read(C53);b=read(C55);r=s['records']
 fa={k:[v for v in rows if S<=set(v)]for k,rows in a['complete_census'].items()}
 fb={k:[v for v in rows if S<=set(v)]for k,rows in b['complete_census'].items()}
 assert fa==fb
 hi=[i for i,v in enumerate(a['complete_census']['hyperideal_rays'])if S<=set(v)]
 touches=[x['index']for x in r if set(hi)&set(x['numerically_cut_ray_indices'])]
 assert not touches
 # Only numerically strict dotted edges enter this graph. The two near-zero
 # tentative label2 fits from the raw selector are deliberately not accepted.
 keep=[i for i in s['no_affine_candidate_indices']if all(m in[2,3,4,5,6,10]for _,m in r[i]['finite_contacts'])]
 G=nx.Graph();G.add_nodes_from(keep)
 G.add_edges_from(x['indices']for x in s['numerical_compatibility']if x['kind']=='dotted'and all(i in G for i in x['indices']))
 best=[];number=0;max_size=0;max_cover=0
 for c in nx.find_cliques(G):
  number+=1;max_size=max(max_size,len(c));u=set().union(*(set(r[i]['numerically_cut_ray_indices'])for i in c));z=len(u)
  if z>max_cover:max_cover=z;best=[]
  if z==max_cover:best.append(sorted(c))
 selected=[844,877];cover=set().union(*(set(r[i]['numerically_cut_ray_indices'])for i in selected))
 thirds=sorted(set(G.neighbors(selected[0]))&set(G.neighbors(selected[1])))
 assert all(set(r[i]['numerically_cut_ray_indices'])<=cover for i in thirds)
 paths=[f'scripts/{TAG}.py',SEL,C53,C55,'results/fresh_nonarithmetic_h3_escape_projection_20261016.json','results/fresh_nonarithmetic_nonright_s6_20261016_caps.json']
 return dict(status='PASS_BOUNDED_BANK_SELECTION_AND_EXACT_UNCHANGED_H3_FACE_DIAGNOSTIC',sources={p:sha(p)for p in paths},white_support=sorted(S),unchanged_actual_face_census=fa,unchanged_actual_face_counts={k:len(v)for k,v in fa.items()},all_numerical_bank_candidates=len(r),bank_candidates_numerically_cutting_any_H3_exterior_ray=touches,finite_label_subset=[2,3,4,5,6,10],graph_candidates=len(keep),numerical_dotted_edges=G.number_of_edges(),maximal_cliques=number,maximum_clique_size=max_size,maximum_old_exterior_ray_coverage=max_cover,maximum_coverage_cliques=sorted(best),selected_candidates=selected,selected_old_ray_coverage=len(cover),compatible_third_candidates=thirds,third_candidates_add_no_old_ray_coverage=True,
 scope='The unchanged actual H3 face follows exactly from equality of the complete53/55 cone active-set censuses and the already certified positive-definite S6. Clique and ray-coverage statistics concern only the saved numerical bank and its tentative contact graph; they are not all-label or all-real obstructions. Two apparent near-zero mutual pair fits in the raw selector are excluded from this dotted-only graph. Only the two selected roots have been exactified here.',
 next_mechanism='Use the residual H3 ideal-boundary region, not old-ray coverage alone, as the closure objective. It still has13 finite vertices and7 exterior directions. A coordinated right-S6 side cluster may close that region while separate nonright-S6 cutters handle the ambient escape. Any proposed compact regional closure must pass a complete four-dimensional cone calculation; no finite sample of ideal directions proves closure.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved selection/face diagnostic')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['unchanged_actual_face_counts'],d['maximum_old_exterior_ray_coverage'])
