#!/usr/bin/env python3
"""Finite-incidence capacity for a three-cap release52 repair, no ideals."""
import argparse,json
from hashlib import sha256
from itertools import product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_b_release_incidence_20261017'
BOUNDARY=ROOT/'results/fresh_nonarithmetic_replace47_close12_boundary_20261017.json'
COSTS=ROOT/'results/fresh_nonarithmetic_b_release_costs_20261017.json'

def run():
    q=json.loads(BOUNDARY.read_text())['after64'];cost=json.loads(COSTS.read_text())
    finite=[r for r in q['rays']if r['kind']=='finite_vertex']
    # Retain only the three explicit geodesic classes, never exterior-only
    # edges of the auxiliary projective cone.
    kinds={'compact_geodesic_edge','one_open_end','two_open_ends'}
    geodesic=[e for e in q['edges']if e['kind']in kinds]
    removed_edges=[e['sides']for e in geodesic if 52 in e['sides']]
    removed_vertices=[r['active']for r in finite if 52 in r['active']]
    assert len(finite)==25 and len(geodesic)==45
    assert len(removed_edges)==6 and len(removed_vertices)==5
    e0=len(geodesic)-len(removed_edges);v0=len(finite)-len(removed_vertices)
    assert(e0,v0,2*e0-3*v0)==(39,20,18)
    records=[]
    for a,b,x,y,z in product(range(7),range(4),range(4),range(7),range(2)):
        if x>a//2 or x+y>a or y+3*z>2*b or(z and b<3):continue
        records.append(dict(old_new_edges=a,new_new_edges=b,
            vertices_one_new=x,vertices_two_new=y,vertices_three_new=z,
            change_nonfinite_end_incidences=2*(a+b)-3*(x+y+z)))
    minimum=min(r['change_nonfinite_end_incidences']for r in records)
    sharp=[r for r in records if r['change_nonfinite_end_incidences']==minimum]
    assert minimum==-3 and len(sharp)==1
    assert tuple(sharp[0][k]for k in['old_new_edges','new_new_edges','vertices_one_new','vertices_two_new','vertices_three_new'])==(6,3,3,3,1)
    return dict(status='PASS_SCOPED_THREE_CAP_FINITE_INCIDENCE_REPAIR_GATE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),BOUNDARY,COSTS]},
        frozen_control=dict(walls=64,finite_vertices_ambient=198,S6_finite=25,S6_geodesic_edges=45,S6_open_ends=15),
        deleted_wall=52,removed_geodesic_edges=removed_edges,removed_S6_finite_vertices=removed_vertices,
        after_deletion=dict(S6_geodesic_edges=e0,S6_finite_vertices=v0,nonfinite_end_incidences=18),
        hypotheses=['All old walls other than52 retained; at most3 new actual S6 facets.',
            'All walls are mutually acute and have a common strict future timelike interior point.',
            'Each new S6 facet meets at most2 retained S6 facets along actual geodesic edges.',
            'No ideal vertices in the resulting S6 chamber when interpreting the bound as an open-end bound.'],
        variables={'a':'old-new geodesic edges','b':'new-new geodesic edges','x':'new finite vertices with1 new facet','y':'new finite vertices with2 new facets','z':'new finite vertices with3 new facets'},
        necessary_constraints=['0<=a<=6','0<=b<=3','0<=x<=floor(a/2)','x+y<=a','y+3z<=2b','0<=z<=1','z>0 implies b=3'],
        relaxed_integer_cases=records,minimum_change=minimum,unique_minimum_case=sharp[0],minimum_final_nonfinite_end_incidences=15,
        retention_proof='At a common strict interior X, projection onto a PD old support has coefficient vector G^{-1}<X,n><=0 because G is Stieltjes. Its scalar product with any acute new normal r is <X,r> minus a nonnegative sum, hence strictly negative. This preserves every retained old finite vertex and produces an actual point on every retained old PD two-side edge. No old vertex or geodesic edge is silently removed by the new caps. Old-only supports after deleting52 introduce no new finite vertices because every PD old triple was already actual in the acute control.',
        incidence_proof='Finite vertices are simple. New old-new edge endpoint incidences give2x+2y<=2a; new-new edge endpoints give y+3z<=2b. With at most2 old neighbors per cap, each cap supports at most1 one-new finite vertex, hence x<=floor(a/2). Three new facets have at most1 common finite vertex, requiring all3 new-new edges. Every geodesic edge has2 ends, so I=2E-3V_f counts all nonfinite end incidences. All old edges/vertices persist, giving deltaI=2(a+b)-3(x+y+z).',
        ideal_vertex_caveat='If ideal vertices are allowed, I includes their incident edge ends. Actual unclosed ends equal I minus the sum of ideal-vertex valences. Therefore this certificate does not assert an unrestricted finite-volume obstruction. The no-ideal hypothesis is explicit; a separate retained-H4 cusp theorem would be needed for a further finite-volume interpretation.',
        continuous_split_subcase=dict(old_anchor_pairs=[[8,12],[12,47]],additional_hypotheses='Exactly those two old finite contacts per new wall; all other native old contacts dotted; all3 prescribed finite endpoint triples PD.',a=4,b=1,x=2,y=1,z=0,change=1,final_nonfinite_end_incidences=19),
        conclusion='A coherent no-ideal release52 repair by at most3 new facets each with at most2 retained intrinsic finite neighbors cannot beat the15 open ends of the frozen control. Successful repairs within this scope need at least4 new facets or a new facet with at least3 retained geodesic neighbors. This is a necessary-incidence gate, not a prohibition of larger coordinated constructions.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run();path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['relaxed_integer_cases']),out['minimum_change'],flush=True)
if __name__=='__main__':main()
