#!/usr/bin/env python3
"""All-label three-finite-corner gate for twenty literal four-wall rings.

The possible 3/5-degree disk would truncate the fourth central corner with
one additional wall.  Its first center is excluded before that second wall
or any other ambient/native/common-center condition is imposed.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, product
import json
from math import acos, floor, pi
from pathlib import Path

import fresh_nonarithmetic_central_ring_contact_gate_20261017 as C

P,U,I,D=C.P,C.U,C.I,C.D
ROOT=C.ROOT
TAG='fresh_nonarithmetic_central_truncation_domain_20261017'
PARENT=ROOT/'results'/f'{C.TAG}.json'


def disk():
    old=list('abcd');new=list('LPRQST')
    triangles=[list(t) for t in ['abL','bcL','cdP','daQ','cLP','dPR','dRQ','aQL',
                                 'TPR','TRQ','TQL','SLP','SLT','SPT']]
    edges=Counter(tuple(sorted(e)) for t in triangles for e in combinations(t,2))
    boundary={tuple(sorted(e)) for e in ['ab','bc','cd','da']}
    assert {e for e,k in edges.items() if k==1}==boundary
    assert all(k==(1 if e in boundary else 2) for e,k in edges.items())
    links={}
    for v in old+new:
        link=defaultdict(set)
        for tri in triangles:
            if v in tri:
                a,b=[i for i in tri if i!=v];link[a].add(b);link[b].add(a)
        seen=set();todo=[next(iter(link))]
        while todo:
            u=todo.pop()
            if u not in seen:seen.add(u);todo.extend(link[u]-seen)
        assert set(link)==seen
        assert sorted(map(len,link.values()))==([1,1]+[2]*(len(link)-2) if v in old else [2]*len(link))
        links[v]={u:sorted(x) for u,x in sorted(link.items())}
    assert (len(links['S']),len(links['T']))==(3,5)
    a=sum((i in old)!=(j in old) for i,j in edges)
    b=sum(i in new and j in new for i,j in edges)
    xyz=Counter(sum(i in new for i in tri) for tri in triangles)
    assert (a,b,xyz[1],xyz[2],xyz[3])==(8,11,4,4,6)
    assert len(old+new)-len(edges)+len(triangles)==1
    assert 2*(a+b)-3*len(triangles)==-4
    return dict(triangles=triangles,links=links,old_cycle=[7,12,52,30],
        ring_order=['L','P','R','Q'],central_degrees=[3,5],
        counts=dict(a=a,b=b,x=4,y=4,z=6),vertices=10,edges=23,triangles_count=14,
        conditional_delta_ends=-4,
        scope='A triangulated disk only. The genuine-boundary, simple, no-ideal and actual geometric realization hypotheses are not proved by topology.')


def literal_rings():
    control=json.loads(C.CONTROL.read_text())['after64']
    front=json.loads(C.FRONTIER.read_text());rings=json.loads(C.RINGS.read_text())
    for index,st in enumerate(rings['stages']):
        if st['best_coherent'] is None:continue
        sp=C.Space(control)
        group=next(g for g in front['changed_seed_pools'] if front['replacement_trials'][g['replacement']]['label_to52']==st['replacement_label'])
        replacement=front['replacement_trials'][group['replacement']]
        n12=sp.add('12new',[sp.nodes[str(i)] for i in [7,8,52]],0,[2,2,st['replacement_label']],replacement['projection'])
        seed=group['seeds'][st['seed']];patch=group['patches'][st['seed']]
        p=patch['P'][st['P']];q=patch['Q'][st['Q']];r=st['best_coherent']['R']
        assert seed['rotation']==0
        L=sp.add('L',[sp.nodes['7'],n12,sp.nodes['52']],seed['profile'],seed['labels'],seed['projection'])
        V=sp.add('P',[sp.nodes['52'],sp.nodes['30'],L],p['profile'],p['labels'],p['projection'])
        W=sp.add('Q',[sp.nodes['30'],sp.nodes['7'],L],q['profile'],q['labels'],q['projection'])
        R=sp.add('R',[sp.nodes['30'],V,W],r['profile'],r['labels'],r['projection'])
        ring=[L,V,R,W];rawring=[seed,p,r,q]
        edge=[p['labels'][2],r['labels'][1],r['labels'][2],q['labels'][2]]
        yield index,sp,ring,rawring,edge


def determinant3(g):
    return U.rnd(g[0][0]*(g[1][1]*g[2][2]-g[1][2]*g[2][1])-
                 g[0][1]*(g[1][0]*g[2][2]-g[1][2]*g[2][0])+
                 g[0][2]*(g[1][0]*g[2][1]-g[1][1]*g[2][0]))


def finite_corners(gram,contacts,remaining):
    rows=[]
    for i,j in remaining:
        g=[[I(2),contacts[i],contacts[j]],
           [contacts[i],gram[i][i],gram[i][j]],
           [contacts[j],gram[j][i],gram[j][j]]]
        minors=[I(2),U.rnd(2*gram[i][i]-contacts[i]*contacts[i]),determinant3(g)]
        rows.append(dict(ring_pair=[i,j],leading_minors=[x.enc() for x in minors],
                         positive_definite=all(x.l>0 for x in minors),
                         negative_minor=any(x.h<0 for x in minors),
                         auxiliary_only=True))
    return rows


def cosine_gap(h):
    # Floating point only proposes neighboring labels.  Acceptance is by the
    # two independent outward rational cosine bounds and strict inequalities.
    estimate=pi/acos(C.mid(h)/2)
    lower=max(6,floor(estimate))
    lo,hi=D.theta(lower),D.theta(lower+1)
    assert lo.h<h.l and h.h<hi.l,('unresolved all-label cosine comparison',h.enc(),estimate,lower,lo.enc(),hi.enc())
    return dict(lower_label=lower,upper_label=lower+1,
        lower_theta=lo.enc(),upper_theta=hi.enc(),
        lower_gap=U.rnd(h-lo).enc(),upper_gap=U.rnd(hi-h).enc())


def run():
    parent=json.loads(PARENT.read_text())
    assert parent['status']=='PASS_EXACT_ALL_LABEL_CENTRAL_CAP_OBSTRUCTION_FOR20_LITERAL_RINGS'
    for name,h in parent['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h
    previous={r['parent_stage']:r for r in parent['records']}
    assert len(previous)==20
    records=[];totals=Counter();gap_labels=Counter();smallest=None
    sqrt3=D.theta(6);box=I(sqrt3.l,Q(2))
    for index,sp,ring,rawring,edge in literal_rings():
        old=previous[index];gram=[[sp.rp(a,b) for b in ring] for a in ring]
        inv,pivots=P.C53.inverse(gram)
        assert old['residual_ring_gram']==C.encmat(gram)
        assert old['ring_edge_labels']==edge and old['recipes']==sp.recipes
        assert old['inverse_pivots']==pivots
        assert len(old['low_contact_exclusions'])==7
        for low in old['low_contact_exclusions']:
            assert low['contact_quadruples']==256
            assert sum(low['sign_counts'].values())==256
            assert Q(low['minimum_absolute_norm_error'])>0
        cases=[]
        for omitted in range(4):
            remaining=[(i,(i+1)%4) for i in range(4) if i!=omitted]
            neighbors={i:{j for e in remaining if i in e for j in e if j!=i} for i in range(4)}
            eligible=[i for i in range(4) if rawring[i]['profile']==0 and
                      all(edge[j]==2 for j in range(4) if j!=omitted and i in (j,(j+1)%4))]
            high_sets=[list(h) for size in range(1,5) for h in combinations(eligible,size)
                       if all(not set(e).issubset(h) for e in remaining)]
            assert max(map(len,high_sets),default=0)<=2
            single=[];double=[];counts=Counter()
            for j in eligible:
                other=[i for i in range(4) if i!=j]
                alphabets=[(2,) if i in neighbors[j] else range(2,6) for i in other]
                frame=[ring[i] for i in other]
                for labels in product(*alphabets):
                    sols,recipe=sp.solutions(frame,0,labels)
                    counts['one_high_contact_frames']+=1
                    branches=[]
                    if not sols:counts['negative_norm_residual']+=1
                    for sign,v in sols:
                        val=sp.solution_pair(frame,recipe,sign,ring[j]);h=U.rnd(-val)
                        rec=dict(branch=sign,product=val.enc(),positive_bond=h.enc())
                        counts['one_high_real_branches']+=1
                        if h.h<sqrt3.l:
                            rec.update(status='BELOW_HIGH_LABEL_RANGE',gap=U.rnd(sqrt3-h).enc())
                            counts['one_high_below_range']+=1
                        elif h.l>=2:
                            rec.update(status='AT_OR_ABOVE_PARALLEL_THRESHOLD',gap=U.rnd(h-2).enc())
                            counts['one_high_above_range']+=1
                        else:
                            assert h.l>sqrt3.h and h.h<2,('unresolved high-range boundary',index,omitted,j,labels,sign,h.enc())
                            gap=cosine_gap(h);gap_labels[gap['lower_label']]+=1
                            contacts={i:-D.theta(m) for i,m in zip(other,labels)};contacts[j]=val
                            corner_data=finite_corners(gram,contacts,remaining)
                            pd=all(r['positive_definite'] for r in corner_data)
                            rec.update(status='STRICT_CONSECUTIVE_COSINE_GAP',cosine_gap=gap,
                                       finite_corner_checks=corner_data,all_three_corners_finite=pd)
                            counts['one_high_cosine_gaps']+=1
                            counts['gap_branches_three_certified_finite_corners' if pd else 'gap_branches_not_all_certified_finite']+=1
                            minimum=min(Q(gap['lower_gap'][0]),Q(gap['upper_gap'][0]))
                            if smallest is None or minimum<Q(smallest['gap_lower_bound']):
                                smallest=dict(parent_stage=index,omitted_corner=omitted,high_vertex=j,
                                    labels=list(labels),branch=sign,lower_label=gap['lower_label'],
                                    gap_lower_bound=str(minimum),all_three_corners_finite=pd)
                        branches.append(rec)
                    single.append(dict(high_vertex=j,other_ring_vertices=other,low_labels=list(labels),recipe=recipe,branches=branches))
            for hh in high_sets:
                if len(hh)!=2:continue
                forced=set().union(*(neighbors[i] for i in hh))
                assert forced==set(range(4))-set(hh)
                j,k=hh;a,b,c=inv[j][j],inv[j][k],inv[k][k]
                du=U.rnd(2*a*box+2*b*box);dv=U.rnd(2*b*box+2*c*box)
                assert (du.l>0 or du.h<0) and (dv.l>0 or dv.h<0),('nonmonotone high rectangle',index,omitted,hh,du.enc(),dv.enc())
                u=sqrt3 if du.l>0 else I(2);v=sqrt3 if dv.l>0 else I(2)
                minimum=U.rnd(a*u*u+2*b*u*v+c*v*v)
                assert minimum.l>2,('two-high norm not excluded',index,omitted,hh,minimum.enc())
                double.append(dict(high_vertices=hh,low_vertices=sorted(forced),
                    quadratic_coefficients=[a.enc(),b.enc(),c.enc()],
                    closed_bond_box=box.enc(),derivative_u=du.enc(),derivative_v=dv.enc(),
                    minimizing_corner=[u.enc(),v.enc()],minimum_norm=minimum.enc(),
                    gap_above_norm_two=U.rnd(minimum-2).enc()))
                counts['two_high_rectangles']+=1
            totals.update(counts);totals['omitted_corner_cases']+=1
            cases.append(dict(omitted_ring_corner=[omitted,(omitted+1)%4],
                retained_corner_edges=[list(e) for e in remaining],eligible_high_vertices=eligible,
                possible_high_vertex_sets=high_sets,counts=dict(counts),one_high_cases=single,two_high_cases=double))
        records.append(dict(parent_stage=index,ring_edge_labels=edge,
            inherited_low_contact_count=7*256,ring_gram=C.encmat(gram),inverse_pivots=pivots,cases=cases))
        totals['rings']+=1;totals['inherited_low_contact_quadruples']+=7*256
        print('TRUNCATION_RING',index,dict(totals),flush=True)
    assert totals['rings']==20 and totals['omitted_corner_cases']==80
    assert totals['inherited_low_contact_quadruples']==35840
    assert totals['one_high_real_branches']==1408
    assert totals['one_high_cosine_gaps']==41
    assert totals['gap_branches_three_certified_finite_corners']==13
    assert totals['two_high_rectangles']==6
    paths=[Path(__file__),Path(C.__file__),PARENT,C.CONTROL,C.FRONTIER,C.RINGS,
           Path(C.E.__file__),Path(C.F.__file__),Path(D.__file__)]
    front=json.loads(C.FRONTIER.read_text())
    return dict(status='PASS_ALL_LABEL_THREE_FINITE_CORNER_CENTRAL_GATE_FOR20_LITERAL_RINGS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        abstract_disk=disk(),counts=dict(totals),consecutive_gap_label_histogram=dict(sorted(gap_labels.items())),
        smallest_certified_gap=smallest,records=records,
        unchanged_control=dict(walls=64,ambient_finite_vertices=198,S6_finite_vertices=25,S6_exterior_vertices=7,S6_open_ends=15,volume='infinite, inherited exact strict escape'),
        changed_or_deleted_old_walls=[12,47,53],inherited_deletion_cost=front['deletion_cost'],
        theorem='For each of the20 exact specified ring recipes and each choice of omitted consecutive ring corner, no norm2 Coxeter wall inducing an actual S6 facet has finite contacts to all four ring walls and makes the other three consecutive ring corners finite. Thus the first center of the3/5 two-center truncation topology is already impossible.',
        scope='All labels and arbitrary real coefficient fields, but only the20 literal recipe rings. Low labels2..5 inherit the exact four-contact norm exclusion after literal recipe/Gram equality. One-high labels are solved from the other three contacts and separated between consecutive exact cosine thresholds; two-high cases are excluded on the whole closed real bond square[sqrt3,2]^2. No claim of ambient coherence for the parent numerical rings, no new chamber, and no obstruction to changed rings, caps missing S6 or different patch topology.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
