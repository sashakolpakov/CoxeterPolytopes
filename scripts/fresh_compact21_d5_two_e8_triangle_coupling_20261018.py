#!/usr/bin/env python3
"""Complete nonadjacent triangle domains at D5+2E8, then golden-field gates.

Finite extension trees give the domain; arbitrary original spectators E8
are forced orthogonal. No quadrilateral or global charge claim.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import fresh_compact29_small_faces_20261015 as P
import fresh_compact29_small_faces_20261015_scout as S

ROOT=P.ROOT;TAG='fresh_compact21_d5_two_e8_triangle_coupling_20261018'
BASE='fresh_compact21_triangle_transfer_scope_20261018'

def psd(g):
    if any(any(x.v[1:])for row in g for x in row):return P.psd_status(g)
    a=[[x.v[0]for x in row]for row in g];zero=False
    for i in range(len(a)):
        if a[i][i]<0:return -1
        if not a[i][i]:
            if any(a[i][j]for j in range(i+1,len(a))):return -1
            zero=True;continue
        for j in range(i+1,len(a)):
            for k in range(j,len(a)):a[k][j]=a[j][k]=a[j][k]-a[j][i]*a[i][k]/a[i][i]
    return 0 if zero else 1

def extensions(old,e,omit):
    nodes=[i for i in range(old)if i!=omit];parts,_=S.components(nodes,e)
    where={j:k for k,part in enumerate(parts)for j in part};out=[];ledger=sha256();trials=0
    for size in range(4):
        for support in combinations(nodes,size):
            if len({where[i]for i in support})!=size:continue
            for ms in ((3,3,3),)if size==3 else product((3,4,5),repeat=size):
                contact=dict(zip(support,ms));ee=e|{(i,old):m for i,m in contact.items()}
                ok=S.finite(nodes+[old],ee);trials+=1
                ledger.update(json.dumps([sorted(contact.items()),ok],separators=(',',':')).encode()+b'\n')
                if ok:out.append(contact)
    return out,dict(trials=trials,accepted=len(out),digest=ledger.hexdigest())

def compact_principals(old,e,contact):
    ee=e|{(i,old):m for i,m in contact.items()};g=P.base_gram(old,ee,old+1)
    parts,_=S.components(list(range(old+1)),ee);comp=next(c for c in parts if old in c)
    other=sorted(set(comp)-{old});counts=Counter()
    for size in range(1,len(other)+1):
        for ix in combinations(other,size):
            ns=list(ix)+[old];counts['subsets']+=1
            cs,_=S.components(ns,ee)
            if len(cs)!=1:counts['disconnected']+=1;continue
            if any(m==5 and i in ns and j in ns for(i,j),m in ee.items()):
                counts['label5_not_affine']+=1;continue
            if S.finite(ns,ee):counts['finite']+=1;continue
            counts['exact_PSD_tests']+=1
            if psd(P.sub(g,ns))==0:return dict(parabolic_subset=ns,counts=dict(counts))
    return dict(parabolic_subset=None,counts=dict(counts))

def field_gate(contact,e,h):
    # Every simple cycle uses c and the unique old path between two contacts.
    # Old labels are3, so its field membership is that of the endpoint product.
    parts,_=S.components(sorted({v for ij in e for v in ij}),e)
    where={j:k for k,part in enumerate(parts)for j in part};cycles=[]
    for i,j in combinations(sorted(contact),2):
        if where[i]!=where[j]:continue
        z=P.COS[contact[i]]*P.COS[contact[j]]
        cycles.append(dict(endpoints=[i,j],unsigned_product=z,in_ground_field=not(z.v[1]or z.v[3])))
    if any(not r['in_ground_field']for r in cycles):return dict(status='cycle_outside_Qsqrt5',cycles=cycles)
    assert not(h.v[1]or h.v[3])
    conjugate=P.K((h.v[0],0,-h.v[2],0))
    return dict(status='golden_conjugate_PSD'if conjugate.sign()<=0 else'conjugate_not_PSD',
        cycles=cycles,conjugate_height_squared=conjugate)

def integral5(x):
    a,b,c,d=x.v;assert not b and not d
    return (2*c).denominator==1 and(a-c).denominator==1

def run():
    baseline=json.loads((ROOT/f'results/{BASE}.json').read_text())
    assert baseline['status']=='PASS_COMPACT21_TRIANGLE_ONLY_ALLOWANCE_SCOPE_DIAGNOSTIC'
    assert any(r['types']==['D5','E8','E8']for r in baseline['exact_dual']['largest_forest_multipliers'])
    cache={};rows=[];roots={};total=Counter()
    for kind,n,extra in [('D',5,0),('E',8,5)]:
        e=S.diagram(kind,n)
        if extra:e|={(i+n,j+n):m for(i,j),m in S.diagram('D',extra).items()}
        old=n+extra;g=P.base_gram(old,e,old);J=P.inverse(g);assert psd(g)==1
        for a,b in combinations(range(n),2):
            if(a,b)in e:continue
            opts,first=extensions(old,e,a);eligible=[];second=sha256();rejected=0
            for opt in opts:
                for ma in range(2,6):
                    contact=opt|({a:ma}if ma>2 else{});ee=e|{(i,old):m for i,m in contact.items()}
                    ok=S.finite([i for i in range(old+1)if i!=b],ee)
                    second.update(json.dumps([sorted(contact.items()),ok],separators=(',',':')).encode()+b'\n')
                    if not ok:rejected+=1;continue
                    h=sum((P.COS[mi]*J[i][j]*P.COS[mj]for i,mi in contact.items()for j,mj in contact.items()),P.K(-2))
                    r=dict(contacts=sorted(contact.items()),height_squared=h)
                    if h.sign()<=0:r['status']='nonpositive_height'
                    else:
                        key=(kind,tuple(sorted(contact.items())))
                        if key not in cache:cache[key]=compact_principals(old,e,contact)
                        r['compactness_check']=cache[key]
                        if cache[key]['parabolic_subset']is not None:r['status']='parabolic_obstruction'
                        else:
                            r['arithmetic_gate']=field_gate(contact,e,h);r['status']=r['arithmetic_gate']['status']
                            if r['status']=='golden_conjugate_PSD':
                                assert all(i<n for i in contact)
                                roots[key]=dict(component=kind+str(n),contacts=sorted(contact.items()),height_squared=h)
                    eligible.append(r);total[r['status']]+=1
            assert len(eligible)+rejected==4*len(opts)
            rows.append(dict(component=kind+str(n),marked=[a,b],first_corner=first,
                second_corner_trials=4*len(opts),second_corner_rejections=rejected,
                second_corner_digest=second.hexdigest(),eligible=eligible))
    assert len(rows)==27 and len(cache)==52
    passed={k:[r['marked']for r in rows if r['component']==k and any(q['status']=='golden_conjugate_PSD'for q in r['eligible'])]for k in('D5','E8')}
    assert passed=={'D5':[[0,2],[1,3],[1,4]],'E8':[[0,2],[3,5],[3,6],[4,6]]}
    ds=[r for r in roots.values()if r['component']=='D5'];es=[r for r in roots.values()if r['component']=='E8']
    assert len(ds)==5 and len(es)==2
    pairs=[]
    for a,b in product(ds,es):
        q=a['height_squared']*b['height_squared'];ok=integral5(q)
        assert ok==(a['contacts']==[(0,5)])
        pairs.append(dict(D5_contacts=a['contacts'],E8_contacts=b['contacts'],squared_mutual_product=q,
            algebraic_integer_in_ground_field=ok,integral_basis_coordinates=[q.v[0]-q.v[2],2*q.v[2]],
            status='remaining_necessary_pair'if ok else'nonintegral_cyclic_product'))
    paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',f'scripts/{S.TAG}.py',
        'scripts/check_compact29_triangle_star_saved_20261018.py']
    for s in(BASE,'fresh_compact29_three_triangle_star_constraint_20261018'):
        paths += [f'scripts/{s}.py',f'results/{s}.json',f'results/{s}.md']
    paths += ['results/fresh_arithmetic_dimension21_hostile_20261016.md']
    return P.enc(dict(status='PASS_COMPACT21_D5_TWO_E8_NONADJACENT_TRIANGLE_COUPLING',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        vertex_forest=['D5','E8','E8'],rank=21,actual_ground_field='Q(sqrt5)',
        D5_edges=sorted(S.diagram('D',5).items()),E8_edges=sorted(S.diagram('E',8).items()),
        domains=rows,quotient_nonadjacent_domains=27,literal_nonadjacent_faces=48,
        status_counts=dict(total),distinct_positive_normals_tested=52,surviving_roots=list(roots.values()),
        possible_triangle_pairs=passed,cross_pairs=pairs,nonintegral_cross_pairs=8,
        near_D5_triangle_excludes_both_E8_triangle_domains=True,
        maximum_designated_triangles=9,minimum_designated_nontriangles=39,
        arithmetic_conditions_are_necessary_only=True,adjacent_marked_pairs_classified=False,
        new_face_charge_saving=False,new_dimension_bound=None))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert a.verify or not path.exists();out=run()
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['status_counts'],out['maximum_designated_triangles'],flush=True)
