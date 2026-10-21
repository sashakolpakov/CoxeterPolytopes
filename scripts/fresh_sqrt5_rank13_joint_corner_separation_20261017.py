#!/usr/bin/env python3
"""Full finite-type coverage audit and two exact omitted compact local quads."""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_adjacent_discount_witness_20261016 as A
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_rank13_joint_corner_separation_20261017'
LP='fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017'
WITNESS='fresh_sqrt5_e6_discount_witness_20261016'
K,Z=A.K,A.Z

def dump(x):return json.dumps(x,separators=(',',':'))
def read(tag):
    d=json.loads((ROOT/f'results/{tag}.json').read_text())
    for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return d

def direct_sum(M,N):
    a,b=len(M),len(N);out=[[Z]*(a+b)for _ in range(a+b)]
    for i in range(a):
        for j in range(a):out[i][j]=M[i][j]
    for i in range(b):
        for j in range(b):out[a+i][a+j]=N[i][j]
    return out

def geometry(rec,spectators,graphs,finite,lookup,weights):
    row=rec['row'];n=row['nodes'];G=[[K(*map(F,x))for x in r]for r in rec['scaled_gram']]
    assert A.inertia(G)==[n-2,1,1]and A.inertia([[x.conj()for x in r]for r in G])==[n-1,0,1]
    ids=rec['basis_indices'];Q=A.sub(G,ids);Qi=A.inverse(Q)
    roots=[A.mv(Qi,[G[j][i]for j in ids])for i in range(n)]
    assert [[A.pair(Q,x,y)for y in roots]for x in roots]==G
    labels=dict(zip(combinations(range(n),2),row['labels']));g=C.graph(range(n),labels)
    spectator_nodes=[]
    for name in spectators:
        h=graphs[name];off=len(g);spectator_nodes+=list(range(off,off+len(h)))
        assert all(m==3 for r in h.values()for m in r.values())
        for i,r in h.items():g[off+i]={off+j:m for j,m in r.items()}
    assert len(g)==15
    r=len(spectator_nodes);S=[[K(2)if i==j else K(-1)if spectator_nodes[j]in g[spectator_nodes[i]]else Z for j in range(r)]for i in range(r)]
    assert A.inertia(S)==[r,0,0]
    form=direct_sum(Q,S);normals=[v+[Z]*r for v in roots]
    for i in range(r):
        v=[Z]*14;v[n-1+i]=K(1);normals.append(v)
    big=direct_sum(G,S)
    assert [[A.pair(form,x,y)for y in normals]for x in normals]==big
    assert A.inertia(big)==[13,1,1]and A.inertia([[x.conj()for x in t]for t in big])==[14,0,1]
    assert A.inertia(form)==[13,1,0]and A.inertia([[x.conj()for x in t]for t in form])==[14,0,0]
    assert all((2*big[i][j]/big[j][j]).integral()for i in range(15)for j in range(15))
    # Check every normalized bond independently of the stored Gram scaling.
    phi=K(F(1,2),F(1,2));squares={2:Z,3:K(1),5:1+phi,10:2+phi}
    for i,j in combinations(range(15),2):
        m=g[i].get(j,2);assert 4*big[i][j]*big[i][j]/(big[i][i]*big[j][j])==squares[m]
        assert big[i][j].sign()<=0
    center=A.mv(A.inverse(form),[K(-1)]*(n-1)+[K(F(-1,1024))]*r)
    products=[A.pair(form,center,v)for v in normals]
    assert A.pair(form,center,center).sign()<0 and all(x.sign()<0 for x in products)
    white=row['white']+spectator_nodes;assert len(white)==11
    corners=[]
    for pair in [(0,1),(1,2),(2,3),(0,3)]:
        active=white+list(pair);B=A.sub(big,active)
        assert A.inertia(B)==[13,0,0]and A.inertia([[x.conj()for x in t]for t in B])==[13,0,0]
        coeff=A.mv(A.inverse(B),[products[i]for i in active])
        v=[center[j]-sum((c*normals[i][j]for c,i in zip(coeff,active)),Z)for j in range(14)]
        ip=[A.pair(form,v,z)for z in normals]
        assert all(not ip[i]for i in active)and all(ip[i].sign()<0 for i in range(15)if i not in active)
        assert A.pair(form,v,v).sign()<0 and A.pair(form,v,center).sign()<0
        nodes=set(active);h={i:{j:m for j,m in g[i].items()if j in nodes}for i in nodes};pieces=C.pieces(h)
        forest=tuple(sorted(finite[C.code(p)]for p in pieces));code=dump(sorted(C.code(p,pair)for p in pieces))
        index=lookup[forest,code]
        corners.append(dict(black=list(pair),forest=list(forest),marked_forest_code=json.loads(code),
            variable=index,payment=str(weights.get(index,F(0))),coordinates=[x.enc()for x in v],
            norm=A.pair(form,v,v).enc(),products=[x.enc()for x in ip]))
    tests=[];accepted=[]
    for pair in combinations(range(4),2):
        ker=A.kernel([A.mv(form,normals[i])for i in white+list(pair)]);assert len(ker)==1
        orientations=[]
        for sign in(-1,1):
            v=[sign*x for x in ker[0]];ip=[A.pair(form,v,z)for z in normals]
            if all(x.sign()<=0 for x in ip):
                assert A.pair(form,v,v).sign()<0 and A.pair(form,v,center).sign()<0
                accepted.append(list(pair));orientations.append(sign)
        tests.append(dict(black=list(pair),feasible_orientations=orientations))
    assert sorted(accepted)==sorted(c['black']for c in corners)
    payment=sum((F(c['payment'])for c in corners),F(0));assert payment==0
    return dict(class_index=rec['class_index'],spectator_forest=spectators,
        full_row=dict(nodes=15,white=white,labels=[g[i].get(j,2)for i,j in combinations(range(15),2)]),
        gram=A.enc(big),basis_form=A.enc(form),root_coordinates=A.enc(normals),
        physical_inertia=[13,1,1],conjugate_inertia=[14,0,1],integral_cartan_entries=225,
        strict_center=[x.enc()for x in center],center_norm=A.pair(form,center,center).enc(),
        corners=corners,complete_face_cone=tests,total_payment=str(payment),required_payment='1',
        corner_variable_indices=[c['variable']for c in corners])

def run():
    lp=read(LP);old=read(WITNESS);models=[(name,g)for name,g in C.models()if len(g)<=13]
    graphs=dict(models);finite={C.code(g):name for name,g in models}
    primal=next(c for c in lp['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    weights={i:F(x)for i,x in primal['weights']};used={i for f in lp['faces']for i,c in f['coefficients']}
    lookup={(tuple(v['forest']),dump(v['marked_forest_code'])):v['index']for v in lp['variables']}
    entries=[(name,len(g),sum((C.base(g,a,b)for a,b in combinations(g,2)),F(0)))for name,g in models]
    def forests(rank,start=0):
        if not rank:yield (),F(0);return
        for j in range(start,len(entries)):
            name,n,b=entries[j]
            if n<=rank:
                for tail,c in forests(rank-n,j):yield tuple(sorted((name,)+tail)),b+c
    all_forests=dict(forests(13));present={tuple(v['forest'])for v in lp['vertices']}
    forest_census=[dict(forest=list(f),baseline=str(b),in_partial_model=f in present)for f,b in sorted(all_forests.items())]
    counts=dict(all_rank13_finite_forests=len(all_forests),present_forests=len(present),
        absent_forests=len(set(all_forests)-present),baseline_at_least12=sum(b>=12 for b in all_forests.values()),
        baseline_above12=sum(b>12 for b in all_forests.values()),
        absent_baseline_at_least12=sum(b>=12 for f,b in all_forests.items()if f not in present),
        absent_baseline_above12=sum(b>12 for f,b in all_forests.items()if f not in present))
    assert counts==dict(all_rank13_finite_forests=6789,present_forests=630,absent_forests=6159,
        baseline_at_least12=3350,baseline_above12=2681,absent_baseline_at_least12=2794,absent_baseline_above12=2174)
    changes=[];orbit_counts=Counter();changed_vertices=set()
    for vertex in lp['vertices']:
        g={}
        for name in vertex['forest']:
            h=graphs[name];off=len(g)
            for i,row in h.items():g[off+i]={off+j:m for j,m in row.items()}
        where={i:p for p in C.pieces(g)for i in p}
        for index,multiplicity in vertex['variables']:
            var=lp['variables'][index];values={C.base(where[a],a,b)if b in where[a]else F(0)for a,b in var['literal_pairs']}
            assert len(values)==1;base=next(iter(values));weight=weights.get(index,F(0))
            orbit_counts['variables']+=1;orbit_counts['used_in_faces']+=index in used;orbit_counts['nonzero']+=bool(weight)
            orbit_counts['discounted']+=weight<base;orbit_counts['discounted_unused']+=weight<base and index not in used
            orbit_counts['increased']+=weight>base
            if weight!=base:
                changes.append(dict(variable=index,vertex=vertex['index'],baseline=str(base),payment=str(weight),
                    in_saved_faces=index in used,multiplicity=multiplicity))
            if weight<base:changed_vertices.add(vertex['index'])
    assert dict(orbit_counts)==dict(variables=21492,used_in_faces=3634,nonzero=839,discounted=5844,discounted_unused=4310,increased=648)
    witnesses=[geometry(next(r for r in old['records']if r['class_index']==ci),spectators,graphs,finite,lookup,weights)
        for ci,spectators in [(4266,['A2','D5']),(7746,['D5'])]]
    paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'scripts/{C.TAG}.py',f'results/{LP}.json',f'results/{WITNESS}.json']
    return dict(status='PASS_EXACT_RANK13_COVERAGE_AUDIT_AND_TWO_ZERO_PAYMENT_LOCAL_QUADS',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},counts=counts,
        forest_census=forest_census,orbit_counts=dict(orbit_counts),discounted_vertex_count=len(changed_vertices),
        changed_orbits=changes,witnesses=witnesses,
        completeness='Finite Coxeter components over the squared-label domain Qsqrt5 are A,B,D,E,F4,H3,H4 and I2(5),I2(6),I2(10), with A2 and B2 representing I2(3),I2(4). Every nondecreasing component multiset of totalrank13 is enumerated. This is finite vertex-type completeness, not face or compact realizability completeness.',
        geometry='Each omitted face is an independently rechecked frozen compact local E6 quadrilateral, orthogonally extended by the displayed spherical spectator forest to rank13. All225 Cartan entries are integral. The actual form has inertia(13,1) physically and is positive definite at the other embedding; all four finite future endpoint rays and all six possible black-pair intersections are checked exactly.',
        scope='The two actual compact local quadrilaterals disprove coverage of the frozen partial payment vector, even on already represented whole vertex types. They are local admissible integral-reflection chambers, not claimed compact finite-volume ambient polytopes. The old partial-domain certificate remains valid on its own input rows. No dimension theorem is claimed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],[(w['class_index'],w['total_payment'])for w in d['witnesses']],flush=True)
