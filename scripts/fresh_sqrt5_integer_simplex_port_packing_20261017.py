#!/usr/bin/env python3
"""Exact integer simplex-port packing on the actual-provenance dual.

Adds all three-port matching inequalities and full-thirteen-port integer
capacities. The geometric proof uses pairwise disjoint maximal-simplex
port blocks. Quad rows receive zero new coefficients.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import fresh_sqrt5_forced_base_dual_repair_20261017 as R
import fresh_sqrt5_actual_provenance_dual_20261017 as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_integer_simplex_port_packing_20261017'
INPUT=A.TAG

def freeze(x):
    return tuple(freeze(v) for v in x) if isinstance(x,(list,tuple)) else x

def code(pieces,marks):
    return tuple(sorted(R.B.C.code(p,marks) for p in pieces))

def run():
    old=R.P.read(INPUT)
    alpha={i:F(a) for i,a in old['face_multipliers']}
    mu={i:F(a) for i,a in old['vertex_multipliers']}
    d,patch,shared=R.domain()
    models=dict(R.B.C.models())
    triples=old['triple_orbit_columns']
    lookup={(c['vertex'],freeze(c['marked_forest_code'])):c['index'] for c in triples}
    assert len(lookup)==len(triples)
    pairlookup={(v['vertex'],freeze(v['marked_forest_code'])):v['index']
                for v in d['variables'][:203565]}
    pair_to_triples=defaultdict(list)
    incidence=[]
    for serial,v in enumerate(sorted(mu)):
        g={}
        for name in d['vertices'][v]['forest']:
            offset=len(g)
            for i,row in models[name].items():
                g[offset+i]={offset+j:m for j,m in row.items()}
        assert len(g)==13
        pieces=R.B.C.pieces(g)
        pair={p:pairlookup[v,code(pieces,p)] for p in combinations(range(13),2)}
        counts=defaultdict(Counter);orbit_sizes=Counter()
        for U in combinations(range(13),3):
            j=lookup[v,code(pieces,U)]
            orbit_sizes[j]+=1
            counts[j].update(pair[p] for p in combinations(U,2))
        assert sum(orbit_sizes.values())==286
        for j,terms in sorted(counts.items()):
            assert orbit_sizes[j]==triples[j]['vertex_cost']
            c=[]
            for p,m in sorted(terms.items()):
                den=d['variables'][p]['multiplicity']
                assert m%den==0
                amount=m//den
                c.append([p,amount]);pair_to_triples[p].append([j,amount])
            assert sum(d['variables'][p]['multiplicity']*m for p,m in c)==3*orbit_sizes[j]
            incidence.append(dict(triple=j,pair_incidence_coefficients=c))
        if serial%200==0:print('INTEGER_PORT_ORBITS',serial,len(mu),flush=True)
    assert len(incidence)==len(triples)
    columns=[dict(index=c['index'],kind='three_port_matching',vertex=c['vertex'],
                  vertex_cost=c['vertex_cost'],marked_forest_code=c['marked_forest_code'])
             for c in triples]
    full={}
    for v in sorted(mu):
        for t in (2,3,4):
            j=len(columns);full[v,t]=j
            columns.append(dict(index=j,kind='full_port_integer',vertex=v,
                                divisor=t,vertex_cost=13//t))
    features={r['face']:r for r in old['positive_face_mixed_features']}
    assert set(features)==set(alpha)
    loads=defaultdict(F);faces=[];corner_count=0
    for fi,a in sorted(alpha.items()):
        face=d['faces'][fi]
        if face['kind']=='canonical_dotted_quad_aggregate':
            continue
        dim={0:1,2:2,8:3,20:4}[face['required']]
        coeff=Counter();ports=Counter()
        for j,m in face['coefficients']:
            var=d['variables'][j]
            if j<203565:
                for tj,n in pair_to_triples.get(j,[]):coeff[tj]+=m*n
            elif var.get('kind')=='maximal_simplex_disjoint_port_charge':
                ports[var['vertex']]+=m
        for j,m in features[fi]['coefficients']:
            coeff[j]-=2*m
        assert all(n>=0 for n in coeff.values()),fi
        # A d-port block meets a triple in at least two elements in
        # C(d,2)(13-d)+C(d,3) ways.
        expected=F(dim*(dim-1),2)*(13-dim)+F(dim*(dim-1)*(dim-2),6)
        assert sum(coeff.values())==(dim+1)*expected
        assert sum(ports.values())==dim*(dim+1)
        corners=[]
        for v,m in sorted(ports.items()):
            assert m%dim==0
            count=m//dim;corner_count+=count;corners.append([v,count])
            for t in (2,3,4):
                coeff[full[v,t]]+=count*(dim//t)
        coeff={j:m for j,m in coeff.items() if m}
        for j,m in coeff.items():loads[j]+=a*m
        faces.append(dict(face=fi,simplex_dimension=dim,vertex_multiplicities=corners,
                          coefficients=[[j,m] for j,m in sorted(coeff.items())]))
    augmented=dict(mu);violations=[]
    for col in columns:
        j,v,cost=col['index'],col['vertex'],col['vertex_cost']
        need=loads[j]/cost
        if need>mu[v]:
            violations.append(dict(column=j,kind=col['kind'],vertex=v,
                 forest=d['vertices'][v]['forest'],old_vertex_mass=str(mu[v]),
                 load=str(loads[j]),required_vertex_mass=str(need)))
        augmented[v]=max(augmented[v],need)
    # All original nonnegative columns still pass, and actual signed
    # provenance loads are unchanged because alpha is unchanged.
    inherited=defaultdict(F)
    for fi,a in alpha.items():
        for j,m in d['faces'][fi]['coefficients']:inherited[j]+=a*m
    for var in d['variables']:
        assert inherited[var['index']]<=var['multiplicity']*augmented.get(var['vertex'],F())
    mixed_loads=defaultdict(F)
    for fi,a in alpha.items():
        for j,m in features[fi]['coefficients']:mixed_loads[j]+=a*m
    for col in triples:
        assert mixed_loads[col['index']]<=col['vertex_cost']*augmented[col['vertex']]
    for col in columns:
        assert loads[col['index']]<=col['vertex_cost']*augmented[col['vertex']]
    signed=R.signed_loads(d,alpha);assert all(x<=0 for x in signed)
    for rec in old['actual_shared_base_constraints']:
        load=sum((m*F(old['provenance_masses'][j]['mass'])
                  for j,m in rec['positive_provenance_coefficients']),F())
        assert load==F(rec['actual_positive_load'])
        assert load<=rec['capacity']*alpha[rec['simplex_face']]
    numerator=sum((a*d['faces'][i]['required'] for i,a in alpha.items()),F())
    assert numerator==F(old['objective_numerator'])
    denominator=sum(augmented.values(),F())
    paths=[f'scripts/{TAG}.py',f'scripts/{R.TAG}.py',f'scripts/{A.TAG}.py',
           f'scripts/{R.B.C.TAG}.py',f'results/{INPUT}.json',
           f'results/{R.N.TAG}.json',f'results/{R.S.TAG}.json',
           f'results/{R.N.I.LP}.json','results/fresh_simplex_disjoint_ports_20261017.md']
    return dict(status='PASS_EXACT_INTEGER_SIMPLEX_PORT_PACKING_DUAL_CHECK',
        sources={p:R.P.sha(p) for p in paths},dimension=13,target=12,
        old_bound=old['bound'],bound=str(numerator/denominator),excludes_target12=numerator/denominator>12,
        objective_numerator=str(numerator),normalizer=str(denominator),
        counts=dict(new_columns=len(columns),triple_columns=len(triples),
            full_port_columns=len(columns)-len(triples),positive_simplex_rows=len(faces),
            reconstructed_corner_multiplicities=corner_count,violating_columns=len(violations),
            violating_forests=len({r['vertex'] for r in violations}),
            increased_vertex_masses=sum(augmented[v]>mu[v] for v in mu)),
        theorem='Disjoint maximal-simplex port blocks C obey sum_C floor(|C intersect U|/t)<=floor(|U|/t). For |U|=3,t=2 this is pair-incidence minus twice complete-triple incidence, and costs1 per literal U. Quad coefficients are0.',
        literal_orbit_pair_triple_incidence=incidence,columns=columns,
        positive_simplex_features=faces,column_loads=[[j,str(loads[j])] for j in range(len(columns))],
        violations=violations,face_multipliers=old['face_multipliers'],
        vertex_multipliers=[[v,str(a)] for v,a in sorted(augmented.items())],
        inherited_signed_loads=list(map(str,signed)),
        scope='All triple subsets and full13-port divisors2,3,4 on the positive-support forests, with exact original/mixed/signed checks and unchanged actual provenance split. No optimizer, full primal export for zero-mass rows, other subset sizes, complete packing polytope or dimension theorem.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();data=json.loads(json.dumps(run()))
    target=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(target.read_text())
    else:target.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],data['bound'],float(F(data['bound'])),flush=True)

