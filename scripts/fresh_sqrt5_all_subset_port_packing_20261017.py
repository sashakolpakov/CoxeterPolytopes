#!/usr/bin/env python3
"""Exact all-subset integer packing separation of the saved rank13 dual.

Every maximal simplex corner has at most four ports. Uniformly averaging its
literal block over its labelled-forest automorphism orbit lets subset zeta
transforms evaluate all8192 subsets, with divisors2,3,4, exactly.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb, lcm
from pathlib import Path
import json
import fresh_sqrt5_labelled_packing_full_model_20261017 as M

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_all_subset_port_packing_20261017'
STARS='fresh_sqrt5_rank13_full_triangle_simplex_separation_20261017'
B=M.R.B
C=B.C
N=13


def freeze(x):
    return tuple(freeze(a) for a in x) if isinstance(x,(list,tuple)) else x


def code(pieces,marks=()):
    return tuple(sorted(C.code(p,marks) for p in pieces))


def restrict(g,omitted):
    return {i:{j:m for j,m in row.items() if j!=omitted}
            for i,row in g.items() if i!=omitted}


def corner_features(d,alpha):
    """Reconstruct actual full blocks, including all positive four-simplexes."""
    old,patch,_=M.R.domain()
    stars=M.R.P.read(STARS)['stars']
    incidence=M.R.P.read(M.R.N.I.TAG)
    models=dict(C.models());finite={C.code(g):name for name,g in models.items()}
    vlookup={tuple(v['forest']):v['index'] for v in old['vertices']}
    pairlookup={(v['vertex'],freeze(v['marked_forest_code'])):v['index']
                for v in old['variables'][:203565]}
    portlookup={(v['vertex'],freeze(v['marked_forest_code'])):v['index']
                for v in old['variables'] if v.get('kind')=='maximal_simplex_disjoint_port_charge'}
    groups={};columns=[];faces=[];dimensions=Counter();corners=0
    for fi,a in sorted(alpha.items()):
        face=d['faces'][fi]
        if face['kind']=='canonical_dotted_quad_aggregate':continue
        src=old['faces'][fi];dim={0:1,2:2,8:3,20:4}[face['required']]
        if 'provenance' in src:
            assert len(src['provenance'])==1
            star=stars[src['provenance'][0]['star']]
            g=C.graph(range(14),{(i,j):m for i,j,m in star['edges']})
            S={r['omitted'] for r in star['spherical_deletions']}
        else:
            core=incidence['cores'][src['core']]
            g={i:{} for i in range(core['nodes'])}
            for i,j,(m,aa,bb) in core['edges']:g[i][j]=g[j][i]=m
            offset=core['nodes']
            for name in src['spectator_forest']:
                for i,row in models[name].items():g[offset+i]={offset+j:m for j,m in row.items()}
                offset+=len(models[name])
            S=set(core['spherical_deletions'])
        assert len(g)==14 and len(S)==dim+1
        cc=Counter();pp=Counter();uu=Counter()
        for omit in sorted(S):
            h=restrict(g,omit);pieces=C.pieces(h)
            assert len(h)==13 and all(m>0 for row in h.values() for m in row.values())
            forest=tuple(sorted(finite[C.code(p)] for p in pieces));v=vlookup[forest]
            block=tuple(sorted(S-{omit}));key=(v,dim,code(pieces,block))
            if key not in groups:
                j=len(columns);groups[key]=j
                columns.append(dict(index=j,vertex=v,block_size=dim,marked_forest_code=key[2]))
            cc[groups[key]]+=1
            pp.update(pairlookup[v,code(pieces,p)] for p in combinations(block,2))
            uu.update(portlookup[v,code(pieces,(p,))] for p in block)
        assert dict(pp)=={j:m for j,m in src['coefficients'] if j<203565},fi
        assert dict(uu)=={j:m for j,m in src['coefficients'] if old['variables'][j].get('kind')=='maximal_simplex_disjoint_port_charge'},fi
        assert sum(cc.values())==dim+1
        faces.append(dict(face=fi,simplex_dimension=dim,block_orbits=[list(x) for x in sorted(cc.items())]))
        dimensions[dim]+=1;corners+=dim+1
    return columns,faces,dict(dimensions),corners


def zeta(values):
    values=list(values)
    for bit in range(N):
        step=1<<bit
        for start in range(0,1<<N,step*2):
            for k in range(start+step,start+2*step):values[k]+=values[k-step]
    return values


def separate(d,alpha,mu,columns,features):
    models=dict(C.models());weights=defaultdict(F);byvertex=defaultdict(list)
    assert {fi for fi,a in alpha.items() if a and d['faces'][fi]['kind']!='canonical_dotted_quad_aggregate'}<={f['face'] for f in features}
    for f in features:
        for j,m in f['block_orbits']:weights[j]+=alpha.get(f['face'],F())*m
    assert all(a>=0 for a in weights.values())
    for c in columns:
        if weights[c['index']]:byvertex[c['vertex']].append(c['index'])
    repaired=dict(mu);records=[];orbits=[];critical=[];totals=Counter();digest=sha256()
    for serial,(v,js) in enumerate(sorted(byvertex.items())):
        g={}
        for name in d['budgets'][v]['forest']:
            off=len(g)
            for i,row in models[name].items():g[off+i]={off+j:m for j,m in row.items()}
        assert len(g)==N;pieces=C.pieces(g)
        wanted={(columns[j]['block_size'],freeze(columns[j]['marked_forest_code'])):j for j in js}
        literal=defaultdict(list)
        for size in sorted({columns[j]['block_size'] for j in js}):
            for block in combinations(range(N),size):
                j=wanted.get((size,code(pieces,block)))
                if j is not None:literal[j].append(sum(1<<p for p in block))
        assert set(literal)==set(js)
        blockweights=defaultdict(F)
        for j in js:
            values=literal[j];assert len(set(values))==len(values)>0
            c=columns[j]
            orbits.append(dict(index=j,vertex=v,block_size=c['block_size'],literal_masks=values))
            for mask in values:blockweights[mask]+=weights[j]/len(values)
        den=lcm(*(a.denominator for a in blockweights.values()))
        ints={mask:int(a*den) for mask,a in blockweights.items()}
        assert sum(ints.values())==den*sum((weights[j] for j in js),F())
        raw={k:[0]*(1<<N) for k in (2,3,4)}
        sizeweights=Counter()
        for block,w in ints.items():
            ports=[i for i in range(N) if block>>i&1];sizeweights[len(ports)]+=w
            for k in (2,3,4):
                for sub in combinations(ports,k):raw[k][sum(1<<p for p in sub)]+=w
        zz={k:zeta(raw[k]) for k in raw}
        for t in (2,3,4):
            bestnum=0;bestcost=1;bestmask=(1<<t)-1;checks=violations=0;size_sums=Counter()
            old=mu.get(v,F())
            for mask in range(1<<N):
                size=mask.bit_count();cost=size//t
                value=(zz[2][mask]-2*zz[3][mask]+4*zz[4][mask] if t==2 else
                       zz[3][mask]-3*zz[4][mask] if t==3 else zz[4][mask])
                assert value>=0
                size_sums[size]+=value
                if not cost:assert value==0;continue
                checks+=1
                if value*old.denominator>old.numerator*den*cost:violations+=1
                if value*bestcost>bestnum*cost:bestnum,bestcost,bestmask=value,cost,mask
            for size in range(N+1):
                expected=0
                for k,w in sizeweights.items():
                    expected+=w*sum(comb(k,j)*comb(N-k,size-j)*(j//t)
                        for j in range(k+1) if 0<=size-j<=N-k)
                assert size_sums[size]==expected,(v,t,size)
            assert bestnum==sum(w*((block&bestmask).bit_count()//t) for block,w in ints.items())
            need=F(bestnum,den*bestcost)
            repaired[v]=max(repaired.get(v,F()),need)
            row=dict(vertex=v,forest=d['budgets'][v]['forest'],divisor=t,
                checked_literal_subsets=checks,violating_literal_subsets=violations,
                best_subset_mask=bestmask,best_subset=[i for i in range(N) if bestmask>>i&1],
                best_subset_size=bestmask.bit_count(),vertex_cost=bestcost,
                maximum_load=str(F(bestnum,den)),required_vertex_mass=str(need),old_vertex_mass=str(old))
            records.append(row);totals['checked_literal_subsets']+=checks
            totals['violating_literal_subsets']+=violations
            # Rationally averaged charge for a maximizing subset. These
            # exact coefficients allow a later explicitly restricted model.
            if need>old:
                coefficients={j:F(sum((b&bestmask).bit_count()//t for b in literal[j]),len(literal[j])) for j in js}
                terms=[]
                for f in features:
                    value=sum((m*coefficients.get(j,F()) for j,m in f['block_orbits']),F())
                    if value:terms.append([f['face'],str(value)])
                assert sum((alpha.get(fi,F())*F(a) for fi,a in terms),F())==F(bestnum,den)
                critical.append(dict(index=len(critical),vertex=v,divisor=t,subset_mask=bestmask,
                    vertex_cost=bestcost,positive_simplex_coefficients=terms))
            digest.update((B.dump(row)+'\n').encode())
        if serial%100==0:print('ALL_SUBSET_FORESTS',serial,'/',len(byvertex),flush=True)
    return repaired,records,orbits,critical,dict(totals),digest.hexdigest()


def run():
    saved=json.loads((ROOT/f'results/{M.TAG}.json').read_text());d,prepared=M.build()
    assert saved['prepared']==json.loads(json.dumps(prepared))
    dual=next(c for c in saved['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    M.verify_dual(d,dual)
    alpha={i:F(a) for i,a in dual['face_multipliers']};mu={v:F(a) for v,a in dual['vertex_multipliers']}
    assert all(a>=0 for a in alpha.values()) and all(a>=0 for a in mu.values())
    columns,features,dimensions,corners=corner_features(d,alpha)
    print('FULL_BLOCK_FEATURES',len(features),len(columns),corners,dimensions,flush=True)
    augmented,records,orbits,critical,totals,trace=separate(d,alpha,mu,columns,features)
    # Alpha and every actual-provenance mass stay fixed. Increasing mu can
    # only improve every inherited positive-cost constraint; signed loads
    # remain exactly those checked above by the full preserved model.
    assert all(augmented[v]>=a for v,a in mu.items())
    numerator=F(dual['objective_numerator']);denominator=sum(augmented.values(),F())
    assert numerator==sum((a*d['faces'][i]['required'] for i,a in alpha.items()),F())
    paths=[f'scripts/{TAG}.py',f'scripts/{M.TAG}.py',f'results/{M.TAG}.json',
        f'results/{STARS}.json',f'results/{M.R.N.TAG}.json',f'results/{M.R.N.I.TAG}.json',
        f'results/{M.R.N.I.LP}.json',f'scripts/{C.TAG}.py',
        'results/fresh_simplex_disjoint_ports_20261017.md',
        'results/fresh_sqrt5_integer_simplex_port_packing_20261017.md']
    return dict(status='PASS_EXACT_ALL_SUBSET_INTEGER_PACKING_SEPARATION',
        sources={p:M.R.P.sha(p) for p in paths},dimension=13,target=12,
        input_model_sha256=prepared['model_sha256'],old_bound=dual['bound'],
        bound=str(numerator/denominator),excludes_target12=numerator/denominator>12,
        objective_numerator=str(numerator),normalizer=str(denominator),
        counts=dict(totals,positive_simplex_faces=len(features),corner_occurrences=corners,
            occupied_block_orbits=len(columns),tested_vertex_forests=len({r['vertex'] for r in records}),
            positive_vertex_forests=len(mu),violating_forests=len({r['vertex'] for r in records if r['violating_literal_subsets']}),
            critical_columns=len(critical)),simplex_dimensions=dimensions,
        block_orbit_columns=columns,literal_block_orbits=sorted(orbits,key=lambda r:r['index']),
        positive_simplex_features=features,subset_maxima=records,critical_columns=critical,
        face_multipliers=dual['face_multipliers'],vertex_multipliers=[[v,str(a)] for v,a in sorted(augmented.items())],
        checked_record_digest=trace,
        theorem='For disjoint maximal-simplex port blocks C, sum_C floor(|C intersect U|/t)<=floor(|U|/t). On compact corners |C|<=4: floor(k/2)=C(k,2)-2C(k,3)+4C(k,4), floor(k/3)=C(k,3)-3C(k,4), floor(k/4)=C(k,4). All8192 literal U and all nonzero divisors2,3,4 are checked after exact orbit averaging.',
        validation='Full saved dual replay; reconstructed full corner blocks reproduce every inherited pair and one-port coefficient; zeta arrays checked nonnegative; every subset-size aggregate checked by an independent hypergeometric counting identity; maximizing loads checked directly from literal block intersections. All old loads remain valid under monotone mu increase.',
        scope='Exact fixed-alpha separation on the optimized saved partial model. All positive simplex rows have complete full-block features; all quad payments from this new family are zero. No optimizer, primal feasibility, complete packing polytope, missing global face coverage or dimension claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(out.read_text())
    else:out.write_text(B.dump(result)+'\n')
    print(result['status'],result['counts'],result['bound'],float(F(result['bound'])),flush=True)
