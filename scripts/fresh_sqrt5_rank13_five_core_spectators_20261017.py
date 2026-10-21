#!/usr/bin/env python3
"""All finite rank13 spectator extensions of five certified compact local quads."""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_adjacent_discount_witness_20261016 as A
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_rank13_five_core_spectators_20261017'
PRE='fresh_sqrt5_rank13_joint_corner_repair1_20261017'
SOURCES=['fresh_sqrt5_e6_discount_witness_20261016','fresh_sqrt5_adjacent_discount_witness_20261016',
    'fresh_sqrt5_multiplicity_geometry_review_20261016','fresh_sqrt5_global_b6_discount_witness_20261016',
    'fresh_sqrt5_adjacent_discount_frontier_20261016']
K,Z=A.K,A.Z;PHI=K(F(1,2),F(1,2));SQUARED={2:Z,3:K(1),4:K(2),5:1+PHI,6:K(3),10:2+PHI}

def dump(x):return json.dumps(x,separators=(',',':'))
def read(tag):
    d=json.loads((ROOT/f'results/{tag}.json').read_text())
    for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return d

def finite_component(name,g):
    n=len(g);scale={0:K(1)};parent={};queue=[0]
    for i in queue:
        for j,m in g[i].items():
            if j not in scale:scale[j]=scale[i]*SQUARED[m];parent[j]=i;queue.append(j)
    assert len(scale)==n and all(x.sign()>0 and x.conj().sign()>0 for x in scale.values())
    G=[[2*scale[i]if i==j else Z for j in range(n)]for i in range(n)]
    for j,i in parent.items():G[i][j]=G[j][i]=-scale[j]
    assert A.inertia(G)==[n,0,0]and A.inertia([[x.conj()for x in r]for r in G])==[n,0,0]
    assert all((2*G[i][j]/G[j][j]).integral()for i in range(n)for j in range(n))
    for i,j in combinations(range(n),2):assert 4*G[i][j]*G[i][j]/(G[i][i]*G[j][j])==SQUARED[g[i].get(j,2)]
    inverse=A.inverse(G);p=sum((x for r in inverse for x in r),Z);assert p.sign()>0
    return dict(name=name,rank=n,edges=[[i,j,m]for i in g for j,m in g[i].items()if i<j],
        scales=[scale[i].enc()for i in range(n)],gram=A.enc(G),positive_center_cost=p.enc(),
        physical_inertia=[n,0,0],conjugate_inertia=[n,0,0],integral_cartan_entries=n*n)

def core_geometry(key,row,encoded,basis):
    n=row['nodes'];G=[[K(*map(F,x))for x in r]for r in encoded]
    assert A.inertia(G)==[n-2,1,1]and A.inertia([[x.conj()for x in r]for r in G])==[n-1,0,1]
    Q=A.sub(G,basis);Qi=A.inverse(Q);N=[A.mv(Qi,[G[j][i]for j in basis])for i in range(n)]
    assert [[A.pair(Q,x,y)for y in N]for x in N]==G
    assert A.inertia(Q)==[n-2,1,0]and A.inertia([[x.conj()for x in r]for r in Q])==[n-1,0,0]
    assert all((2*G[i][j]/G[j][j]).integral()for i in range(n)for j in range(n))
    for(i,j),m in zip(combinations(range(n),2),row['labels']):
        assert G[i][j].sign()<=0 and 4*G[i][j]*G[i][j]/(G[i][i]*G[j][j])==SQUARED[m]
    H=A.mv(Qi,[K(-1)]*(n-1));hp=[A.pair(Q,H,x)for x in N];hq=A.pair(Q,H,H)
    assert hq.sign()<0 and all(x.sign()<0 for x in hp)
    white=row['white'];vertices=[]
    for pair in [(0,1),(1,2),(2,3),(0,3)]:
        ids=white+list(pair);B=A.sub(G,ids)
        assert A.inertia(B)==[n-2,0,0]and A.inertia([[x.conj()for x in r]for r in B])==[n-2,0,0]
        coefficients=A.mv(A.inverse(B),[hp[i]for i in ids])
        v=[H[j]-sum((c*N[i][j]for c,i in zip(coefficients,ids)),Z)for j in range(n-1)]
        ips=[A.pair(Q,v,x)for x in N]
        assert all(not ips[i]for i in ids)and all(ips[i].sign()<0 for i in range(n)if i not in ids)
        assert A.pair(Q,v,v).sign()<0 and A.pair(Q,v,H).sign()<0
        vertices.append(dict(black=list(pair),coordinates=[x.enc()for x in v],norm=A.pair(Q,v,v).enc(),products=[x.enc()for x in ips]))
    tests=[];accepted=[]
    for pair in combinations(range(4),2):
        kernel=A.kernel([A.mv(Q,N[i])for i in white+list(pair)]);assert len(kernel)==1;signs=[]
        for sign in(-1,1):
            v=[sign*x for x in kernel[0]];ip=[A.pair(Q,v,x)for x in N]
            if all(x.sign()<=0 for x in ip):
                assert A.pair(Q,v,v).sign()<0 and A.pair(Q,v,H).sign()<0
                signs.append(sign);accepted.append(list(pair))
        tests.append(dict(black=list(pair),feasible_orientations=signs))
    assert sorted(accepted)==sorted(v['black']for v in vertices)
    return dict(key=key,row=row,gram=encoded,basis_indices=basis,basis_form=A.enc(Q),root_coordinates=A.enc(N),
        physical_inertia=[n-2,1,1],conjugate_inertia=[n-1,0,1],integral_cartan_entries=n*n,
        strict_center=[x.enc()for x in H],center_norm=hq.enc(),center_products=[x.enc()for x in hp],
        vertices=vertices,complete_face_cone=tests,spectator_rank=15-n)

def run():
    data={t:read(t)for t in SOURCES};pre=read(PRE)
    primal=next(c for c in pre['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    weights={i:F(x)for i,x in primal['weights']}
    lookup={(tuple(v['forest']),dump(v['marked_forest_code'])):v['index']for v in pre['variables']}
    rows={r['class_index']:r['row']for r in data[SOURCES[4]]['records']}
    inputs=[]
    for rec in data[SOURCES[0]]['records']:inputs.append((str(rec['class_index']),rec['row'],rec['scaled_gram'],rec['basis_indices']))
    rec=data[SOURCES[1]]['geometry'];inputs.append(('4202',rows[4202],rec['scaled_gram'],rec['basis_indices']))
    rec=data[SOURCES[2]]['geometry'];inputs.append(('7661',rows[7661],rec['gram'],rec['basis_indices']))
    rec=data[SOURCES[3]];inputs.append(('B6_first39',rec['row'],rec['scaled_gram'],rec['basis_indices']))
    cores=[core_geometry(*x)for x in inputs]
    models=[(name,g)for name,g in C.models()if len(g)<=13];finite={C.code(g):name for name,g in models};graphs=dict(models)
    components=[finite_component(name,g)for name,g in models if len(g)<=7];component_lookup={r['name']:r for r in components}
    ordered=[(name,g)for name,g in models if len(g)<=7]
    def spectators(rank,first=0):
        if not rank:yield ();return
        for j in range(first,len(ordered)):
            name,g=ordered[j]
            if len(g)<=rank:
                for tail in spectators(rank-len(g),j):yield tuple(sorted((name,)+tail))
    contexts=[];counts=Counter();all_forests=set();known_failing=[]
    for ci,core in enumerate(cores):
        row=core['row'];n=row['nodes'];base=C.graph(range(n),dict(zip(combinations(range(n),2),row['labels'])))
        hq=K(*map(F,core['center_norm']));core_count=0
        for names in sorted(spectators(core['spectator_rank'])):
            g={i:dict(r)for i,r in base.items()};extra=[];p=Z
            for name in names:
                h=graphs[name];off=len(g);extra+=list(range(off,off+len(h)))
                for i,r in h.items():g[off+i]={off+j:m for j,m in r.items()}
                p+=K(*map(F,component_lookup[name]['positive_center_cost']))
            assert len(g)==15 and p.sign()>0
            epsilon=F(1)
            while(hq+p*epsilon*epsilon).sign()>=0:epsilon/=2
            white=row['white']+extra;assert len(white)==11
            corners=[]
            for pair in[(0,1),(1,2),(2,3),(0,3)]:
                ids=set(white+list(pair));h={i:{j:m for j,m in g[i].items()if j in ids}for i in ids};pieces=C.pieces(h)
                forest=tuple(sorted(finite[C.code(t)]for t in pieces));code=dump(sorted(C.code(t,pair)for t in pieces));all_forests.add(forest)
                index=lookup.get((forest,code))
                corners.append(dict(forest=list(forest),marked_forest_code=json.loads(code),black=list(pair),
                    predecessor_variable=index,predecessor_payment=None if index is None else str(weights.get(index,F(0)))))
            known=all(c['predecessor_variable']is not None for c in corners)
            payment=sum((F(c['predecessor_payment'])for c in corners),F(0))if known else None
            index=len(contexts)
            contexts.append(dict(index=index,core=ci,core_key=core['key'],spectator_forest=list(names),
                center_epsilon=str(epsilon),extended_center_norm=(hq+p*epsilon*epsilon).enc(),
                physical_inertia=[13,1,1],conjugate_inertia=[14,0,1],integral_cartan_entries=225,
                corners=corners,predecessor_known=known,predecessor_payment=None if payment is None else str(payment)))
            counts['contexts']+=1;counts['known_to_predecessor']+=known;counts['requires_new_vertex_types']+=not known
            counts['known_unpaid']+=known and payment<1
            if known and payment<1:known_failing.append(index)
            core_count+=1
        counts['core_'+core['key']]=core_count
    counts['distinct_corner_forests']=len(all_forests)
    counts['new_corner_forests']=len(all_forests-{tuple(v['forest'])for v in pre['vertices']})
    paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'scripts/{C.TAG}.py']+[f'results/{t}.json'for t in SOURCES+[PRE]]
    return dict(status='PASS_COMPLETE_FINITE_SPECTATOR_FAMILIES_FOR_FIVE_LOCAL_QUADS',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},counts=dict(counts),
        components=components,cores=cores,contexts=contexts,known_failing_context_indices=known_failing,
        complete_context_ranks=sorted({c['spectator_rank']for c in cores}),
        proof='At a finite tree root set s0=1, s_child=s_parent*t_squared, diagonal2s and parent-child product−s_child. Every scale is totally positive and the Cartan coefficients are−t_squared and−1. Exact two-place positive definiteness and every bond are checked per component. Each spectator forest is their orthogonal sum. The core Lorentz form and conjugate positive form extend orthogonally; the old face corner rays acquire zero spectator coordinates. Their common face cone is exactly the old four-ray cone because all spectator walls are active. A strictly feasible center is obtained by adding−epsilon*S_inverse*1 in the spectator subspace; the saved epsilon makes its norm strictly negative. All225 Cartan entries consist of checked core/component blocks and zeros.',
        completeness='Every finite component type allowed by squared labels2,3,4,5,6,10 and every component multiset of the exact remaining rank5 or7 is included for each of the five fixed cores. This is complete for these orthogonal spectator families only.',
        scope='Exact local arithmetic compact-quadrilateral families and partial-rule separation, with absent whole vertex types identified. No compact finite-volume ambient extension, complete rank13 face domain, or dimension bound is claimed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
