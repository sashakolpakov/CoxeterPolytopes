#!/usr/bin/env python3
"""Independent full height30 replay and exact silver chamber symmetry audit.

Only the old standalone Fraction verifier supplies field/linear arithmetic.
Positive-height blocks use a direct fourfold integer-covector product, with no
producer recursion, Schur pruning, or new-primary arithmetic imports.
"""
import argparse,json,sys
from collections import Counter
from fractions import Fraction as F
from functools import cmp_to_key,lru_cache
from hashlib import sha256
from itertools import combinations,groupby,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import fresh_silver_f4_independent_verify as I
BASE='fresh_silver_full_chamber_20261015';TAG=BASE+'_review'
OLD='results/fresh_silver_f4_h9_unit_h30_20261013.json'
Z,O,S=I.Z,I.O,I.P

def imul(x,y):return (x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def integer_pair(x):assert I.integral(x);return tuple(map(int,x))

@lru_cache(None)
def direct_blocks(target,q):
    alphabet=tuple(x for x in I.bounded((I.scale(target,2),I.scale(target,2)))if I.sign(x)>=0 and I.integral(I.div(I.scale(x,2),q)))
    alpha=list(map(integer_pair,alphabet));J=[[integer_pair(x)for x in row]for row in I.FI];n=len(alpha)
    # All16^4 possible tuples at the largest level are inspected. Tables
    # accelerate only exact integer field operations, never prune a tuple.
    terms={}
    for i,j in combinations(range(4),2):terms[i,j]=[[tuple(2*z for z in imul(J[i][j],imul(a,b)))for b in alpha]for a in alpha]
    diagonal=[[imul(J[i][i],imul(a,a))for a in alpha]for i in range(4)]
    out=[]
    for ids in product(range(n),repeat=4):
        a=sum(diagonal[i][ids[i]][0]for i in range(4));b=sum(diagonal[i][ids[i]][1]for i in range(4))
        for i,j in combinations(range(4),2):x,y=terms[i,j][ids[i]][ids[j]];a+=x;b+=y
        energy=(F(a),F(b));diff=I.sub(target,energy)
        if I.sign(diff)<0 or I.sign(I.conj(diff))<0:continue
        d=tuple(alphabet[k]for k in ids);v=I.mv(I.FI,d)
        assert all(map(I.integral,v))and I.pair(v,I.F4,v)==energy
        out.append((energy,v))
    return tuple(out),dict(alphabet=n,tuples=n**4,admissible_block_vectors=len(out))

def field_compare(x,y):
    for a,b in zip(x,y):
        z=I.sign(I.sub(a,b))
        if z:return z
    return 0

def full_prefix(old):
    G=I.host(I.UNIT,2,True);pool=(O,I.ETA,I.scale(O,2));Hcap=I.parse(old['height']);assert Hcap==I.scale(O,30)
    assert I.matrix(old['gram'])==G and tuple(map(I.parse,old['norm_pool']))==pool
    assert I.determinant(I.F4)==O and I.determinant(G)==I.neg(I.UNIT)and I.integral(I.inv(I.UNIT))
    assert I.mul(I.ETA,I.ETA)==I.scale(I.mul(I.UNIT,I.UNIT),2)and I.ETA==I.mul(S,I.UNIT)
    f4=I.finite_f4();zero=I.zero_census(pool,2,True)
    assert zero=={(I.parse(z['norm']),tuple(map(I.parse,z['vector'])))for z in old['zero_height_dominant_census']}and len(zero)==6
    controls=[]
    for q in pool:
        caps=(I.div(I.mul(Hcap,q),I.UNIT),I.neg(I.div(q,I.UNIT)))
        for x in I.bounded(caps):
            if I.sign(x)>0 and I.integral(I.div(I.scale(I.mul(I.UNIT,x),2),q)):
                controls.append((I.div(I.mul(I.UNIT,I.mul(x,x)),q),q,x))
    controls.sort(key=cmp_to_key(field_compare));assert len(controls)==17
    old_controls=[(I.parse(r['height']),I.parse(r['norm']),I.parse(r['control']))for r in old['levels']]
    assert controls==old_controls
    walls=[tuple(I.neg(O)if i==j else Z for i in range(10))for j in range(1,10)];keys={I.oriented_key(v)for v in walls};rows=[];tuplecount=0
    for h,packet in groupby(controls,key=lambda z:z[0]):
        add=[]
        for _,q,x in packet:
            target=I.add(q,I.mul(I.UNIT,I.mul(x,x)));blocks,diag=direct_blocks(target,q);tuplecount+=diag['tuples']
            scalars=[z for z in I.bounded((target,target))if I.sign(z)>=0 and I.integral(I.div(I.scale(z,2),q))]
            bysquare={I.mul(z,z):z for z in scalars};assert len(bysquare)==len(scalars)
            candidates=[]
            for first,second in product(blocks,repeat=2):
                remaining=I.sub(I.sub(target,first[0]),second[0])
                if remaining not in bysquare:continue
                r=(x,*first[1],*second[1],bysquare[remaining]);cov=I.mv(G,r)
                assert I.dot(r,cov)==q and all(map(I.integral,r))and all(I.integral(I.div(I.scale(c,2),q))for c in cov)
                candidates.append((r,cov))
            kept=duplicates=0
            for r,cov in candidates:
                if any(I.sign(I.dot(w,cov))>0 for w in walls):continue
                key=I.oriented_key(r)
                if key in keys:duplicates+=1;continue
                keys.add(key);add.append(r);kept+=1
            oldrow=old['levels'][len(rows)]
            assert len(candidates)==oldrow['candidates']and kept==oldrow['roots']and duplicates==oldrow['parallel_duplicates']
            rows.append(dict(height=I.enc(h),norm=I.enc(q),control=I.enc(x),**diag,candidates=len(candidates),selected=kept,parallel_duplicates=duplicates))
        assert all(I.sign(I.pair(r,G,t))<=0 for r,t in combinations(add,2));walls.extend(add)
        saved=[r for r in I.matrix(old['all_walls'])if I.sign(I.sub(I.div(I.mul(I.UNIT,I.mul(r[0],r[0])),I.pair(r,G,r)),h))<=0]
        assert {I.oriented_key(r)for r in walls}=={I.oriented_key(r)for r in saved}
    assert len(walls)==92 and {I.oriented_key(r)for r in walls}=={I.oriented_key(r)for r in I.matrix(old['all_walls'])}
    return dict(zero_height_roots=len(zero),finite_F4=f4,control_levels=len(controls),direct_covector_tuples=tuplecount,levels=rows,walls=92)

def positive(M):
    A=[list(row)for row in M];pivots=[]
    while A:
        p=A[0][0];assert I.sign(p)>0;pivots.append(p)
        A=[[I.sub(A[i][j],I.div(I.mul(A[i][0],A[0][j]),p))for j in range(1,len(A))]for i in range(1,len(A))]
    return pivots

def witness(old,primary):
    G=I.matrix(old['gram']);N=I.matrix(old['all_walls']);H=(O,)+(Z,)*9;cov=[I.mv(G,n)for n in N]
    norm=[I.dot(n,c)for n,c in zip(N,cov)];center=tuple(map(I.parse,old['strict_interior']))
    assert I.sign(I.pair(center,G,center))<0 and I.sign(I.pair(center,G,H))<0
    assert all(I.sign(I.dot(center,c))<0 for c in cov)
    positive(tuple(tuple(I.conj(x)for x in row)for row in G));positive(tuple(row[1:]for row in G[1:]))
    angles=Counter();allowed={Z:2,O:3,I.scale(O,2):4,I.scale(O,3):6,I.ETA:8}
    for i,j in combinations(range(92),2):
        p=I.dot(N[i],cov[j]);assert I.sign(p)<=0
        q=I.div(I.scale(I.mul(p,p),4),I.mul(norm[i],norm[j]))
        if I.sign(I.sub(q,I.scale(O,4)))<0:assert q in allowed;angles['m'+str(allowed[q])]+=1
        else:assert I.sign(I.sub(q,I.scale(O,4)))>0;angles['dotted']+=1
    assert dict(angles)==primary['finite_angle_counts']
    for i,n in enumerate(N):
        coefficient=I.div(I.dot(center,cov[i]),norm[i]);v=tuple(I.sub(x,I.mul(coefficient,y))for x,y in zip(center,n))
        assert v==tuple(map(I.parse,primary['actual_facet_points'][i]))
        assert I.sign(I.pair(v,G,v))<0 and I.sign(I.pair(v,G,H))<0
        assert all(I.sign(I.dot(v,c))==(0 if i==j else -1)for j,c in enumerate(cov))
    corners=[]
    for name in('source_corner','target_corner'):
        raw=primary[name];ids=raw['ordered_indices'];B=tuple(tuple(I.dot(N[i],cov[j])for j in ids)for i in ids);pivots=positive(B);J=I.inverse(B)
        assert all(I.sign(z)>=0 for row in J for z in row)
        v=I.kernel_one([cov[i]for i in ids]);assert v
        if I.sign(I.pair(v,G,H))>0:v=tuple(map(I.neg,v))
        saved=tuple(map(I.parse,raw['point']));factor=I.div(next(a for a in saved if a!=Z),next(a for a in v if a!=Z));assert I.sign(factor)>0
        assert saved==tuple(I.mul(factor,a)for a in v)
        vv=I.pair(v,G,v);hv=I.pair(H,G,v);assert I.sign(vv)<0 and I.sign(hv)<0
        values=[I.dot(v,c)for c in cov];assert [i for i,a in enumerate(values)if a==Z]==sorted(ids)and all(I.sign(a)<=0 for a in values)
        a=I.div(hv,vv);b=I.mv(J,[I.dot(H,cov[i])for i in ids]);assert I.sign(a)>0 and all(I.sign(z)<=0 for z in b)
        reconstructed=tuple(I.add(I.mul(a,v[k]),I.total(I.mul(b[j],N[ids[j]][k])for j in range(9)))for k in range(10));assert reconstructed==H
        assert b==tuple(map(I.parse,raw['basepoint_decomposition']['normal_coefficients']))
        corners.append(dict(ordered_indices=ids,point=list(map(I.enc,v)),norm=I.enc(vv),base_vertex_coefficient=I.enc(a),normal_coefficients=list(map(I.enc,b)),pivots=list(map(I.enc,pivots))))
    # Literal connected finite-diagram type, independent of producer naming.
    ids=primary['source_corner']['ordered_indices'];adj={i:[]for i in ids};labels={}
    for i,j in combinations(ids,2):
        p=I.dot(N[i],cov[j])
        if p!=Z:
            q=I.div(I.scale(I.mul(p,p),4),I.mul(norm[i],norm[j]));assert q in allowed
            adj[i].append(j);adj[j].append(i);labels[i,j]=allowed[q]
    remaining=set(ids);types=[]
    while remaining:
        todo=[min(remaining)];component=set()
        while todo:
            i=todo.pop()
            if i in component:continue
            component.add(i);todo.extend(adj[i])
        remaining-=component
        if len(component)==1:types.append('A1')
        elif len(component)==2:
            i,j=sorted(component);assert labels[i,j]==4;types.append('I2(4)')
        else:
            assert len(component)==4 and sorted(len(adj[i])for i in component)==[1,1,2,2]
            edges=[(i,j,m)for(i,j),m in labels.items()if i in component and j in component]
            assert sorted(m for i,j,m in edges)==[3,3,4]
            i,j,_=next(e for e in edges if e[2]==4);assert min(len(adj[i]),len(adj[j]))==1;types.append('B4')
    assert sorted(types)==['A1','B4','I2(4)','I2(4)']
    T=I.matrix(primary['matrix']);Tinverse=I.inverse(T);assert Tinverse==I.matrix(primary['inverse_matrix'])
    assert all(I.integral(x)for row in T+Tinverse for x in row)
    assert I.mm(T,Tinverse)==I.mm(Tinverse,T)==I.eye(10)and I.mm(I.tr(T),I.mm(G,T))==G
    src,dst=(primary[k]['ordered_indices']for k in('source_corner','target_corner'))
    assert all(I.mv(T,N[a])==N[b]for a,b in zip(src,dst))
    v,w=(tuple(map(I.parse,primary[k]['point']))for k in('source_corner','target_corner'));scale=I.parse(primary['positive_corner_scale'])
    assert I.sign(scale)>0 and I.mv(T,v)==tuple(I.mul(scale,x)for x in w)
    tr=I.total(T[i][i]for i in range(10));assert tr==(F(8),F(3))and I.sign(I.sub(tr,I.scale(O,10)))>0
    det=I.determinant(T);assert det==I.parse(primary['determinant'])and det in(O,I.neg(O))
    return dict(pair_counts=dict(angles),actual_facets=92,corners=corners,corner_type=sorted(types),integral_matrix_and_inverse=True,outward_images=9,future_corner_scale=I.enc(scale),determinant=I.enc(det),trace=I.enc(tr),trace_minus_dimension=I.enc(I.sub(tr,I.scale(O,10))))

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    old=json.loads((ROOT/OLD).read_text());primary=json.loads((ROOT/f'results/{BASE}.json').read_text());enum=full_prefix(old);proof=witness(old,primary)
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',OLD,'scripts/fresh_silver_f4_h9_20261013.py','scripts/fresh_silver_f4_independent_verify.py','results/fresh_silver_f4_independent_report_20261013.md']
    visited=set()
    def check(p):
        if p in visited:return
        visited.add(p)
        if p.endswith('.json'):
            for q,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():assert sha(q)==h,(p,q);check(q)
    for p in paths:check(p)
    return dict(status='PASS_INDEPENDENT_COMPLETE_SILVER_PREFIX_AND_INFINITE_ORDER_CHAMBER_SYMMETRY',sources={p:sha(p)for p in paths},recursive_manifest_files=len(visited),
       independent_complete_height30_prefix=enum,independent_witness=proof,
       proof_scope='Primitive norm classes1,eta,2 and complete zero sector audited. Independent direct covector enumeration repeats every bounded positive-height control. Complete outward corners persist for all later simple walls; integral future isometry maps their local sectors and hence the full arrangement chamber to itself. Infinite-order discrete stabilizer forces infinite volume. Literal lattice/full integral reflection group only; no other commensurable-lattice or dimension exclusion.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert data==json.loads(out.read_text());print('PASS saved independent silver review')
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['independent_complete_height30_prefix']['direct_covector_tuples'])
