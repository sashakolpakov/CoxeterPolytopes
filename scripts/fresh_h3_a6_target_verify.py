#!/usr/bin/env python3
"""Independent Fraction/integer replay of targeted H3+A6 cuts and obstructions."""
import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
from scripts.fresh_h3_blocks_h9_independent_verify import (
    H3,A6,Z,O,P,parse,enc,add,neg,sub,scale,mul,div,conj,sign,integral,total,dot,mv,tr,mm,pair,matrix,
    inverse,diagonal_block,read,ROOT,bounded_integers,valid_coxeter_value,escape,
)

def imul(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def isum(xs):
    xs=tuple(xs);return sum(x[0] for x in xs),sum(x[1] for x in xs)
def idot(x,y):return isum(imul(a,b) for a,b in zip(x,y))
def integer(x):
    assert integral(x);return int(x[0]),int(x[1])
def angle(p,qq):
    if sign(p)>0:return False
    return valid_coxeter_value(scale(imul(p,p),F(4,qq)))

def setup(name):
    d=read('results/fresh_h3_blocks_h9_a6_'+name+'_h30_20261013.json');alpha=parse(d['alpha'])
    B=diagonal_block(((neg(alpha),),),H3,A6);assert B==matrix(d['gram'])
    hi,ai=inverse(H3),inverse(A6)
    e=(O,*(Z for _ in range(9)))
    h=(Z,*(scale(hi[i][0],2) for i in range(3)),*(Z for _ in range(6)))
    k=(Z,*(hi[i][1] for i in range(3)),*(Z for _ in range(6)))
    a=(Z,Z,Z,Z,*(add(ai[i][2],ai[i][3]) for i in range(6)))
    b=(Z,Z,Z,Z,*(scale(ai[i][2],7) for i in range(6)))
    assert all(integral(x) for r in (e,h,k,a,b) for x in r)
    # Integral inverses of coordinate minors prove that these bases give
    # the full integral intersection with their respective field spans.
    for columns,indices in (((e,h,a,b),(0,3,4,9)),((e,h,k,a),(0,1,3,4))):
        minor=tuple(tuple(col[i] for col in columns) for i in indices)
        assert all(integral(x) for row in inverse(minor) for x in row)
    return d,alpha,B,(e,h,k,a,b)

def verify_host(name):
    base,alpha,B,basis=setup(name);d=read('results/fresh_h3_a6_target_'+name+'_neighbor_20261013.json')
    c=read('results/fresh_h3_a6_target_'+name+'_cones_20261013.json');old=matrix(base['all_walls']);inside=tuple(map(parse,base['strict_interior']))
    e,h,k,a,b=basis;newbasis=(e,h,k,a)
    assert matrix(d['basis'])==newbasis and matrix(d['gram'])==B
    assert d['bounds']=={'s_primary_max':100,'t_primary_max':100,'u_primary_max':100}
    expected_pool=[1,2,4] if name=='phi' else [2,4];assert d['norm_pool']==expected_pool
    # Covectors (2s,t,0) and (0,0,u,u,0,0), together with the
    # unit time coefficient at7, force q divisible by7 to be imprimitive.
    assert pair(h,B,h)==scale(mul(P,P),6) and pair(k,B,k)==scale(mul(P,P),2)
    assert pair(h,B,k)==scale(mul(mul(P,P),P),2) and pair(a,B,a)==scale(O,6)
    target=tuple(map(parse,d['target']));targetnorm=pair(target,B,target)
    assert targetnorm==parse(d['target_norm']) and sign(targetnorm)>0
    assert all(sign(pair(target,B,r))<=0 for r in old)
    assert sign(pair(inside,B,inside))<0 and sign(inside[0])>0
    roots=[];qs=[];parameters=[];cut_counts={q:0 for q in expected_pool};counts={q:0 for q in expected_pool}
    Bint=tuple(tuple(map(integer,r)) for r in B)
    oldint=tuple(tuple(map(integer,r)) for r in old)
    oldcov=tuple(tuple(idot(row,r) for row in Bint) for r in oldint)
    oldqs=[]
    for r,n in zip(oldint,oldcov):
        q=idot(r,n);assert q[1]==0 and q[0]>0;oldqs.append(q[0])
    for record in d['cuts']:
        r=tuple(map(parse,record['root']));p=tuple(map(parse,record['parameters']));q=record['norm']
        x,s,t,u=p;assert q in expected_pool and len(p)==4
        assert all(integral(z) for z in (*r,*p)) and sign(x)>0 and sign(s)>=0 and sign(t)>0 and sign(u)>=0
        assert all(sign(sub(scale(O,100),z))>=0 for z in (s,t,u))
        assert r==tuple(total(mul(z,col[j]) for z,col in zip(p,newbasis)) for j in range(10))
        ri=tuple(map(integer,r));ni=tuple(idot(row,ri) for row in Bint)
        assert idot(ri,ni)==(q,0)
        assert all((2*v)%q==0 for z in ni for v in z)
        if q==4:assert any(v%2 for z in ri for v in z)
        assert sign(dot(inside,ni))<0
        assert all(angle(idot(ri,n),q*oq) for n,oq in zip(oldcov,oldqs))
        height=div(mul(alpha,mul(x,x)),scale(O,q));assert height==parse(record['height'])
        targetpair=dot(target,ni);assert targetpair==parse(record['target_pairing'])
        cut_counts[q]+=int(sign(targetpair)>0);counts[q]+=1
        roots.append(r);qs.append(q);parameters.append(p)
    assert len(set(roots))==len(roots)==c['dictionary_size'] and c['source_walls']==len(old)
    # Independent rectangular embedding census of the parameter alphabets.
    for row in d['censuses']:
        q=row['norm'];ss=tuple(x for x in bounded_integers(scale(O,10000),scale(O,F(q,2))) if sign(x)>=0)
        ts=tuple(x for x in bounded_integers(scale(O,10000),scale(O,F(3*q,2))) if sign(x)>0)
        us=tuple(x for x in bounded_integers(scale(O,10000),scale(O,F(q,6))) if sign(x)>=0)
        if q==4:
            ts=tuple(x for x in ts if integral(scale(x,F(1,2))));us=tuple(x for x in us if integral(scale(x,F(1,2))))
        assert (len(ss),len(ts),len(us))==(row['s_alphabet'],row['t_positive_alphabet'],row['u_alphabet'])
        assert row['triples']==len(ss)*len(ts)*len(us) and counts[q]==row['coxeter_compatible'] and cut_counts[q]==row['cut_old_target']
    # Build the restricted Gram independently from the full host matrix.
    small=tuple(tuple(integer(pair(v,B,w)) for w in newbasis) for v in newbasis)
    params=tuple(tuple(map(integer,p)) for p in parameters)
    cov=tuple(tuple(idot(row,p) for row in small) for p in params)
    conflicts=set()
    for i,p in enumerate(params):
        for j in range(i):
            if not angle(idot(p,cov[j]),qs[i]*qs[j]):conflicts.add((j,i))
    assert conflicts=={tuple(edge) for edge in c['conflict_edges']}
    assert len(conflicts)==len(c['conflict_edges'])
    neighbors=[set() for _ in roots]
    for i,j in conflicts:neighbors[i].add(j);neighbors[j].add(i)
    escaped=[]
    for case in c['cases']:
        ids=case['selection'];assert len(ids)==len(set(ids)) and all(0<=i<len(roots) for i in ids)
        chosen=set(ids);assert all(not neighbors[i]&chosen for i in ids)
        assert all(neighbors[i]&chosen for i in range(len(roots)) if i not in chosen)
        walls=(*old,*(roots[i] for i in ids));assert len(walls)==case['wall_count']
        v=escape(case['hyperideal_witness'],walls,B);escaped.append({'walls':len(walls),'norm':enc(pair(v,B,v))})
    common=escape(c['common_dictionary_escape'],(*old,*roots),B)
    return {'alpha':name,'accepted_dictionary_roots':len(roots),'old_wall_count':len(old),'cuts_initial_target':sum(cut_counts.values()),
        'conflict_edges':len(conflicts),'all_new_old_and_new_new_angles_checked':True,'maximal_compatible_cases':escaped,
        'common_dictionary_escape_norm':enc(pair(common,B,common))}

def obstruction_payload():
    out=[]
    for name in ('phi','2phi'):
        base,alpha,B,basis=setup(name);e,h,k,a,b=basis;old=matrix(base['all_walls'])
        closed=(e,h,a,b);closedgram=tuple(tuple(pair(v,B,w) for w in closed) for v in closed)
        assert all(integral(scale(c,F(1,3))) for row in closedgram[1:] for c in row[1:])
        assert closedgram[1][1]==scale(mul(P,P),6) and closedgram[2][2]==scale(O,6)
        assert closedgram[2][3]==scale(O,21) and closedgram[3][3]==scale(O,84)
        # Enumerate O/(3)=F9. The polynomial z²-z-1 is irreducible mod3.
        assert all((j*j-j-1)%3 for j in range(3))
        pool=[1,2,4,7,14,28] if name=='phi' else [2,4,14,28]
        residue_solutions=[]
        for q in pool:
            for x in product(range(3),repeat=2):
                lhs=neg(imul(integer(alpha),imul(x,x)))
                if lhs[0]%3==q%3 and lhs[1]%3==0:residue_solutions.append((q,x))
        assert not residue_solutions
        d=read('results/fresh_h3_a6_target_'+name+'_neighbor_20261013.json');target=tuple(map(parse,d['target']))
        geometric=(add(target[0],scale(O,F(1,100))),*target[1:])
        assert sign(pair(geometric,B,geometric))>0 and sign(pair(geometric,B,target))>0
        assert all(sign(pair(geometric,B,r))<=0 for r in old) and all(sign(pair(geometric,B,r))<0 for r in old[9:])
        c=read('results/fresh_h3_a6_target_'+name+'_cones_20261013.json')
        if name=='phi':
            witness=c['common_dictionary_escape'];ids=(0,1,21,33)
            coef=(scale(P,F(1,2)),O,sub(P,scale(O,F(3,2))),sub(scale(O,F(5,2)),scale(P,F(3,2))))
            source='common_dictionary_escape'
        else:
            witness=c['cases'][2]['hyperideal_witness'];ids=(5,12);coef=(sub(scale(O,4),scale(P,2)),sub(P,O));source='cases[2]'
        v=tuple(map(parse,witness['vector']));qv=pair(v,B,v)
        assert qv==parse(witness['norm']) and sign(qv)>0 and sign(v[0])>0
        assert all(sign(pair(v,B,r))<=0 for r in old)
        newbasis=(e,h,k,a);project=lambda r:tuple(pair(z,B,r) for z in newbasis)
        projected=project(v)
        assert all(sign(z)>0 for z in coef)
        assert projected==tuple(total(mul(z,project(old[i])[j]) for z,i in zip(coef,ids)) for j in range(4))
        # The identity covers the whole finite dictionary as a corollary.
        assert all(sign(pair(v,B,tuple(map(parse,r['root']))))<=0 for r in d['cuts'])
        out.append({'alpha':name,'original_source':'results/fresh_h3_blocks_h9_a6_'+name+'_h30_20261013.json',
            'root_free_middle_slice_basis':[[enc(c) for c in r] for r in closed],
            'root_free_slice_integral_inverse_minor_coordinates':[0,3,4,9],
            'root_free_middle_slice_gram':[[enc(c) for c in r] for r in closedgram],
            'primitive_candidate_norm_pool':pool,'mod3_root_solutions':residue_solutions,
            'geometric_cutter':[enc(z) for z in geometric],'geometric_cutter_norm':enc(pair(geometric,B,geometric)),
            'geometric_cutter_target_pairing':enc(pair(geometric,B,target)),
            'neighbor_slice_basis':[[enc(c) for c in r] for r in newbasis],
            'neighbor_slice_integral_inverse_minor_coordinates':[0,1,3,4],
            'escape_source':'results/fresh_h3_a6_target_'+name+'_cones_20261013.json','escape_source_member':source,
            'escape':[enc(z) for z in v],'escape_norm':enc(qv),'projection':[enc(z) for z in projected],
            'original_wall_indices':ids,'strictly_positive_coefficients':[enc(z) for z in coef]})
    return {'status':'EXACT_MOD3_AND_UNBOUNDED_NEIGHBOR_SLICE_OBSTRUCTIONS','hosts':out,
        'scope':'No primitive reflection normal in span(time,omega1(H3),omega3(A6),omega4(A6)). Separately, no full-dimensional finite-volume completion using only old-compatible real normals in span(time,omega1(H3),omega2(H3),omega3(A6)+omega4(A6)). Neither statement proves full-host nonreflectivity.'}

def verify():
    out={'status':'INDEPENDENT_H3_A6_TARGET_REPLAY_PASS','hosts':[verify_host(n) for n in ('phi','2phi')],
        'scope':'Independent full accepted-root, angle, conflict-graph, cone and escape replay; independent parameter-alphabet census. Positive-height candidate enumeration is source-reviewed, not independently repeated in full.'}
    obstruction=obstruction_payload()
    (ROOT/'results/fresh_h3_a6_target_obstructions_20261013.json').write_text(json.dumps(obstruction,indent=2,sort_keys=True)+'\n')
    (ROOT/'results/fresh_h3_a6_target_independent_audit_20261013.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    return out

if __name__=='__main__':
    print(json.dumps(verify(),sort_keys=True))
