#!/usr/bin/env python3
"""Independent Fraction-only replay of released-H8 H9 cap certificates.

No project imports, NumPy, SciPy or symbolic-algebra package. The cap
proposal searches are not rerun; every accepted normal, field identity,
conflict and strict escape used in the finite-bank obstruction is replayed.
"""
from fractions import Fraction as F
from itertools import combinations, product
from functools import lru_cache
from pathlib import Path
import argparse
import json

Z=(F(0),F(0));O=(F(1),F(0));P=(F(0),F(1))
def parse(x):return F(x['a']),F(x['b'])
def enc(x):return {'a':str(x[0]),'b':str(x[1])}
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def sub(x,y):return add(x,neg(y))
def scale(x,k):return x[0]*k,x[1]*k
def mul(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def inv(x):
    n=x[0]*x[0]+x[0]*x[1]-x[1]*x[1]
    return (x[0]+x[1])/n,-x[1]/n
def div(x,y):return mul(x,inv(y))
def conj(x):return x[0]+x[1],-x[1]
def sign(x):
    a=2*x[0]+x[1];b=x[1]
    if not b:return (a>0)-(a<0)
    if not a or (a>0)==(b>0):return (b>0)-(b<0)
    return ((a>0)-(a<0))*((a*a>5*b*b)-(a*a<5*b*b))
def integral(x):return all(c.denominator==1 for c in x)
def total(values):
    out=Z
    for x in values:out=add(out,x)
    return out
def dot(x,y):return total(mul(a,b) for a,b in zip(x,y,strict=True))
def mv(a,v):return tuple(dot(row,v) for row in a)
def tr(a):return tuple(zip(*a))
def mm(a,b):return tuple(tuple(dot(row,col) for col in tr(b)) for row in a)
def pair(x,a,y):return dot(x,mv(a,y))
def matrix(x):return tuple(tuple(map(parse,row)) for row in x)
def eye(n):return tuple(tuple(O if i==j else Z for j in range(n)) for i in range(n))

def rref(a):
    a=[list(row) for row in a];cols=len(a[0]);pivots=[]
    for j in range(cols):
        k=len(pivots)
        pivot=next((i for i in range(k,len(a)) if a[i][j]!=Z),None)
        if pivot is None:continue
        a[k],a[pivot]=a[pivot],a[k]
        q=a[k][j];a[k]=[div(x,q) for x in a[k]]
        for i in range(len(a)):
            if i==k:continue
            q=a[i][j];a[i]=[sub(x,mul(q,y)) for x,y in zip(a[i],a[k])]
        pivots.append(j)
        if len(pivots)==len(a):break
    return a,pivots

def inverse(a):
    n=len(a);rows,pivots=rref(tuple(tuple(row)+ident for row,ident in zip(a,eye(n))))
    assert pivots==list(range(n))
    return tuple(tuple(row[n:]) for row in rows)

def determinant(a):
    a=[list(row) for row in a];out=O
    for j in range(len(a)):
        k=next((i for i in range(j,len(a)) if a[i][j]!=Z),None)
        if k is None:return Z
        if k!=j:a[k],a[j]=a[j],a[k];out=neg(out)
        pivot=a[j][j];out=mul(out,pivot)
        for i in range(j+1,len(a)):
            q=div(a[i][j],pivot)
            a[i]=[sub(x,mul(q,y)) for x,y in zip(a[i],a[j])]
    return out

def kernel_one(a):
    rows,pivots=rref(a);free=[i for i in range(len(a[0])) if i not in pivots]
    if len(free)!=1:return None
    out=[Z]*len(a[0]);out[free[0]]=O
    for row,pivot in zip(rows,pivots):out[pivot]=neg(row[free[0]])
    return tuple(out)

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/enumerations/fresh_h8_release_h9_caps_expanded_20261013.json'
COMMON=ROOT/'data/enumerations/fresh_h8_release_h9_volume_expanded_20261013.json'
FREE=ROOT/'data/enumerations/fresh_h8_release_h9_free_caps_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h8_release_h9_independent_verify_20261013.json'
EZ=(Z,Z)
def ea(x,y):return add(x[0],y[0]),add(x[1],y[1])
def en(x):return neg(x[0]),neg(x[1])
def es(x,y):return ea(x,en(y))
def em(x,y,d):return add(mul(x[0],y[0]),mul(d,mul(x[1],y[1]))),add(mul(x[0],y[1]),mul(x[1],y[0]))
def sg(x,d):
    a,b=x;sa,sb=sign(a),sign(b)
    if not sb:return sa
    if not sa or sa==sb:return sb
    return sa*sign(sub(mul(a,a),mul(d,mul(b,b))))
def ep(cs,x,d):
    p=EZ
    for c in cs:p=ea(em(p,x,d),(scale(O,c),Z))
    return p
def pe(x):return parse(x['rational']),parse(x['radical'])
def ve(x):return tuple(map(parse,x))
def vv(x):return [enc(t)for t in x]

def eligible(p,d):
    if sg(p,d)>0:return False
    q=em(p,p,d)
    if sg(es(q,(scale(O,4),Z)),d)>=0:return True
    if any(sg(es(q,(scale(O,k),Z)),d)==0 for k in (0,1,2,3)):return True
    t=es(q,(scale(O,2),Z))
    return any(sg(ep(cs,t,d),d)==0 for cs in ((1,1,-1),(1,-1,-1),(1,0,-2),(1,0,-3),
        (1,-1,-4,4,1),(1,0,-4,0,2),(1,0,-5,0,5),(1,0,-4,0,1),(1,1,-4,-4,1)))

def bsg(x,d,e):
    a=x[:2];b=x[2:];sa,sb=sg(a,d),sg(b,d)
    if not sb:return sa
    if not sa or sa==sb:return sb
    q=em(b,b,d);return sa*sg(es(em(a,a,d),(mul(e,q[0]),mul(e,q[1]))),d)
def bm(x,y,d,e):
    z=[Z]*4
    for i,a in enumerate(x):
        for j,b in enumerate(y):
            p=mul(a,b)
            if i&j&1:p=mul(p,d)
            if i&j&2:p=mul(p,e)
            z[i^j]=add(z[i^j],p)
    return tuple(z)

def connected(g,ids):
    seen={ids[0]};todo=[ids[0]]
    while todo:
        i=todo.pop()
        for j in ids:
            if j not in seen and sign(g[i][j])<0:seen.add(j);todo.append(j)
    return len(seen)==len(ids)

def check_strict_escape(y,I,center,old,eligible_ids,d,ca,cb,extra):
    assert sign(pair(y,I,y))>0 and sign(pair(y,I,center))<0
    assert all(sign(x)>0 for x in y)
    assert all(sign(dot(old[k],y))>0 for k in eligible_ids)
    assert sg((dot(ca,y),dot(cb,y)),d)>0
    assert all(sg((dot(a,y),dot(b,y)),rad)>0 for rad,a,b in extra)

def affine_and_rank_gate(original,Q,roots):
    rows=[]
    for label,k in ((3,scale(P,2)),(4,sub(P,O))):
        q=tuple(tuple((*r,Z))for r in Q)+((*(Z for _ in range(9)),k),)
        rr=tuple(tuple((*r,O if i==0 else Z))for i,r in enumerate(roots))
        rr=(tuple(div(x,P)for x in rr[0]),*rr[1:]);g=tuple(tuple(pair(r,q,s)for s in rr)for r in rr)
        aff=[]
        for size in range(2,12):
            for ids in combinations(range(11),size):
                if not connected(g,ids):continue
                if any(sign(sub(mul(g[i][i],g[j][j]),mul(g[i][j],g[i][j])))<0 for i,j in combinations(ids,2)):continue
                block=tuple(tuple(g[i][j]for j in ids)for i in ids)
                v=kernel_one(block)
                if v is None:continue
                if sign(v[0])<0:v=tuple(neg(x)for x in v)
                if not all(sign(x)>0 for x in v):continue
                assert all(x==Z for x in mv(block,v))
                aff.append({'indices':list(ids),'positive_null_coefficients':vv(v)})
        for omit in range(11):
            selected=[i for i in range(11)if i!=omit]
            rk=len(rref(tuple(rr[i]for i in selected))[1])
            obstructions=[r for r in aff if omit not in r['indices']]
            rows.append({'lift_label':label,'deleted_wall':omit,'rank':rk,'retained_affine_subdiagrams':obstructions,
                         'passes_rank_and_no_affine_gate':rk==10 and not obstructions})
    passed=[(r['lift_label'],r['deleted_wall'])for r in rows if r['passes_rank_and_no_affine_gate']]
    assert passed==[(3,4),(3,5),(3,6),(4,4)]
    return rows

def verify():
    raw=json.loads(BASE.read_text());common=json.loads(COMMON.read_text());free=json.loads(FREE.read_text())
    original=json.loads((ROOT/'data/known_examples/bugaenko_compact_h8.json').read_text())
    oldQ=matrix(original['lattice_gram']);oldroots=matrix(original['simple_roots'])
    assert oldQ[0][0]==scale(P,-2) and all(oldQ[0][i]==oldQ[i][0]==Z for i in range(1,9))
    assert all(sign(determinant(tuple(tuple(oldQ[i][j]for j in range(1,n+1))for i in range(1,n+1))))>0 for n in range(1,9))
    assert matrix(original['wall_gram'])==tuple(tuple(pair(r,oldQ,s)for s in oldroots)for r in oldroots)
    Q=tuple(tuple((*r,Z))for r in oldQ)+((*(Z for _ in range(9)),scale(P,2)),)
    lifted=tuple(tuple((*r,O if i==0 else Z))for i,r in enumerate(oldroots))
    lifted=(tuple(div(x,P)for x in lifted[0]),*lifted[1:]);ids=[i for i in range(11)if i!=4];R=tuple(lifted[i]for i in ids)
    assert raw['retained_original_walls']==ids and matrix(raw['host_gram'])==Q and matrix(raw['seed_normals'])==R
    assert len(rref(R)[1])==10
    G=tuple(tuple(pair(r,Q,s)for s in R)for r in R);I=inverse(G)
    assert G==matrix(raw['seed_gram']) and I==matrix(raw['seed_gram_inverse'])
    assert all(G[i][i]==scale(O,2)for i in range(10))
    assert I[0][0]==scale(P,F(1,2))
    polar_cosine=div(O,mul(scale(O,2),I[0][0]));assert polar_cosine==sub(P,O)
    assert sign(sub(polar_cosine,scale(O,F(1,2))))>0
    assert sign(sub(scale(add(O,P),F(1,4)),polar_cosine))>0
    old=[];oldp=[]
    for cap in raw['caps']:
        r=ve(cap['normal']);p=ve(cap['seed_pairings']);c=ve(cap['seed_coefficients'])
        assert pair(r,Q,r)==scale(O,2) and tuple(pair(r,Q,s)for s in R)==p
        assert c==mv(I,p) and r==tuple(total(mul(a,s[j])for a,s in zip(c,R))for j in range(10))
        assert all(sign(t)<=0 and (t in (Z,neg(O),neg(P)) or sign(sub(mul(t,t),scale(O,4)))>=0)for t in p)
        old.append(c);oldp.append(p)
    assert len(old)==434 and len(set(old))==434
    center=ve(common['strict_center_y']);center_coeff=ve(common['positive_center_seed_coefficients'])
    assert all(sign(x)>0 for x in center_coeff)
    expected=tuple(neg(x)for x in mv(G,center_coeff));expected=tuple(div(x,total(expected))for x in expected)
    assert center==expected and all(sign(x)>0 for x in center) and sign(pair(center,I,center))<0
    assert all(sign(dot(c,center))>0 for c in old)
    extra=[]
    for tag,d in (('2',scale(O,2)),('3',scale(O,3)),('2phi',add(scale(O,2),P))):
        data=json.loads((ROOT/f'data/enumerations/fresh_h8_release_h9_quadratic_caps_{tag}_20261013.json').read_text())
        assert parse(data['extension_radicand'])==d
        for row in data['caps']:
            a,b=ve(row['normal_rational']),ve(row['normal_radical']);p,z=ve(row['seed_pairings_rational']),ve(row['seed_pairings_radical'])
            ca,cb=ve(row['seed_coefficients_rational']),ve(row['seed_coefficients_radical'])
            assert add(pair(a,Q,a),mul(d,pair(b,Q,b)))==scale(O,2) and pair(a,Q,b)==Z
            assert tuple(pair(a,Q,s)for s in R)==p and tuple(pair(b,Q,s)for s in R)==z
            assert ca==mv(I,p) and cb==mv(I,z)
            for x,y in zip(p,z):
                assert sg((x,y),d)<=0
                v=em((x,y),(x,y),d)
                assert v[1]==Z and (v[0]in(Z,O,scale(O,2),scale(O,3),add(O,P),add(scale(O,2),P)) or sign(sub(v[0],scale(O,4)))>=0)
            extra.append((d,ca,cb))
    assert len(extra)==2
    common_y=ve(common['strict_escape_y']);assert pair(common_y,I,common_y)==parse(common['strict_escape_norm'])
    assert sign(pair(common_y,I,common_y))>0 and sign(pair(common_y,I,center))<0
    assert all(sign(x)>0 for x in common_y) and all(sign(dot(c,common_y))>0 for c in old)
    assert all(sg((dot(a,common_y),dot(b,common_y)),d)>0 for d,a,b in extra)
    cases=[];proofs=[]
    for index,row in enumerate(free['cases']):
        assert row['index']==index;fixed=ve(row['fixed_pairing_magnitudes']);d=parse(row['discriminant']);x=pe(row['free_magnitude']);j=row['free_pairing_seed_index']
        assert row['free_pairing_original_wall']==ids[j] and fixed[j]==Z and all(t in(Z,O,P)for t in fixed)
        aa=I[j][j];bb=scale(dot(I[j],fixed),2);cc=sub(pair(fixed,I,fixed),scale(O,2))
        assert (aa,bb,cc)==ve(row['norm_equation_coefficients_descending'])
        assert d==sub(mul(bb,bb),scale(mul(aa,cc),4)) and sign(d)>0
        assert sg(es(x,(scale(O,2),Z)),d)>0
        assert sg(ea(ea(em((aa,Z),em(x,x,d),d),em((bb,Z),x,d)),(cc,Z)),d)==0
        p=tuple(neg(add(t,x[0]if k==j else Z))for k,t in enumerate(fixed));z=tuple(neg(x[1])if k==j else Z for k in range(10))
        a,b=ve(row['normal_rational']),ve(row['normal_radical']);ca,cb=mv(I,p),mv(I,z)
        assert ca==ve(row['seed_coefficients_rational']) and cb==ve(row['seed_coefficients_radical'])
        assert tuple(pair(a,Q,r)for r in R)==p and tuple(pair(b,Q,r)for r in R)==z
        assert sg((sub(add(pair(a,Q,a),mul(d,pair(b,Q,b))),scale(O,2)),scale(pair(a,Q,b),2)),d)==0
        assert sg((dot(ca,center),dot(cb,center)),d)>0
        eligible_ids=[k for k,ps in enumerate(oldp)if eligible((dot(ca,ps),dot(cb,ps)),d)]
        assert eligible_ids==row['potentially_compatible_old_cap_indices']
        proof=row['volume_gate'];y=ve(proof['y']);q=pair(y,I,y)
        assert proof['status']=='EXACT_STRICT_SPACELIKE_ESCAPE' and q==parse(proof['norm']) and sign(q)>0
        assert all(sign(t)>0 for t in y) and sign(pair(y,I,center))<0
        assert pair(y,I,center)==parse(proof['center_pairing'])
        assert all(sign(dot(old[k],y))>0 for k in eligible_ids)
        capvalue=(dot(ca,y),dot(cb,y));assert sg(capvalue,d)>0 and capvalue==pe(proof['free_cap_constraint_value'])
        assert all(sg((dot(ae,y),dot(be,y)),de)>0 for de,ae,be in extra)
        assert [pe(v)for v in proof['extra_cap_constraint_values']]==[(dot(ae,y),dot(be,y))for de,ae,be in extra]
        check_strict_escape(y,I,center,old,eligible_ids,d,ca,cb,extra)
        cases.append((d,p,z,ca,cb));proofs.append({'free_cap':index,'eligible_old_caps':len(eligible_ids),'strict_escape_norm':enc(q)})
    assert len(cases)==16
    pair_rows=[];positive_count=0;small_angle_count=0
    for i,j in combinations(range(16),2):
        d,p,z,ca,cb=cases[i];e,pp,zz,cca,ccb=cases[j]
        v=(dot(p,cca),dot(z,cca),dot(p,ccb),dot(z,ccb))
        assert bsg((sub(v[0],scale(O,2)),*v[1:]),d,e)!=0
        if bsg(v,d,e)>0:reason='positive_outward_pairing';positive_count+=1
        else:
            assert bsg(v,d,e)<0 and bsg((add(v[0],P),*v[1:]),d,e)>0
            square=bm(v,v,d,e)
            assert bsg((sub(square[0],O),*square[1:]),d,e)!=0
            assert bsg((sub(square[0],scale(O,2)),*square[1:]),d,e)!=0
            reason='finite_nonright_angle_below_pi5_cosine_but_neither_label3_nor4';small_angle_count+=1
        pair_rows.append({'pair':[i,j],'reason':reason,'biquadratic_pairing_coefficients':vv(v)})
    assert pair_rows==free['free_cap_conflicts'] and free['all_free_cap_pairs_conflict'] is True
    assert free['all_eligible_supersystems_have_exact_escapes'] is True
    # Negative controls exercise the mathematical gates without rerunning
    # any LP or depending on a copied expected status string.
    assert sign(pair(center,I,center))<0
    try:check_strict_escape(center,I,center,old,eligible_ids,cases[-1][0],cases[-1][3],cases[-1][4],extra)
    except AssertionError:timelike_rejected=True
    else:raise AssertionError('timelike fake escape accepted')
    try:assert sg(es(EZ,(scale(O,2),Z)),cases[3][0])>0
    except AssertionError:fake_cap_rejected=True
    else:raise AssertionError('zero dotted-distance fake accepted')
    gates=affine_and_rank_gate(original,oldQ,oldroots)
    return {'status':'PASS','standard_library_only':True,'exact_single_deletion_gates':gates,
            'polar_cosine_squared':enc(polar_cosine),'base_caps':434,'extra_quadratic_caps':2,'free_distance_caps':16,
            'free_cap_positive_pairing_conflicts':positive_count,'free_cap_small_angle_conflicts':small_angle_count,
            'strict_escape_cases':proofs,'timelike_fake_escape_rejected':timelike_rejected,'zero_dotted_distance_rejected':fake_cap_rejected,
            'conclusion':'No finite-volume Coxeter polytope uses these ten retained seed walls and any subset of the 452 archived cap normals.',
            'scope':'Exact finite-bank obstruction and exact rank/affine gates for deleting one wall AFTER either of the two fixed lifts. Does not classify all lifts after wall deletion or any wider cap family.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();out=verify()
    if args.verify:assert out==json.loads(OUTPUT.read_text())
    else:OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'PASS','caps':452,'free_pairs':120,'escape_cases':16,'positive_pairing_conflicts':out['free_cap_positive_pairing_conflicts'],'small_angle_conflicts':out['free_cap_small_angle_conflicts'],'rank_affine_gate_survivors':[(r['lift_label'],r['deleted_wall'])for r in out['exact_single_deletion_gates']if r['passes_rank_and_no_affine_gate']]},indent=2))

if __name__=='__main__':main()
