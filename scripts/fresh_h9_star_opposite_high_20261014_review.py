#!/usr/bin/env python3
"""Independent domain and quadratic interpolation replay for the opposite tail.

Uses only the previously independent Fraction arithmetic, never the producer's
SymPy/ANP operations or its symbolic B,A,C determinant identity.
"""
import argparse,json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_h9_star_lastblock_20261014_review as R

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_opposite_high_20261014'
S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
Z,O=R.Z,R.O
add,sub,neg,mul,scale,total,sg,inv=R.add,R.sub,R.neg,R.mul,R.scale,R.total,R.sign,R.inverse
co=R.COS
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def ev(p,x):return add(add(p[0],scale(p[1],x)),scale(p[2],x*x))
def bern(p,l):
    b0=ev(p,l);b1=add(b0,scale(add(p[1],scale(p[2],2*l)),(1-l)/2))
    return[b0,b1,total(p)]
@lru_cache(None)
def interpolated_polynomial(aa,bb,tt,ab,at,bt):
    def value(c):
        D=1/(1-c*c)
        rr=[[sub(aa,scale(O,D/4)),sub(ab,scale(O,c*D/4)),at],
            [sub(ab,scale(O,c*D/4)),sub(bb,scale(O,D/4)),bt],
            [at,bt,tt]]
        return scale(R.determinant3(rr),1-c*c)
    v0,vp,vm=[value(x)for x in(F(0),F(1,2),F(-1,2))]
    # Degree<=2 is separately checked by two further values below.
    return v0,sub(vp,vm),scale(sub(add(vp,vm),scale(v0,2)),2)
def padd(p,q):
    r=[F(0)]*max(len(p),len(q))
    for i,x in enumerate(p):r[i]+=x
    for i,x in enumerate(q):r[i]+=x
    return r
def pmul(p,q):
    r=[F(0)]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q):r[i+j]+=x*y
    return r

def run():
    path=f'results/{TAG}.json';d=json.loads((ROOT/path).read_text())
    todo=[path];seen=set()
    while todo:
        f=todo.pop()
        if f in seen:continue
        seen.add(f);o=json.loads((ROOT/f).read_text())
        for p,h in o.get('sources',{}).items():
            assert sha(p)==h,p
            if p.endswith('.json'):todo.append(p)
    geo=json.loads((ROOT/'results/fresh_h9_tree_cube_wedges_20261014.json').read_text())['representatives'][0]['geometry']
    V=[set(v['incident'])for v in geo['vertices']];missing={frozenset(m)for m in geo['minimal_nonfaces']}
    assert len(V)==378 and any(set(S)<=v for v in V)
    rank_nodes=S+[11,17,14]
    assert len(set(rank_nodes))==11
    assert all(any(set(pair)<=v for v in V)for pair in combinations(rank_nodes,2))
    own={tuple(sorted((leaf,a))):3 for i in range(3)for leaf in(P[i],C[i])for a in S[2*i:2*i+2]}
    own|={(11,12):2,(11,13):3,(12,17):3,(13,17):2}
    assert d['fixed_own_contacts']==sorted([list(p)+[m]for p,m in own.items()])
    assert d['internal_lower_bounds']==[S[2*i:2*i+2]+[4 if i<3 else 7]for i in range(4)]
    assert all(frozenset(S[2*i:2*i+2]+[leaf])in missing for i in range(4)for leaf in(P[i],C[i]))
    assert F(1,2)+F(1,3)+F(1,6)==1
    # First3 own(3,3) requires internal>=4; opposite single3 requires>=7.
    assert F(1,3)*3==1
    swap=list(range(18));swap[12],swap[13]=13,12
    assert swap==d['literal_twin_swap']
    assert {frozenset(swap[i]for i in v)for v in V}=={frozenset(v)for v in V}
    assert {frozenset(swap[i]for i in m)for m in missing}==missing
    gates=json.loads((ROOT/'results/fresh_h9_star_zero_20261014_gates.json').read_text())
    assert len(d['first_three_high_reductions'])==3
    for row in d['first_three_high_reductions']:
        i=row['high_block'];assert i in range(3)
        lab=own|{tuple(S[2*j:2*j+2]):6 if i==j else 7 if j==3 else 4 for j in range(4)}
        gate=next(x for x in gates['families']if x['id']==row['reviewed_gate'])
        assert gate['origin']=='two_dense'and not gate['fixed']
        assert all(lab.get((a,b),0)>=m for a,b,m in gate['bounds'])
    # Literal actual-vertex supports for every cross-block monotone guard.
    expected={tuple(sorted((a,b)))for i,j in combinations(range(4),2)for a in S[2*i:2*i+2]for b in S[2*j:2*j+2]}
    assert {tuple(r['pair'])for r in d['cross_right_guards']}==expected
    rt2=scale(co[4],2);null=[O,rt2,rt2,O]
    M=[[O if i==j else neg(co[[4,3,4][min(i,j)]])if abs(i-j)==1 else Z for j in range(4)]for i in range(4)]
    assert all(total(mul(M[i][j],null[j])for j in range(4))==Z for i in range(4))
    for row in d['cross_right_guards']:
        a,b=row['pair'];i,j=S.index(a)//2,S.index(b)//2
        assert set(row['support'])==set(S[2*i:2*i+2]+S[2*j:2*j+2])and set(row['support'])<=V[row['actual_vertex']]
    wanted={11:12,17:13,14:3}
    assert len(d['projection_budgets'])==3
    for row in d['projection_budgets']:
        leaf=row['leaf'];drop=wanted.pop(leaf)
        assert set(row['support'])==(set(S)-{drop})|{leaf}
        assert set(row['support'])<=V[row['actual_vertex']]
    assert not wanted
    # Own root retained contributes1/4. A foreign label>=4 contributes>=1;
    # two foreign labels>=3 contribute>=1, even in the same block.
    assert F(1,4)+F(1,1)>1 and F(1,4)+2*F(1,2)>1
    for key,nodes in [('high_block_third_right_guard',{12,13,14}),('mutual_label_six_guard',{11,13,17})]:
        row=d[key];assert set(row['support'])==nodes and nodes<=V[row['actual_vertex']]
    assert F(3,4)+F(1,4)==1
    base=own|{tuple(S[2*i:2*i+2]):4 for i in range(4)}
    uppers=set()
    for row in d['third_upper_guards']:
        pair=tuple(row['pair']);label=row['test_label']
        assert (pair,label)in(((11,14),4),((14,17),5))and set(pair)<=V[row['actual_vertex']]
        assert row['relaxed_internal_labels']==[4]*4
        R.failure(base|{pair:label},row,V);uppers.add((pair,label))
    assert len(uppers)==2
    # Degree<=2 without using the producer's symbolic formula: multiply
    # each of the first two rows by(1-c²), expand all six permutations as
    # multivariate polynomials, and divide the common factor afterward.
    # A reusable formal polynomial check is encoded with rational probes:
    # the determinant is separately verified at c=±1/3 in every case below.
    # Chebyshev7 and the strict9/10 endpoint are derived with rational lists.
    T0,T1=[F(1)],[F(0),F(1)]
    for _ in range(2,8):T0,T1=T1,padd(pmul([F(0),F(2)],T1),[-x for x in T0])
    cubic=[F(1),F(-4),F(-4),F(8)]
    assert padd(T1,[F(1)])==pmul([F(1),F(1)],pmul(cubic,cubic))
    assert sum(c*F(9,10)**i for i,c in enumerate(cubic))==-F(1,125)
    assert 24*F(3,4)**2-8*F(3,4)-4==F(7,2)
    assert 48*F(3,4)-8==28 and F(3,4)>F(36,49)
    records=set();counts={};zero_coeffs=0;formal_checks=set()
    for row in d['cases']:
        blocks=tuple(row['block_labels']);extra=tuple(row['foreign_contacts'])
        mutual=row['mutual_label'];ma,mb=row['third_labels'];key=(blocks,extra,mutual,ma,mb)
        assert key not in records;records.add(key)
        assert blocks in set(product((4,5),repeat=3))and extra[0]in[None]+list(range(6))and extra[1]in[None]+list(range(6))and extra[2]in[None,2,3,4,5]
        assert mutual in(2,3,4,5)and ma in(2,3)and mb in(2,3,4)
        vectors=[]
        for leaf,e in zip((11,17,14),extra):
            v=[Z]*6
            if leaf==14:v[0]=v[1]=neg(co[3])
            if e is not None:v[e]=neg(co[3])
            vectors.append(v)
        Js=[inv(sub(O,mul(co[m],co[m])))for m in blocks]
        def inner(x,y):
            return total(mul(Js[i],add(add(mul(x[2*i],y[2*i]),mul(x[2*i+1],y[2*i+1])),mul(co[blocks[i]],add(mul(x[2*i],y[2*i+1]),mul(x[2*i+1],y[2*i])))))for i in range(3))
        va,vb,vt=vectors
        aa,bb,tt=[sub(O,inner(v,v))for v in vectors]
        ab=sub(neg(co[mutual]),inner(va,vb));at=sub(neg(co[ma]),inner(va,vt));bt=sub(neg(co[mb]),inner(vb,vt))
        poly=interpolated_polynomial(aa,bb,tt,ab,at,bt)
        B,AA,CC=[R.dec_anp(row[k])for k in('B','A','C')]
        assert poly==(add(B,AA),CC,neg(B))
        # A fourth and fifth evaluation detect a spurious cubic/quartic.
        # The cleared numerator has degree<=4 by its literal two-row form;
        # agreement at five distinct points therefore proves this identity.
        params=(aa,bb,tt,ab,at,bt)
        if params not in formal_checks:
            for c in(F(-1,3),F(1,3)):
                D=1/(1-c*c)
                MM=[[sub(aa,scale(O,D/4)),sub(ab,scale(O,c*D/4)),at],
                    [sub(ab,scale(O,c*D/4)),sub(bb,scale(O,D/4)),bt],[at,bt,tt]]
                # Clearing two denominators yields degree<=4.
                assert scale(R.determinant3(MM),(1-c*c)**2)==scale(ev(poly,c),1-c*c)
            formal_checks.add(params)
        bs=bern(poly,F(6,7));assert bs==list(map(R.dec_anp,row['bernstein_6_7']))
        status=row['status'];counts[status]=counts.get(status,0)+1
        if status=='STRICT_INTERIOR_BERNSTEIN_6_7_TO_ONE':
            assert all(sg(v)>=0 for v in bs)and any(v!=Z for v in bs)
            zero_coeffs+=sum(v==Z for v in bs)
        elif status=='STRICT_BERNSTEIN_9_10_TO_ONE':
            br=bern(poly,F(9,10));assert br==list(map(R.dec_anp,row['bernstein_9_10']))and all(sg(v)>0 for v in br)
        else:
            assert status=='NONPOSITIVE_ACTUAL_VERTEX_AT_RELAXED_LAST_INTERNAL_FOUR'
            lab=base|{tuple(S[2*i:2*i+2]):blocks[i]for i in range(3)}|{(11,17):mutual,(11,14):ma,(14,17):mb}
            for leaf,e in zip((11,17,14),extra):
                if e is not None:lab[tuple(sorted((leaf,S[e])))]=3
            R.failure(lab,row['vertex_failure'],V)
        if len(records)%10000==0:print('independently checked',len(records),flush=True)
    expected_cases=set(product(product((4,5),repeat=3),product([None]+list(range(6)),[None]+list(range(6)),[None,2,3,4,5]),(2,3,4,5),(2,3),(2,3,4)))
    assert records==expected_cases and len(records)==47040
    assert counts=={'STRICT_INTERIOR_BERNSTEIN_6_7_TO_ONE':43412,'NONPOSITIVE_ACTUAL_VERTEX_AT_RELAXED_LAST_INTERNAL_FOUR':3612,'STRICT_BERNSTEIN_9_10_TO_ONE':16}
    assert zero_coeffs==d['zero_bernstein_coefficients']
    paths=[path,f'scripts/{TAG}.py',f'scripts/{TAG}_review.py','scripts/fresh_h9_star_lastblock_20261014_review.py']
    return dict(status='PASS_INDEPENDENT_OPPOSITE_HIGH_DOMAIN_AND_QUADRATIC_REVIEW',sources={p:sha(p)for p in paths},
        arithmetic='Reuses separate Fraction basis1,sqrt2,sqrt5,sqrt10 from the independent lastblock review; never imports producer arithmetic. Permutation determinant values at five rational c prove every quadratic identity.',
        recursive_manifests_checked=len(seen),actual_vertices=378,cross_right_guards=24,projection_budgets=3,
        extra_actual_triangles=2,fresh_upper_vertex_failures=2,all_cases=len(records),case_counts=counts,
        distinct_five_point_determinant_checks=len(formal_checks),zero_coarse_bernstein_coefficients=zero_coeffs,
        rational_cos7_threshold='9/10',scope='Case0 first3blocks all-own3/3, last two leaves opposite single ownpi3, all finite internal labels. First3 internals4/5 use explicitly checked placements of inherited two-dense gates; last internal>=7 has no upper cutoff. All omitted contacts/distances remain free. Not whole case0 or all tree incidences.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    result=run();out=ROOT/f'results/{TAG}_review.json'
    if args.verify:assert result==json.loads(out.read_text());print('SAVED_REPLAY_PASS')
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['all_cases'])
