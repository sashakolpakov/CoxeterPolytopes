#!/usr/bin/env python3
"""Standalone physical and numerical replay of the old7 escape cap.

No project imports, numerical root solving or historical native-bank replay.
The four saved numerical cones are independently reconstructed below.
The same Schur identities are implemented separately, with directed rational
rounding. This is not an alternative interval algorithm or scaffold proof.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations, product, combinations
from functools import lru_cache
from math import isqrt
from pathlib import Path
import json
import math
import numpy as np
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_exact64_escape7_cap_20261018'
PRIOR = 'fresh_oblique_high_s_polar_interval_20261018'
PRIMARY = 'fresh_nonarithmetic_release49_62_replacement_20261017'
HIGH = 'fresh_oblique_high_n_right_r_join_20261018'
WHITE, OLD = [0, 1, 2, 3, 5, 6], [8, 10, 52, 63, 64]
SCALE = 10**70


def iv(x):
    return tuple(map(F, x)) if isinstance(x, (list, tuple)) else (F(x), F(x))


def roundout(v):
    a, b = v
    assert a <= b
    return F(a.numerator*SCALE//a.denominator, SCALE), F(-((-b.numerator*SCALE)//b.denominator), SCALE)


def add(a, b): return roundout((a[0]+b[0], a[1]+b[1]))
def neg(a): return -a[1], -a[0]
def sub(a, b): return add(a, neg(b))
def mul(a, b):
    z = [x*y for x in a for y in b]
    return roundout((min(z), max(z)))
def div(a, b):
    assert b[0]*b[1] > 0
    return mul(a, (1/b[1], 1/b[0]))
def total(values):
    out = iv(0)
    for v in values: out = add(out, v)
    return out
def sqrt(a):
    assert a[0] >= 0
    lo = isqrt(a[0].numerator*SCALE*SCALE//a[0].denominator)
    hi = isqrt(a[1].numerator*SCALE*SCALE//a[1].denominator)
    if F(hi, SCALE)**2 < a[1]: hi += 1
    return F(lo, SCALE), F(hi, SCALE)
def dot(x, g, y): return total(mul(mul(a, g[i][j]), b) for i, a in enumerate(x) for j, b in enumerate(y))
def enc(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, dict): return {str(k): enc(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [enc(v) for v in x]
    return x
def determinant(g):
    n, terms = len(g), []
    for p in permutations(range(n)):
        value = iv((-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1, n)))
        for i in range(n): value = mul(value, g[i][p[i]])
        terms.append(value)
    return total(terms)
def inverse(g):
    n, d = len(g), determinant(g)
    assert d[0]*d[1] > 0
    return [[div(mul(iv((-1)**(i+j)), determinant([[g[a][b] for b in range(n) if b != i]
        for a in range(n) if a != j])), d) for j in range(n)] for i in range(n)], d
def gauss(g):
    n = len(g)
    a = [r[:]+[iv(int(i == j)) for j in range(n)] for i, r in enumerate(g)]
    pivots = []
    for i in range(n):
        p = a[i][i]
        assert p[0]*p[1] > 0
        pivots.append(p)
        a[i] = [div(x, p) for x in a[i]]
        for j in range(n):
            if i != j:
                q = a[j][i]
                a[j] = [sub(x, mul(q, y)) for x, y in zip(a[j], a[i])]
    return [r[n:] for r in a], pivots
def theta(m):
    if m in (2, 3): return iv(m-2)
    if m == 4: return sqrt(iv(2))
    if m == 5: return div(add(iv(1), sqrt(iv(5))), iv(2))
    if m == 6: return sqrt(iv(3))
    assert isinstance(m,int) and m>=7
    def atan(x):
        z = sum((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(100))
        return z, z+x**201/F(201)
    pi = sub(mul(iv(16), atan(F(1, 5))), mul(iv(4), atan(F(1, 239))))
    x = div(pi, iv(m))
    xx, term, value = mul(x, x), iv(1), iv(1)
    for k in range(1, 60):
        term = div(mul(neg(xx), term), iv((2*k-1)*(2*k)))
        value = add(value, term)
    e = xx[1]*max(abs(t) for t in term)/F(119*120)
    return mul(iv(2), (value[0]-e, value[1]+e))


def row_times(v,g):return [total(mul(v[k],g[k][j])for k in range(len(v)))for j in range(len(g))]
def linear(coeff,vecs):return [total(mul(c,v[k])for c,v in zip(coeff,vecs))for k in range(len(vecs[0]))]
def vector_sub(a,b):return [sub(x,y)for x,y in zip(a,b)]


def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())

def refine(z,i,j):
    exact=None;label=None
    if i==j:exact=iv(2);label='norm2'
    elif -1<z[0]<=0<=z[1]<1:exact=iv(0);label=2
    elif -2<z[0]<=z[1]<0:
        m=round(math.pi/math.acos(-float(sum(z)/2)/2))
        if m>=3 and z[0]>-theta(m+1)[0] and z[1]<-theta(m-1)[1]:exact=neg(theta(m));label=m
    if exact is None:return z,None
    val=max(z[0],exact[0]),min(z[1],exact[1]);assert val[0]<=val[1]
    return val,dict(pair=[i,j],raw=z,identity=label,identity_interval=exact,intersection=val)

V=iv
FRAME=[7,30,47];NEW=76
A=SimpleNamespace(TAG='fresh_remote60_exact_original_extension_20261018',PRIMARY=PRIMARY,
    P=SimpleNamespace(inverse=inverse,determinant=determinant),I=SimpleNamespace(inverse=gauss),
    refine=refine,F=F,row_times=row_times,linear=linear,vector_sub=vector_sub)

def physical():
    d=read(A.TAG);p=read(A.PRIMARY)
    for x in(d,p):
        for path,h in x['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    z=d['physical_input'];metric=[list(map(V,r))for r in z['metric']];h=list(map(V,z['center']));h2=V(z['center_square'])
    roots={int(i):list(map(V,v))for i,v in z['roots'].items()};roots[75]=list(map(V,d['normal']));ids=sorted(roots)
    hc={int(i):V(v)for i,v in z['center_products'].items()};hc[75]=V(d['center_product'])
    wi=[list(map(V,r))for r in d['white_inverse']];wc={int(i):list(map(V,v))for i,v in d['white_covectors'].items()};wc[75]=[V(0)]*6
    hw=[hc[w]for w in WHITE];hp={i:sub(hc[i],dot(hw,wi,wc[i]))for i in ids};hp2=sub(h2,dot(hw,wi,hw))
    native={};identities=[]
    for i in FRAME:
        native[i]={}
        for j in ids:
            raw=dot(roots[i],metric,roots[j]);v,pr=A.refine(raw,i,j);native[i][j]=v
            if pr:identities.append(pr)
    assert next(r['identity']for r in identities if r['pair']==[7,30])==5
    assert next(r['identity']for r in identities if r['pair']==[7,47])==3
    assert native[30][47][1]<-2
    g={i:{j:sub(native[i][j],dot(wc[i],wi,wc[j]))for j in ids}for i in FRAME}
    gram=[[g[i][j]for j in FRAME]for i in FRAME];inv,det=A.P.inverse(gram)
    assert gram[0][0][0]>0 and A.P.determinant([r[:2]for r in gram[:2]])[0]>0 and det[1]<0
    hf=[hp[i]for i in FRAME];q=sub(hp2,dot(hf,inv,hf));assert q[0]>0
    # All-right member chosen; the other two label variants are numerical.
    rhs=[V(0)]*3;rem=V(2);fac=neg(sqrt(div(rem,q)))
    products={j:mul(fac,sub(hp[j],dot(hf,inv,[g[i][j]for i in FRAME])))for j in ids}
    for j in WHITE+FRAME:products[j]=V(0)
    cp=mul(fac,q);assert cp[1]<0
    pj={j:A.vector_sub(roots[j],A.linear(A.row_times(wc[j],wi),[roots[w]for w in WHITE]))for j in ids}
    ph=A.vector_sub(h,A.linear(A.row_times(hw,wi),[roots[w]for w in WHITE]))
    res=A.vector_sub(ph,A.linear(A.row_times(hf,inv),[pj[i]for i in FRAME]));v=[mul(fac,x)for x in res]
    direct={j:dot(v,metric,roots[j])for j in ids};schur=products.copy()
    products={j:(max(x[0],direct[j][0]),min(x[1],direct[j][1]))for j,x in products.items()}
    assert all(x[0]<=x[1]for x in products.values())
    dotted=set(ids)-set(WHITE+FRAME);assert len(dotted)==55 and all(products[i][1]<-2 for i in dotted),[(i,enc(products[i]))for i in dotted if products[i][1]>=-2]
    norm=dot(v,metric,v);assert norm[0]<=2<=norm[1]
    oldw={i:sub(cp,mul(div(hc[i],V(2)),products[i]))for i in ids}
    neww={i:sub(hc[i],mul(div(cp,V(2)),products[i]))for i in ids}
    nsq=sub(h2,div(mul(cp,cp),V(2)));assert nsq[1]<0 and all(x[1]<0 for x in oldw.values())and all(x[1]<0 for x in neww.values())
    for i in FRAME:g[i][NEW]=products[i]
    g[NEW]={**products,NEW:V(2)};hp[NEW]=cp;pj[NEW]=v
    vertices=[]
    for support in([7,30,NEW],[7,47,NEW]):
        gv=[[g[i][j]for j in support]for i in support];vi,piv=A.I.inverse(gv);assert all(x[0]>0 for x in piv)
        cv=[hp[i]for i in support];coeff=A.row_times(cv,vi);point=A.vector_sub(ph,A.linear(coeff,[pj[i]for i in support]))
        sq=sub(hp2,dot(cv,vi,cv));assert sq[1]<0
        pp={}
        for j in ids+[NEW]:
            val=sub(hp[j],dot(cv,vi,[g[i][j]for i in support]))
            if j in WHITE+support:val=V(0)
            else:assert val[1]<0,(support,j,enc(val))
            pp[j]=val
        vertices.append(dict(support=WHITE+support,section_support=support,Gram=gv,inverse=vi,pivots=piv,point=point,square=sq,products=pp))
    trials=[]
    for subset in combinations(sorted(WHITE+FRAME),8):
        support=sorted([NEW,*subset])
        if 30 in subset and 47 in subset:trials.append(dict(support=support,dotted_pair=[30,47]))
        else:
            assert support in[sorted(r['support'])for r in vertices]
            trials.append(dict(support=support,finite=True))
    assert len(trials)==9 and sum('finite'in r for r in trials)==2
    oldescape=p['strict_escape'];oldprod=dot(v,metric,list(map(V,oldescape['point'])))
    assert oldprod[0]>0
    # A fresh exact ray at the retained exterior30/59/75 region, followed by
    # a strict center perturbation. All gates below are physical inequalities.
    ef=[30,59,75];eg={}
    for i in ef:
        eg[i]={}
        for j in ids:
            raw=dot(roots[i],metric,roots[j]);val,proof=A.refine(raw,i,j)
            eg[i][j]=sub(val,dot(wc[i],wi,wc[j]))
    ei,edet=A.P.inverse([[eg[i][j]for j in ef]for i in ef]);assert edet[1]<0
    ec=[hp[i]for i in ef];er=A.vector_sub(ph,A.linear(A.row_times(ec,ei),[pj[i]for i in ef]));er=[neg(x)for x in er]
    escape=None
    for power in range(1,50):
        epsilon=V(A.F(1,2**power));point=[add(x,mul(epsilon,y))for x,y in zip(er,h)]
        sq=dot(point,metric,point);future=dot(point,metric,h)
        pp={i:dot(point,metric,r)for i,r in{**roots,NEW:v}.items()}
        if sq[0]>0 and future[1]<0 and all(x[1]<0 for x in pp.values()):
            escape=dict(frame=ef,power=power,point=point,square=sq,future=future,products=pp);break
    return dict(old_walls=ids,new_wall=NEW,deleted_original_walls=[],frame=FRAME,profile=0,labels=[2,2,2],
        identity_refinements=identities,frame_Gram=gram,frame_inverse=inv,frame_determinant=det,
        center_complement_square=q,factor=fac,normal=v,direct_norm=norm,center_product=cp,
        schur_products=schur,direct_products=direct,native_products=products,
        old_facet_witness_products=oldw,new_facet_witness_products=neww,new_facet_witness_square=nsq,
        new_finite_vertices=vertices,new_finite_support_trials=trials,previous_escape_new_product=oldprod,
        remaining_strict_escape=escape,
        summary=dict(total_walls=65,old_finite_vertices_retained=201,new_finite_vertices=2,complete_finite_vertices=203,
            strict_dotted_contacts=55,right_contacts=9,infinite_volume_proved=escape is not None))


def numerical_check(n):
    from itertools import combinations
    assert n['frame']==FRAME and n['profile']==0 and n['declared_labels']==[2,3,4]
    snap=n['snapshot'];metric=np.array(snap['metric']);center=np.array(snap['center'])
    assert np.array_equal(metric,np.diag([-1,1,1,1]))
    old={int(i):np.array(v)for i,v in snap['old_normals'].items()}
    prior=read('fresh_remote60_compact_face_bridge_20261018');base=prior['numerical_snapshot']
    assert len(old)==58 and set(old)==set(map(int,base['old_normals']))|{75}
    assert all(np.array_equal(old[int(i)],np.array(v))for i,v in base['old_normals'].items())
    ub=next(b for b in prior['remote_bank']['recipes'][1]['branches']if b['native_admissible'])
    assert np.array_equal(old[75],np.array(ub['normal']))
    assert np.array_equal(metric,np.array(base['metric']))and np.array_equal(center,np.array(base['center']))
    assert len(n['recipes'])==3
    caps=[]
    for row,m in zip(n['recipes'],[2,3,4]):
        assert row['label']==m and len(row['branches'])==2
        frame=np.array([old[i]for i in FRAME]);g=frame@metric@frame.T
        assert np.linalg.det(g)<0
        rhs=np.array([0,0,-float(sum(theta(m))/2)])
        found=[]
        for b,sign in zip(row['branches'],[-1,1]):
            assert b['branch']==sign
            v=np.array(b['normal']);assert abs(v@metric@v-2)<1e-9
            assert np.max(abs(frame@metric@v-rhs))<1e-9
            pp={i:float(v@metric@r)for i,r in old.items()}
            assert len(b['products'])==58 and [i for i,_ in b['products']]==sorted(old)
            assert max(abs(q-pp[i])for i,q in b['products'])<1e-8
            cp=float(v@metric@center);assert abs(cp-b['center_product'])<1e-9
            if sign==-1:
                assert cp<0 and not b['blockers']and b['native']
                assert all(pp[i]<-2 for i in old if i not in FRAME)
                caps.append(v)
            else:assert cp>0 and not b['native']
            found.append(v)
        assert np.linalg.norm(found[0]-found[1])>1
        # Rank3 equations plus a nonzero norm quadratic admit at most two
        # solutions; the distinct saved solutions exhaust the two sheets.
        assert np.linalg.matrix_rank(frame)==3
    assert len(n['cones'])==4
    for entry,cap,label in zip(n['cones'],[None]+caps,[None,2,3,4]):
        assert entry['label']==label
        full=old.copy()
        if cap is not None:full[NEW]=cap
        ids=sorted(full);cov=np.array([full[i]@metric for i in ids]);cov=cov/np.linalg.norm(cov,axis=1)[:,None]
        dual=cov.sum(axis=0);accepted={}
        for tri in combinations(range(len(ids)),3):
            a=cov[list(tri)]
            v=np.array([(-1)**j*np.linalg.det(np.delete(a,j,axis=1))for j in range(4)])
            norm=np.linalg.norm(v)
            if norm<1e-13:continue
            v=v/norm
            if dual@v>0:v=-v
            pp=cov@v
            if np.max(pp)>2e-9:continue
            active=tuple(ids[j]for j,x in enumerate(pp)if abs(x)<2e-9)
            assert len(active)==3
            accepted[active]=v
        raw=entry['cone'];assert not raw['warnings']
        assert {tuple(sorted(r['active']))for r in raw['rays']}==set(accepted)
        for r in raw['rays']:
            v=accepted[tuple(sorted(r['active']))]
            assert np.linalg.norm(v-np.array(r['coordinates']))<1e-7
            sq=float(v@metric@v);assert abs(sq)>1e-12
            assert r['kind']==('finite'if sq<0 else'exterior')and abs(sq-r['square'])<1e-8
        from collections import Counter,defaultdict
        assert raw['counts']==dict(Counter(r['kind']for r in raw['rays']))
        pairs=defaultdict(list)
        for r in raw['rays']:
            for pair in combinations(sorted(r['active']),2):pairs[pair].append(r['index'])
        assert all(len(v)==2 for v in pairs.values())and len(pairs)==len(raw['edges'])
        ends=0
        for e in raw['edges']:
            key=tuple(sorted(e['sides']));assert sorted(e['rays'])==sorted(pairs[key])
            a,b=[raw['rays'][i]for i in e['rays']];nf=sum(r['kind']=='finite'for r in(a,b))
            ab=float(np.array(a['coordinates'])@metric@np.array(b['coordinates']))
            kind='compact'if nf==2 else'one_open'if nf==1 else'two_open'if ab<-1e-12 and a['square']*b['square']-ab*ab<-1e-15 else'outside'
            assert e['kind']==kind;ends+=1 if kind=='one_open'else 2 if kind=='two_open'else 0
        assert ends==raw['open_ends']==10
        for face in(7,60):
            ledger=[r for r in raw['rays']if face in r['active']]
            assert entry[f'section{face}']==ledger
            if face==60:assert len(ledger)==3 and all(r['kind']=='finite'for r in ledger)
            elif cap is not None:assert len(ledger)==11 and all(r['kind']=='finite'for r in ledger)
        assert raw['counts']==dict(finite=28 if cap is None else 30,exterior=6)


def check():
    d=read(TAG)
    for path,h in d['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    assert d['status']=='PASS_EXACT65_ESCAPE7_CAP'
    seed=read(A.TAG)
    labels={tuple(r['pair']):r['identity']for r in seed['native_isolation']}
    special={(0,1):5,(1,2):3,(2,3):3,(5,6):3}
    assert all(labels[i,j]==special.get((i,j),2)for i,j in combinations(WHITE,2))
    assert seed['white_covectors']['7']==enc([iv(0)]*5+[iv(-1)])
    assert all(seed['white_covectors'][str(i)]==enc([iv(0)]*6)for i in(30,47))
    expected={2:[['A1','H4'],['A1','A4']],3:[['A1','H4'],['A5']],4:[['A1','H4'],['B5']]}
    assert all(r['duties']==expected[r['label']]for r in d['numerical']['recipes'])
    numerical_check(d['numerical']);assert d['physical']==enc(physical())
    assert d['physical']['summary']==dict(total_walls=65,old_finite_vertices_retained=201,new_finite_vertices=2,
        complete_finite_vertices=203,strict_dotted_contacts=55,right_contacts=9,infinite_volume_proved=True)
    assert d['scope']==dict(full_two_corner_duty_domain=True,all_right_member_exact=True,other_two_members_numerical=True,
        compact7_section_numerical=True,full_ambient_facet7_compact=False,finite_volume_proved=False)
    print('PASS exact65-wall cap/203 finite; six native roots; four independently enumerated numerical cones; strict remaining escape')


if __name__=='__main__':check()

