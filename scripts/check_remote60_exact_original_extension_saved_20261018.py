#!/usr/bin/env python3
"""Standalone physical interval replay of the exact remote U extension.

No project imports, numerical root solving, native-bank or cone reconstruction.
The same Schur identities are implemented separately, with directed rational
rounding. This is not an alternative interval algorithm or scaffold proof.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations, product
from functools import lru_cache
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_remote60_exact_original_extension_20261018'
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


def check():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    p=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    assert d['status']=='PASS_EXACT_REMOTE60_ORIGINAL_EXTENSION'
    assert d['frame']==[60,58,59]and d['profile']==0 and d['labels']==[2,2,3]
    assert d['fixed_white_walls']==WHITE and d['new_wall']==75 and d['deleted_original_walls']==[]
    assert p['total_wall_count']==63 and p['complete_finite_vertex_count']==199 and p['new_wall']==66
    metric=[list(map(iv,r))for r in p['computational_gram']];center=list(map(iv,p['common_center']))
    roots={r['wall']:[div(mul(iv(2),sub(x,y)),iv(p['center_products'][str(r['wall'])]))for x,y in zip(center,map(iv,r['point']))]for r in p['facet_witnesses']}
    roots[66]=list(map(iv,p['new_normal']));ids=sorted(roots);frame=[60,58,59];new=75
    assert len(ids)==63 and d['old_walls']==ids
    hc={i:iv(p['center_products'][str(i)])for i in ids};h2=iv(p['center_norm'])
    assert h2[1]<0 and all(z[1]<0 for z in hc.values())
    assert d['physical_input']==enc(dict(metric=metric,center=center,roots=roots,center_products=hc,center_square=h2))
    native={i:{j:dot(roots[i],metric,roots[j])for j in ids}for i in WHITE+frame}
    seen=set()
    for row in d['native_isolation']:
        i,j=row['pair'];assert i in WHITE+frame and j in ids and(i,j)not in seen;seen.add((i,j))
        z=native[i][j];assert row['raw']==enc(z);m=row['identity']
        if m=='norm2':assert i==j;exact=iv(2)
        elif m==2:assert i!=j and -1<z[0]<=0<=z[1]<1;exact=iv(0)
        else:
            assert i!=j and isinstance(m,int)and m>=3
            assert -theta(m+1)[0]<z[0]<=z[1]<-theta(m-1)[1]
            exact=neg(theta(m))
        intersection=max(z[0],exact[0]),min(z[1],exact[1]);assert intersection[0]<=intersection[1]
        assert row['identity_interval']==enc(exact)and row['intersection']==enc(intersection)
        native[i][j]=intersection
    assert all((i,i)in seen for i in WHITE+frame)
    assert next(r['identity']for r in d['native_isolation']if r['pair']==[58,59])==9
    wi,wp=gauss([[native[i][j]for j in WHITE]for i in WHITE]);assert all(z[0]>0 for z in wp)
    wc={i:[native[w][i]for w in WHITE]for i in ids};hw=[hc[w]for w in WHITE]
    assert all(wc[i][k]==iv(0)for i in(58,59)for k in range(6))
    hp={i:sub(hc[i],dot(hw,wi,wc[i]))for i in ids};hp2=sub(h2,dot(hw,wi,hw));assert hp2[1]<0
    projected={i:{j:sub(native[i][j],dot(wc[i],wi,wc[j]))for j in ids}for i in frame}
    gf=[[projected[i][j]for j in frame]for i in frame];inv,piv=gauss(gf);assert all(z[0]>0 for z in piv)
    rhs=list(map(iv,[0,0,-1]));hf=[hp[i]for i in frame]
    rem=sub(iv(2),dot(rhs,inv,rhs));q=sub(hp2,dot(hf,inv,hf));assert rem[1]<0 and q[1]<0
    fac=sqrt(div(rem,q));coef=row_times(rhs,inv);hcoef=row_times(hf,inv)
    prods={}
    for j in ids:
        cj=[projected[i][j]for i in frame]
        prods[j]=add(dot(rhs,inv,cj),mul(fac,sub(hp[j],dot(hf,inv,cj))))
    for j in WHITE:prods[j]=iv(0)
    for j,z in zip(frame,rhs):prods[j]=z
    cp=add(dot(rhs,inv,hf),mul(fac,q));op=sub(dot(rhs,inv,hf),mul(fac,q))
    assert cp[1]<0 and op[0]>0
    pj={j:vector_sub(roots[j],linear(row_times(wc[j],wi),[roots[w]for w in WHITE]))for j in ids}
    ph=vector_sub(center,linear(row_times(hw,wi),[roots[w]for w in WHITE]))
    residual=vector_sub(ph,linear(hcoef,[pj[j]for j in frame]))
    u=[add(a,mul(fac,b))for a,b in zip(linear(coef,[pj[j]for j in frame]),residual)]
    norm=dot(u,metric,u);assert norm[0]<=2<=norm[1]
    direct={j:dot(u,metric,roots[j])for j in ids}
    schur=prods.copy();prods={j:(max(z[0],direct[j][0]),min(z[1],direct[j][1]))for j,z in prods.items()}
    assert all(z[0]<=z[1]for z in prods.values())
    dotted=set(ids)-set(WHITE+frame);assert len(dotted)==54 and all(prods[j][1]<-2 for j in dotted)
    assert [prods[i]for i in WHITE+frame]==[iv(0)]*8+[iv(-1)]
    expected=dict(white_inverse=wi,white_pivots=wp,white_covectors=wc,projected_center_square=hp2,
        projected_frame_products=projected,frame_Gram=gf,frame_inverse=inv,frame_pivots=piv,rhs=rhs,
        remainder=rem,center_complement_square=q,positive_factor=fac,normal=u,direct_norm=norm,
        direct_products=direct,schur_products=schur,native_products=prods,center_product=cp,opposite_sheet_center_product=op)
    assert all(d[k]==enc(v)for k,v in expected.items())
    fp=vector_sub(center,[mul(div(cp,iv(2)),z)for z in u]);fs=sub(h2,div(mul(cp,cp),iv(2)))
    fprod={i:sub(hc[i],mul(div(cp,iv(2)),prods[i]))for i in ids}
    oldprod={i:sub(cp,mul(div(hc[i],iv(2)),prods[i]))for i in ids}
    assert fs[1]<0 and all(z[1]<0 for z in fprod.values())and all(z[1]<0 for z in oldprod.values())
    assert d['facet_witness']==enc(dict(point=fp,square=fs,products=fprod))
    assert d['retained_old_facet_witness_products']==enc(oldprod)
    allroots={**pj,new:u};hp[new]=cp;g={}
    for i in frame+[new]:
        g[i]={j:iv(2)if i==j==new else prods[j]if i==new else prods[i]if j==new else projected[i][j]for j in ids+[new]}
    assert d['projected_center_products']==enc(hp)
    assert len(d['section_vertices'])==3
    for record,pair in zip(d['section_vertices'],([58,59],[58,new],[59,new])):
        support=[60]+pair;gv=[[g[i][j]for j in support]for i in support]
        vi,vp=gauss(gv);assert all(z[0]>0 for z in vp)
        cv=[hp[i]for i in support];vc=row_times(cv,vi);point=vector_sub(ph,linear(vc,[allroots[i]for i in support]))
        square=sub(hp2,dot(cv,vi,cv));assert square[1]<0
        pp={}
        for j in ids+[new]:
            z=sub(hp[j],dot(cv,vi,[g[i][j]for i in support]))
            if j in WHITE+support:z=iv(0)
            else:assert z[1]<0
            pp[j]=z
        assert record==enc(dict(section_support=support,ambient_support=WHITE+support,Gram=gv,pivots=vp,inverse=vi,
            coefficients=vc,point=point,square=square,products=pp))
    tg=[[g[i][j]for j in(58,59,new)]for i in(58,59,new)];td=determinant(tg);assert td[1]<0
    assert d['triangle_Gram']==enc(tg)and d['triangle_determinant']==enc(td)
    from itertools import combinations
    trials=[]
    for subset in combinations(sorted(WHITE+frame),8):
        support=sorted([new,*subset])
        if 58 in subset and 59 in subset:trials.append(dict(support=support,indefinite_principal=[58,59,new]))
        else:
            assert support in[sorted(r['ambient_support'])for r in d['section_vertices']if new in r['ambient_support']]
            trials.append(dict(support=support,certified_finite=True))
    assert len(trials)==9 and sum('certified_finite'in r for r in trials)==2 and d['new_finite_vertex_trials']==trials
    escape=p['strict_escape'];eprod=dot(u,metric,list(map(iv,escape['point'])))
    assert iv(escape['norm'])[0]>0 and iv(escape['future'])[1]<0 and eprod[1]<0
    assert len(escape['products'])==63 and all(iv(z)[1]<0 for z in escape['products'].values())
    assert d['retained_strict_escape']==enc(dict(inherited=escape,new_product=eprod))
    assert d['summary']==dict(original_walls=63,total_walls=64,strict_dotted_original_contacts=54,
        right_original_contacts=8,simple_original_contacts=1,retained_old_finite_vertices=199,
        added_certified_finite_vertices=2,complete_finite_vertices=201,compact_section_face=60,section_face_angles=[9,2,3])
    assert d['scope']==dict(exact_physical_extension=True,numerical_chain_required=False,complete_new_finite_vertex_census=True,
        ambient_facet60_compact=False,infinite_volume_proved=True,finite_volume_proved=False,new_dimension_example_proved=False)
    print('PASS exact64-wall extension;63 native gates;64 strict facet witnesses;201 finite vertices;compact60 section;strict escape')


if __name__=='__main__':check()

