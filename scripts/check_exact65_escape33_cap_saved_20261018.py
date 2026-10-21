#!/usr/bin/env python3
"""Standalone physical and numerical replay of the changed old33 escape cap.

No project imports, numerical root solving or historical native-bank replay.
The28 saved numerical cones are independently reconstructed below.
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
TAG = 'fresh_exact65_escape33_cap_20261018'
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
FRAME=[30,33,66];NEW=77
F_TAG='fresh_exact64_escape7_cap_20261018'
A=SimpleNamespace(TAG='fresh_remote60_exact_original_extension_20261018',PRIMARY=PRIMARY,
    P=SimpleNamespace(inverse=inverse,determinant=determinant),I=SimpleNamespace(inverse=gauss),
    refine=refine,F=F,row_times=row_times,linear=linear,vector_sub=vector_sub)

def physical():
    d=read(A.TAG);p=read(A.PRIMARY);prior=read(F_TAG)
    for x in(d,p,prior):
        for path,h in x['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    z=d['physical_input'];metric=[list(map(V,r))for r in z['metric']];h=list(map(V,z['center']));h2=V(z['center_square'])
    roots={int(i):list(map(V,v))for i,v in z['roots'].items()};roots[75]=list(map(V,d['normal']));roots[76]=list(map(V,prior['physical']['normal']));ids=sorted(roots)
    hc={int(i):V(v)for i,v in z['center_products'].items()};hc[75]=V(d['center_product']);hc[76]=V(prior['physical']['center_product'])
    wi=[list(map(V,r))for r in d['white_inverse']];wc={int(i):list(map(V,v))for i,v in d['white_covectors'].items()};wc[75]=[V(0)]*6;wc[76]=[V(0)]*6
    hw=[hc[w]for w in WHITE];hp={i:sub(hc[i],dot(hw,wi,wc[i]))for i in ids};hp2=sub(h2,dot(hw,wi,hw))
    native={};identities=[]
    for i in FRAME:
        native[i]={}
        for j in ids:
            raw=dot(roots[i],metric,roots[j]);v,pr=A.refine(raw,i,j);native[i][j]=v
            if pr:identities.append(pr)
    assert next(r['identity']for r in identities if r['pair']==[30,33])==2
    assert next(r['identity']for r in identities if r['pair']==[33,66])==4
    assert native[30][66][1]<-2
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
    dotted=set(ids)-set(WHITE+FRAME);assert len(dotted)==56 and all(products[i][1]<-2 for i in dotted),[(i,enc(products[i]))for i in dotted if products[i][1]>=-2]
    norm=dot(v,metric,v);assert norm[0]<=2<=norm[1]
    oldw={i:sub(cp,mul(div(hc[i],V(2)),products[i]))for i in ids}
    neww={i:sub(hc[i],mul(div(cp,V(2)),products[i]))for i in ids}
    nsq=sub(h2,div(mul(cp,cp),V(2)));assert nsq[1]<0 and all(x[1]<0 for x in oldw.values())and all(x[1]<0 for x in neww.values())
    for i in FRAME:g[i][NEW]=products[i]
    g[NEW]={**products,NEW:V(2)};hp[NEW]=cp;pj[NEW]=v
    vertices=[]
    for support in([30,33,NEW],[33,66,NEW]):
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
        if 30 in subset and 66 in subset:trials.append(dict(support=support,dotted_pair=[30,66]))
        else:
            assert support in[sorted(r['support'])for r in vertices]
            trials.append(dict(support=support,finite=True))
    assert len(trials)==9 and sum('finite'in r for r in trials)==2
    # The prior strict escape need not be the region this changed anchor cuts.
    oldescape=prior['physical']['remaining_strict_escape'];oldprod=dot(v,metric,list(map(V,oldescape['point'])))
    # A fresh exact ray at the retained exterior30/52/76 region, followed by
    # a strict center perturbation. All gates below are physical inequalities.
    ef=[30,52,76];eg={}
    for i in ef:
        eg[i]={}
        for j in ids:
            raw=dot(roots[i],metric,roots[j]);val,proof=A.refine(raw,i,j)
            eg[i][j]=sub(val,dot(wc[i],wi,wc[j]))
    ei,edet=A.P.inverse([[eg[i][j]for j in ef]for i in ef]);assert edet[1]<0
    ec=[hp[i]for i in ef];eq=sub(hp2,dot(ec,ei,ec));assert eq[0]>0
    er=A.vector_sub(ph,A.linear(A.row_times(ec,ei),[pj[i]for i in ef]));er=[neg(x)for x in er]
    escape=None
    for power in range(1,50):
        epsilon=V(A.F(1,2**power));point=[add(x,mul(epsilon,y))for x,y in zip(er,h)]
        direct_sq=dot(point,metric,point);direct_future=dot(point,metric,h)
        sq=sub(add(eq,mul(mul(epsilon,epsilon),h2)),mul(V(2),mul(epsilon,eq)))
        future=add(neg(eq),mul(epsilon,h2))
        direct_pp={i:dot(point,metric,r)for i,r in{**roots,NEW:v}.items()}
        pp={}
        for i in ids+[NEW]:
            cols=[eg[j][i]for j in ef]if i!=NEW else[products[j]for j in ef]
            val=add(sub(dot(ec,ei,cols),hp[i]),mul(epsilon,hc[i]if i!=NEW else cp))
            if i in WHITE+ef:val=mul(epsilon,hc[i])
            pp[i]=(max(val[0],direct_pp[i][0]),min(val[1],direct_pp[i][1]));assert pp[i][0]<=pp[i][1]
        sq=(max(sq[0],direct_sq[0]),min(sq[1],direct_sq[1]));assert sq[0]<=sq[1]
        future=(max(future[0],direct_future[0]),min(future[1],direct_future[1]));assert future[0]<=future[1]
        if sq[0]>0 and future[1]<0 and all(x[1]<0 for x in pp.values()):
            escape=dict(frame=ef,power=power,point=point,square=sq,future=future,products=pp,complement_square=eq,direct_square=direct_sq,direct_future=direct_future,direct_products=direct_pp);break
    return dict(old_walls=ids,new_wall=NEW,deleted_original_walls=[],frame=FRAME,profile=0,labels=[2,2,2],
        identity_refinements=identities,frame_Gram=gram,frame_inverse=inv,frame_determinant=det,
        center_complement_square=q,factor=fac,normal=v,direct_norm=norm,center_product=cp,
        schur_products=schur,direct_products=direct,native_products=products,
        old_facet_witness_products=oldw,new_facet_witness_products=neww,new_facet_witness_square=nsq,
        new_finite_vertices=vertices,new_finite_support_trials=trials,previous_escape_new_product=oldprod,
        remaining_strict_escape=escape,
        summary=dict(total_walls=66,old_finite_vertices_retained=203,new_finite_vertices=2,complete_finite_vertices=205,
            strict_dotted_contacts=56,right_contacts=9,infinite_volume_proved=escape is not None))


@lru_cache(None)
def graph_key(flat):
    n=math.isqrt(len(flat));g=[flat[i*n:(i+1)*n]for i in range(n)]
    return min(tuple(g[p[i]][p[j]]for i in range(n)for j in range(i+1,n))for p in permutations(range(n)))

@lru_cache(None)
def finite_catalogue():
    out={}
    def add_graph(name,n,edges):
        g=[[2]*n for _ in range(n)]
        for i in range(n):g[i][i]=1
        for i,j,m in edges:g[i][j]=g[j][i]=m
        out[n,graph_key(tuple(x for r in g for x in r))]=name
    for n in range(1,6):add_graph('A'+str(n),n,[(i,i+1,3)for i in range(n-1)])
    for n in range(2,6):add_graph('B'+str(n),n,[(i,i+1,4 if i==0 else 3)for i in range(n-1)])
    for n in(4,5):add_graph('D'+str(n),n,[(i,i+1,3)for i in range(n-3)]+[(n-3,n-2,3),(n-3,n-1,3)])
    add_graph('F4',4,[(0,1,3),(1,2,4),(2,3,3)])
    for n in(3,4):add_graph('H'+str(n),n,[(i,i+1,5 if i==0 else 3)for i in range(n-1)])
    add_graph('I2(5)',2,[(0,1,5)])
    return out

def finite_duty(pi,b,m,n):
    g=[[2]*5 for _ in range(5)]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    if pi:g[0 if pi<=3 else 1][4]=g[4][0 if pi<=3 else 1]=3+(pi-1)%3
    for i,j,l in[(2,3,b),(2,4,m),(3,4,n)]:g[i][j]=g[j][i]=l
    unseen=set(range(5));names=[]
    while unseen:
        todo=[min(unseen)];comp=set()
        while todo:
            i=todo.pop()
            if i in comp:continue
            comp.add(i);todo.extend(j for j in unseen if j!=i and g[i][j]!=2)
        unseen-=comp;vs=sorted(comp);key=len(vs),graph_key(tuple(g[i][j]for i in vs for j in vs))
        name=finite_catalogue().get(key)
        if name is None:return None
        names.append(name)
    return sorted(names)

def native(x):
    if x<=-2+2e-7:return True
    if abs(x)<2e-7:return True
    if x>0 or x<=-2:return False
    m=round(math.pi/math.acos(-x/2))
    return m>=3 and abs(x+2*math.cos(math.pi/m))<2e-7

def numerical_check(n):
    from collections import Counter
    snap=n['snapshot'];metric=np.array(snap['metric']);center=np.array(snap['center'])
    assert np.array_equal(metric,np.diag([-1,1,1,1]))
    old={int(i):np.array(v)for i,v in snap['old_normals'].items()}
    prior=read(F_TAG);base=prior['numerical']['snapshot']
    assert len(old)==59 and set(old)==set(map(int,base['old_normals']))|{76}
    assert all(np.array_equal(old[int(i)],np.array(v))for i,v in base['old_normals'].items())
    assert np.array_equal(old[76],np.array(prior['numerical']['recipes'][0]['branches'][0]['normal']))
    assert np.array_equal(metric,np.array(base['metric']))and np.array_equal(center,np.array(base['center']))
    raw=read('fresh_remote60_compact_face_bridge_20261018')['numerical_snapshot'];profiles=snap['profiles']
    assert profiles==raw['profiles']and len(profiles)==7
    original=sorted(map(int,raw['old_normals']));index={i:k for k,i in enumerate(original)}
    admitted=[]
    for dom,frame,bonds,count in zip(n['domains'],[[30,59,75],FRAME],[[5,3],[2,4]],[17,38],strict=True):
        assert dom['frame']==frame and dom['bonds']==bonds and dom['profiles']==list(range(7))and dom['labels']==[2,3,4,5]
        expected={(pi,ms)for pi,ms in product(range(7),product(range(2,6),repeat=3))}
        seen=set()
        for row in dom['rejected_finite_duties']+dom['recipes']:
            pi=row['profile'];ms=tuple(row['labels']);key=pi,ms;assert key in expected and key not in seen;seen.add(key)
            duties=[finite_duty(pi,b,*xy)for b,xy in zip(bonds,(ms[:2],ms[1:]))]
            # Catalogue aliases for I2(5) are accepted only as the same type.
            actual=[None if x is None else sorted({'I2(3)':'A2','I2(4)':'B2'}.get(y,y)for y in x)for x in row['duties']]
            assert actual==duties
            assert (any(x is None for x in duties))==('branches'not in row)
            if 'branches'not in row:continue
            assert all(x is not None for x in actual)
            assert len(row['branches'])==2 and [b['branch']for b in row['branches']]==[-1,1]
            ff=np.array([old[i]for i in frame]);rhs=np.array([-float(sum(theta(m))/2)for m in ms]);vs=[]
            assert np.linalg.matrix_rank(ff)==3
            for b in row['branches']:
                v=np.array(b['normal']);vs.append(v);pr=profiles[pi]
                assert abs(v@metric@v-pr['square'])<1e-8 and np.max(abs(ff@metric@v-rhs))<1e-8
                pp={i:float(v@metric@r)+(pr['offsets'][index[i]]if i in index else 0)for i,r in old.items()}
                assert [i for i,_ in b['products']]==sorted(old)and len(b['products'])==59
                assert max(abs(q-pp[i])for i,q in b['products'])<1e-7
                cp=float(v@metric@center)+pr['center_offset'];assert abs(cp-b['center_product'])<1e-8
                bad=sorted(i for i,x in pp.items()if not native(x))
                assert bad==[i for i,_ in b['blockers']]
                assert b['native']==(not bad and cp<0)
                if b['native']:
                    assert frame==FRAME
                    admitted.append((pi,list(ms),b['branch'],v))
            assert np.linalg.norm(vs[0]-vs[1])>1e-5
        assert seen==expected and len(dom['recipes'])==count
    assert len(admitted)==len(n['accepted'])==28
    assert n['summary']==dict(eligible_recipes=[17,38],norm_roots=[34,76],native=[0,28],deleted_original_walls=[])
    for entry,(pi,ms,sign,cap)in zip(n['accepted'],admitted,strict=True):
        assert(entry['profile'],entry['labels'],entry['branch'])==(pi,ms,sign)
        assert np.array_equal(np.array(entry['normal']),cap)
        full={**old,NEW:cap}
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
        for face in(7,33,60,75,76,77):
            ledger=[r for r in raw['rays']if face in r['active']]
            assert entry['faces'][str(face)]==ledger
            if face==60:assert len(ledger)==3 and all(r['kind']=='finite'for r in ledger)
            elif face==7:assert len(ledger)==11 and all(r['kind']=='finite'for r in ledger)
            elif face==33:assert len(ledger)==5 and all(r['kind']=='finite'for r in ledger)
        assert raw['counts']==dict(finite=32,exterior=6)


def check():
    d=read(TAG)
    for path,h in d['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    assert d['status']=='PASS_EXACT66_ESCAPE33_CAP'
    seed=read(A.TAG)
    assert all(seed['white_covectors'][str(i)]==enc([iv(0)]*6)for i in FRAME+[59])
    numerical_check(d['numerical']);assert d['physical']==enc(physical())
    assert d['physical']['summary']==dict(total_walls=66,old_finite_vertices_retained=203,new_finite_vertices=2,
        complete_finite_vertices=205,strict_dotted_contacts=56,right_contacts=9,infinite_volume_proved=True)
    assert d['scope']==dict(first_two_corner_duty_domain_complete=True,second_low_label_domain_complete=True,
        second_high_labels_unsearched=True,all_right_member_exact=True,other_members_numerical=True,
        compact33_section_numerical=True,full_ambient_facet33_compact=False,finite_volume_proved=False)
    print('PASS exact66-wall cap/205 finite;110 norm roots;28 independently enumerated numerical cones; strict remaining escape')

if __name__=='__main__':check()
