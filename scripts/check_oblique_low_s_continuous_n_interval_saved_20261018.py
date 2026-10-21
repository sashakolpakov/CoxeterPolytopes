#!/usr/bin/env python3
"""Standalone physical interval replay for low-S continuous-N event exclusion.

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
TAG = 'fresh_oblique_low_s_continuous_n_interval_20261018'
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


def canonical(g):
    return min(tuple(g[p[i]][p[j]]for i in range(len(g))for j in range(i))for p in permutations(range(len(g))))


@lru_cache(None)
def finite_corner(p0,p1,p2,m01,m02,m12):
    # Complete spherical connected types of order at most5; the four fixed
    # H4 whites form a separate positive component and are omitted here.
    g=[[2]*5 for _ in range(5)];g[0][1]=g[1][0]=3
    for i,p in enumerate((p0,p1,p2),2):
        if p:g[(p-1)//3][i]=g[i][(p-1)//3]=3+(p-1)%3
    for(i,j),m in zip(((2,3),(2,4),(3,4)),(m01,m02,m12)):g[i][j]=g[j][i]=m
    left=set(range(5))
    while left:
        seen=set();todo=[min(left)]
        while todo:
            x=todo.pop()
            if x in seen:continue
            seen.add(x);todo.extend(y for y in left if g[x][y]>2 and y not in seen)
        left-=seen;n=len(seen)
        if n<=2:continue
        ids=sorted(seen);gg=[[g[i][j]for j in ids]for i in ids]
        if sum(gg[i][j]>2 for i in range(n)for j in range(i))!=n-1:return False
        templates=[]
        edges=[(i,i+1,3)for i in range(n-1)];templates.append(edges)
        templates.append([(i,j,4 if i==0 else m)for i,j,m in edges])
        if n in(3,4):templates.append([(i,j,5 if i==0 else m)for i,j,m in edges])
        if n==4:templates.extend([[(0,1,3),(1,2,4),(2,3,3)],[(0,1,3),(0,2,3),(0,3,3)]])
        if n==5:templates.append([(0,1,3),(1,2,3),(2,3,3),(2,4,3)])
        codes=[]
        for es in templates:
            t=[[2]*n for _ in range(n)]
            for i,j,m in es:t[i][j]=t[j][i]=m
            codes.append(canonical(t))
        if canonical(gg)not in codes:return False
    return True


def check():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    prior=json.loads((ROOT/f'results/{PRIOR}.json').read_text())
    num=json.loads((ROOT/'results/fresh_oblique_low_s_continuous_n_events_20261018.json').read_text())
    pri=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    assert d['status']=='PASS_PHYSICAL_LOW_S_LOW_R_CONTINUOUS_HIGH_N_EXCLUSION'and d['interval_scale']==SCALE
    assert len(d['records'])==len(d['parents'])==len(num['records'])==2 and d['parents']==prior['parents']
    ids=WHITE+OLD;metric=[list(map(iv,r))for r in pri['computational_gram']];center=list(map(iv,pri['common_center']));normals={}
    for row in pri['facet_witnesses']:
        if row['wall']in ids:
            hp=iv(pri['center_products'][str(row['wall'])]);assert hp[1]<0
            normals[row['wall']]=[div(mul(iv(2),sub(x,y)),hp)for x,y in zip(center,map(iv,row['point']))]
    native={a:{b:dot(normals[a],metric,normals[b])for b in ids}for a in ids}
    wi,pivots=gauss([[native[a][b]for b in WHITE]for a in WHITE]);assert all(p[0]>0 for p in pivots)
    wc0={a:[native[a][w]for w in WHITE]for a in OLD}
    g0={a:{b:sub(native[a][b],dot(wc0[a],wi,wc0[b]))for b in OLD}for a in OLD}
    assert prior['primary_projection']==enc(dict(walls=ids,physical_Gram=native,white_inverse=wi,white_pivots=pivots,old_white_covectors=wc0,projected_Gram=g0))
    # Original Coxeter norm/finite-label identities only; dotted products
    # retain their original enclosures. No new numerical root is rounded.
    seen=set()
    for proof in d['refined_primary']['identities']:
        i,j=proof['pair'];assert (i,j)not in seen and i in ids and j in ids;seen.add((i,j))
        z=native[i][j];assert proof['old']==enc(z);label=proof['identity']
        if label=='norm2':assert i==j;exact=iv(2)
        elif label==2:assert i!=j and -1<z[0]<=0<=z[1]<1;exact=iv(0)
        else:
            assert i!=j and isinstance(label,int)and label>=3 and -theta(label+1)[0]<z[0]<=z[1]<-theta(label-1)[1]
            exact=neg(theta(label))
        intersection=max(z[0],exact[0]),min(z[1],exact[1]);assert intersection[0]<=intersection[1]
        assert proof['identity_interval']==enc(exact)and proof['intersection']==enc(intersection)
        native[i][j]=native[j][i]=intersection
    assert all((i,i)in seen for i in ids)
    wi,pivots=gauss([[native[a][b]for b in WHITE]for a in WHITE]);assert all(p[0]>0 for p in pivots)
    wc0={a:[native[a][w]for w in WHITE]for a in OLD};g0={a:{b:sub(native[a][b],dot(wc0[a],wi,wc0[b]))for b in OLD}for a in OLD}
    expect=dict(physical_Gram=native,white_inverse=wi,white_pivots=pivots,old_white_covectors=wc0,projected_Gram=g0)
    assert all(d['refined_primary'][k]==enc(v)for k,v in expect.items())
    profiles=[[iv(0)]*6]
    for pos in(4,5):
        for m in(3,4,5):
            v=[iv(0)]*6;v[pos]=neg(theta(m));profiles.append(v)
    assert d['profiles']==prior['profiles']==enc(profiles)
    assert d['parameter_domain']==enc([theta(6)[0],F(2)])

    def extend(g,wc,row,old=False):
        wall,frame,pi,labels=(row[k]for k in('wall','frame','profile','labels'))
        gram=[[g[a][b]for b in frame]for a in frame];inv,det=inverse(gram)
        assert det[1]<0 and gram[0][0][0]>0
        rhs=[sub(neg(theta(m)),dot(profiles[pi],wi,wc[i]))for i,m in zip(frame,labels)]
        square=sub(iv(2),dot(profiles[pi],wi,profiles[pi]));rem=sub(square,dot(rhs,inv,rhs))
        w,sign=row['witness'],row['witness_sign'];assert w in g and w not in frame and sign in(-1,1)
        cv=[g[i][w]for i in frame];q=sub(g[w][w],dot(cv,inv,cv));assert q[0]>0 and rem[0]>0
        linear=dot(rhs,inv,cv);radical=sqrt(mul(rem,q));branches=[sub(linear,radical),add(linear,radical)]
        factor=mul(iv(sign),sqrt(div(rem,q)));out={}
        for j in g:
            cj=[g[i][j]for i in frame]
            out[j]=rhs[frame.index(j)]if j in frame else add(dot(rhs,inv,cj),mul(factor,sub(g[w][j],dot(cv,inv,cj))))
        expected=dict(frame_Gram=gram,inverse=inv,determinant=det,rhs=rhs,projected_square=square,remainder=rem,
                      witness_residual=q,branch_witness_products=branches,new_products=out)
        assert all(row[k]==enc(v)for k,v in expected.items())
        if not old:
            distances=[abs(row['numerical_witness_product']-float(sum(z)/2))for z in branches]
            assert distances==row['branch_comparison']and distances[(sign+1)//2]<distances[(1-sign)//2]/100
        for j,z in out.items():g[j][wall]=z
        g[wall]={**out,wall:square};wc[wall]=profiles[pi]

    def reject(g,wc,w,proof,allowed):
        j=proof['wall'];assert j in allowed and j not in(63,64)
        value=add(g[w][j],dot(wc[w],wi,wc[j]));assert proof['product']==enc(value)
        kind=proof['kind']
        if kind=='positive':assert value[0]>0
        elif kind=='between_right_and_simple':assert -1<value[0]<=value[1]<0
        else:
            assert kind=='between_consecutive_native_labels';a,b=proof['labels'];assert b==a+1 and a>=3
            lo,hi=theta(a),theta(b)
            assert proof['lower_magnitude']==enc(lo)and proof['upper_magnitude']==enc(hi)
            assert -hi[0]<value[0]<=value[1]<-lo[1]

    count=dict(S=0,SR=0,R=0,RR=0,E=0,roots=0,outside=0,gaps=0);margins=[]
    for row,parent,nr in zip(d['records'],d['parents'],num['records']):
        assert row['index']==parent['index']==nr['index']and row['instruction_key']==parent['instruction_key']==nr['instruction_key']and row['origins']==nr['origins']
        g={i:r.copy()for i,r in g0.items()};wc=wc0.copy()
        assert len(parent['steps'])==len(row['refined_parent_steps'])==5
        for oldstep,step in zip(parent['steps'],row['refined_parent_steps']):
            for k in('wall','frame','profile','labels','witness','witness_sign'):assert step[k]==oldstep[k]
            extend(g,wc,step,True)
        del g[70],wc[70]
        for r in g.values():del r[70]
        retained=[i for i in g if i not in(63,64)]
        lprofile=int(nr['prefix_profiles']['68']);assert lprofile in(1,4)
        for i,j in((52,68),(68,71),(69,71)):
            native_contact=add(g[i][j],dot(wc[i],wi,wc[j]));assert -1<native_contact[0]<=0<=native_contact[1]<1
        srecipes=[(p,list(ms))for p in range(7)for ms in product((2,3,4,5),repeat=3)
            if finite_corner(0,lprofile,p,2,ms[0],ms[1])and finite_corner(lprofile,0,p,2,ms[1],ms[2])]
        assert[(q['profile'],q['labels'])for q in nr['S_bank']['recipes']]==srecipes and len(srecipes)==42
        want=[(sq,sb)for sq in nr['S_bank']['recipes']for sb in sq['branches']]
        assert len(want)==len(row['S'])==84
        ssigns={}
        for sr,(sq,sb)in zip(row['S'],want):
            count['S']+=1
            assert(sr['recipe'],sr['branch'],sr['profile'],sr['labels'])==(sq['index'],sb['branch'],sq['profile'],sq['labels'])
            gs={i:r.copy()for i,r in g.items()};ws=wc.copy();ext=sr['extension']
            assert[ext[k]for k in('wall','frame','profile','labels')]==[72,[52,68,71],sq['profile'],sq['labels']]
            extend(gs,ws,ext);ssigns.setdefault(sq['index'],[]).append((ext['witness'],ext['witness_sign']))
            if sr['blocker']:
                reject(gs,ws,72,sr['blocker'],retained);count['SR']+=1;assert sr['R']==[];continue
            si=sr['numerical_cap_index'];assert si==sb['cap'];bank=nr['R_banks'][si]
            rrecipes=[(p,[2,m,r])for p,m,r in product(range(7),(2,3,4,5),(2,3,4,5))
                if finite_corner(0,0,p,2,2,m)and finite_corner(0,sq['profile'],p,sq['labels'][2],m,r)]
            assert[(q['profile'],q['labels'])for q in bank['recipes']]==rrecipes
            rwant=[(rq,rb)for rq in bank['recipes']for rb in rq['branches']];assert len(rwant)==len(sr['R'])
            rsigns={}
            for rr,(rq,rb)in zip(sr['R'],rwant):
                count['R']+=1
                assert(rr['recipe'],rr['branch'],rr['profile'],rr['labels'])==(rq['index'],rb['branch'],rq['profile'],rq['labels'])
                gr={i:r.copy()for i,r in gs.items()};wr=ws.copy();ext=rr['extension']
                assert[ext[k]for k in('wall','frame','profile','labels')]==[70,[69,71,72],rq['profile'],rq['labels']]
                extend(gr,wr,ext);rsigns.setdefault(rq['index'],[]).append((ext['witness'],ext['witness_sign']))
                if rr['blocker']:
                    reject(gr,wr,70,rr['blocker'],retained+[72]);count['RR']+=1;assert rr['events']==[];continue
                assert len(rr['events'])==2
                for a,ev in enumerate(rr['events']):
                    count['E']+=1;assert ev['a']==a and ev['N_frame']==[10,67,69]
                    frame=ev['N_frame'];inv,det=inverse([[gr[i][j]for j in frame]for i in frame])
                    r0=[iv(-a),iv(0),iv(0)];r1=[iv(0),iv(0),iv(-1)];cv=[gr[i][64]for i in frame];rv=[gr[i][70]for i in frame]
                    sig=sub(gr[64][64],dot(cv,inv,cv));assert sig[0]>0
                    aa,bb=dot(r0,inv,rv),dot(r1,inv,rv);cc=div(sub(gr[64][70],dot(cv,inv,rv)),sqrt(sig))
                    rho=[sub(iv(2),dot(r0,inv,r0)),neg(mul(iv(2),dot(r0,inv,r1))),neg(dot(r1,inv,r1))]
                    p=[sub(mul(aa,aa),mul(mul(cc,cc),rho[0])),sub(mul(iv(2),mul(aa,bb)),mul(mul(cc,cc),rho[1])),sub(mul(bb,bb),mul(mul(cc,cc),rho[2]))]
                    disc=sub(mul(p[1],p[1]),mul(iv(4),mul(p[2],p[0])));assert p[2][0]>0 and disc[0]>0
                    expect=dict(N_inverse=inv,N_determinant=det,old64_residual_square=sig,coefficients=[aa,bb,cc],norm_remainder_coefficients=rho,polynomial=p,discriminant=disc)
                    assert all(ev[k]==enc(v)for k,v in expect.items())
                    roots=[div(sub(neg(p[1]),sqrt(disc)),mul(iv(2),p[2])),div(add(neg(p[1]),sqrt(disc)),mul(iv(2),p[2]))]
                    assert len(ev['roots'])==2
                    for rt,decision in zip(roots,ev['roots']):
                        count['roots']+=1;assert decision['root']==enc(rt)
                        if decision['kind']=='below_high_domain':assert rt[1]<theta(6)[0];count['outside']+=1
                        elif decision['kind']=='above_finite_domain':assert rt[0]>2;count['outside']+=1
                        else:
                            assert decision['kind']=='strict_native_label_gap';lo,hi=decision['labels'];assert hi==lo+1 and(lo,hi)in((8,9),(10,11))
                            lb,ub=theta(lo),theta(hi);assert decision['lower_magnitude']==enc(lb)and decision['upper_magnitude']==enc(ub)
                            gaps=[rt[0]-lb[1],ub[0]-rt[1]];assert min(gaps)>0
                            assert[decision['lower_gap'],decision['upper_gap']]==enc(gaps);margins.extend(gaps);count['gaps']+=1
            assert all(len({w for w,sgn in v})==1 and sorted(sgn for w,sgn in v)==[-1,1]for v in rsigns.values())and set(rsigns)==set(range(len(bank['recipes'])))
        assert all(len({w for w,sgn in v})==1 and sorted(sgn for w,sgn in v)==[-1,1]for v in ssigns.values())and set(ssigns)==set(range(42))
    assert count==dict(S=168,SR=92,R=2476,RR=2416,E=120,roots=240,outside=236,gaps=4)
    expected=dict(parents=2,parent_origins=4,S_recipes=84,S_branches=168,S_strict_rejections=92,S_overdomain_survivors=76,
        R_recipes=1238,R_branches=2476,R_strict_rejections=2416,R_overdomain_survivors=60,event_quadratics=120,event_roots=240,
        outside_domain=236,native_label_gaps=4,minimum_native_gap=min(margins),unresolved=0,accepted=0,new_cones=0)
    assert d['summary']==enc(expected)
    assert d['deleted_walls']==num['deleted_walls']and d['exact_deletion_debit']==num['exact_deletion_debit']and d['topology']==num['topology']
    print('PASS_PHYSICAL_CONTINUOUS_N_EVENT_REPLAY',count,'minimum_gap',float(min(margins)),flush=True)


if __name__=='__main__':check()

