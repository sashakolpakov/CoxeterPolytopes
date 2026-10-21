#!/usr/bin/env python3
"""Complete exact local census of one-omitted rank-three cusp sections.

This does not certify the global ridges of the resulting H20 sections.
Positive roots of rational conics replace a bounded angular grid.
"""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import permutations
import hashlib
import json
from pathlib import Path
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_rank3_omitted_sections_20261014'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK = 'results/fresh_h21_source_corner_walk_20261013.json'
x, y, t = S.symbols('x y t')

def sha(p):
    return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()

def coeff(p, v=t):
    return [str(a) for a in S.Poly(p, v).all_coeffs()]

def interval_eval(p, lo, hi):
    a = b = S.Rational(0)
    for c in S.Poly(p, t).all_coeffs():
        vals = [a*lo, a*hi, b*lo, b*hi]
        a, b = min(vals)+c, max(vals)+c
    return a, b

def signed(p, minimal, interval):
    p = S.rem(p, minimal, t)
    if p == 0:
        return 0
    lo, hi = interval
    for digits in (15, 30, 60, 120, 240):
        a, b = interval_eval(p, lo, hi)
        if a > 0:
            return 1
        if b < 0:
            return -1
        lo, hi = S.Poly(minimal, t).refine_root(lo, hi, eps=S.Rational(1,10**digits))
    raise AssertionError(('unresolved exact sign', p, minimal, interval))

@lru_cache(None)
def conic_census(btuple, u, d):
    B = S.Matrix(3, 3, btuple)
    J = B.inv()
    X = S.Matrix([a*b[0] for a,b in zip((1,x,y),u)])
    Y = S.Matrix([a*b[1] for a,b in zip((1,x,y),u)])
    f = S.expand((X.T*J*Y)[0])
    g = S.expand((X.T*J*X-d*Y.T*J*Y)[0])
    rx = S.Poly(S.resultant(f,g,y),x)
    ry = S.Poly(S.resultant(f,g,x),y)
    assert not rx.is_zero and not ry.is_zero
    # x,y are strictly positive, so remove zero roots before radicalizing.
    px = rx.sqf_part().as_expr()
    py = ry.sqf_part().as_expr()
    if px.subs(x,0) == 0:
        px = S.cancel(px/x)
    if py.subs(y,0) == 0:
        py = S.cancel(py/y)
    base = [f,g,px,py]
    if S.groebner(base,x,y).is_zero_dimensional is False:
        assert list(S.groebner(base,x,y)) == [1]
        return {'conics':[str(f),str(g)], 'empty':True}, []
    for c in range(1,9):
        gb = list(S.groebner(base+[t-x-c*y],x,y,t))
        if len(gb) == 3 and gb[0].diff(x) == 1 and gb[1].diff(y) == 1:
            xp,yp,h = S.expand(x-gb[0]),S.expand(y-gb[1]),gb[2]
            if not (xp.free_symbols | yp.free_symbols | h.free_symbols)-{t}:
                break
    else:
        raise AssertionError(('no separating primitive element',B,u,d,gb))
    assert S.Poly(h,t).degree() <= 4
    certificates = {'conics':[str(f),str(g)], 'primitive_element':f'x+{c}*y',
                    'x_polynomial':coeff(xp), 'y_polynomial':coeff(yp),
                    'eliminant':coeff(h), 'empty':False}
    accepted=[]
    for p, multiplicity in S.factor_list(h,t)[1]:
        assert multiplicity == 1
        p = S.Poly(p,t).monic().as_expr()
        red=lambda z: S.rem(S.expand(z),p,t)
        xx,yy=red(xp),red(yp)
        assert red(f.subs({x:xx,y:yy}, simultaneous=True)) == 0
        assert red(g.subs({x:xx,y:yy}, simultaneous=True)) == 0
        for (lo,hi), mult in S.Poly(p,t).intervals(eps=S.Rational(1,10**15)):
            assert mult==1
            interval=(lo,hi)
            sg=lambda z:signed(z,p,interval)
            if sg(xx)<=0 or sg(yy)<=0:
                continue
            XX=X.subs({x:xx,y:yy}, simultaneous=True)
            YY=Y.subs({x:xx,y:yy}, simultaneous=True)
            q=red((XX.T*J*XX)[0])
            assert sg(q)>0
            R=(q*B-XX*XX.T-d*YY*YY.T).applyfunc(red)
            assert all(red(R[i,j]*R[k,l]-R[i,l]*R[k,j])==0
                       for i in range(3) for j in range(3)
                       for k in range(3) for l in range(3))
            index=next(i for i in range(3) if sg(R[i,i])>0)
            a=R[:,index]
            k=(J*a).applyfunc(red)
            if sg(k[0])<0:
                a,k=-a,-k
            assert all(sg(v)>0 for v in k)
            norm=red((a.T*J*a)[0])
            assert sg(norm)>0
            assert all(red(norm*R[i,j]-q*a[i]*a[j])==0
                       for i in range(3) for j in range(3))
            accepted.append({'minimal':p,'interval':interval,'x':xx,'y':yy,
                             'a':a,'k':k,'norm':norm})
    return certificates, accepted

def run():
    raw=json.loads((ROOT/RAW).read_text())['coxeter_input']
    Q=-S.Matrix(raw['basis_gram'])
    W=S.Matrix(raw['wall_coordinates_as_columns_transposed'])
    G=W*Q*W.T
    walk=json.loads((ROOT/WALK).read_text())
    assert walk['sources'][RAW]==sha(RAW) and not walk['pending_orbit_ids']
    targets=[]
    for kind,d,us in [('A2',3,((2,0),(-1,1),(-1,-1))),
                      ('C2',1,((1,0),(0,1),(-1,-1))),
                      ('G2',3,((1,0),(0,1),(-1,-1)))]:
        seen=set()
        for u in permutations(us):
            norms=[a*a+d*b*b for a,b in u]
            signature=tuple(S.Rational((u[i][0]*u[j][0]+d*u[i][1]*u[j][1])**2,
                                       norms[i]*norms[j]) for i,j in ((0,1),(0,2),(1,2)))
            if signature in seen:
                continue
            seen.add(signature)
            targets.append((kind,d,u,signature))
    assert len(targets)==10
    rows=[]; systems=[]; factors=[]; counts=Counter()
    for cusp in walk['corner_orbits']:
        if cusp['type']!='ideal':
            continue
        for factor in cusp['affine_components']:
            ids=factor['walls']
            if len(ids)!=4:
                continue
            factors.append({'cusp':cusp['id'],'walls':ids})
            marks=list(map(S.Rational,factor['positive_kernel']))
            D=W[ids,:].T*S.Matrix(marks)
            assert (D.T*Q*D)[0]==0
            for omit in range(4):
                rr=[i for i in range(4) if i!=omit]
                retained=[ids[i] for i in rr]
                B=G.extract(retained,retained)
                assert B.is_positive_definite
                for kind,d,u,signature in targets:
                    cert,candidates=conic_census(tuple(B),u,d)
                    sid=len(systems)
                    systems.append({'cusp':cusp['id'],'factor':ids,'omitted_wall':ids[omit],
                                    'retained_walls':retained,'target':kind,
                                    'squared_cosines':list(map(str,signature)),
                                    'planar_vectors':u,'y_square_factor':d,
                                    'positive_directions':len(candidates),**cert})
                    for a in candidates:
                        p=a['minimal']; interval=a['interval']
                        red=lambda z:S.rem(S.expand(z),p,t)
                        sg=lambda z:signed(z,p,interval)
                        e=(W[retained,:].T*a['k']).applyfunc(red)
                        assert red((e.T*Q*e)[0]-a['norm'])==0
                        assert red((e.T*Q*D)[0])==0
                        quot=[red(a['k'][j]/marks[i]) for j,i in enumerate(rr)]
                        low=0
                        for j in range(1,3):
                            if sg(quot[j]-quot[low])<0:
                                low=j
                        rows.append({'system':sid,'minimal_polynomial':coeff(p),
                                     'root_interval':list(map(str,interval)),
                                     'x_polynomial':coeff(a['x']),'y_polynomial':coeff(a['y']),
                                     'retained_pairings':[coeff(v) for v in a['a']],
                                     'positive_kernel':[coeff(v) for v in a['k']],
                                     'normal_constant':[coeff(v) for v in e],
                                     'normal_square':coeff(a['norm']),
                                     'null_shift':list(map(str,D)),
                                     'parameter_lower_magnitude':coeff(quot[low]),
                                     'lower_bound_witness_wall':retained[low]})
                        counts[kind]+=1
                        counts['field_degree_'+str(S.degree(p,t))]+=1
    assert len(factors)==3 and len(systems)==120
    return {'status':'PASS_EXACT_RANK_THREE_LOCAL_DIRECTION_CENSUS',
            'sources':{p:sha(p) for p in [RAW,WALK,f'scripts/{TAG}.py']},
            'factors':factors,'system_count':len(systems),'direction_count':len(rows),
            'counts':dict(counts),'systems':systems,'directions':rows,
            'scope':'All irreducible affine triangular links obtained by omitting exactly one facet of one rank-three affine factor of the literal source cusps. Directions carry exact real algebraic representations. Global H20 Coxeter conditions are not checked. The local parameter interval is [-L,0), with degeneracies at its lower endpoint requiring separate treatment.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert json.loads(p.read_text())==json.loads(json.dumps(out))
    else:
        p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print(out['direction_count'],out['counts'])
