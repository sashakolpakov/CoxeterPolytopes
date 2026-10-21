#!/usr/bin/env python3
"""Independent saved-certificate audit using explicit xy-saturation.

No imports from the producer.  SymPy's exact rational Groebner/root-isolation
implementation is a shared external dependency, not separate CAS arithmetic.
"""
import argparse
from collections import Counter
from itertools import permutations
import hashlib
import json
from pathlib import Path
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_rank3_omitted_sections_20261014'
SOURCE = f'results/{TAG}.json'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK = 'results/fresh_h21_source_corner_walk_20261013.json'
x,y,t,z = S.symbols('x y t z')

def sha(p):
    return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

def poly(c):
    return S.Poly.from_list(list(map(S.Rational,c)),t).as_expr()

def enclosure(p,lo,hi):
    a=b=S.Rational(0)
    for c in S.Poly(p,t).all_coeffs():
        v=(a*lo,a*hi,b*lo,b*hi)
        a,b=min(v)+c,max(v)+c
    return a,b

def sign(p,m,iv):
    p=S.rem(S.expand(p),m,t)
    if p==0:return 0
    lo,hi=iv
    for precision in (20,40,80,160,320):
        a,b=enclosure(p,lo,hi)
        if a>0:return 1
        if b<0:return -1
        lo,hi=S.Poly(m,t).refine_root(lo,hi,eps=S.Rational(1,10**precision))
    raise AssertionError(('unresolved sign',p,m,iv))

def run():
    saved=json.loads((ROOT/SOURCE).read_text())
    for p,h in saved['sources'].items():assert sha(p)==h
    raw=json.loads((ROOT/RAW).read_text())['coxeter_input']
    Q=-S.Matrix(raw['basis_gram'])
    W=S.Matrix(raw['wall_coordinates_as_columns_transposed'])
    G=W*Q*W.T
    walk=json.loads((ROOT/WALK).read_text())
    assert not walk['pending_orbit_ids'] and walk['sources'][RAW]==sha(RAW)
    factors={(c['id'],tuple(a['walls'])):a for c in walk['corner_orbits']
             if c['type']=='ideal' for a in c['affine_components'] if len(a['walls'])==4}
    assert len(factors)==3
    target_signatures={}
    for kind,d,rays in [('A2',3,((2,0),(-1,1),(-1,-1))),
                        ('C2',1,((1,0),(0,1),(-1,-1))),
                        ('G2',3,((1,0),(0,1),(-1,-1)))]:
        signatures=set()
        for u in permutations(rays):
            n=[a*a+d*b*b for a,b in u]
            signatures.add(tuple(S.Rational((u[i][0]*u[j][0]+d*u[i][1]*u[j][1])**2,n[i]*n[j])
                                 for i,j in ((0,1),(0,2),(1,2))))
        target_signatures[kind]=signatures
    assert sum(map(len,target_signatures.values()))==10
    seen=set();positive=Counter();field_counts=Counter();saturation_empty=0
    audit_rows=[]
    for sid,r in enumerate(saved['systems']):
        factor=factors[(r['cusp'],tuple(r['factor']))]
        ids=r['retained_walls'];omit=r['omitted_wall']
        assert len(ids)==3 and sorted(ids+[omit])==sorted(r['factor'])
        signature=tuple(map(S.Rational,r['squared_cosines']))
        assert signature in target_signatures[r['target']]
        key=(r['cusp'],tuple(r['factor']),omit,r['target'],signature)
        assert key not in seen;seen.add(key)
        B=G.extract(ids,ids);J=B.inv()
        assert all(B[:k,:k].det()>0 for k in (1,2,3))
        u=r['planar_vectors'];d=r['y_square_factor']
        X=S.Matrix([q*a[0] for q,a in zip((1,x,y),u)])
        Y=S.Matrix([q*a[1] for q,a in zip((1,x,y),u)])
        f=S.expand((X.T*J*Y)[0]);g=S.expand((X.T*J*X-d*Y.T*J*Y)[0])
        assert [f,g]==[S.sympify(p,locals={'x':x,'y':y})for p in r['conics']]
        rows=[a for a in saved['directions']if a['system']==sid]
        if r['empty']:
            # Direct ideal saturation, without either producer resultant.
            assert list(S.groebner([f,g,z*x*y-1],z,x,y))==[1]
            assert not rows and r['positive_directions']==0
            saturation_empty+=1;continue
        c=S.sympify(r['primitive_element'],locals={'x':x,'y':y}).coeff(y)
        assert c>0
        xp,yp,h=map(poly,(r['x_polynomial'],r['y_polynomial'],r['eliminant']))
        # Independent elimination starts from the ORIGINAL conics and an
        # inverse-xy variable, rather than radicalizing their resultants.
        gb=list(S.groebner([f,g,z*x*y-1,t-x-c*y],z,x,y,t))
        without_z=[a for a in gb if z not in a.free_symbols]
        assert len(without_z)==3
        assert without_z[0].diff(x)==1 and without_z[1].diff(y)==1
        sx=S.expand(x-without_z[0]);sy=S.expand(y-without_z[1]);sh=without_z[2]
        assert not (sx.free_symbols|sy.free_symbols|sh.free_symbols)-{t}
        # Some original intersections are tangent (nonreduced).  In shape
        # position their radical is obtained just by squarefreeing the
        # univariate eliminant and reducing the coordinate polynomials.
        radical=S.Poly(sh,t).sqf_part().monic().as_expr()
        assert radical==h and S.rem(sx-xp,h,t)==S.rem(sy-yp,h,t)==0
        assert S.gcd(S.Poly(h,t),S.Poly(h,t).diff()).degree()==0
        accepted=[]
        for m,mult in S.factor_list(h,t)[1]:
            assert mult==1
            m=S.Poly(m,t).monic().as_expr()
            assert S.Poly(m,t).is_irreducible
            xx,yy=S.rem(xp,m,t),S.rem(yp,m,t)
            for iv,rootmult in S.Poly(m,t).intervals(eps=S.Rational(1,10**12)):
                assert rootmult==1
                if sign(xx,m,iv)>0 and sign(yy,m,iv)>0:accepted.append((m,iv))
        assert len(accepted)==len(rows)==r['positive_directions']
        matched=set()
        for row in rows:
            m=poly(row['minimal_polynomial']);iv=tuple(map(S.Rational,row['root_interval']))
            assert S.Poly(m,t).count_roots(*iv)==1
            options=[k for k,(mm,jv)in enumerate(accepted)
                     if mm==m and max(jv[0],iv[0])<=min(jv[1],iv[1])
                     and S.Poly(m,t).count_roots(max(jv[0],iv[0]),min(jv[1],iv[1]))==1]
            assert len(options)==1 and options[0]not in matched;matched.add(options[0])
            red=lambda a:S.rem(S.expand(a),m,t)
            sg=lambda a:sign(a,m,iv)
            xx,yy=poly(row['x_polynomial']),poly(row['y_polynomial'])
            assert red(xx-xp)==red(yy-yp)==0 and sg(xx)>0 and sg(yy)>0
            XX=X.subs({x:xx,y:yy},simultaneous=True);YY=Y.subs({x:xx,y:yy},simultaneous=True)
            q=red((XX.T*J*XX)[0]);assert sg(q)>0
            a=S.Matrix(list(map(poly,row['retained_pairings'])))
            k=S.Matrix(list(map(poly,row['positive_kernel'])))
            assert all(red(v)==0 for v in B*k-a) and all(sg(v)>0 for v in k)
            norm=poly(row['normal_square']);assert sg(norm)>0 and red((k.T*B*k)[0]-norm)==0
            # Projected Gram identity proves the exact affine-triangle shape.
            assert all(red(norm*(q*B-XX*XX.T-d*YY*YY.T)[i,j]-q*a[i]*a[j])==0
                       for i in range(3)for j in range(3))
            e=S.Matrix(list(map(poly,row['normal_constant'])))
            assert all(red(v)==0 for v in e-W[ids,:].T*k)
            marks=list(map(S.Rational,factor['positive_kernel']))
            D=W[r['factor'],:].T*S.Matrix(marks)
            assert D==S.Matrix(list(map(S.Rational,row['null_shift'])))
            assert (D.T*Q*D)[0]==0 and all((W[i,:]*Q*D)[0]==0 for i in r['factor'])
            assert red((e.T*Q*e)[0]-norm)==red((e.T*Q*D)[0])==0
            L=poly(row['parameter_lower_magnitude']);assert sg(L)>0
            quot=[red(k[j]/marks[r['factor'].index(i)])for j,i in enumerate(ids)]
            assert all(sg(v-L)>=0 for v in quot)
            assert red(quot[ids.index(row['lower_bound_witness_wall'])]-L)==0
            positive[r['target']]+=1;field_counts[str(S.degree(m,t))]+=1
        audit_rows.append({'system':sid,'saturated_eliminant_degree':S.degree(h,t),
                           'positive_real_solutions':len(accepted)})
    assert len(seen)==120 and sum(positive.values())==91
    for (cid,fs)in factors:
        for omit in fs:
            for kind,signatures in target_signatures.items():
                assert all((cid,fs,omit,kind,ss)in seen for ss in signatures)
    assert dict(positive)=={k:saved['counts'][k]for k in positive}
    paths=[SOURCE,RAW,WALK,f'scripts/{TAG}.py',str(Path(__file__).resolve().relative_to(ROOT))]
    return {'status':'PASS_INDEPENDENT_SATURATED_CONIC_AND_DIRECTION_REVIEW',
            'sources':{p:sha(p)for p in paths},'systems':len(seen),
            'saturated_empty_systems':saturation_empty,'positive_directions':dict(positive),
            'field_degree_counts':dict(field_counts),'nonempty_system_audits':audit_rows,
            'scope':'Independent original-conic xy-saturation, complete real-root coverage, '
                    'all 91 saved direction/parameter identities and source-factor/target coverage. '
                    'No producer imports. Shared SymPy rational CAS; no global section or ridge claim.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run()
    result=json.loads(json.dumps(result,default=lambda a:int(a)))
    path=ROOT/'results'/f'{TAG}_review.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status']);print(result['positive_directions'],result['field_degree_counts'])
