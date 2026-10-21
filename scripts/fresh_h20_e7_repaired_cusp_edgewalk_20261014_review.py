#!/usr/bin/env python3
"""Independent all-parameter review of four cusp-edge cap obstructions.

Rebuilds each physical dual basis and its finite-extension census. Uses the
separate sparse squarefree arithmetic, not the producer's tower arithmetic,
numerical intervals, or profile enumeration.
"""
import argparse
import itertools
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_e7_replacement_escape_20261014_review as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e7_repaired_cusp_edgewalk_20261014'
BASE = 'results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
WALK = 'results/fresh_h20_e7_repaired_cusp_walk_20261014.json'
PRIMARY = f'results/{TAG}_exact.json'
RAD = [1,2,79,158,11,22,869,1738,3,6,237,474,33,66,2607,5214]
RAD32 = RAD + [5*r for r in RAD]

def sha(p):
    return sha256((ROOT/p).read_bytes()).hexdigest()

def read(p):
    d = json.loads((ROOT/p).read_text())
    for s,h in d.get('sources',{}).items():
        assert sha(s) == h, s
    return d

def run():
    d,w,p = map(read,[BASE,WALK,PRIMARY])
    def base(a):
        assert len(a) in [16,32]
        return {r:F(x) for r,x in zip(RAD32,a) if F(x)}
    D = base(d['field']['discriminant'])
    assert R.sign(D) > 0
    zero = ({},{})
    def num(q): return (R.c(q),{})
    def add(a,b): return (R.add(a[0],b[0]),R.add(a[1],b[1]))
    def scale(a,q): return (R.scale(a[0],q),R.scale(a[1],q))
    def mul(a,b):
        return (R.add(R.mul(a[0],b[0]),R.mul(D,R.mul(a[1],b[1]))),
                R.add(R.mul(a[0],b[1]),R.mul(a[1],b[0])))
    def sign(a):
        sa,sb = R.sign(a[0]),R.sign(a[1])
        if not sb: return sa
        if not sa or sa == sb: return sb
        return sa*R.sign(R.sub(R.mul(a[0],a[0]),R.mul(D,R.mul(a[1],a[1]))))
    def total(xs):
        z = zero
        for x in xs: z = add(z,x)
        return z
    def parse(a): return (base(a[0]),base(a[1]))
    def vec(a): return tuple(map(parse,a))
    def rat(a):
        assert not a[1]
        return R.rational(a[0])
    Q = [list(map(F,row)) for row in d['form']]
    def cov(v):
        return tuple(total(scale(v[j],q) for j,q in enumerate(row) if q) for row in Q)
    def rawdot(v,c): return total(mul(x,y) for x,y in zip(v,c))
    def dot(v,z): return rawdot(v,cov(z))
    def lin(vs,cs):
        return tuple(total(mul(v[k],c) for v,c in zip(vs,cs)) for k in range(21))
    one,two = num(1),num(2)
    sq2,sq3,sq5 = [({i:F(1)},{}) for i in [2,3,5]]
    magnitude = {3:one,4:sq2,5:scale(add(one,sq5),F(1,2)),6:sq3}
    N = list(map(vec,d['normals']))
    CN = list(map(cov,N))
    H,E,nu = map(vec,[d['strict_center'],d['strict_escape']['point'],w['cusp_point']])
    assert len(N) == 95 and all(dot(n,n) == num(F(q)) for n,q in zip(N,d['norms']))
    assert nu == lin([N[9],N[88]],[one,one])
    assert dot(nu,nu) == zero and sign(dot(H,H)) < 0 and sign(dot(E,E)) > 0
    assert sign(dot(H,E)) < 0
    h,en = dot(H,nu),dot(E,nu)
    assert sign(h) < 0 and sign(en) < 0
    nn = [rawdot(nu,c) for c in CN]
    assert all(sign(x) <= 0 for x in nn)
    assert all(sign(rawdot(H,c)) < 0 and sign(rawdot(E,c)) < 0 for c in CN)
    finite = {(r['edge'],tuple(map(tuple,r['contacts']))):r for r in p['finite_profile_certificates']}
    tails = {(r['edge'],r['contact_wall']):r for r in p['I2_tail_certificates']}
    assert len(finite) == len(p['finite_profile_certificates']) == 418
    assert len(tails) == len(p['I2_tail_certificates']) == 15
    seen_finite,seen_tails = set(),set()
    counts = Counter()
    records = []
    for edge in p['edges']:
        ei = edge['edge']
        germ = w['edge_germs'][ei]
        ids = germ['support']
        assert germ['type'] == 'open_opposite_end' and ids == edge['support']
        assert germ['omitted'] == edge['omitted'] and len(ids) == 19
        assert all(dot(N[i],N[i]) == two for i in ids)
        scales = [sq2 if i in [17,79] else one for i in ids]
        roots = [lin([N[i]],[s]) for i,s in zip(ids,scales)]
        M = [[rat(dot(a,b)) for b in roots] for a in roots]
        R.pd(M)
        I = R.inverse(M)
        dual = [lin(roots,[scale(scales[i],q) for q in I[i]]) for i in range(19)]
        assert all(dot(a,N[j]) == num(i==k) for i,a in enumerate(dual) for k,j in enumerate(ids))
        assert all(dot(a,nu) == zero for a in dual)
        W = lin([H]+dual,[one]+[scale(dot(H,N[i]),-1) for i in ids])
        w2 = dot(W,W)
        assert sign(w2) < 0 and dot(W,nu) == h
        assert all(dot(W,N[i]) == zero for i in ids)
        assert h == parse(edge['h']) and w2 == parse(edge['W_norm'])
        wn = [rawdot(W,c) for c in CN]
        we = dot(W,E)
        zn = [add(scale(mul(h,x),-2),mul(w2,y)) for x,y in zip(wn,nn)]
        ze = add(scale(mul(h,we),-2),mul(w2,en))
        assert ze == parse(edge['escape_null_pairing']) and sign(ze) < 0
        assert all(sign(zn[i]) < 0 for i in range(95) if i not in ids)
        # The inverse normalized Gram is entrywise nonnegative. These exact
        # inequalities give a complete profile census without classification:
        # four contacts cost >=2; two in one component cost >=2; any label
        # >=7 must be the sole contact and attach to an isolated A1.
        GI = [[mul(mul(scales[i],scales[j]),num(I[i][j])) for j in range(19)] for i in range(19)]
        assert all(sign(x) >= 0 for row in GI for x in row)
        diag = [rat(GI[i][i]) for i in range(19)]
        assert min(diag) >= F(1,2)
        components = []
        remaining = set(range(19))
        while remaining:
            cc = {min(remaining)}
            todo = list(cc)
            for i in todo:
                for j in sorted(remaining-cc):
                    if M[i][j]: cc.add(j); todo.append(j)
            components.append(sorted(cc))
            remaining -= cc
        which = {i:k for k,cc in enumerate(components) for i in cc}
        isolated = [cc[0] for cc in components if len(cc)==1]
        for cc in components:
            for i,j in itertools.combinations(cc,2):
                assert sign(add(num(diag[i]+diag[j]-2),scale(GI[i][j],2))) >= 0
            if len(cc)>1: assert min(diag[i] for i in cc) >= F(2,3)
        assert all(diag[i] == F(1,2) for i in isolated)
        # At three contacts, any label >=4 costs >= (2+1+1)/2 =2.
        profiles = {()}
        for size in [1,2,3]:
            for positions in itertools.combinations(range(19),size):
                if len({which[i] for i in positions}) < size: continue
                labels = itertools.product([3,4,5,6],repeat=size) if size<3 else [(3,3,3)]
                for mm in labels:
                    cost = total(scale(mul(magnitude[m],magnitude[m]),diag[i]) for i,m in zip(positions,mm))
                    if sign(add(two,scale(cost,-1))) > 0:
                        profiles.add(tuple((ids[i],m) for i,m in zip(positions,mm)))
        assert profiles == {key[1] for key in finite if key[0]==ei}
        assert {ids[i] for i in isolated} == {key[1] for key in tails if key[0]==ei}
        dual_pairings = [[rawdot(u,c) for c in CN]+[dot(u,E)] for u in dual]
        def numerator(pn,q,a,j):
            vn,ww = (en,we) if j==95 else (nn[j],wn[j])
            return total([scale(mul(h,pn),2*a),scale(mul(h,ww),2*a*a),mul(add(q,scale(w2,-a*a)),vn)])
        for contacts in sorted(profiles):
            key = (ei,contacts)
            row = finite[key]
            seen_finite.add(key)
            co = [(ids.index(i),scale(magnitude[m],-1)) for i,m in contacts]
            cost = total(mul(mul(x,y),GI[i][j]) for i,x in co for j,y in co)
            q = add(two,scale(cost,-1))
            assert sign(q)>0 and q==parse(row['residual_norm'])
            a = F(row['separating_parameter'])
            j = row['blocking_wall']
            assert a>0 and j not in ids and 0<=j<95
            pe = total(mul(x,dual_pairings[i][95]) for i,x in co)
            pn = total(mul(x,dual_pairings[i][j]) for i,x in co)
            ke,kn = numerator(pe,q,a,95),numerator(pn,q,a,j)
            assert ke==parse(row['escape_numerator']) and kn==parse(row['wall_numerator'])
            assert sign(ke)>0 and sign(kn)<0
            counts['finite_profiles'] += 1
            counts['bridge_blockers'] += (j==87)
        for i in isolated:
            key = (ei,ids[i])
            row = tails[key]
            seen_tails.add(key)
            a,j = F(row['separating_parameter']),row['blocking_wall']
            assert a>0 and j not in ids and 0<=j<95
            def polynomial(k):
                vn,ww = (en,we) if k==95 else (nn[k],wn[k])
                return (scale(vn,F(-1,2)),scale(mul(h,dual_pairings[i][k]),-2*a),
                        add(scale(mul(h,ww),2*a*a),mul(add(two,scale(w2,-a*a)),vn)))
            def at(poly,c):
                A,B,C=poly
                return total([mul(A,mul(c,c)),mul(B,c),C])
            wall,escape = polynomial(j),polynomial(95)
            assert wall==tuple(map(parse,row['wall_polynomial']))
            assert escape==tuple(map(parse,row['escape_polynomial']))
            wl,wr = at(wall,sq3),at(wall,two)
            assert sign(wall[0])>=0 and sign(wl)<0 and sign(wr)<0
            assert [wl,wr]==list(map(parse,row['wall_endpoint_values']))
            A,B,C = escape
            el,er = at(escape,sq3),at(escape,two)
            dl,dr = add(scale(mul(A,sq3),2),B),add(scale(A,4),B)
            assert sign(A)>0 and sign(el)>0 and sign(er)>0
            assert [el,er]==list(map(parse,row['escape_endpoint_values']))
            assert [dl,dr]==list(map(parse,row['escape_derivatives']))
            if sign(dl)>=0: case,certificate='left',el
            elif sign(dr)<=0: case,certificate='right',er
            else:
                case,certificate='interior',add(scale(mul(A,C),4),scale(mul(B,B),-1))
                assert sign(certificate)>0
            assert case==row['escape_minimum_case'] and certificate==parse(row['escape_minimum_certificate'])
            counts['all_label_I2_tails'] += 1
            counts['tail_minimum_'+case] += 1
        records.append(dict(edge=ei,finite_profiles=len(profiles),all_label_I2_tails=len(isolated),
                            independent_dual_basis=True,complete_profile_budget=True))
        print('independent four-support review',ei,len(profiles),len(isolated),flush=True)
    assert seen_finite==set(finite) and seen_tails==set(tails) and len(records)==4
    paths=[BASE,WALK,PRIMARY,f'scripts/{TAG}_review.py','scripts/fresh_h20_e7_replacement_escape_20261014_review.py']
    return dict(status='PASS_INDEPENDENT_FOUR_SUPPORT_ALL_PARAMETER_ESCAPE_OBSTRUCTION',
        sources={s:sha(s) for s in paths},counts=dict(counts),edges=records,
        proof='For every spherical extension, r=p+aW+(q-a^2 W^2)/(2ah)*nu with q=2-p^2>0. Acuteness with N9 and N88 forces r.nu<0, hence a>0. For each old nonshared wall and E, r.x=p.x+B_x*a+C_x*q/a is strictly decreasing: B_x=Z.x/(-2h)<0 and C_x=nu.x/(2h)>=0. At the recorded positive a0 one old wall has r.N>0 and r.E<0. Thus every a violates either acuteness or escape cutting. For isolated I2 tails, exact quadratic bounds establish the same separator for every c=2cos(pi/m) in [sqrt3,2], covering all m>=7.',
        scope='Only the four listed cusp supports and spherical finite-endpoint extensions retaining all95 walls. This proves no such cap cuts the specified strict escape. It does not exclude non-cutting caps, other supports, ideal endpoint extensions, changed wall sets, or the H20 construction route. No finite-volume example is produced.')

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--verify',action='store_true')
    args=ap.parse_args()
    result=run()
    target=ROOT/f'results/{TAG}_review.json'
    if args.verify: assert result==json.loads(target.read_text())
    else: target.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'])
