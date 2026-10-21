#!/usr/bin/env python3
"""Standalone physical interval replay for continuous high N and high S.

No project imports, numerical root solving, native-bank or cone reconstruction.
The same Schur identities are implemented separately, with directed rational
rounding. This is not an alternative interval algorithm or scaffold proof.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations, product
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_oblique_continuous_high_n_high_s_interval_20261018'
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
    if m == 8: return sqrt(add(iv(2), sqrt(iv(2))))
    if m == 10: return sqrt(div(add(iv(5), sqrt(iv(5))), iv(2)))
    if m == 12: return div(add(sqrt(iv(2)), sqrt(iv(6))), iv(2))
    assert m == 7
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


def check():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    pri=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    prior=json.loads((ROOT/f'results/{PRIOR}.json').read_text())
    old=json.loads((ROOT/f'results/{HIGH}.json').read_text())
    assert d['status']=='PASS_CONTINUOUS_HIGH_N_ALL_NINE_HIGH_S_FAMILIES_INTERVAL_OBSTRUCTION'
    assert d['interval_scale']==SCALE and len(d['records'])==66
    assert d['parents']==prior['parents'] and len(d['parents'])==len(old['records'])==2
    ids = WHITE+OLD
    metric, center = [list(map(iv, r)) for r in pri['computational_gram']], list(map(iv, pri['common_center']))
    roots = {}
    for row in pri['facet_witnesses']:
        if row['wall'] in ids:
            hp = iv(pri['center_products'][str(row['wall'])])
            assert hp[1] < 0
            roots[row['wall']] = [div(mul(iv(2), sub(x, y)), hp) for x, y in zip(center, map(iv, row['point']))]
    native = {a: {b: dot(roots[a], metric, roots[b]) for b in ids} for a in ids}
    wi, pivots = gauss([[native[a][b] for b in WHITE] for a in WHITE])
    assert all(p[0] > 0 for p in pivots)
    wc0 = {a: [native[a][w] for w in WHITE] for a in OLD}
    g0 = {a: {b: sub(native[a][b], dot(wc0[a], wi, wc0[b])) for b in OLD} for a in OLD}
    assert prior['primary_projection'] == enc(dict(walls=ids, physical_Gram=native, white_inverse=wi,
        white_pivots=pivots, old_white_covectors=wc0, projected_Gram=g0))
    profiles = [[iv(0)]*6]
    for pos in (4, 5):
        for m in (3, 4, 5):
            v = [iv(0)]*6
            v[pos] = neg(theta(m))
            profiles.append(v)
    assert enc(profiles) == prior['profiles']

    def extend(g, wc, row):
        wall, frame, pi, labels = (row[k] for k in ('wall', 'frame', 'profile', 'labels'))
        gram = [[g[a][b] for b in frame] for a in frame]
        inv, det = inverse(gram)
        assert det[1] < 0 and gram[0][0][0] > 0
        rhs = [sub(neg(theta(m)), dot(profiles[pi], wi, wc[i])) for i, m in zip(frame, labels)]
        square = sub(iv(2), dot(profiles[pi], wi, profiles[pi]))
        rem = sub(square, dot(rhs, inv, rhs))
        w, sign = row['witness'], row['witness_sign']
        assert w in g and w not in frame and sign in (-1, 1)
        cv = [g[i][w] for i in frame]
        q = sub(g[w][w], dot(cv, inv, cv))
        assert q[0] > 0 and rem[0] > 0
        linear, radical = dot(rhs, inv, cv), sqrt(mul(rem, q))
        branches = [sub(linear, radical), add(linear, radical)]
        factor = mul(iv(sign), sqrt(div(rem, q)))
        out = {}
        for j in g:
            cj = [g[i][j] for i in frame]
            out[j] = rhs[frame.index(j)] if j in frame else add(dot(rhs, inv, cj),
                mul(factor, sub(g[w][j], dot(cv, inv, cj))))
        expected = dict(frame_Gram=gram, inverse=inv, determinant=det, rhs=rhs,
            projected_square=square, remainder=rem, witness_residual=q,
            branch_witness_products=branches, new_products=out)
        assert all(row[k] == enc(v) for k, v in expected.items())
        distances = [abs(row['numerical_witness_product']-float(sum(v)/2)) for v in branches]
        assert distances == row['branch_comparison'] and distances[(sign+1)//2] < 2e-6 and distances[(1-sign)//2] > 1e-4
        assert row['all_product_comparison_error'] < 2e-6
        for j, v in out.items(): g[j][wall] = v
        g[wall], wc[wall] = {**out, wall: square}, profiles[pi]


    families=[[0,[k,2,6],[2]]for k in range(2,6)]+[[0,[6,2,k],[0]]for k in range(3,6)]+[[0,[6,2,6],[0,2]]]
    assert d['nonright_high_S_families']==families
    assert d['N_families']==[[2,2,'m>=6'],[3,2,'m>=6']]
    domain=(theta(6)[0],F(2));assert d['N_parameter_domain']==enc(domain)
    remainders,positives,bounds,derivatives=[],[],[],[]
    n8values,nonright_bounds,polar_dets,polar_margins=[],[],[],[]
    rowkeys=[]
    def verify_box(g,wc,row):
        a=row['N10_magnitude'];n=iv(row['N']['N_domain'])
        assert a in(0,1) and domain[0]<=n[0]<n[1]<=domain[1]
        frame=[10,67,69];ni,nd=inverse([[g[i][j]for j in frame]for i in frame])
        cv=[g[i][64]for i in frame];sig=sub(g[64][64],dot(cv,ni,cv))
        rhs=[iv(-a),iv(0),neg(n)];nrho=sub(iv(2),dot(rhs,ni,rhs))
        assert nd[1]<0 and sig[0]>0 and nrho[0]>0
        factor=sqrt(div(nrho,sig));out={}
        for j in g:
            cj=[g[i][j]for i in frame];linear=dot(rhs,ni,cj)
            radical=mul(factor,sub(g[64][j],dot(cv,ni,cj)))
            out[j]=rhs[frame.index(j)]if j in frame else add(linear,radical)
            if j==8:opposite8=sub(linear,radical);assert opposite8[0]>0
        expected=dict(frame=frame,inverse=ni,determinant=nd,witness=64,witness_residual=sig,
            N_domain=n,N10_magnitude=a,rhs=rhs,norm_remainder=nrho,
            opposite_N_sheet_old8=opposite8,new_products=out)
        assert row['N']==enc(expected);n8values.append(opposite8[0])
        gn={i:r.copy()for i,r in g.items()};wn={**wc,73:profiles[0]}
        for j in gn:gn[j][73]=out[j]
        gn[73]={**out,73:iv(2)}
        del gn[70],wn[70]
        for rr in gn.values():del rr[70]
        assert 70 not in gn and all(70 not in rr for rr in gn.values())
        # S profile zero makes its projected product the full product.
        sf=[52,68,71];si,sd=inverse([[gn[i][j]for j in sf]for i in sf])
        sr=[iv(-1),iv(0),iv(0)];sc=[iv(0),iv(0),iv(-1)];s8=[gn[i][8]for i in sf]
        alpha=neg(dot(sr,si,sr));q8=sub(gn[8][8],dot(s8,si,s8))
        assert sd[1]<0 and alpha[0]>0 and q8[0]>0
        z=(div(theta(6),sqrt(add(iv(2),mul(iv(3),alpha))))[0],div(iv(2),sqrt(add(iv(2),mul(iv(4),alpha))))[1])
        assert add(mul(z,dot(sr,si,s8)),div(q8,sqrt(q8)))[0]>0
        rec=row['R_cube']
        # No numerical R selection: all28 recipes and both exact signs.
        frame, witness = [69, 71, 73], 10
        gram = [[gn[a][b] for b in frame] for a in frame]
        inv, det = inverse(gram)
        cv = [gn[i][witness] for i in frame]
        sigma = sub(gn[witness][witness], dot(cv, inv, cv))
        assert det[1] < 0 and gram[0][0][0] > 0 and sigma[0] > 0
        expected = dict(R_frame=frame, frame_Gram=gram, inverse=inv, determinant=det,
            witness=witness, witness_covector=cv, witness_residual=sigma)
        assert all(rec[k] == enc(v) for k, v in expected.items())
        assert len(rec['recipes']) == 28
        assert [(r['profile'], r['RM_label']) for r in rec['recipes']] == list(product(range(7), range(2, 6)))
        for index, (profile, label) in enumerate(product(range(7), range(2, 6))):
            saved_recipe = rec['recipes'][index]
            labels = [2, label, 2]
            rhs = [sub(neg(theta(m)), dot(profiles[profile], wi, wn[i])) for i, m in zip(frame, labels)]
            square = sub(iv(2), dot(profiles[profile], wi, profiles[profile]))
            rho = sub(square, dot(rhs, inv, rhs))
            assert rho[0] > 0
            factor = sqrt(div(rho, sigma))
            sheets = {}
            for sign in (-1, 1):
                contacts = {}
                for j in gn:
                    cj = [gn[i][j] for i in frame]
                    contacts[j] = rhs[frame.index(j)] if j in frame else add(dot(rhs, inv, cj),
                        mul(mul(iv(sign), factor), sub(gn[witness][j], dot(cv, inv, cj))))
                sheets[sign] = contacts
            positive8 = add(sheets[1][8], dot(profiles[profile], wi, wn[8]))
            assert positive8[0] > 0
            cc = [sheets[-1][i] for i in sf]
            a = dot(sr, si, cc)
            b = div(sub(sheets[-1][8], dot(s8, si, cc)), sqrt(q8))
            derivative = sub(a, mul(mul(alpha, z), b))
            assert derivative[1] < 0
            def evaluate(t):
                return sub(mul(t, a), mul(sqrt(add(iv(2), mul(alpha, mul(t, t)))), b))
            at6, limit = evaluate(theta(6)), evaluate(iv(2))
            product_range = limit[0], at6[1]
            assert product_range[1] < -2
            continuous = dict(S_frame=sf, S_frame_inverse=si, S_frame_determinant=sd,
                alpha=alpha, residual8_square=q8, A_dot_R=a, B_dot_R=b,
                theta_domain=[theta(6)[0], F(2)], transformed_z_domain=z,
                derivative_range=derivative, at_label6=at6, at_limit=limit,
                full_product_range=product_range, monotonicity='strictly_decreasing')
            expected = dict(index=index, profile=profile, RM_label=label, R_labels=labels,
                rhs=rhs, projected_square=square, norm_remainder=rho,
                positive_R_sheet_old8=positive8,
                negative_R_sheet_products={i: sheets[-1][i] for i in (8, 52, 68, 71)},
                continuous=continuous)
            assert saved_recipe == enc(expected)
            remainders.append(rho[0]); positives.append(positive8[0])
            bounds.append(product_range[1]); derivatives.append(derivative[1])
        aa,ac,cc=dot(sr,si,sr),dot(sr,si,sc),dot(sc,si,sc)
        a8,c8=dot(sr,si,s8),dot(sc,si,s8)
        basis=dict(AA=aa,AC=ac,CC=cc,A8=a8,C8=c8,residual8_square=q8)
        assert row['S_basis']==enc(basis) and len(row['nonright_families'])==8
        coeff=[]
        for r in rec['recipes']:
            contacts=[iv(r['negative_R_sheet_products'][str(i)])for i in sf]
            coeff.append((r['index'],dot(sr,si,contacts),dot(sc,si,contacts),
                div(sub(iv(r['negative_R_sheet_products']['8']),dot(s8,si,contacts)),sqrt(q8))))
        for family,saved in zip(families,row['nonright_families']):
            _,labels,high=family
            th=domain if 0 in high else theta(labels[0]);ta=domain if 2 in high else theta(labels[2])
            rem=sub(iv(2),add(add(mul(aa,mul(th,th)),mul(iv(2),mul(ac,mul(th,ta)))),mul(cc,mul(ta,ta))))
            assert rem[0]>0
            positive8=add(add(mul(th,a8),mul(ta,c8)),mul(sqrt(rem),sqrt(q8)));assert positive8[0]>0
            joins=[]
            for i,ar,cr,br in coeff:
                value=sub(add(mul(th,ar),mul(ta,cr)),mul(sqrt(rem),br));assert value[1]<-2
                joins.append(dict(R_index=i,full_SR_range=value));nonright_bounds.append(value[1])
            expected=dict(family=family,theta_domain=th,tau_domain=ta,norm_remainder_range=rem,
                positive_S_sheet_old8=positive8,negative_S_R_joins=joins)
            assert saved==enc(expected)
        # Full 3-frame avoids division by the limiting singular 2-frame X,N.
        ap={j:dot(sr,si,[gn[i][j]for i in sf])for j in gn}
        bp={j:div(sub(gn[8][j],dot(s8,si,[gn[i][j]for i in sf])),sqrt(q8))for j in gn}
        vv=sub(iv(1),mul(alpha,mul(z,z)));vx=sub(mul(z,ap[69]),bp[69]);vn=sub(mul(z,ap[73]),bp[73])
        gram=[[gn[69][69],gn[69][73],vx],[gn[73][69],gn[73][73],vn],[vx,vn,vv]]
        det=determinant(gram);minor=sub(mul(gn[69][69],vv),mul(vx,vx))
        h=sub(mul(sub(iv(3),mul(iv(2),gn[71][71])),neg(det)),mul(iv(2),mul(mul(gn[71][73],gn[71][73]),minor)))
        assert gn[69][71]==iv(0) and gn[71][71]==iv(2) and sr[2]==iv(0)
        assert det[1]<0 and h[0]>0
        expected=dict(S_frame=sf,S_inverse=si,S_determinant=sd,alpha=alpha,residual8_square=q8,
            transformed_z_domain=z,V_square=vv,VX=vx,VN=vn,full_Gram=gram,determinant=det,
            M_covector=[iv(0),gn[71][73],iv(0)],middle_cofactor=minor,strict_obstruction=h)
        assert row['high_RM']==enc(expected)
        polar_dets.append(det[1]);polar_margins.append(h[0])

    for parent,source in zip(prior['parents'],old['records']):
        assert parent['index']==source['index'] and parent['instruction_key']==source['instruction_key']
        chain=json.loads(source['instruction_key']);assert len(chain)==len(parent['steps'])==5
        assert [x[0]for x in chain]==[67,68,69,70,71]
        g={i:r.copy()for i,r in g0.items()};wc=wc0.copy()
        for step,instruction in zip(parent['steps'],chain):
            assert [step[k]for k in('wall','frame','profile','labels','saved_numerical_branch')]==instruction
            extend(g,wc,step)
        for a in(0,1):
            rows=[r for r in d['records']if r['parent']==parent['index']and r['N10_magnitude']==a]
            parts=sorted(iv(r['N']['N_domain'])for r in rows)
            assert parts and parts[0][0]==domain[0]and parts[-1][1]==domain[1]
            assert all(x[1]==y[0]for x,y in zip(parts,parts[1:]))
            for row in rows:
                n=iv(row['N']['N_domain']);depth=row['depth']
                assert isinstance(depth,int)and 0<=depth<=8
                step=(domain[1]-domain[0])/2**depth
                assert n[1]-n[0]==step and ((n[0]-domain[0])/step).denominator==1
                rowkeys.append(row['index']);verify_box(g,wc,row)
    assert sorted(rowkeys)==list(range(66))
    assert len(bounds)==len(positives)==66*28 and len(nonright_bounds)==66*8*28
    summary=dict(parent_classes=2,parent_origins=4,continuous_N_families=4,attempted_boxes=128,covering_boxes=66,
        low_R_recipe_boxes=len(bounds),low_R_sheet_boxes=2*len(bounds),SM_right_continuous_joins=len(bounds),
        nonright_SM_continuous_joins=len(nonright_bounds),high_RM_polar_boxes=66,
        minimum_opposite_N_old8=min(n8values),minimum_positive_R_old8=min(positives),
        maximum_SM_right_SR=max(bounds),maximum_nonright_SM_SR=max(nonright_bounds),
        maximum_high_RM_determinant=max(polar_dets),minimum_high_RM_polynomial=min(polar_margins),
        unresolved=0,accepted=0,new_cones=0)
    assert d['summary']==enc(summary)
    assert d['exact_deletion_debit']==old['exact_deletion_debit']and d['deleted_walls']==old['deleted_walls']
    assert d['exact_deletion_debit']['lost_finite_count']==19
    assert d['exact_deletion_debit']['compact_faces_by_dimension']=={'2':12,'3':1}
    assert d['topology']==old['topology'] and 8-2*(3+1)==0
    print('PASS_STANDALONE_CONTINUOUS_HIGH_N_HIGH_S',len(rowkeys),len(bounds),len(nonright_bounds),flush=True)


if __name__=='__main__':check()
