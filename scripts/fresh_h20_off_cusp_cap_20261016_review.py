#!/usr/bin/env python3
"""Independent physical review of the rational off-cusp141st wall.

No producer arithmetic imports. The old parameter's existence/isolation and
old chamber acuteness are inherited from the frozen independent chain review.
Every new contact, center, finite-corner exclusion and remaining escape is replayed.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from fresh_h20_bridge_edge_closure_20261014_review import Field, parse
import fresh_h20_bridge_faces_20261014_review as B
import fresh_h20_s61_multicap_20261014_exact_review as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_off_cusp_cap_20261016_review'
TARGET = 'results/fresh_h20_affine_a17_completion_20261016.json'
CAP_REVIEW = 'results/fresh_h20_complementary_cap_20261016_review.json'
INPUT = 'results/fresh_h20_complementary_cap_20261016.json'
PRIOR = 'results/fresh_h20_dotted_arc_cap_20261016.json'
PRIOR_REVIEW = 'results/fresh_h20_dotted_arc_cap_20261016_review.json'
BASE = 'results/fresh_h20_s61_multicap_20261014_exact.json'
CHAIN_REVIEW = 'results/fresh_h20_s61_multicap_20261014_exact_review.json'
BANK = 'results/fresh_h20_nonright_s61_20261014_bank_gate.json'
FILES = R.FILES[:5]
pt, plus, minus, neg, times, div, root = R.pt, R.plus, R.minus, R.neg, R.times, R.div, R.root
total, inner, lor = R.total, R.inner, R.lor


def saved(x):
    den = 1 << 192
    return [str(F(R.floor(x[0] * den), den)), str(F(R.ceil(x[1] * den), den))]


def vec(v):
    return [saved(x) for x in v]


def overlap(a, raw):
    b = tuple(map(F, raw))
    assert R.overlap(a, b)


def integrity(paths):
    queue, seen, count = list(paths), set(), 0
    while queue:
        p = queue.pop()
        if not p.endswith('.json') or p in seen:
            continue
        seen.add(p)
        for q, digest in json.loads((ROOT / p).read_text()).get('sources', {}).items():
            assert sha256((ROOT / q).read_bytes()).hexdigest() == digest, (p, q)
            count += 1
            if q.endswith('.json'):
                queue.append(q)
    return dict(manifests_checked=len(seen), source_references_checked=count)


def run():
    old, bridge, edge, polar, cur = [json.loads((ROOT / p).read_text()) for p in FILES]
    d, base, reviewed, prior, prior_review = [json.loads((ROOT / p).read_text()) for p in (INPUT, BASE, CHAIN_REVIEW, PRIOR, PRIOR_REVIEW)]
    assert d['status'] == 'PASS_EXACT_135_WALL_COMPLEMENTARY_E0_CUTTER'
    assert prior_review['sources'][PRIOR] == sha256((ROOT / PRIOR).read_bytes()).hexdigest()
    assert reviewed['status'] == 'PASS_INDEPENDENT_133_WALL_RIGHT_CHAIN_AND_COMPACT_S61_FACE'
    assert base['known_total_vertices_at_least'] == 168
    ff = Field([parse(old['field']['discriminant'], 0), parse(bridge['field']['new_positive_sqrt_radical'], 1), parse(edge['field']['new_positive_sqrt_radical'], 2)])
    add, mul, sc, ft, const = ff.add, ff.mul, ff.scale, ff.total, ff.const
    zero = ff.zeros[3]
    Q = [list(map(F, row)) for row in old['form']]
    assert B.inertia(Q) == [20, 1, 0]
    form_rows = [[(j, q) for j, q in enumerate(row) if q] for row in Q]
    cache = {}
    def vp(raw, level):
        return tuple(ff.lift(parse(x, level), level) for x in raw)
    def dot(v, w):
        key = id(w)
        if key not in cache:
            cache[key] = (w, tuple(ft(sc(w[j], q) for j, q in row) for row in form_rows))
        return ft(mul(x, y) for x, y in zip(v, cache[key][1]))
    def lin(vectors, coefficients):
        return tuple(ft(sc(v[j], c) for v, c in zip(vectors, coefficients)) for j in range(21))
    def rat(x):
        for k in (3, 2, 1):
            assert x[1] == ff.zeros[k - 1]
            x = x[0]
        assert not any(x[1:])
        return x[0]
    @lru_cache(None)
    def iv(x):
        return R.out(ff.interval(x, 3, R.BITS))
    N = [vp(v, 1) for v in old['normals']] + [vp(bridge['new_normal'], 2)]
    N += [vp(v, 3) for v in edge['new_normals']]
    N += [tuple(const(x) for x in v) for v in polar['new_normals']]
    N += [vp(cur['cap_R']['normal'], 1)]
    norms = [const(x) for x in old['norms']] + [const(2)] + [parse(edge['new_normal_norm'], 3)] * 2
    norms += [const(x) for x in polar['new_norms']] + [ff.lift(parse(cur['cap_R']['norm'], 1), 1)]
    assert len(N) == len(norms) == 101
    normden = [root(div(iv(q), pt(2))) for q in norms]
    S = base['support']
    assert S == edge['recipe']['support'] and len(S) == 18 and 31 in S and 5 in S
    G = [[rat(dot(N[i], N[j])) for j in S] for i in S]
    J = B.inv(G)
    support_pivots = B.pd(G)
    assert G == [list(map(F, row)) for row in base['support_gram']]
    assert J == [list(map(F, row)) for row in base['support_inverse']]
    index5 = S.index(5)
    assert J[index5][index5] == F(30,11)
    omega = lin([N[i] for i in S], J[index5])
    assert dot(omega,omega) == const(F(30,11))
    assert [dot(omega,N[i]) for i in S] == [const(i==5) for i in S]
    assert dot(omega,N[0]) == const(F(-6,11))
    assert norms[0] == const(2)
    assert dot(omega,N[95]) == zero
    frame = [vp(base['frame'][key], 2) for key in ('X', 'r95', 'W')]
    squares = [ff.lift(parse(base['frame']['A'], 2), 2), const(2), ff.lift(parse(base['frame']['W_norm'], 2), 2)]
    assert [[dot(v, w) for w in frame] for v in frame] == [[sc(squares[0], -1), zero, zero], [zero, squares[1], zero], [zero, zero, squares[2]]]
    assert all(dot(v, N[i]) == zero for v in frame for i in S)
    assert all(dot(v, omega) == zero for v in frame)
    frame_den = [root(iv(x)) for x in squares]
    hb = [vp(cur['hex_basis'][key], 1) for key in ('C', 'U', 'V')]
    hq = [ff.lift(parse(x, 1), 1) for x in cur['hex_basis']['squared_lengths']]
    hex_coeff = [[ff.lift(parse(x, 0), 0) for x in row] for row in cur['hexagon']['frame_coefficients']]
    def old_pairs(v):
        pairs = [div(iv(dot(v, n)), den) for n, den in zip(N, normden)]
        for co in hex_coeff:
            pairs.append(total(div(iv(mul(c, dot(v, w))), root(iv(q))) for c, w, q in zip(co, hb, hq)))
        assert len(pairs) == 107
        return pairs
    # Recompute all old107 covectors from literal normals, using the prior
    # independent sparse-radical field engine, not saved producer intervals.
    by_axis = [[div(x, den) for x in old_pairs(v)] for v, den in zip(frame, frame_den)]
    P = [list(row) for row in zip(*by_axis)]
    for i in S:
        assert all(dot(v, N[i]) == zero for v in frame)
        P[i] = [pt(0)] * 3
    for row, raw in zip(P, base['frame']['old107_unit_normal_pair_bounds']):
        for value, value_raw in zip(row, raw):
            overlap(value, value_raw)
    a = tuple(map(F, base['parameter']['isolating_interval']))
    assert list(map(str, a)) == reviewed['parameter_interval']
    coeff, rotation = R.chain(a)
    assert len(coeff) == 26
    for row, raw in zip(coeff, base['chain']['frame_coefficient_bounds']):
        for value, value_raw in zip(row, raw):
            overlap(value, value_raw)
    P += [[neg(row[0]), row[1], row[2]] for row in coeff]
    wp = old_pairs(omega) + [pt(0)] * 26
    for i in S:
        overlap(wp[i], [int(i == 5)] * 2)
        wp[i] = pt(i == 5)
    print('PASS independent field frame,107 covectors,26 chain normals and omega5', flush=True)

    # Anchor95 has zero e-coordinate contact exactly; anchor0 has pi/6.
    assert frame[1] == N[95] and dot(N[95],N[95]) == const(2)
    P[95] = [pt(0),root(pt(2)),pt(0)]
    wp[0] = pt(F(-6,11)); wp[95] = pt(0)
    p,q = P[0][0],P[0][2]
    assert p[1]<0 and q[1]<0
    h = minus(pt(F(-6,11)),root(pt(3)))
    alpha,beta = div(h,p),neg(div(q,p))
    qa = minus(pt(1),times(beta,beta))
    qb = neg(times(pt(2),times(alpha,beta)))
    qc = minus(pt(F(8,11)),times(alpha,alpha))
    disc = minus(times(qb,qb),times(pt(4),times(qa,qc)))
    assert qa[0]>0 and disc[0]>0
    z = div(plus(neg(qb),root(disc)),times(pt(2),qa))
    zminus = div(minus(neg(qb),root(disc)),times(pt(2),qa))
    t = plus(alpha,times(beta,z))
    y = [t,pt(0),z]
    assert t[0]>0 and z[0]>0 and zminus[1]<0
    for value,raw in zip(y,d['definition']['frame_bounds']):overlap(value,raw)
    # Symbolic equations certify exact norm/finite contacts, beyond intervals.
    import sympy as sy
    aa,bb,zz,pp,qq,hh=sy.symbols('aa bb zz pp qq hh')
    polynomial=(1-bb**2)*zz**2-2*aa*bb*zz+sy.Rational(8,11)-aa**2
    assert sy.expand((aa+bb*zz)**2-zz**2-sy.Rational(8,11)+polynomial)==0
    assert sy.simplify(pp*(hh/pp-qq*zz/pp)+qq*zz-hh)==0
    assert R.overlap(lor(y,y),pt(F(-8,11)))
    assert F(-8,11)+F(30,11)==2
    finite={i:2 for i in S}|{5:3,0:6,95:2}
    pairs=[];counts=Counter()
    for i,(cov,w) in enumerate(zip(P,wp)):
        value=minus(inner(y,cov),w)
        if i in finite:
            m=finite[i]
            expected=pt(0) if m==2 else pt(-1) if m==3 else neg(root(pt(3)))
            assert R.overlap(value,expected)
            relation='right' if m==2 else f'pi/{m}'
        else:
            assert value[1]<-2,(i,saved(value))
            relation='dotted'
        assert relation==d['pairs'][i]['relation']
        overlap(value,d['pairs'][i]['pair_bounds'])
        pairs.append(dict(wall=i,relation=relation,bounds=saved(value)))
        counts[relation]+=1
    # Reconstruct the preceding arc cap frame from its independently reviewed recipe.
    sqrt2=root(pt(2));h2=minus(sqrt2,pt(1));a2=times(a,a)
    y133=[div(plus(times(h2,a2),times(pt(2),sqrt2)),times(pt(2),a)),
          times(div(h2,pt(2)),coeff[22][1]),times(div(h2,pt(2)),coeff[22][2])]
    c133=root(plus(pt(5),div(pt(4),a2)))
    assert J[index5][S.index(31)]==0
    mutual=lor(y,y133)
    assert mutual[1]<-2
    overlap(mutual,d['mutual_arc_cap_pair_bounds'])
    pairs.append(dict(wall=133,relation='dotted',bounds=saved(mutual)));counts['dotted']+=1
    assert counts==dict(right=18,**{'pi/3':1,'pi/6':1,'dotted':114})
    assert sorted(finite)==d['finite_neighbors']
    print('PASS all134 independent contacts and symbolic finite equalities',flush=True)

    target=json.loads((ROOT/TARGET).read_text())
    capreview=json.loads((ROOT/CAP_REVIEW).read_text())
    assert capreview['sources'][INPUT]==sha256((ROOT/INPUT).read_bytes()).hexdigest()
    assert target['status']=='PASS_EXACT_140_WALL_AFFINE_G2_A17_CUSP_COMPLETION'
    def ilinear(vectors,weights):
        return [total(times(c,v[j]) for v,c in zip(vectors,weights)) for j in range(21)]
    def ivadd(a,b):return [plus(x,z) for x,z in zip(a,b)]
    def ivscale(v,c):return [times(x,c) for x in v]
    def ipair(a,b):
        return total(times(pt(q),times(a[i],b[j])) for i,row in enumerate(form_rows) for j,q in row)
    unitframe=[[div(iv(x),den) for x in v] for v,den in zip(frame,frame_den)]
    physical=[[div(iv(x),den) for x in v] for v,den in zip(N,normden)]
    for co in hex_coeff:
        physical.append(ilinear([[div(iv(x),root(iv(q))) for x in w] for w,q in zip(hb,hq)],[iv(x) for x in co]))
    assert len(physical)==107
    physical += [ilinear(unitframe,row) for row in coeff]
    dual=[ilinear([physical[i] for i in S],list(map(pt,row))) for row in J]
    physical += [ivadd(ilinear(unitframe,y133),ivscale(dual[S.index(31)],neg(c133))),
                 ivadd(ilinear(unitframe,y),ivscale(dual[index5],pt(-1)))]
    assert len(physical)==135
    Hraw=vp(cur['strict_center']['point'],1)
    Hq=dot(Hraw,Hraw);assert ff.sign(Hq)<0
    Hunit=[div(iv(x),root(iv(sc(Hq,-1)))) for x in Hraw]
    H=ivadd(ilinear(unitframe,[pt(1),pt(F(-1,2)),pt(0)]),ivscale(Hunit,pt(F(base['strict_center']['epsilon']))))
    assert ipair(H,H)[1]<0 and all(ipair(H,n)[1]<0 for n in physical)
    nu=ilinear([physical[134],physical[0],physical[1]],[root(pt(3)),pt(2),pt(1)])
    nh=ipair(nu,H);assert nh[1]<0
    paths=[[3,4],[6,7,8,9,10],[12,13,14,15],[17],[31]]
    white=sum(paths,[]);uids=white+[0,1]
    GU=[[rat(dot(N[i],N[j])) for j in uids] for i in uids]
    assert all(norms[i]==const(2) for i in uids)
    assert all(dot(N[i],N[j])==zero for i in white for j in (0,1))
    assert dot(N[0],N[1])==const(-1)
    import sympy as sy
    def smatrix(m):return sy.Matrix([[sy.Rational(v.numerator,v.denominator) for v in row] for row in m])
    SG=smatrix(GU);assert SG.is_positive_definite
    JU=[[F(x) for x in row] for row in SG.inv().tolist()]
    U=[physical[i] for i in uids]
    def project_off_U(v):
        products=[ipair(v,w) for w in U]
        coefficients=[neg(total(times(pt(JU[i][j]),products[j]) for j in range(15))) for i in range(15)]
        return ivadd(v,ilinear(U,coefficients))
    hbar=project_off_U(H);assert ipair(hbar,hbar)[1]<0
    k=ivscale(hbar,div(pt(-1),nh))
    eta=ivadd(k,ivscale(nu,div(ipair(k,k),pt(2))))
    assert R.overlap(ipair(nu,eta),pt(-1)) and R.overlap(ipair(eta,eta),pt(0))
    # The affineG2 identities imply nu is exactly null and orthogonal to U.
    # Projection and null-companion formulas then give exact orthogonality;
    # interval overlaps below are sanity checks, not proofs of equality.
    free=[];pivot_checks=[]
    for coordinate in (0,3,6,12):
        v=project_off_U([pt(int(i==coordinate)) for i in range(21)])
        v=ivadd(v,ilinear([nu,eta],[ipair(v,eta),ipair(v,nu)]))
        for w in free:v=ivadd(v,ivscale(w,neg(ipair(v,w))))
        square=ipair(v,v);assert square[0]>0
        free.append(ivscale(v,div(pt(1),root(square))))
        pivot_checks.append(dict(coordinate=coordinate,norm_bounds=saved(square)))
    assert [r['coordinate'] for r in target['definition']['free_pivots']]==[0,3,6,12]
    for v,raw in zip(free,target['definition']['free_frame_bounds']):
        for x,bounds in zip(v,raw):overlap(x,bounds)
    for i,v in enumerate(free):
        for w in U+[nu,eta]:assert R.overlap(ipair(v,w),pt(0))
        for j,w in enumerate(free):assert R.overlap(ipair(v,w),pt(i==j))
    JW=SG[:13,:13].inv()
    contact=sy.zeros(5,13);ends=[]
    for i in range(5):
        ends.append([paths[i][-1],paths[(i+1)%5][0]])
        for wall in ends[i]:contact[i,white.index(wall)]=-1
    projection=contact*JW
    schur=2*sy.eye(5)-projection*contact.T
    weights=[sy.Rational(1,len(p)+1) for p in paths]
    assert schur*sy.ones(5,1)==sy.zeros(5,1) and schur.rank()==4
    expected=sy.zeros(5)
    for i,weight in enumerate(weights):
        j=(i-1)%5
        expected[i,i]+=weight;expected[j,j]+=weight
        expected[i,j]-=weight;expected[j,i]-=weight
    assert schur==expected
    assert list(sy.ones(1,5)*projection)==[-1]*13
    L,diag=schur[:4,:4].LDLdecomposition(hermitian=False)
    diagonal=[diag[i,i] for i in range(4)]
    assert diagonal==[sy.Rational(1,2),sy.Rational(14,45),sy.Rational(4,7),sy.Rational(9,16)]
    assert L*diag*L.T==schur[:4,:4]
    assert [str(v) for v in diagonal]==target['definition']['Schur_LDL_D']
    residuals=[ilinear(free,[times(pt(F(L[i,j])),root(pt(F(diagonal[j])))) for j in range(4)]) for i in range(4)]
    residuals.append(ilinear(residuals,[pt(-1)]*4))
    unshifted=[ivadd(ilinear([physical[j] for j in white],[pt(F(projection[i,j])) for j in range(13)]),residuals[i]) for i in range(5)]
    heights=[4,4,3,4,4];assert heights==target['definition']['heights']
    shifted=[ivadd(v,ivscale(nu,pt(height))) for v,height in zip(unshifted,heights)]
    incident=set([134,0,1]+white)
    null_pairs=[ipair(nu,n) for n in physical]
    for i,g in enumerate(null_pairs):
        if i in incident:assert R.overlap(g,pt(0))
        else:assert g[1]<0,(i,saved(g))
    # Exact incident zeros follow from the old affine triangle and the
    # literal old-support orthogonality checked above, including r134.S.
    pair_records=[];counts=Counter()
    for i,v in enumerate(shifted):
        for j,n in enumerate(physical):
            g=ipair(v,n)
            if j in incident:
                expected=-1 if j in ends[i] else 0
                assert R.overlap(g,pt(expected));relation='pi/3' if expected else 'right'
            else:
                assert g[1]<-2,(i,j,saved(g));relation='dotted'
            row=target['new_old_pairs'][135*i+j]
            assert row['new_wall']==135+i and row['old_wall']==j and row['relation']==relation
            overlap(g,row['pair_bounds']);counts[relation]+=1
            pair_records.append(dict(new_wall=135+i,old_wall=j,relation=relation,bounds=saved(g)))
        for j,w in enumerate(shifted):assert R.overlap(ipair(v,w),pt(2*(i==j)))
    assert counts==dict(right=70,**{'pi/3':10,'dotted':595})
    assert len(pair_records)+10==685
    full=physical+shifted;center=[ipair(v,H) for v in full];hnorm=ipair(H,H)
    assert all(g[1]<0 for g in center)
    facets=[minus(hnorm,div(times(g,g),pt(2))) for g in center]
    assert all(g[1]<0 for g in facets)
    print('PASS independent140facets and685 newpair gates',flush=True)
    cycle=[]
    for i,path in enumerate(paths):cycle+=path+[135+i]
    star=sorted([134,0,1]+cycle)
    assert len(star)==21 and target['new_proper_cusp']['incident_walls']==star
    assert target['new_proper_cusp']['A17_null_sum_multiplier']==sum(heights)==19
    # Both positive null dependencies are literal construction identities:
    # sqrt3*r134+2*n0+n1=nu; sum(A17roots)=19nu.
    nullsum=ilinear([full[j] for j in cycle],[pt(1)]*18)
    assert all(R.overlap(x,times(pt(19),z)) for x,z in zip(nullsum,nu))
    exclusions=[]
    for omitted in star:
        active=set(star)-{omitted}
        component='affineG2' if {134,0,1}<=active else 'affineA17'
        assert component=='affineG2' or set(cycle)<=active
        exclusions.append(dict(omitted=omitted,full_affine_component=component))
    # All finite neighbors of any newwall lie in this21-node star, so these
    #21 subsets cover every possible new20-wall finite corner.
    ez=list(map(F,d['strict_escape']['frame_coefficients']))
    eb=list(map(F,d['strict_escape']['dual_coefficients']))
    E=ivadd(ilinear(unitframe,list(map(pt,ez))),ilinear(dual,list(map(pt,eb))))
    Enew=[ipair(E,v) for v in shifted];assert all(g[1]<0 for g in Enew)
    assert capreview['strict_escape_norm']=='1842000049/4000000000'
    cusps=[]
    for raw in [old['completed_cusp']['point']]+[r['point'] for r in old['preserved_old_bounded_cusps']]:
        cv=[iv(x) for x in vp(raw,0)]
        vals=[ipair(v,cv) for v in shifted];assert all(g[1]<0 for g in vals);cusps.append(vec(vals))
    newfile='results/fresh_h20_off_cusp_cap_20261016.json'
    new=json.loads((ROOT/newfile).read_text())
    assert new['status']=='PASS_EXACT_141_WALL_OFF_CUSP_E1_CUTTER'
    coefficients=[F(8),F(-3,2),F(-77,10)]
    d17=F(14,5);i17=S.index(17);i31=S.index(31)
    r=ivadd(ilinear(unitframe,list(map(pt,coefficients))),ivscale(dual[i17],pt(-d17)))
    q=-coefficients[0]**2+coefficients[1]**2+coefficients[2]**2+d17*d17*J[i17][i17]
    assert q==F(73,50) and new['definition']['norm']==str(q)
    assert new['definition']['frame_coefficients']==list(map(str,coefficients))
    assert new['definition']['dual_coefficient']==str(d17)
    for value,raw in zip(r,new['definition']['physical_normal_bounds']):overlap(value,raw)
    finite=set(S)-{17};relations=[]
    for i,n in enumerate(full):
        g=ipair(r,n)
        if i in finite:
            relation='right'
            assert R.overlap(g,pt(0))
            # This zero is an exact frame/support-dual identity.
        else:
            relation='dotted'
            assert g[1]<0 and g[1]*g[1]>2*q,(i,saved(g))
        assert new['pairs'][i]['wall']==i and new['pairs'][i]['relation']==relation
        overlap(g,new['pairs'][i]['pair_bounds'])
        relations.append(dict(wall=i,relation=relation,bounds=saved(g)))
    assert len(finite)==17 and new['right_neighbors']==sorted(finite)
    rh=ipair(r,H);assert rh[1]<0
    new_facet_square=minus(ipair(H,H),div(times(rh,rh),pt(q)))
    assert new_facet_square[1]<0
    cut=-coefficients[0]*ez[0]+coefficients[1]*ez[1]+coefficients[2]*ez[2]-d17*sum(J[i17][j]*eb[j] for j in range(18))
    assert cut==F(961,1000) and ipair(r,E)[0]>0
    eframe=[F(1),F(-1,100),F(-49,50)]
    ecoeff=[F(-1,100000)]*18;ecoeff[i17]=F(-3,10);ecoeff[i31]=F(-2,5)
    E2=ivadd(ilinear(unitframe,list(map(pt,eframe))),ilinear(dual,list(map(pt,ecoeff))))
    eq=-eframe[0]**2+eframe[1]**2+eframe[2]**2+sum(ecoeff[i]*J[i][j]*ecoeff[j] for i in range(18) for j in range(18))
    assert eq==F(427500061,5000000000)>0
    assert new['strict_escape']['norm']==str(eq)
    assert new['strict_escape']['frame_coefficients']==list(map(str,eframe))
    assert new['strict_escape']['dual_coefficients']==list(map(str,ecoeff))
    e2pairs=[ipair(E2,n) for n in full+[r]]
    assert all(g[1]<0 for g in e2pairs) and ipair(E2,H)[1]<0
    last=-coefficients[0]*eframe[0]+coefficients[1]*eframe[1]+coefficients[2]*eframe[2]-d17*sum(J[i17][j]*ecoeff[j] for j in range(18))
    assert last==F(-19,1000)
    cusp_pairs=[ipair(r,nu)]
    for raw in [old['completed_cusp']['point']]+[row['point'] for row in old['preserved_old_bounded_cusps']]:
        cusp_pairs.append(ipair(r,[iv(x) for x in vp(raw,0)]))
    assert len(cusp_pairs)==4 and all(g[1]<0 for g in cusp_pairs)
    paths=[f'scripts/{TAG}.py',newfile,
           'scripts/fresh_h20_affine_a17_completion_20261016_review.py',
           'results/fresh_h20_affine_a17_completion_20261016_review.json',
           TARGET,CAP_REVIEW,INPUT,BASE,CHAIN_REVIEW,PRIOR,PRIOR_REVIEW]+FILES+[
           'scripts/fresh_h20_bridge_edge_closure_20261014_review.py',
           'scripts/fresh_h20_bridge_faces_20261014_review.py',
           'scripts/fresh_h20_s61_multicap_20261014_exact_review.py']
    return dict(status='PASS_INDEPENDENT_141_WALL_RATIONAL_OFF_CUSP_CUTTER',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                norm=str(q),relations=relations,right_neighbors=sorted(finite),
                actual_facets=141,new_center_pair=saved(rh),new_facet_projection_square=saved(new_facet_square),
                E1_cut=str(cut),four_cusp_new_pairs=vec(cusp_pairs),
                complete_new_finite_vertices=0,known_finite_vertices_at_least=170,
                escape=dict(norm=str(eq),all141_pairs=vec(e2pairs),new_pair=str(last),center_pair=saved(ipair(E2,H))),
                proof='Independent140normal reconstruction and exact rational frame/dual recipe; norm and right contacts follow algebraically, all123dotted contacts by strict rational inequalities. Strict center proves all141facets actual. Only17finite neighbors excludes any new20-root finite corner; oldfinitevertices persist by acuteness. All4propercusps preserved. E1cut exactly, but strictE2 remains.',
                scope='Coherent141-wall infinite-volume Coxeter chamber; no complete ambientcone or new dimension. Old chain parameter and oldgeometry remain reviewed dependencies; no producer arithmetic imported.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
