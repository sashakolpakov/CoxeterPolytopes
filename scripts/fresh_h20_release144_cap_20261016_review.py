#!/usr/bin/env python3
"""Independent reconstruction of the release144 E6 cutter.

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
TAG = 'fresh_h20_release144_cap_20261016_review'
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
    clusterfile='results/fresh_h20_residual_h4_corner_20261016.json'
    mixedfile='results/fresh_h20_mixed_isolated_cap_20261016.json'
    cluster=json.loads((ROOT/clusterfile).read_text())
    mixed=json.loads((ROOT/mixedfile).read_text())
    assert cluster['status']=='PASS_EXACT_145_WALL_RESIDUAL_H4_FINITE_CORNER_CLUSTER'
    recipes=[
        [F(17),F(-3,2),F(-67,4),F(-2),F(-3)],
        [F(13),F(0),F(-64,5),F(-21,10),F(-3)],
        [F(2727,100),F(-773,300),F(-267,10),F(-5),F(-5)],
        [F(276851,39846),F(-8739,13282),F(-135452,19923),F(-1),F(-3,2)]]
    assert cluster['definition']['four_recipes']==[list(map(str,v)) for v in recipes]
    # The two isolated duals have exact squares 1/2, are mutually
    # orthogonal, and are orthogonal to the other sixteen support roots.
    assert all(J[i17][j]==F(j==i17,2) and J[i31][j]==F(j==i31,2) for j in range(18))
    metric=[F(-1),F(1),F(1),F(1,2),F(1,2)]
    def dot5(v,w):return sum(a*b*c for a,b,c in zip(metric,v,w))
    def lift5(v):
        return ivadd(ilinear(unitframe,list(map(pt,v[:3]))),
                     ilinear([dual[i17],dual[i31]],list(map(pt,v[3:]))))
    q4=[dot5(v,v) for v in recipes]
    G4=[[dot5(v,w) for w in recipes] for v in recipes]
    assert all(q>0 for q in q4)
    assert G4==[[q4[i] if i==j else F(0) for j in range(4)] for i in range(4)]
    assert list(map(str,q4))==cluster['definition']['four_norms']
    allnormals=full+[r]
    squared=[F(2)]*140+[q]
    S16=set(S)-{17,31}
    records=[];counts=Counter()
    for k,recipe in enumerate(recipes):
        n=lift5(recipe);norm=q4[k]
        rights=S16|set(range(141,141+k))
        if k==1:rights.add(95)
        for j,oldn in enumerate(allnormals):
            value=ipair(n,oldn)
            if j in rights:
                assert R.overlap(value,pt(0));relation='right'
                # Support-dual orthogonality, G4, and the zero e
                # coefficient of wall142 give these exact identities.
            else:
                assert value[1]<0 and value[1]*value[1]>norm*squared[j],(141+k,j,saved(value))
                relation='dotted'
            if k==0:
                row=mixed['pairs'][j]
                assert row['wall']==j and row['relation']==relation
            else:
                idx=sum(range(142,141+k))+j
                row=cluster['pairs'][idx]
                assert row['walls']==[j,141+k] and row['relation']==relation
            overlap(value,row['pair_bounds'])
            records.append(dict(walls=[j,141+k],relation=relation,bounds=saved(value)))
            counts[relation]+=1
        assert R.overlap(ipair(n,n),pt(norm))
        allnormals.append(n);squared.append(norm)
    assert len(records)==570 and len(allnormals)==145
    assert counts==dict(right=71,dotted=499)
    allcenter=[ipair(n,H) for n in allnormals]
    assert all(x[1]<0 for x in allcenter)
    facets=[minus(ipair(H,H),div(times(x,x),pt(nq))) for x,nq in zip(allcenter,squared)]
    assert all(x[1]<0 for x in facets)
    # Solve the four independent rational orthogonality equations afresh
    # with time coordinate one, reversing their order.
    equations=[[metric[j]*v[j] for j in range(1,5)] for v in reversed(recipes)]
    rhs=[v[0] for v in reversed(recipes)]
    J4=B.inv(equations)
    corner=[F(1)]+[sum(a*b for a,b in zip(row,rhs)) for row in J4]
    assert all(dot5(corner,v)==0 for v in recipes)
    assert list(map(str,corner))==cluster['new_finite_corner']['chart']
    cornerq=dot5(corner,corner);assert cornerq<0
    assert str(cornerq)==cluster['new_finite_corner']['norm']
    cv=lift5(corner);incident=sorted(S16|set(range(141,145)))
    cvals=[ipair(cv,n) for n in allnormals]
    for j,value in enumerate(cvals):
        if j in incident:assert R.overlap(value,pt(0))
        else:assert value[1]<0,(j,saved(value))
    assert ipair(cv,H)[1]<0
    g16=[[G[S.index(i)][S.index(j)] for j in sorted(S16)] for i in sorted(S16)]
    pivots=B.pd(g16)
    assert len(pivots)==16
    # Complete finite-neighbor lists, in both directions, are determined
    # by the strict pair checks above. Only the displayed corner survives.
    neighbors={i:S16|(set(range(141,145))-{i})|({95} if i==142 else set()) for i in range(141,145)}
    positive=[];rejected=[]
    for i,ns in neighbors.items():
        for co in combinations(sorted(ns),19):
            support=sorted(set(co)|{i})
            if 95 in support:
                witness=next(j for j in (141,143,144) if j in support)
                z=ipair(allnormals[95],allnormals[witness])
                assert z[1]<0 and z[1]*z[1]>squared[95]*squared[witness]
                rejected.append(dict(support=support,dotted_pair=[95,witness]))
            else:
                assert support==incident
                positive.append(support)
    assert len(positive)==4 and len(rejected)==19
    # Retain all four proper cusp links: every new wall is strict there.
    cusp_vectors=[nu]+[[iv(x) for x in vp(raw,0)] for raw in
        [old['completed_cusp']['point']]+[row['point'] for row in old['preserved_old_bounded_cusps']]]
    cusp_checks=[[ipair(n,c) for c in cusp_vectors] for n in allnormals[141:]]
    assert all(x[1]<0 for row in cusp_checks for x in row)
    # Independent whole-rectangle proof: wall pairings are affine and
    # the exact norm is constant plus (a^2+b^2)/2.
    ecoeff4=[F(-1,100000)]*18;ecoeff4[i17]=F(-3,10);ecoeff4[i31]=F(-1,20)
    def esc(aa,bb):
        cs=ecoeff4[:];cs[i17]=-aa;cs[i31]=-bb
        v=ivadd(ilinear(unitframe,list(map(pt,eframe))),ilinear(dual,list(map(pt,cs))))
        norm=-eframe[0]**2+eframe[1]**2+eframe[2]**2+sum(cs[i]*J[i][j]*cs[j] for i in range(18) for j in range(18))
        pairs=[ipair(v,n) for n in allnormals]
        assert norm>0 and all(x[1]<0 for x in pairs) and ipair(v,H)[1]<0
        return dict(a=str(aa),b=str(bb),norm=str(norm),all145_pairs=vec(pairs),center_pair=saved(ipair(v,H)))
    E4=esc(F(3,10),F(1,20))
    assert E4['norm']==cluster['strict_escape']['norm']=='33750061/5000000000'
    rectangle=[esc(aa,bb) for aa in (F(29,100),F(31,100)) for bb in (F(3,100),F(3,50))]
    minimum=min(F(row['norm']) for row in rectangle)
    assert minimum==F(15000061,5000000000)
    assert str(minimum)==cluster['residual_escape_rectangle']['norm_minimum']
    partnerfile='results/fresh_h20_h4_parallel_partner_20261016.json'
    partner=json.loads((ROOT/partnerfile).read_text())
    assert partner['status']=='PASS_EXACT_146_WALL_FINITE_PARALLEL_E4_CUTTER'
    kk=-cornerq;qs=q4[0];qw=q4[3]
    nu2=ivadd(ivscale(cv,div(pt(1),root(pt(kk)))),
              ivscale(allnormals[144],neg(div(pt(1),root(pt(qw))))))
    Avec=ivscale(allnormals[141],neg(root(pt(2/qs))))
    Bvec=ivscale(nu2,root(pt(2)))
    apair,bpair=ipair(Avec,allnormals[140]),ipair(Bvec,allnormals[140])
    target=neg(root(pt(squared[140]/2)))
    assert apair[0]>0 and bpair[1]<0
    zz=div(minus(target,apair),bpair);assert zz[0]>1
    rr=ivadd(Avec,ivscale(Bvec,zz))
    assert R.overlap(ipair(rr,rr),pt(2))
    # Exact orthogonal chart identities give A^2=2,B^2=0,A.B=0,
    # hence r^2=2; the defining quotient gives the pi/3 contact.
    assert R.overlap(ipair(nu2,nu2),pt(0)) and ipair(nu2,H)[1]<0
    rights=S16|{142,143};newpairs=[];relations=Counter()
    for j,n in enumerate(allnormals):
        value=ipair(rr,n)
        if j in rights:
            relation='right';assert R.overlap(value,pt(0))
        elif j==140:
            relation='pi/3';assert R.overlap(value,target)
        elif j==141:
            relation='parallel';assert R.overlap(value,neg(root(pt(2*qs))))
        else:
            relation='dotted'
            assert value[1]<0 and value[1]*value[1]>2*squared[j],(j,saved(value))
        row=partner['pairs'][j]
        assert row['wall']==j and row['relation']==relation
        overlap(value,row['pair_bounds'])
        newpairs.append(dict(wall=j,relation=relation,bounds=saved(value)))
        relations[relation]+=1
    assert relations==dict(right=18,**{'pi/3':1,'parallel':1,'dotted':125})
    for value,raw in zip(rr,partner['definition']['physical_normal_bounds']):overlap(value,raw)
    all146=allnormals+[rr];norms146=squared+[F(2)]
    center146=[ipair(n,H) for n in all146]
    assert all(v[1]<0 for v in center146)
    facets146=[minus(ipair(H,H),div(times(v,v),pt(qn))) for v,qn in zip(center146,norms146)]
    assert all(v[1]<0 for v in facets146)
    ideal_support=sorted(S16|{141,142,143,145})
    assert len(ideal_support)==20 and ideal_support==partner['new_ideal_vertex']['incident']
    nullpairs=[ipair(nu2,n) for n in all146]
    for j,v in enumerate(nullpairs):
        if j in ideal_support:assert R.overlap(v,pt(0))
        else:assert v[1]<0,(j,saved(v))
    # The18 positive spectator normals plus independent s,r give normal
    # rank20. Their affineA1 positive relation leaves Gramrank19 and this
    # unique nullray. It is an actual ideal vertex, not yet a bounded cusp.
    parallel_sum=ivadd(ivscale(rr,div(pt(1),root(pt(2)))),
                       ivscale(allnormals[141],div(pt(1),root(pt(qs)))))
    assert all(R.overlap(x,y) for x,y in zip(parallel_sum,ivscale(nu2,zz)))
    p16=sorted(S16)
    components=[];unused=set(p16)
    while unused:
        component={min(unused)};unused-=component
        frontier=list(component)
        while frontier:
            i=frontier.pop()
            adjacent={j for j in unused if G[S.index(i)][S.index(j)]!=0}
            component|=adjacent;unused-=adjacent;frontier+=sorted(adjacent)
        components.append(sorted(component))
    assert sorted(map(len,components))==[2,4,10]
    for comp in components:
        mat=[[G[S.index(i)][S.index(j)] for j in comp] for i in comp]
        assert all(mat[i][j] in (F(0),F(-1)) for i in range(len(comp)) for j in range(i))
        degrees=[sum(x==-1 for j,x in enumerate(row) if j!=i) for i,row in enumerate(mat)]
        assert sorted(degrees)==[1,1]+[2]*(len(comp)-2)
    finite_support=sorted(rights|{140})
    assert len(finite_support)==19 and finite_support==partner['finite_old_neighbors']
    obstructing=ipair(allnormals[140],allnormals[142])
    assert obstructing[1]<0 and obstructing[1]**2>squared[140]*squared[142]
    e4point=ivadd(ilinear(unitframe,list(map(pt,eframe))),ilinear(dual,list(map(pt,ecoeff4))))
    cut=ipair(rr,e4point);assert cut[0]>0
    overlap(cut,partner['E4_cut_pair_bounds'])
    ecoeff5=ecoeff4[:];ecoeff5[i31]=F(-1,25)
    E5=ivadd(ilinear(unitframe,list(map(pt,eframe))),ilinear(dual,list(map(pt,ecoeff5))))
    eq5=-eframe[0]**2+eframe[1]**2+eframe[2]**2+sum(ecoeff5[i]*J[i][j]*ecoeff5[j] for i in range(18) for j in range(18))
    assert eq5==F(31500061,5000000000) and str(eq5)==partner['strict_escape']['norm']
    e5pairs=[ipair(E5,n) for n in all146]
    assert all(x[1]<0 for x in e5pairs) and ipair(E5,H)[1]<0
    oldcuspchecks=[ipair(rr,c) for c in cusp_vectors]
    assert len(oldcuspchecks)==4 and all(v[1]<0 for v in oldcuspchecks)
    completionfile='results/fresh_h20_second_cusp_completion_20261016.json'
    closure=json.loads((ROOT/completionfile).read_text())
    assert closure['status']=='PASS_EXACT_151_WALL_SECOND_PROPER_CUSP_COMPLETION'
    spherical=[list(range(1,11)),list(range(12,16)),[61,87],[142],[143]]
    assert spherical==closure['definition']['spherical_components']
    assert set(sum(spherical,[]))==S16|{142,143}
    heights=[3,3,3,5,6]
    assert heights==closure['definition']['heights']
    newnorms=[F(2),F(2),F(2),q4[1],q4[2]]
    added=[]
    for comp,height,nq in zip(spherical,heights,newnorms):
        opposite=ilinear([all146[j] for j in comp],[pt(-1)]*len(comp))
        value=ivadd(opposite,ivscale(nu2,pt(height)))
        assert R.overlap(ipair(value,value),pt(nq))
        # Exact A_n path identities (or the literal A1 normal) give
        # opposite^2=nq. Orthogonality to nu2 makes the shift norm-neutral.
        if len(comp)>1:
            for i,a in enumerate(comp):
                for j,b in enumerate(comp):
                    expected=2 if i==j else -1 if abs(i-j)==1 else 0
                    assert G[S.index(a)][S.index(b)]==expected
        added.append(value)
    newrelations=[];counts151=Counter()
    oldincident=set(ideal_support)
    for i,(v,nq,comp) in enumerate(zip(added,newnorms,spherical)):
        for j,n in enumerate(all146+added[:i]):
            value=ipair(v,n)
            if j>=146:
                relation='right';expected=pt(0)
            elif j in oldincident:
                if len(comp)==1 and j==comp[0]:
                    relation='parallel';expected=pt(-nq)
                elif len(comp)>1 and j in (comp[0],comp[-1]):
                    relation='pi/3';expected=pt(-1)
                else:
                    relation='right';expected=pt(0)
            else:
                relation='dotted';expected=None
            if expected is not None:assert R.overlap(value,expected)
            else:assert value[1]<0 and value[1]*value[1]>nq*norms146[j],(146+i,j,saved(value))
            offset=sum(146+k for k in range(i))+j
            raw=closure['new_pairs'][offset]
            assert raw['walls']==[j,146+i] and raw['relation']==relation
            overlap(value,raw['pair_bounds'])
            newrelations.append(dict(walls=[j,146+i],relation=relation,bounds=saved(value)))
            counts151[relation]+=1
    assert len(newrelations)==740
    assert counts151==dict(right=102,**{'pi/3':6,'parallel':2,'dotted':630})
    full151=all146+added;norms151=norms146+newnorms
    center151=[ipair(v,H) for v in full151]
    assert all(v[1]<0 for v in center151)
    facets151=[minus(ipair(H,H),div(times(v,v),pt(nq))) for v,nq in zip(center151,norms151)]
    assert all(v[1]<0 for v in facets151)
    cusp_components=[comp+[146+i] for i,comp in enumerate(spherical)]+[[141,145]]
    cusp_support=sorted(sum(cusp_components,[]));assert len(cusp_support)==len(set(cusp_support))==25
    assert cusp_support==closure['new_proper_cusp']['incident']
    index={wall:i for i,wall in enumerate(cusp_support)}
    gram25=[[F(0)]*25 for _ in range(25)]
    for comp in cusp_components:
        for i,a in enumerate(comp):
            for j,b in enumerate(comp):
                value=2 if i==j else -2 if len(comp)==2 else -1 if (i-j)%len(comp) in (1,len(comp)-1) else 0
                gram25[index[a]][index[b]]=F(value)
    M=sy.Matrix([[sy.Rational(x.numerator,x.denominator) for x in row] for row in gram25])
    # Distinct finite A-components become literal affine cycles. Their
    # full nullspaces, rather than numerical eigensigns, certify the rank.
    kernels=[];component_ranks=[]
    for comp in cusp_components:
        ids=[index[j] for j in comp];sub=M.extract(ids,ids)
        assert sub*sy.ones(len(comp),1)==sy.zeros(len(comp),1)
        assert sub.rank()==len(comp)-1
        assert sub[:len(comp)-1,:len(comp)-1].is_positive_definite
        kernels.append([1 if j in ids else 0 for j in range(25)])
        component_ranks.append(len(comp)-1)
    assert sum(component_ranks)==19 and M.rank()==19 and len(M.nullspace())==6
    assert [[str(x) for x in row] for row in gram25]==closure['new_proper_cusp']['normalized_gram']
    for i,a in enumerate(cusp_support):
        for j,b in enumerate(cusp_support):
            normalized=div(times(pt(2),ipair(full151[a],full151[b])),root(pt(norms151[a]*norms151[b])))
            assert R.overlap(normalized,pt(gram25[i][j]))
    # Each kernel maps to the same future null ray by exact recipe sums.
    # These are algebraic identities; interval overlaps only check arithmetic.
    multipliers=[]
    for i,comp in enumerate(cusp_components):
        coefficients=[root(pt(2/norms151[j])) for j in comp]
        nullsum=ilinear([full151[j] for j in comp],coefficients)
        mult=times(pt(heights[i]),root(pt(2/newnorms[i]))) if i<5 else times(root(pt(2)),zz)
        assert mult[0]>0
        assert all(R.overlap(a,b) for a,b in zip(nullsum,ivscale(nu2,mult)))
        multipliers.append(saved(mult))
    null151=[ipair(nu2,v) for v in full151]
    for i,v in enumerate(null151):
        if i in cusp_support:assert R.overlap(v,pt(0))
        else:assert v[1]<0
    retained_cusps=[[ipair(v,c) for v in added] for c in cusp_vectors]
    assert all(v[1]<0 for row in retained_cusps for v in row)
    e151=[ipair(E5,n) for n in full151]
    assert all(v[1]<0 for v in e151) and ipair(E5,H)[1]<0
    assert str(eq5)==closure['strict_escape']['norm']
    # Every new root is dotted to all126 walls outside the cusp. Any
    # finite20-root support containing it would lie inside this rank19
    # positive-semidefinite star, which cannot contain a PD20 block.
    assert len(set(range(146))-oldincident)==126

    capfile='results/fresh_h20_second_cusp_pi8_cap_20261016.json'
    cap=json.loads((ROOT/capfile).read_text())
    assert cap['status']=='PASS_EXACT_ALTERNATIVE_152_WALL_PI8_CUSP_EDGE_CAP'
    def interval_inverse(matrix):
        size=len(matrix)
        work=[list(row)+[pt(int(i==j)) for j in range(size)] for i,row in enumerate(matrix)]
        pivots=[]
        for k in range(size):
            pivot=next(i for i in range(k,size) if work[i][k][0]>0 or work[i][k][1]<0)
            work[k],work[pivot]=work[pivot],work[k]
            p=work[k][k];pivots.append(saved(p))
            work[k]=[div(x,p) for x in work[k]]
            for i in range(size):
                if i==k:continue
                t=work[i][k]
                work[i]=[minus(x,times(t,y)) for x,y in zip(work[i],work[k])]
        return [row[size:] for row in work],pivots
    def matvec(matrix,v):
        return [total(times(a,b) for a,b in zip(row,v)) for row in matrix]
    U=[full151[j] for j in (95,142,145)]
    GU=[[ipair(v,w) for w in U] for v in U]
    assert GU[0][2][1]<-2
    JU,upivots=interval_inverse(GU)
    def project(v):
        return ilinear(U,matvec(JU,[ipair(v,u) for u in U]))
    vv=ivadd(full151[150],ivscale(project(full151[150]),pt(-1)))
    vq=ipair(vv,vv);assert vq[0]>0
    raw=dual[i31]
    kval=ivadd(raw,ivscale(project(raw),pt(-1)))
    kval=ivadd(kval,ivscale(vv,neg(div(ipair(kval,vv),vq))))
    kq=ipair(kval,kval);assert kq[0]>0
    cosine=div(root(plus(pt(2),root(pt(2)))),pt(2))
    target8=neg(times(root(pt(2*norms151[150])),cosine))
    projection=ivscale(vv,div(target8,vq))
    rho=div(minus(pt(2),div(times(target8,target8),vq)),kq)
    assert rho[0]>0
    rootrho=root(rho)
    rights152=S16|{95,142,145}
    trials=[];accepted=[]
    for sign in (-1,1):
        value=ivadd(projection,ivscale(kval,times(pt(sign),rootrho)))
        allpairs=[ipair(value,v) for v in full151]
        bad=[];records=[]
        for j,pair in enumerate(allpairs):
            if j in rights152:
                relation='right';good=R.overlap(pair,pt(0))
            elif j==150:
                relation='pi/8';good=R.overlap(pair,target8)
            else:
                relation='dotted';good=pair[1]<0 and pair[1]**2>2*norms151[j]
            if not good:bad.append(j)
            records.append(dict(wall=j,relation=relation,bounds=saved(pair)))
        hp=ipair(value,H);cut=ipair(value,E5)
        trials.append(dict(sign=sign,bad_walls=bad,center=saved(hp),E5=saved(cut)))
        if not bad and hp[1]<0<cut[0]:accepted.append((sign,value,records,hp,cut))
    assert len(accepted)==1
    sign,new,records,hnew,cut=accepted[0]
    assert sign==cap['definition']['sign']
    assert R.overlap(ipair(new,new),pt(2))
    # The projection identities, not interval overlap, establish every
    # prescribed zero, the target8 equality and squared norm2 exactly.
    assert Counter(row['relation'] for row in records)==dict(right=19,**{'pi/8':1,'dotted':131})
    for row,raw in zip(records,cap['new_pairs']):
        assert (row['wall'],row['relation'])==(raw['wall'],raw['relation'])
        overlap(tuple(map(F,row['bounds'])),raw['pair_bounds'])
    for x,raw in zip(new,cap['definition']['physical_normal_bounds']):overlap(x,raw)
    full152=full151+[new];norms152=norms151+[F(2)]
    hp152=[ipair(H,v) for v in full152]
    assert all(v[1]<0 for v in hp152)
    facets152=[minus(ipair(H,H),div(times(v,v),pt(nq))) for v,nq in zip(hp152,norms152)]
    assert all(v[1]<0 for v in facets152)
    cuspnew=[ipair(new,v) for v in cusp_vectors+[nu2]]
    assert len(cuspnew)==5 and all(v[1]<0 for v in cuspnew)
    # Separate generic four-by-four inversion for the new corner,
    # after removing the already verified positive S16 component.
    small=[142,145,150,151]
    GS=[[pt(norms152[a]) if a==b else target8 if {a,b}=={150,151} else pt(0) for b in small] for a in small]
    check=[row[:] for row in GS];smallpivots=[]
    for k in range(4):
        pivot=check[k][k];assert pivot[0]>0;smallpivots.append(saved(pivot))
        for i in range(k+1,4):
            for j in range(i,4):
                check[j][i]=check[i][j]=minus(check[i][j],div(times(check[i][k],check[k][j]),pivot))
    for i,a in enumerate(small):
        for j,b in enumerate(small):assert R.overlap(ipair(full152[a],full152[b]),GS[i][j])
    JS,jspivots=interval_inverse(GS)
    hsmall=[ipair(H,full152[j]) for j in small]
    point=ivadd(H,ivscale(ilinear([full152[j] for j in small],matvec(JS,hsmall)),pt(-1)))
    point=ivadd(point,ivscale(ilinear([dual[S.index(j)] for j in sorted(S16)],
                  [ipair(H,full152[j]) for j in sorted(S16)]),pt(-1)))
    pq=ipair(point,point);ph=ipair(point,H)
    assert pq[1]<0 and ph[1]<0
    incidence=sorted(S16|set(small));assert len(incidence)==20
    pointpairs=[ipair(point,v) for v in full152]
    for j,v in enumerate(pointpairs):
        if j in incidence:assert R.overlap(v,pt(0))
        else:assert v[1]<0,(j,saved(v))
    assert incidence==cap['new_finite_corner']['incident']
    finite_neighbors=sorted(rights152|{150});assert len(finite_neighbors)==20
    rejections=[]
    for omitted in finite_neighbors:
        if omitted==95:continue
        support=set(finite_neighbors)-{omitted}
        witness=None
        for a,b in combinations(sorted(support),2):
            value=ipair(full151[a],full151[b])
            if value[1]<0 and value[1]**2>norms151[a]*norms151[b]:
                witness=[a,b];break
        assert witness is not None,omitted
        rejections.append(dict(omitted=omitted,dotted_pair=witness))
    assert len(rejections)==19
    shared=sorted(set(incidence)&set(cusp_support))
    assert shared==sorted(S16|{142,145,150}) and len(shared)==19
    assert shared==cap['new_cusp_edge']['shared_support']
    # The shared block is S16 plus3orthogonal positive roots. The
    # feasible finite and ideal endpoints are distinct and every other
    # inequality is strict inside their positive open segment.
    for j in range(152):
        if j in shared:continue
        nullvalue=ipair(nu2,full152[j])
        assert pointpairs[j][1]<0 or nullvalue[1]<0,j
    coeff6=ecoeff5[:];coeff6[i17]=F(-7,25);coeff6[i31]=F(-1,20)
    E6=ivadd(ilinear(unitframe,list(map(pt,eframe))),ilinear(dual,list(map(pt,coeff6))))
    q6=-eframe[0]**2+eframe[1]**2+eframe[2]**2+sum(coeff6[i]*J[i][j]*coeff6[j] for i in range(18) for j in range(18))
    assert q6>0 and str(q6)==cap['strict_escape']['norm']
    e6pairs=[ipair(E6,v) for v in full152]
    assert all(v[1]<0 for v in e6pairs) and ipair(E6,H)[1]<0

    trianglefile='results/fresh_h20_pi8_triangle_face_20261016.json'
    triangle=json.loads((ROOT/trianglefile).read_text())
    assert triangle['status']=='PASS_EXACT_153_WALL_COMPACT_246_TRIANGLE_FACE'
    new152=ilinear([full152[142],full152[145],point],
        [neg(div(pt(1),root(pt(norms152[142])))),neg(div(root(pt(3)),pt(2))),
         div(pt(1),root(times(pt(-2),pq)))])
    assert R.overlap(ipair(new152,new152),pt(2))
    # Orthogonal recipe summands have squares1,3/2,-1/2. This proves
    # norm2 exactly; all prescribed contact identities follow as well.
    assert F(1)+F(3,2)-F(1,2)==2
    targets={142:neg(root(pt(norms152[142]))),145:neg(root(pt(3)))}
    face_support=sorted(S16|{150,151});assert len(face_support)==18
    newpairs=[]
    for j,oldnormal in enumerate(full152):
        v=ipair(new152,oldnormal)
        if j in face_support:
            rel='right';assert R.overlap(v,pt(0))
        elif j in targets:
            rel={142:'pi/4',145:'pi/6'}[j];assert R.overlap(v,targets[j])
        else:
            rel='dotted';assert v[1]<0 and v[1]**2>2*norms152[j],(j,saved(v))
        oldrow=triangle['new_pairs'][j]
        assert (j,rel)==(oldrow['wall'],oldrow['relation'])
        overlap(v,oldrow['pair_bounds'])
        newpairs.append(dict(wall=j,relation=rel,bounds=saved(v)))
    counts=Counter(row['relation'] for row in newpairs)
    assert counts==dict(right=18,**{'pi/4':1,'pi/6':1,'dotted':132})
    for value,raw in zip(new152,triangle['definition']['physical_normal_bounds']):overlap(value,raw)
    full153=full152+[new152];norms153=norms152+[F(2)]
    hp153=[ipair(H,v) for v in full153];assert all(v[1]<0 for v in hp153)
    facets153=[minus(ipair(H,H),div(times(v,v),pt(nq))) for v,nq in zip(hp153,norms153)]
    assert all(v[1]<0 for v in facets153)
    cusps153=[ipair(new152,v) for v in cusp_vectors+[nu2]]
    assert len(cusps153)==5 and all(v[1]<0 for v in cusps153)
    # Each corner is reconstructed independently by a generic4x4 Gram
    # inverse, not by the producer's H18/two-root projection formula.
    corners=[]
    for neighbor in (142,145):
        smallids=[150,151,neighbor,152]
        G4=[[pt(norms153[a]) if a==b else target8 if {a,b}=={150,151}
             else targets[neighbor] if {a,b}=={neighbor,152} else pt(0)
             for b in smallids] for a in smallids]
        schur=[row[:] for row in G4];pivots=[]
        for k in range(4):
            t=schur[k][k];assert t[0]>0;pivots.append(saved(t))
            for i in range(k+1,4):
                for j in range(i,4):
                    schur[j][i]=schur[i][j]=minus(schur[i][j],div(times(schur[i][k],schur[k][j]),t))
        inverse,_=interval_inverse(G4)
        hsmall=[ipair(H,full153[j]) for j in smallids]
        vertex=ivadd(H,ivscale(ilinear([full153[j] for j in smallids],matvec(inverse,hsmall)),pt(-1)))
        vertex=ivadd(vertex,ivscale(ilinear([dual[S.index(j)] for j in sorted(S16)],
                     [ipair(H,full153[j]) for j in sorted(S16)]),pt(-1)))
        norm=ipair(vertex,vertex);time=ipair(vertex,H)
        assert norm[1]<0 and time[1]<0
        ids=sorted(set(face_support)|{neighbor,152});assert len(ids)==20
        vertexpairs=[ipair(vertex,v) for v in full153]
        for j,v in enumerate(vertexpairs):
            if j in ids:assert R.overlap(v,pt(0))
            else:assert v[1]<0,(neighbor,j,saved(v))
        targetcorner=next(row for row in triangle['new_finite_corners'] if row['side_pair']==[neighbor,152])
        assert ids==targetcorner['incident']
        for value,raw in zip(vertex,targetcorner['point_bounds']):overlap(value,raw)
        corners.append(dict(incident=ids,side_pair=[neighbor,152],point=vec(vertex),
            norm=saved(norm),center_pair=saved(time),all153_pairs=vec(vertexpairs),small_PD_pivots=pivots))
    neighbors=sorted(set(face_support)|{142,145});assert len(neighbors)==20
    rejections=[]
    # Literal normalized3Gram has determinant8-2*2-2*3=-2.
    # Its proper2minors are4,2,1, so it is the compactLanner246triangle.
    determinant=F(8)-2*F(2)-2*F(3);assert determinant==-2
    for omitted in neighbors:
        if omitted in (142,145):continue
        assert {142,145}<=set(neighbors)-{omitted}
        rejections.append(dict(omitted=omitted,indefinite_subset=[142,145,152],determinant='-2'))
    assert len(rejections)==18
    oldpointnew=ipair(point,new152);assert oldpointnew[1]<0
    # S18 is positive: S16 plus the exact I2(8) pair. The3side normals
    # are perpendicular to it and independent by the negative Gramdet.
    # All3finite corners satisfy all153inequalities, so the whole compact
    # triangle and no smaller truncation is the actual codimension18face.
    pairdet=minus(pt(2*norms153[150]),times(target8,target8))
    assert pairdet[0]>0
    oldpointpairs=pointpairs+[oldpointnew]
    assert all(oldpointpairs[j][1]<0 for j in range(153) if j not in incidence)
    alltrianglevertices=[dict(incident=incidence,point=vec(point),norm=saved(pq),all153_pairs=vec(oldpointpairs))]+corners
    assert all(set(face_support)<=set(row['incident']) for row in alltrianglevertices)
    assert {tuple(sorted(set(row['incident'])-set(face_support))) for row in alltrianglevertices}=={(142,145),(142,152),(145,152)}
    e153=[ipair(E6,v) for v in full153];assert all(v[1]<0 for v in e153)
    assert str(q6)==triangle['strict_escape']['norm']

    axisfile='results/fresh_h20_dominant142_cap_20261016.json'
    axis=json.loads((ROOT/axisfile).read_text())
    assert axis['status']=='PASS_EXACT_154_WALL_DOMINANT_AXIS_PI5_CAP'
    qu=norms153[142];assert qu==F(309,200)
    gamma=div(plus(pt(1),root(pt(5))),pt(4))
    # Use the general quadratic formula and a separate positive radical
    # for s, rather than the producer's simplified formula for c and s.
    qa=pt(25-qu);qb=times(pt(-2*qu),gamma)
    qc=neg(plus(pt(25),times(pt(qu),times(gamma,gamma))))
    disc=minus(times(qb,qb),times(pt(4),times(qa,qc)))
    assert qa[0]>0 and disc[0]>0
    boost=div(plus(neg(qb),root(disc)),times(pt(2),qa))
    assert boost[0]>1
    shift=root(minus(times(boost,boost),pt(1)));assert shift[0]>0
    anchor_shift=div(times(root(pt(qu)),plus(boost,gamma)),pt(5))
    assert R.overlap(shift,anchor_shift)
    # Quadratic identity and positivity prove s=sqrt(q)*(c+gamma)/5
    # exactly; hence the old149 contact is precisely pi/5.
    new153=ilinear([full153[142],cv],
                  [neg(times(root(pt(2/qu)),boost)),times(root(pt(2/(-cornerq))),shift)])
    assert R.overlap(ipair(new153,new153),pt(2))
    target5=neg(times(root(pt(2*qu)),gamma))
    rightaxis=S16|{141,143,144}
    axis_pairs=[]
    for j,n in enumerate(full153):
        value=ipair(new153,n)
        if j in rightaxis:
            relation='right';assert R.overlap(value,pt(0))
        elif j==149:
            relation='pi/5';assert R.overlap(value,target5)
        else:
            relation='dotted';assert value[1]<0 and value[1]**2>2*norms153[j],(j,saved(value))
        raw=axis['new_pairs'][j]
        assert (j,relation)==(raw['wall'],raw['relation'])
        overlap(value,raw['pair_bounds'])
        axis_pairs.append(dict(wall=j,relation=relation,bounds=saved(value)))
    counts154=Counter(row['relation'] for row in axis_pairs)
    assert counts154==dict(right=19,**{'pi/5':1,'dotted':133})
    for value,raw in zip(new153,axis['definition']['physical_normal_bounds']):overlap(value,raw)
    full154=full153+[new153];norms154=norms153+[F(2)]
    hp154=[ipair(H,v) for v in full154];assert all(v[1]<0 for v in hp154)
    facets154=[minus(ipair(H,H),div(times(v,v),pt(nq))) for v,nq in zip(hp154,norms154)]
    assert all(v[1]<0 for v in facets154)
    cusps154=[ipair(new153,v) for v in cusp_vectors+[nu2]]
    assert len(cusps154)==5 and all(v[1]<0 for v in cusps154)
    axis_corners=[]
    for neighbor in (144,149):
        smallids=[141,143,neighbor,153]
        G4=[[pt(norms154[a]) if a==b else target5 if {a,b}=={149,153} else pt(0)
             for b in smallids] for a in smallids]
        schur=[row[:] for row in G4];pivots=[]
        for k in range(4):
            t=schur[k][k];assert t[0]>0;pivots.append(saved(t))
            for i in range(k+1,4):
                for j in range(i,4):
                    schur[j][i]=schur[i][j]=minus(schur[i][j],div(times(schur[i][k],schur[k][j]),t))
        inverse,_=interval_inverse(G4)
        hsmall=[ipair(H,full154[j]) for j in smallids]
        v=ivadd(H,ivscale(ilinear([full154[j] for j in smallids],matvec(inverse,hsmall)),pt(-1)))
        v=ivadd(v,ivscale(ilinear([dual[S.index(j)] for j in sorted(S16)],
                   [ipair(H,full154[j]) for j in sorted(S16)]),pt(-1)))
        vn,vh=ipair(v,v),ipair(v,H);assert vn[1]<0 and vh[1]<0
        support=sorted(S16|set(smallids));assert len(support)==20
        vpairs=[ipair(v,n) for n in full154]
        for j,value in enumerate(vpairs):
            if j in support:assert R.overlap(value,pt(0))
            else:assert value[1]<0,(neighbor,j,saved(value))
        raw=next(row for row in axis['new_finite_corners'] if row['side_pair']==[neighbor,153])
        assert support==raw['incident']
        for value,rvalue in zip(v,raw['point_bounds']):overlap(value,rvalue)
        axis_corners.append(dict(incident=support,side_pair=[neighbor,153],point=vec(v),
            norm=saved(vn),center_pair=saved(vh),all154_pairs=vec(vpairs),small_positive_pivots=pivots))
    neighbors154=sorted(rightaxis|{149});assert len(neighbors154)==20
    obstructing=ipair(full153[144],full153[149])
    assert obstructing[1]<0 and obstructing[1]**2>norms153[144]*norms153[149]
    rejections154=[]
    for omit in neighbors154:
        if omit in (144,149):continue
        assert {144,149}<=set(neighbors154)-{omit}
        rejections154.append(dict(omitted=omit,dotted_pair=[144,149]))
    assert len(rejections154)==18
    triangle_persistence=[]
    for record in alltrianglevertices:
        v=[tuple(map(F,x)) for x in record['point']]
        value=ipair(v,new153);assert value[1]<0
        triangle_persistence.append(saved(value))
    e154=[ipair(E6,v) for v in full154];assert all(v[1]<0 for v in e154)
    assert str(q6)==axis['strict_escape']['norm']
    cutoff=div(shift,boost);assert 0<cutoff[0]<=cutoff[1]<1
    overlap(cutoff,axis['definition']['axis_cutoff_bounds'])

    pairfile='results/fresh_h20_dominant_axis_cluster_20261016.json'
    pairdoc=json.loads((ROOT/pairfile).read_text())
    assert pairdoc['status']=='PASS_EXACT_156_WALL_COORDINATED_DOMINANT_AXIS_CLUSTER'
    dependency_check=integrity([pairfile,'results/fresh_h20_dominant142_cap_20261016_review.json'])
    # These exact positive radicals give norm2 and product -1 with wall153.
    # They are reconstructed from the independently derived previous boost.
    sstar=div(pt(1),times(pt(2),shift));cstar=root(plus(pt(1),times(sstar,sstar)))
    assert sstar[0]>1
    overlap(sstar,pairdoc['definition']['s_star_bounds'])
    overlap(cstar,pairdoc['definition']['c_star_bounds'])
    newpair=[]
    for a in (141,143):
        n=ilinear([full154[a],cv],
            [neg(times(root(pt(2/norms154[a])),cstar)),times(root(pt(2/(-cornerq))),sstar)])
        assert R.overlap(ipair(n,n),pt(2))
        newpair.append(n)
    full156=full154+newpair;norms156=norms154+[F(2),F(2)]
    mutual=ipair(newpair[0],newpair[1]);assert mutual[1]<-2
    assert R.overlap(mutual,neg(times(pt(2),times(sstar,sstar))))
    overlap(mutual,pairdoc['definition']['mutual_pair_bounds'])
    contactrows=[];corners156=[];census156=[]
    primary_pairs={tuple(x['walls']):x for x in pairdoc['new_pairs']}
    for ni,(a,n) in enumerate(zip((141,143),newpair),154):
        rights=S16|({141,142,143,144}-{a});assert len(rights)==19
        for j in range(ni):
            value=ipair(n,full156[j])
            if j in rights:
                rel='right';assert R.overlap(value,pt(0))
            elif j==153:
                rel='pi/3';assert R.overlap(value,pt(-1))
            else:
                rel='dotted';assert value[1]<0 and value[1]**2>2*norms156[j],(ni,j,saved(value))
            raw=primary_pairs[j,ni];assert rel==raw['relation'];overlap(value,raw['pair_bounds'])
            contactrows.append(dict(walls=[j,ni],relation=rel,bounds=saved(value)))
        raw=next(x for x in pairdoc['definition']['new_roots'] if x['new_wall']==ni)
        for x,y in zip(n,raw['physical_normal_bounds']):overlap(x,y)
        other_axes=sorted({141,142,143,144}-{a,142})
        for neighbor in (142,153):
            ids=other_axes+[neighbor,ni]
            assert len(set(ids))==4
            small=[[pt(norms156[x]) if x==y else pt(-1) if {x,y}=={153,ni} else pt(0)
                    for y in ids] for x in ids]
            schur=[row[:] for row in small];pivots=[]
            for k in range(4):
                q=schur[k][k];assert q[0]>0;pivots.append(saved(q))
                for i in range(k+1,4):
                    for j in range(i,4):
                        schur[i][j]=schur[j][i]=minus(schur[i][j],div(times(schur[i][k],schur[k][j]),q))
            inv,_=interval_inverse(small)
            hs=[ipair(H,full156[j]) for j in ids]
            v=ivadd(H,ivscale(ilinear([full156[j] for j in ids],matvec(inv,hs)),pt(-1)))
            v=ivadd(v,ivscale(ilinear([dual[S.index(j)] for j in sorted(S16)],
                [ipair(H,full156[j]) for j in sorted(S16)]),pt(-1)))
            support=sorted(S16|set(ids));assert len(support)==20
            vn=ipair(v,v);vh=ipair(v,H);assert vn[1]<0 and vh[1]<0
            pairs=[ipair(v,w) for w in full156]
            for j,val in enumerate(pairs):
                if j in support:assert R.overlap(val,pt(0))
                else:assert val[1]<0,(ni,neighbor,j,saved(val))
            raw=next(x for x in pairdoc['new_finite_corners'] if x['side_pair']==[neighbor,ni])
            assert support==raw['incident']
            for x,y in zip(v,raw['point_bounds']):overlap(x,y)
            corners156.append(dict(incident=support,side_pair=[neighbor,ni],point=vec(v),
                norm=saved(vn),center_pair=saved(vh),all156_pairs=vec(pairs),small_gram_positive_pivots=pivots))
        neighbors=sorted(rights|{153});assert len(neighbors)==20
        obstruction=ipair(full156[142],full156[153])
        assert obstruction[1]<0 and obstruction[1]**2>2*norms156[142]
        rejected=[]
        for omitted in neighbors:
            if omitted in (142,153):continue
            assert {142,153}<=set(neighbors)-{omitted}
            rejected.append(dict(omitted=omitted,dotted_pair=[142,153]))
        assert len(rejected)==18
        census156.append(dict(new_wall=ni,finite_neighbors=neighbors,rejected_supports=rejected))
    counts156=Counter(x['relation'] for x in contactrows)
    assert counts156==dict(right=38,**{'pi/3':2,'dotted':269}) and len(contactrows)==309
    hp156=[ipair(H,n) for n in full156];assert all(x[1]<0 for x in hp156)
    facet156=[minus(ipair(H,H),div(times(x,x),pt(q))) for x,q in zip(hp156,norms156)]
    assert all(x[1]<0 for x in facet156)
    cusp156=[[ipair(n,v) for v in cusp_vectors+[nu2]] for n in newpair]
    assert all(x[1]<0 for row in cusp156 for x in row)
    triangle156=[]
    for record in alltrianglevertices:
        v=[tuple(map(F,x)) for x in record['point']]
        row=[ipair(v,n) for n in newpair];assert all(x[1]<0 for x in row)
        triangle156.append(vec(row))
    escape156=[ipair(E6,n) for n in full156];assert all(x[1]<0 for x in escape156)
    paircutoff=div(sstar,cstar);assert 0<paircutoff[0]<=paircutoff[1]<1
    overlap(paircutoff,pairdoc['definition']['axis_cutoff_bounds'])

    releasefile='results/fresh_h20_release144_cap_20261016.json'
    release=json.loads((ROOT/releasefile).read_text())
    assert release['status']=='PASS_EXACT_157_WALL_RELEASE144_E6_CUTTER'
    integrity([releasefile,'results/fresh_h20_dominant_axis_cluster_20261016_review.json'])
    frame157=[95,142,145,151]
    ga=ipair(full156[95],full156[145]);assert ga[1]<-2
    G157=[[pt(norms156[a]) if a==b else ga if {a,b}=={95,145} else pt(0)
           for b in frame157] for a in frame157]
    inverse157,pivot157=interval_inverse(G157)
    target142=neg(root(pt(norms156[142]/2)))
    c157=[pt(0),target142,pt(-1),pt(0)]
    particular157=ilinear([full156[j] for j in frame157],matvec(inverse157,c157))
    # Use the common interior projected off S16 as a different kernel seed.
    # The four-frame Gram has signature(3,1), so its one-dimensional
    # complement inside S16-perpendicular is positive definite.
    Hres=ivadd(H,ivscale(ilinear([dual[S.index(j)] for j in sorted(S16)],
        [ipair(H,full156[j]) for j in sorted(S16)]),pt(-1)))
    hs157=[ipair(Hres,full156[j]) for j in frame157]
    k157=ivadd(Hres,ivscale(ilinear([full156[j] for j in frame157],matvec(inverse157,hs157)),pt(-1)))
    qk157=ipair(k157,k157);assert qk157[0]>0
    pp157=inner(c157,matvec(inverse157,c157))
    delta157=div(minus(pt(2),pp157),qk157);assert delta157[0]>0
    branches157=[];accepted157=[]
    for sign in (-1,1):
        r=ivadd(particular157,ivscale(k157,times(pt(sign),root(delta157))))
        hp=ipair(r,H);ec=ipair(r,E6)
        vals=[ipair(r,n) for n in full156]
        bad=[j for j,x in enumerate(vals) if j not in S16|set(frame157) and not(x[1]<0 and x[1]**2>2*norms156[j])]
        branches157.append(dict(sign=sign,center_pair=saved(hp),E6_pair=saved(ec),failed_dotted_gates=bad))
        if not bad and hp[1]<0<ec[0]:accepted157.append((sign,r,vals))
    assert len(accepted157)==1
    sign157,new156,vals157=accepted157[0]
    assert R.overlap(ipair(new156,new156),pt(2))
    for x,y in zip(new156,release['definition']['physical_normal_bounds']):overlap(x,y)
    contacts157=[]
    for j,value in enumerate(vals157):
        if j in S16|{95,151}:
            rel='right';assert R.overlap(value,pt(0))
        elif j in (142,145):
            rel='pi/3';assert R.overlap(value,target142 if j==142 else pt(-1))
        else:rel='dotted'
        raw=release['new_pairs'][j];assert raw['wall']==j and raw['relation']==rel;overlap(value,raw['pair_bounds'])
        contacts157.append(dict(wall=j,relation=rel,bounds=saved(value)))
    counts157=Counter(x['relation'] for x in contacts157)
    assert counts157==dict(right=18,**{'pi/3':2,'dotted':136})
    full157=full156+[new156];norms157=norms156+[F(2)]
    hp157=[ipair(H,n) for n in full157];assert all(x[1]<0 for x in hp157)
    facets157=[minus(ipair(H,H),div(times(x,x),pt(q))) for x,q in zip(hp157,norms157)]
    assert all(x[1]<0 for x in facets157)
    cusps157=[ipair(new156,v) for v in cusp_vectors+[nu2]];assert all(x[1]<0 for x in cusps157)
    triangle157=[]
    for rec in alltrianglevertices:
        v=[tuple(map(F,x)) for x in rec['point']];val=ipair(v,new156);assert val[1]<0;triangle157.append(saved(val))
    corners157=[]
    for side in (95,145):
        ids=[151,142,side,156]
        G4=[[pt(norms157[a]) if a==b else target142 if {a,b}=={142,156}
             else pt(-1) if {a,b}=={145,156} else pt(0) for b in ids] for a in ids]
        schur=[row[:] for row in G4];pivots=[]
        for k in range(4):
            q=schur[k][k];assert q[0]>0;pivots.append(saved(q))
            for i in range(k+1,4):
                for j in range(i,4):schur[i][j]=schur[j][i]=minus(schur[i][j],div(times(schur[i][k],schur[k][j]),q))
        inv,_=interval_inverse(G4)
        hs=[ipair(H,full157[j]) for j in ids]
        v=ivadd(Hres,ivscale(ilinear([full157[j] for j in ids],matvec(inv,hs)),pt(-1)))
        vn,vh=ipair(v,v),ipair(v,H);assert vn[1]<0 and vh[1]<0
        support=sorted(S16|set(ids));assert len(support)==20
        vp=[ipair(v,n) for n in full157]
        for j,value in enumerate(vp):
            if j in support:assert R.overlap(value,pt(0))
            else:assert value[1]<0,(side,j,saved(value))
        raw=next(x for x in release['new_finite_corners'] if x['old_side']==side)
        assert support==raw['incident']
        for x,y in zip(v,raw['point_bounds']):overlap(x,y)
        corners157.append(dict(incident=support,old_side=side,point=vec(v),norm=saved(vn),center_pair=saved(vh),
            all157_pairs=vec(vp),positive_four_gram_pivots=pivots))
    neighbors157=sorted(S16|set(frame157));assert len(neighbors157)==20
    rejected157=[]
    for omit in neighbors157:
        if omit in (95,145):continue
        assert {95,145}<=set(neighbors157)-{omit}
        rejected157.append(dict(omitted=omit,dotted_pair=[95,145]))
    assert len(rejected157)==18
    z157=list(map(F,['7/10','2/5','7/10','2/5']));eps157=F(1,10**10)
    assert release['strict_escape']['z']==list(map(str,z157)) and release['strict_escape']['epsilon']==str(eps157)
    E7=ivadd(cv,ilinear([full157[j] for j in (141,142,143,144)],
        [neg(times(pt(z),root(pt((-cornerq)/norms157[j])))) for j,z in zip((141,142,143,144),z157)]))
    E7=ivadd(E7,ilinear([dual[S.index(j)] for j in sorted(S16)],[pt(-eps157)]*16))
    q7=(-cornerq)*(sum(x*x for x in z157)-1)+eps157**2*sum(J[S.index(i)][S.index(j)] for i in S16 for j in S16)
    assert q7>0 and str(q7)==release['strict_escape']['norm']
    assert R.overlap(ipair(E7,E7),pt(q7)) and ipair(E7,H)[1]<0
    escape157=[ipair(E7,n) for n in full157];assert all(x[1]<0 for x in escape157)
    for x,y in zip(E7,release['strict_escape']['physical_point_bounds']):overlap(x,y)
    paths=[f'scripts/{TAG}.py',releasefile,pairfile,axisfile,trianglefile,capfile,
        'scripts/fresh_h20_dominant_axis_cluster_20261016_review.py',
        'results/fresh_h20_dominant_axis_cluster_20261016_review.json',
        'scripts/fresh_h20_dominant142_cap_20261016_review.py',
        'results/fresh_h20_dominant142_cap_20261016_review.json',
        'scripts/fresh_h20_pi8_triangle_face_20261016_review.py',
        'results/fresh_h20_pi8_triangle_face_20261016_review.json',
        'scripts/fresh_h20_second_cusp_pi8_cap_20261016_review.py',
        'results/fresh_h20_second_cusp_pi8_cap_20261016_review.json',
        completionfile,partnerfile,clusterfile,mixedfile,newfile,
        'scripts/fresh_h20_bridge_edge_closure_20261014_review.py',
        'scripts/fresh_h20_bridge_faces_20261014_review.py',
        'scripts/fresh_h20_s61_multicap_20261014_exact_review.py']
    return dict(status='PASS_INDEPENDENT_157_WALL_RELEASE144_E6_CUTTER',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        kernel_seed='Common future H projected off S16, then off the four contact roots',
        frame=frame157,inverse_pivots=pivot157,kernel_square=saved(qk157),particular_square=saved(pp157),
        positive_radical_square=saved(delta157),branch_tests=branches157,selected_sign=sign157,
        normal=vec(new156),norm_exact=2,all156_contacts=contacts157,relations=dict(counts157),
        actual_facets=157,all157_center_pairs=vec(hp157),facet_projection_squares=vec(facets157),
        five_proper_cusp_pairs=vec(cusps157),compact_triangle_new_pairs=triangle157,
        E6_cut=saved(ipair(new156,E6)),complete_new_finite_vertices=2,known_finite_vertices_at_least=182,
        new_finite_corners=corners157,finite_neighbors=neighbors157,rejected_supports=rejected157,
        strict_escape=dict(norm=str(q7),point=vec(E7),all157_pairs=vec(escape157),center_pair=saved(ipair(E7,H))),
        proof='Independent156geometry. A generic four-frame inverse and different kernel seed construct the unique accepted norm2 radical root. Exactprojection identities give18rights and2pi3contacts; strict intervals certify136dots, E6cut, commoncenter and5cusps. Two independentfour-root corner projections and18dotted rejections exhaustnewfinitevertices. Exactrational norm formula and157strict signs verify renewedE7.',
        scope='157actualfacets with at least182finitevertices, still infinitevolume. No new construction dimension or complete ambient finite/exterior census.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
