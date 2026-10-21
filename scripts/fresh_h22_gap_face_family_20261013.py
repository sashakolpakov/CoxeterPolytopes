#!/usr/bin/env python3
"""Exact complete cap-orbit family, deficient cusp and uniform strict escape."""
import argparse,json,hashlib
from collections import Counter
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT='results/fresh_h22_gap_face_family_20261013.json'
ARRAY='data/enumerations/fresh_h22_gap_face_orbit_20261013.npz'
PRIOR='results/fresh_h22_cusp_orbits_20261013_verify.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
SYMS='data/enumerations/fresh_h20_global_interior_sections_20261013.json'
S=[1,2,3,9,10,14,21,22,23,28,32,34,39,45,46,48,50,63,198]

def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def frac(x):return F(int(x.p),int(x.q))

def run(saved=None,rebuild=False):
    prior=read(PRIOR)
    assert all(sha(p)==h for p,h in prior['sources'].items())
    raw=read(RAW)['coxeter_input']
    Q=-np.array(raw['basis_gram'],dtype=np.int64)
    W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64)
    h=W.sum(axis=0);h//=np.gcd.reduce(h)
    assert h@Q@h==-14
    gens=[np.array(g['matrix'],dtype=np.int64)for g in read(SYMS)['ambient_symmetry_generators']]
    norm=np.diag(W@Q@W.T);assert set(norm)=={2,4}
    mat=sp.Matrix((W[S+[5,12,26]]@Q).tolist())
    assert mat.det()!=0
    p=mat.inv()*sp.Matrix([0]*19+[-1]*3)
    seed=np.array([int(6*x)for x in p],dtype=np.int64)
    assert all(F(int(x),6)==frac(y)for x,y in zip(seed,p))
    assert seed@Q@seed==-372 and h@Q@seed==-216
    if rebuild:
        rows=[seed];known={seed.tobytes()}
        for row in rows:
            for g in gens:
                v=g@row;k=v.tobytes()
                if k not in known:known.add(k);rows.append(v)
            assert len(rows)<=241920
        np.savez_compressed(ROOT/ARRAY,six_old_projections=np.array(rows,dtype=np.int64))
    R=np.load(ROOT/ARRAY)['six_old_projections']
    assert R.dtype==np.int64 and R.shape==(241920,22)and np.array_equal(R[0],seed)
    U=np.load(ROOT/'data/enumerations/fresh_h22_symmetry_orbit_20261013.npz')['twice_old_projections']
    olddata=read('data/enumerations/fresh_h22_cusp_orbits_20261013.json')
    assert all(sha(path)==digest for path,digest in olddata['sources'].items())
    oldcert=read('results/fresh_h22_symmetry_orbit_verify_20261013.json')
    assert all(sha(path)==digest for path,digest in oldcert['source_sha256'].items())
    data=np.load(ROOT/olddata['array']);D=data['old_null_projections'];four=data['four_heights']
    maximum=max(int(abs(a).max())for a in(Q,W,R,U,D,*gens))
    assert 22**2*maximum**3<2**63
    lookup={r.tobytes():i for i,r in enumerate(R)};assert len(lookup)==len(R)
    reached=np.zeros(len(R),dtype=bool);reached[0]=True
    for g,meta in zip(gens,read(SYMS)['ambient_symmetry_generators']):
        assert np.array_equal(g.T@Q@g,Q)and np.array_equal(g@h,h)
        assert np.array_equal(W@g.T,W[meta['wall_permutation']])
        perm=np.array([lookup[v.tobytes()]for v in R@g.T])
        assert len(np.unique(perm))==len(R)
        forward=perm>np.arange(len(R));reached[perm[forward]]=True
    assert np.all(reached)
    pp=R@(Q@seed);assert pp[0]==-372 and max(pp[1:])==-420
    assert all(x==-420 or x<=-612 for x in pp[1:])
    parallel=np.flatnonzero(pp==-420).tolist();assert len(parallel)==1
    pairs=W@(Q@seed);assert max(pairs)==0
    kinds=Counter()
    for v,n in zip(pairs,norm):
        product=F(int(v),6)
        if not product:kinds['right']+=1
        elif product==-1 and n==2:kinds['pi_over_k']+=1
        else:
            assert product<0 and product**2>F(2,3)*int(n)
            kinds['uniformly_dotted']+=1
    # For q_k=1/(2cos²(pi/k)), k>=6: 1/2<q_k<=2/3 and
    # H_k²=31/3+q_k<=11. Distinct same-sign orbit pairings are
    # <=-35/3+H_k²=q_k-4/3<=-q_k. At k6 the nearest pair is parallel.
    maxu=int(max(U@(Q@seed)));assert maxu==-96
    dp=D@(Q@seed);assert np.all(2*dp<=-15*four)
    assert max(F(int(v),6)*4/int(t)for v,t in zip(dp,four))==-5
    # Cross-family dotted bounds use sqrt(55/2)<6, sqrt11<10/3,
    # sqrt(q_k)<1 and sqrt(2q_k)<2. Orbit transfer checks every new cap.
    assert F(55,2)<36 and F(11)<F(100,9)
    assert -8+6<-1 and -5+F(10,3)<-1 and -36+4*F(10,3)<-2
    partner=R[parallel[0]]
    common=np.flatnonzero((pairs==0)&(W@(Q@partner)==0)).tolist()
    assert len(common)==18
    gram=sp.Matrix((W[common]@Q@W[common].T).tolist())
    assert gram.rank()==18
    # Complete incidence at delta=(p+p',2sqrt11): all old roots need
    # both summand pairings zero; all other cap families are strict.
    null_rows=sp.Matrix(np.vstack([W[common]@Q,seed@Q,partner@Q]).tolist())
    ns=null_rows.nullspace();assert len(ns)==2
    tangent=[frac(x)for x in ns[0]]
    tangent_np=np.array(tangent,dtype=object)
    tangent_norm=sum(tangent_np*(Q.astype(object)@tangent_np));assert tangent_norm>0
    delta=[F(int(a+b),6)for a,b in zip(seed,partner)]
    assert sum(np.array(delta,dtype=object)*(Q.astype(object)@np.array(delta,dtype=object)))==-44
    eps=F(1,1000);weight=F(1,100)
    if saved is None:
        x=[a+eps*b+weight*eps**2*int(c)for a,b,c in zip(delta,tangent,h)]
    else:
        x=list(map(F,saved['uniform_strict_escape']['old_coordinates']))
    den=lcm(*(v.denominator for v in x));num=np.array([int(v*den)for v in x],dtype=object)
    cov=Q.astype(object)@num;qp=F(int(num@cov),den**2)+44;hp=int(h.astype(object)@cov)
    assert qp>0 and hp<0 and hp*hp>704*den*den
    assert max(W.astype(object)@cov)<0
    bound=22*maximum*max(abs(int(v))for v in cov)
    assert bound+132*den<2**63
    cov=np.array(cov,dtype=np.int64)
    up=U@cov;assert max(up)<0 and int(max(up))**2>440*den*den
    dp=D@cov
    null_margins={}
    for t in(1,4):
        v=int(max(dp[four==t]));assert v<0 and 4*v*v>11*t*t*den*den
        null_margins[str(F(t,4))]=str(F(v,den))
    rp=R@cov;assert max(rp)+132*den<0
    # The escape works at maximal H_k=sqrt11 and maximal old center
    # height4, hence for all k>=6 and all old center labels m>=3.
    oldescape=prior['strict_escape'];oldnum=list(map(int,oldescape['numerator']))
    oldcov=Q.astype(object)@np.array(oldnum[:22],dtype=object)
    target_pair=int(seed.astype(object)@oldcov)
    assert target_pair<0 and 390*oldnum[-1]**2>target_pair**2
    result={'status':'PASS_COMPLETE_H22_CAP_ORBIT_FAMILY_AND_UNIFORM_INFINITE_VOLUME',
        'sources':{s:sha(s)for s in(PRIOR,RAW,SYMS,ARRAY)},
        'old_projection':list(map(str,map(frac,p))), 'six_old_projection':seed.tolist(),
        'old_projection_norm':'-31/3','old_center_pairing':'-36',
        'formula':'(g*p,+/-sqrt(31/3+q_k)), q_k=1/(2*cos(pi/k)^2), integers k>=6',
        'new_wall_norm':'q_k','orbit_directions':len(R),'new_signed_walls':2*len(R),
        'total_facets':prior['facets']+2*len(R),'old_pair_counts':dict(kinds),
        'max_twice_prior_timelike_pairing':str(F(maxu,6)),
        'max_prior_null_pairing_divided_by_height':'-5',
        'same_sign_nearest_pair_at_k6':'-2/3',
        'each_cap_parallel_partners_at_k6':1,
        'new_cusp':{'partner_index':parallel[0],'six_partner':partner.tolist(),
            'formula':'(p+p_prime,2*sqrt11)','common_old_walls':common,
            'old_rank':18,'parallel_rank':1,'horosphere_dimension':21,
            'quotient_rank':19,'unbounded_directions_at_least':2},
        'tangent':list(map(str,tangent)),'tangent_norm':str(tangent_norm),
        'uniform_strict_escape':{'old_coordinates':list(map(str,x)),
            'added_coordinate':'2*sqrt11','norm':str(qp),
            'old_center_pairing':str(F(hp,den)),
            'maximum_new_cap_pairing_at_k6':str(F(int(max(rp)),6*den)+22),
            'prior_null_maxima':null_margins},
        'previous_exact_escape_cut_for_all_k':True,
        'scope':'Exact complete241920-element symmetryorbit. All k>=6 and all old center labels m>=3 are compatible and infinite-volume. Deficient cusp statement is at k6. No obstruction to further caps or changed orbits is asserted.'}
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--rebuild-orbit',action='store_true');a=p.parse_args()
    assert not(a.verify and a.rebuild_orbit)
    saved=read(OUT)if a.verify else None;out=run(saved,a.rebuild_orbit)
    if a.verify:assert saved==out
    else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['total_facets'],out['uniform_strict_escape']['norm'])
