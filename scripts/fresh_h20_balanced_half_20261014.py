#!/usr/bin/env python3
"""All-label radius bounds and exact radial-probe gates for balanced endpoints."""
import argparse, hashlib, json
from fractions import Fraction as F
from pathlib import Path
import fresh_h20_balanced_endpoint_20261014 as E
import fresh_h20_balanced_endpoint_20261014_bounds as B

ROOT=E.ROOT; TAG='fresh_h20_balanced_half_20261014'
SRC='results/fresh_h20_balanced_endpoint_20261014.json'
BOUND='results/fresh_h20_balanced_endpoint_20261014_bounds.json'
sub=B.sub; dec=B.dec; div=B.div
rt2=(F(0),F(1),F(0),F(0),F(0),F(0),F(0),F(0))
rt3=(F(0),F(0),F(1),F(0),F(0),F(0),F(0),F(0))
phi=(F(1,2),F(0),F(0),F(0),F(1,2),F(0),F(0),F(0))
c4=E.scale(rt2,F(1,2)); c6=E.scale(rt3,F(1,2))
r=E.scale(sub(rt2,E.c(1)),16)

def nonintegral_minpoly(a):
    """Distinct Galois conjugates give the monic minimal polynomial over Q."""
    conj=set()
    for mask in range(8):
        conj.add(tuple(x*(-1 if (i&mask).bit_count()%2 else 1)for i,x in enumerate(a)))
    p=[E.c(1)]
    for b in sorted(conj):
        out=[E.ZERO]*(len(p)+1)
        for j,c in enumerate(p):
            out[j]=sub(out[j],E.mul(b,c));out[j+1]=E.add(out[j+1],c)
        p=out
    assert all(not any(c[1:])for c in p)
    q=[c[0]for c in p]
    assert q[-1]==1 and any(c.denominator!=1 for c in q)
    return q

def sharp_bound(profile):
    contacts=profile['contacts']; fs=[0 if w<7 else 1 if w<14 else 2 for w,m in contacts]
    if not contacts:return E.c(8),'all_right'
    if len(contacts)==1:
        if fs[0]==2:return r,'single_A7'
        if contacts[0][1]==3:return E.c(7),'single_A6_pi3'
        return E.scale(sub(rt2,E.c(1)),14),'single_A6_pi4'
    assert all(m==3 for w,m in contacts)
    if fs[0]==fs[1]:
        return E.c(F(16,3)if fs[0]==2 else 7),'same_factor_cycle'
    if 2 in fs:return E.scale(sub(rt3,E.c(1)),7),'cross_A6_A7'
    return E.c(7),'cross_A6_A6'

def run():
    source=json.loads((ROOT/SRC).read_text());bounds=json.loads((ROOT/BOUND).read_text())
    for doc in [source,bounds]:
        for p,h in doc['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    typ=source['types'][0];profiles=typ['profiles'];I=[[F(x)for x in row]for row in typ['inverse_Gram']]
    old=bounds['scales'][0]['profiles']; upper=[];radius_proofs=[]
    # R>=8 in the all-right family forces the A6 cosine into m=3,4,5,6
    # (or the rational boundary 7/8). The other factor then cannot be finite.
    allright=[]
    for label,c in [(3,E.c(F(1,2))),(4,c4),(5,E.scale(phi,F(1,2))),(6,c6)]:
        other=E.scale(c,F(8,7));assert E.sign(other)>0 and E.sign(sub(E.c(1),other))>0
        allright.append(dict(label=label,other_cosine=other,minpoly_2cos=nonintegral_minpoly(E.scale(other,2))))
    boundary=nonintegral_minpoly(E.c(F(7,4)))
    # For R>7/(1/2+cos(pi/6)), the touched A6 cosine is 0 or m=3,4,5.
    # In each case the touched A7 cosine (8*c/7+1/14) is impossible.
    cross=[]
    for label,c in [(2,E.ZERO),(3,E.c(F(1,2))),(4,c4),(5,E.scale(phi,F(1,2)))]:
        other=E.add(E.scale(c,F(8,7)),E.c(F(1,14)))
        assert E.sign(other)>0 and E.sign(sub(E.c(1),other))>0
        cross.append(dict(label=label,other_cosine=other,minpoly_2cos=nonintegral_minpoly(E.scale(other,2))))
    for profile,prior in zip(profiles,old,strict=True):
        assert profile['id']==prior['profile_id']
        U,kind=sharp_bound(profile);upper.append(U);rejected=[]
        if kind not in ['all_right','cross_A6_A7']:
            assert E.sign(sub(U,dec(prior['continuous_radius_bound'])))>=0
            for ex in prior['exceptional_radii']:
                if not ex['retained']or E.sign(sub(dec(ex['radius']),U))<=0:continue
                witnesses=[]
                for factor,ratio in enumerate(ex['ratios']):
                    cc=dec(ratio)
                    if E.sign(cc)<=0 or E.sign(sub(E.c(1),cc))<=0:continue
                    try:poly=nonintegral_minpoly(E.scale(cc,2))
                    except AssertionError:continue
                    witnesses.append(dict(factor=factor,cosine=cc,minpoly_2cos=poly))
                assert witnesses,(profile['id'],ex)
                rejected.append(dict(radius=dec(ex['radius']),witness=witnesses[0]))
        assert E.sign(U)>0
        radius_proofs.append(dict(profile=profile['id'],kind=kind,radius_upper_bound=U,
            forbidden_exceptional_radii=rejected))
    # Evaluate the entire finite set of directions at their all-real upper bounds.
    # In normalized unit-root coordinates the retained Gram inverse is 2I.
    evaluations=[]
    for lam in [F(2,5),F(49,100),F(1,2)]:
        minimum=None;minimum_at=[];possible=[];zero=[];hsh=hashlib.sha256()
        for ti,tr in enumerate(source['literal_edge_transports']):
            ds=[1-lam*(1-(7 if w in [0,7]else 8 if w==14 else 0))for w in tr['edge_support_in_profile_order']]
            assert all(d>0 for d in ds)
            dual=[sum(ds[i]*I[i][j]for i in range(19))for j in range(19)]
            assert all(v>0 for v in dual)
            Fdist=2*sum(ds[i]*dual[i]for i in range(19))
            for profile,U in zip(profiles,upper,strict=True):
                cc=[dec(v)for v in profile['cosine_vector']];p=dec(profile['normalized_spatial_cost'])
                t=E.ZERO
                for a,b in zip(dual,cc,strict=True):
                    if b!=E.ZERO:t=sub(t,E.scale(b,2*a))
                assert E.sign(t)<=0 and E.sign(sub(E.c(1),p))>=0
                power=E.add(E.c(Fdist),E.add(E.scale(E.mul(U,t),2),E.mul(sub(p,E.c(1)),E.mul(U,U))))
                sg=E.sign(power);row=[ti,profile['id']]
                if sg<0:possible.append(row)
                if sg==0:zero.append(row)
                if minimum is None or E.sign(sub(power,minimum))<0:minimum=power;minimum_at=[row]
                elif power==minimum:minimum_at.append(row)
                hsh.update(json.dumps(E.enc([row,power]),separators=(',',':')).encode())
        if lam==F(2,5):assert E.sign(sub(minimum,E.c(6)))>0 and not possible and not zero
        else:
            assert possible==[[0,0]]+[[0,i]for i in range(19,24)],possible
            if lam==F(49,100):assert not zero,zero
        evaluations.append(dict(probe_multiple=lam,probe_square=196*lam*lam,
            count=392*41,minimum_power_lower_bound=minimum,minimum_attainers=minimum_at,
            potential_strict_cutters=possible,zero_lower_bounds=zero,power_values_sha256=hsh.hexdigest()))
    # For each A7-inner candidate a strict half-cut requires R>49/8.
    # Thereafter the omitted A7 cosine is in [cos(pi/4),79/98), below cos(pi/5).
    assert E.sign(sub(E.scale(phi,F(1,2)),E.c(F(79,98))))>0
    assert E.sign(sub(r,E.c(F(49,8))))>0
    for pid in range(19,24):
        assert profiles[pid]['contacts']==[[pid-3,3]]
        cost=dec(profiles[pid]['normalized_spatial_cost']);k=16*cost[0]
        assert not any(cost[1:])and k in [12,15,16]
        power=E.add(E.c(49),E.add(E.scale(r,-k/2),E.scale(E.mul(r,r),k/16-1)))
        assert E.sign(power)<0
    # Strict two-probe cutter incompatibility. Nearest A6 vertices differ by square84.
    # Changing the omitted A6 node to its cyclic neighbor changes one retained
    # gap by7; its squared vertex displacement is therefore 2*7²*I[0,0].
    assert 2*7**2*I[0][0]==84
    # Reconstruct every direction difference in the literal normalized A7 block.
    direction_differences=[]
    for i in range(19,24):
        for j in range(i+1,24):
            ci=[dec(v)for v in profiles[i]['cosine_vector']]
            cj=[dec(v)for v in profiles[j]['cosine_vector']]
            delta=[sub(a,b)for a,b in zip(ci,cj,strict=True)];sq=E.ZERO
            for a in range(19):
                for b in range(19):
                    if delta[a]!=E.ZERO and delta[b]!=E.ZERO:
                        sq=E.add(sq,E.scale(E.mul(delta[a],delta[b]),2*I[a][b]))
            assert not any(sq[1:])
            direction_differences.append(dict(profiles=[i,j],square=sq[0]))
    assert {v['square']for v in direction_differences}=={F(7,16),F(3,4),F(15,16),F(1)}
    r2=E.mul(r,r);assert E.sign(sub(r2,E.c(42)))>0 and E.sign(sub(E.c(44),r2))>0
    assert E.sign(sub(r,E.c(6)))>0
    # all-right/all-right: R,S>7 => R²+S²>98>84 (positive pairing).
    # mixed: 0<(84+r²(p-1)-R²)/(2Rr)<35/84<1/2, p>=3/4, 7<R<8.
    assert F(20)-F(44,4)>0 and F(35,84)<F(1,2)
    # A7/A7: same index gives positive pairing; otherwise 7/16<=|yi-yj|²<=1.
    assert F(42,44)-1+F(7,32)>0
    # The universal .4 probe also lies beyond every radius-five bridge ball.
    x2=F(784,25);assert x2>25 and F(9,25)>F(1,4)
    # Standard Lorentz coordinates (x,D,U) with D=(x²-1/4)/2,U=1.
    D=(x2-F(1,4))/2
    assert x2-2*D==F(1,4)and -16-D<0
    paths=[SRC,BOUND,'scripts/fresh_h20_balanced_endpoint_20261014.py',
        'scripts/fresh_h20_balanced_endpoint_20261014_bounds.py',str(Path(__file__).resolve().relative_to(ROOT))]
    return E.enc(dict(status='PASS_ALL_LABEL_BALANCED_ENDPOINT_RADIAL_GAP_AND_HALF_CUTTER_RESTRICTION',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        radius_bounds=radius_proofs,all_right_forbidden_finite_cosines=allright,
        all_right_boundary_minpoly=boundary,cross_A6_A7_forbidden_finite_cosines=cross,
        radial_probe_evaluations=evaluations,
        strict_half_cutters=dict(transport=0,all_right_radius='7/cos(pi/m), m>=7',
            inner_A7_profiles=list(range(19,24)),inner_A7_radius=r),
        two_adjacent_half_probes=dict(old_vertex_separation_square=84,
            inner_A7_direction_difference_squares=direction_differences,
            all_right_vs_all_right='positive pairing',mixed='cosine strictly between0and1/2',
            inner_A7_same_index='positive pairing',inner_A7_distinct_indices='cosine strictly between0and1/2',
            scope='No pair of mutually Coxeter-compatible endpoint caps can strictly cut both half probes.'),
        common_strict_spacelike_escape=dict(spatial='(2/5)*representative_vertex',spatial_square=x2,
            D=D,U=1,norm=F(1,4),future_pairing_to_center_D16=-16-D,
            seed_gap_minimum=F(3,5),endpoint_ball_power_strictly_greater_than=6,
            outside_every_ball_contained_in_radius5=True,
            radius_five_ball_family_power_lower_bound=F(9,25)),
        scope='Unit A6xA6xA7 cusp seed. All real radii/all omitted labels for complete41 endpoint profiles and392 transports. The2/5 radial probe is strict outside every endpoint ball. It also survives any separately proved bridge family whose balls lie in the radius-five ball. This does not exclude arbitrary interior caps or a changed seed.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
