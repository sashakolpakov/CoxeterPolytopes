#!/usr/bin/env python3
"""Complete necessary-contact test inside adjacent spacelike two-ray faces.

The original three cones are frozen inputs. Only relative interiors of
their adjacent spacelike-ray faces are covered. Passive contacts are
unrestricted; all candidates fail before further native/conjugate gates.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
import fresh_raw23_union_polar_completion_preflight_20261018 as P
H=P.H;Q=P.Q;ROOT=P.ROOT
TAG='fresh_raw23_two_ray_completion_preflight_20261018'
OLD=P.TAG
PINS={f'scripts/{OLD}.py':'801198f128ee0e8d3ff9587fbd65d656edaedc0d413e4b4665c094c290010ee8',
 f'results/{OLD}.json':'eee8c3836a16c6342b2882c3b65d84999bb7fbc7363413fbc59626c739ff8361',
 f'results/{OLD}.md':'42bd1f4baeecd9c8188a1590adf1e689c5f9582ae5df378a75d94dc5116da973',
 'scripts/check_raw23_same_family_saved_20261018.py':'5b1ab9873b2e069f18966a7052856d8a9df8a8c1b0eef4794522a812a54616a7'}
SQUARES=P.FINITE_SQUARES[1:]

def add(a,b):return[x+y for x,y in zip(a,b)]
def scale(a,q):return[x*q for x in a]

def bounds(h,contacts,scales):
    inv=H.inverse(h);out=[]
    for i,p in enumerate(contacts):
        if p==[Q(),Q()]:continue
        z=[sum((inv[u][v]*p[v]for v in range(2)),Q())for u in range(2)]
        boundary=[2*p[k]*p[k]/(h[k][k]*scales[i])for k in range(2)]
        interior=all(x.sign()>0 for x in z)
        if interior:upper=2*sum((p[k]*z[k]for k in range(2)),Q())/scales[i]
        else:upper=boundary[0]if(boundary[0]-boundary[1]).sign()>=0 else boundary[1]
        assert all((upper-x).sign()>=0 for x in boundary)
        finite=(upper-4).sign()<0 or(upper-4)==0 and not interior
        out.append(dict(wall=i,contact=p,stationary_coefficients=z,boundary_bounds=boundary,
            upper_square=upper,upper_attained_in_interior=interior,strictly_below_four_in_interior=finite))
    return out

def candidates(h,contacts,scales,anchors):
    i,j=anchors;u,v=contacts[i],contacts[j]
    det=u[0]*v[1]-u[1]*v[0];assert det!=0
    out=[]
    for d,e in product(SQUARES,repeat=2):
        rad1,rad2=d*scales[i],e*scales[j];f=H.Field([rad1,rad2])
        aa=scale(add(scale(f.sqrt(rad1),v[1]),scale(f.sqrt(rad2),-u[1])),Q(1)/det)
        bb=scale(add(scale(f.sqrt(rad1),-v[0]),scale(f.sqrt(rad2),u[0])),Q(1)/det)
        assert add(scale(aa,u[0]),scale(bb,u[1]))==f.sqrt(rad1)
        assert add(scale(aa,v[0]),scale(bb,v[1]))==f.sqrt(rad2)
        rec=dict(anchor_squares=[d,e],field_radicands=f.basis,coefficients=[aa,bb],decision=None)
        if f.sign(aa)<=0 or f.sign(bb)<=0:
            rec['decision']='outside_relative_interior';out.append(rec);continue
        norm=add(add(scale(f.mul(aa,aa),h[0][0]),scale(f.mul(bb,bb),h[1][1])),scale(f.mul(aa,bb),2*h[0][1]))
        rec['projection_norm']=norm;assert f.sign(norm)>0
        if f.sign(add(norm,f.scalar(-2)))>0:rec['decision']='physical_norm_exceeds_two'
        else:rec['decision']='requires_further_known_wall_gates'
        out.append(rec)
    assert len(out)==25
    return out


def run():
    for p,h in PINS.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    old=H.read(OLD);H.bind(old);cases=[];hist=Counter();tests=0;survivors=[]
    for ci,c in enumerate(old['cases']):
        cone=c['initial_cone'];metric=H.mat(c['metric']);g=H.mat(c['original_scaled_gram'])
        rows=H.mat(cone['normal_rows']);sc=[Q(*q)for q in c['original_squared_scales']]
        rays=cone['rays']
        outside=[r for r in rays if r['norm_sign']>0];decisions=[];faces=[]
        assert H.positive([[Q(x.a,-x.b)for x in row]for row in metric])
        for x,y in combinations(outside,2):
            common=sorted(set(x['zero_facets'])&set(y['zero_facets']))
            rank=len(P.rref([rows[k]for k in common])[1])if common else 0
            decisions.append(dict(rays=[x['index'],y['index']],common_walls=common,common_rank=rank,adjacent=rank==8))
            if rank!=8:continue
            a,b=([Q(*q)for q in r['coordinates']]for r in(x,y))
            h=[[P.inner(a,metric,a),P.inner(a,metric,b)],[P.inner(a,metric,b),P.inner(b,metric,b)]]
            assert H.positive(h)
            contacts=[[-Q(*r['products'][k])for r in(x,y)]for k in range(len(g))]
            assert all(q.sign()>=0 for p in contacts for q in p)
            bb=bounds(h,contacts,sc);eligible=[z['wall']for z in bb if z['strictly_below_four_in_interior']]
            anchors=next((i,j)for i,j in combinations(eligible,2)if contacts[i][0]*contacts[j][1]-contacts[i][1]*contacts[j][0]!=0)
            cc=candidates(h,contacts,sc,anchors)
            hist.update(z['decision']for z in cc);tests+=len(cc)
            for k,z in enumerate(cc):
                if z['decision']=='requires_further_known_wall_gates':survivors.append(dict(case=ci,rays=[x['index'],y['index']],candidate=k,record=z))
            faces.append(dict(rays=[x['index'],y['index']],plane_gram=h,contact_magnitudes=contacts,
                sharp_cone_bounds=bb,anchors=anchors,candidates=cc))
        assert len(decisions)==496 and len(faces)==112
        cases.append(dict(kinds=c['kinds'],adjacency_decisions=decisions,faces=faces))
        print('CASE',c['kinds'],'faces',len(faces),'retained',sum(z['case']==ci for z in survivors),flush=True)
    assert tests==8400 and not survivors
    assert hist==Counter(outside_relative_interior=4864,physical_norm_exceeds_two=3536)
    paths=set(PINS)|{f'scripts/{TAG}.py'}
    return P.enc(dict(status='PASS_RAW23_ADJACENT_SPACELIKE_TWO_RAY_PREFLIGHT',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},cases=cases,
        complete_pair_domain=1488,adjacent_spacelike_faces=336,native_anchor_candidates=tests,
        decision_counts=dict(hist),retained_candidates=survivors,
        relative_interiors_only=True,boundaries_inherited_from_frozen_polar=True,
        timelike_generator_faces_tested=False,higher_dimensional_faces_tested=False,
        arbitrary_passive_lifts_relaxed=True,unknown_passive_contacts_assumed_zero=False,
        candidates_reaching_other_native_or_conjugate_gates=0,
        geometric_realization_claimed=False,pair_exclusion_claimed=False,exported_columns=[],optimizer_invocations=0,
        scope='Complete necessary finite-anchor/physical-norm check on relative interiors of all336 adjacent-spacelike two-ray faces of the original three cones. Positive-cone Cauchy bounds force25 exact anchor pairs per face. All8400 candidates fail interior positivity or physical projection norm before any further native or conjugate gate. Passive components cannot reduce that norm. Other projected faces, unknown contacts and full completion remain open.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();d=run()
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['decision_counts'],'retained',len(d['retained_candidates']),flush=True)
