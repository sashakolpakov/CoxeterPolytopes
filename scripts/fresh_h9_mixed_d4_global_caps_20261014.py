#!/usr/bin/env python3
"""A cheap exact instance of the general strict-escape dotted-cap lemma.

This is search triage, not a finite-volume construction: the new wall is
dotted to every old wall, and its entire H8 hyperplane is an actual facet.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from fresh_h9_mixed_d4_coupling_20261014_review import E,Z,dec,dv,dm,enc,pair

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_mixed_d4_global_caps_20261014'
BASE='results/fresh_h9_mixed_d4_coupling_20261014.json'
POINT=['3.445009490103295','4.274481843159041','3.0011044891310377',
       '1.7539516778449502','0.5329849937376314','0.15617993073424302',
       '-1.048155530053127','-0.7473950044921694','-0.656929489224263',
       '-0.665829155523975']
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    old=json.loads((ROOT/BASE).read_text())
    for p,h in old['sources'].items():assert sha(p)==h,p
    Q=dm(old['form']);N=dm(old['final_normals']);H=dv(old['center'])
    assert len(N)==12 and len(Q)==10
    esc=[E(F(x))for x in POINT]
    qH=pair(H,H,Q);qE=pair(esc,esc,Q);EH=pair(esc,H,Q)
    oldq=[pair(n,n,Q)for n in N]
    olde=[pair(esc,n,Q)for n in N]
    oldh=[pair(H,n,Q)for n in N]
    assert qH.sign()<0 and qE.sign()>0 and EH.sign()<0
    assert all(x.sign()>0 for x in oldq)
    assert all(x.sign()<0 for x in olde+oldh)
    lo,hi=F(0),F(1);trace=[]
    for it in range(80):
        t=(lo+hi)/2
        v=[(1-t)*a+t*b for a,b in zip(H,esc)];q=pair(v,v,Q)
        trace.append(dict(parameter=str(t),norm_sign=q.sign()))
        if q.sign()<=0:lo=t
        else:
            hi=t;g=[pair(v,n,Q)for n in N];cut=pair(v,esc,Q);vh=pair(v,H,Q)
            margins=[x*x-q*y for x,y in zip(g,oldq)]
            if cut.sign()>0 and vh.sign()<0 and all(x.sign()<0 for x in g)and all(x.sign()>0 for x in margins):break
    else:raise AssertionError('No exact near-null parameter within bounded bisection')
    assert 0<t<1 and q.sign()>0
    # Positively normalizing v by sqrt(2/q) gives a norm-two Coxeter cap.
    # Its projected old walls are all future timelike in v-perp.
    p=[a-vh*b/q for a,b in zip(H,v)]
    assert not pair(p,v,Q) and pair(p,p,Q).sign()<0
    pn=[pair(p,n,Q)for n in N];assert all(x.sign()<0 for x in pn)
    projected_norms=[y-x*x/q for x,y in zip(g,oldq)]
    assert all(x.sign()<0 for x in projected_norms)
    # All old facets survive at their strict center projections.
    facet_signs=[]
    for n,qn,hn in zip(N,oldq,oldh):
        f=[a-hn*b/qn for a,b in zip(H,n)]
        val=pair(f,v,Q);assert val.sign()<0;facet_signs.append(val)
    paths=[BASE,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',f'scripts/{TAG}.py']
    return enc(dict(status='PASS_EXACT_STRICT_ESCAPE_DOTTED_CAP_LEMMA_EXAMPLE',
        sources={p:sha(p)for p in paths},old_wall_count=12,
        literal_escape=POINT,escape_norm=qE,escape_center_pair=EH,escape_old_pairs=olde,
        interpolation_parameter=str(t),bisection_trace=trace,
        raw_cap=v,raw_cap_norm=q,normalized_cap_recipe='sqrt(2/raw_cap_norm)*raw_cap',
        raw_cap_escape_pair=cut,raw_cap_center_pair=vh,raw_cap_old_pairs=g,
        strict_dotted_squared_margins=margins,new_facet_point=p,new_facet_norm=pair(p,p,Q),
        new_facet_old_pairs=pn,old_facet_new_pairings=facet_signs,
        old_projected_norms=projected_norms,
        consequences=dict(new_finite_contacts=0,new_facet='entire H8',
            finite_volume=False,additional_meeting_facets_necessary_at_least=9),
        scope='An exact example of a general escape-cutting dotted-wall lemma. It creates no finite corner and is not a proposed finite-volume completion. No global no-single-cap exclusion can hold for a strict spacelike escape.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],'t',d['interpolation_parameter'],'steps',len(d['bisection_trace']),flush=True)
