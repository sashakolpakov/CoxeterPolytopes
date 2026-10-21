#!/usr/bin/env python3
"""Independent constraint assembly and dual inequality replay.

Uses independently reconstructed physical coefficients and a different
interval engine; no optimization library or producer arithmetic imports.
"""
import argparse
import json
import math
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sympy as sy
import fresh_h20_s61_multicap_20261014_exact_review as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_null_link_escape_gate_20261016_review'
GEO='results/fresh_h20_null_link_escape_gate_20261016_geometry_review.json'
PRO='results/fresh_h20_null_link_profiles_20261016.json'
PRO_REVIEW='results/fresh_h20_null_link_profiles_20261016_review.json'
PRIMARY='results/fresh_h20_null_link_escape_gate_20261016.json'
pt,plus,minus,neg,times,div,total=R.pt,R.plus,R.minus,R.neg,R.times,R.div,R.total


def dec(raw):
    return tuple(map(F,raw))


def midpoint(a):
    return float((a[0]+a[1])/2)


def upper_sqrt(q):
    assert q>=0
    if not q:
        return F(0)
    den=1<<256
    a=math.isqrt((q*den*den).__floor__())
    if F(a*a,den*den)<q:
        a+=1
    assert F(a*a,den*den)>=q
    return F(a,den)


def norm_upper(values):
    return upper_sqrt(sum(max(a*a,b*b) for a,b in values))


def run():
    geo,pro,review,primary=[json.loads((ROOT/p).read_text()) for p in (GEO,PRO,PRO_REVIEW,PRIMARY)]
    assert geo['status']=='PASS_INDEPENDENT_NULL_LINK_DUAL_COEFFICIENTS'
    assert review['status']=='PASS_INDEPENDENT_COMPLETE_125_PLUS_TWO_INCIDENT_DOMAIN'
    assert review['sources'][PRO]==sha256((ROOT/PRO).read_bytes()).hexdigest()
    assert primary['status']=='PASS_EXACT_ALL_THROUGH_NULL_COXETER_ESCAPE_GATE'
    assert primary['unresolved']==[] and primary['extra_center_constraint'] is False
    assert primary['profile_domains']==127 and primary['cases']==2032
    white=geo['white']
    inverse=sy.Matrix([[sy.Rational(x) for x in row] for row in geo['positive_Gram']]).inv()
    J=[[F(x) for x in row] for row in inverse.tolist()]
    sq={d:R.root(pt(d)) for d in (2,3,5,6,10,15,30)}
    sq[1]=pt(1)
    magnitude={3:pt(1),4:sq[2],5:div(plus(pt(1),sq[5]),pt(2)),6:sq[3]}
    fixed=[dict(free=list(map(dec,row['free'])),null=dec(row['null']),white=list(map(dec,row['white'])))
           for row in geo['normals']+[geo['escape']]]
    incident=set(geo['incident'])
    nonincident=[i for i in range(135) if i not in incident]
    assert len(nonincident)==119
    assert all(fixed[i]['null'][1]<0 for i in nonincident)
    profiles=pro['profiles']+pro['continuous_tails']
    records=[]
    for number,profile in enumerate(profiles):
        tail='wall' in profile
        coeff=[pt(0)]*13
        if tail:
            radius_upper=sq[2][1]
            tail_index=white.index(profile['wall'])
        else:
            magnitudes=[pt(0)]*13
            for wall,label in zip(profile['walls'],profile['labels']):
                magnitudes[white.index(wall)]=magnitude[label]
            coeff=[neg(total(times(pt(J[j][i]),magnitudes[i]) for i in range(13))) for j in range(13)]
            if profile['cost']==['2']+['0']*7:
                radius_upper=F(0)
            else:
                cost=total(times(pt(F(q)),sq[d]) for q,d in zip(profile['cost'],pro['squarefree_basis']))
                assert cost[1]<2
                radius_upper=upper_sqrt(2-cost[0])
        def form(row):
            a=row['free'][:]
            if tail:
                a.append(neg(div(row['white'][tail_index],sq[2])))
            a.append(row['null'])
            constant=total(times(c,w) for c,w in zip(coeff,row['white']))
            return a,constant
        forms=[form(row) for row in fixed]
        c,e=forms[-1]
        for mask in range(16):
            saved=primary['records'][16*number+mask]
            assert saved['profile']==number and saved['mask']==mask
            witness=saved['certificate']
            if number==mask==0:
                assert witness['kind']=='zero_line_exclusion'
                zero=geo['zero_line']
                assert dec(zero['pair95'])[1]<0<dec(zero['pair107'])[0]
                assert dec(zero['r134_third_frame_coordinate'])[0]>0
                a,b=dec(zero['frame0_second_covector'])
                assert a>0 or b<0
                records.append(dict(profile=number,mask=mask,kind='zero_line_exclusion'))
                continue
            rows=[forms[i][0] for i in nonincident]
            bounds=[neg(forms[i][1]) for i in nonincident]
            labels=[['acute',i] for i in nonincident]
            for bit,wall in enumerate((2,5,61,87)):
                a,p=forms[wall]
                if mask & (1<<bit):
                    rows.append(a);bounds.append(minus(pt(-1),p));labels.append(['nonright',wall])
                else:
                    rows.append(list(map(neg,a)));bounds.append(p);labels.append(['right_reverse',wall])
            if tail:
                a=[pt(0)]*len(c);a[-2]=pt(1)
                rows.extend([a,list(map(neg,a))])
                bounds.extend([sq[2],neg(R.root(pt(F(3,2))))])
                labels.extend([['tail_upper'],['tail_lower']])
            # Reconstruct the same harmless rational row scaling. These
            # choices need not approximate an exact geometric coefficient.
            scales=[]
            for a,b in zip(rows,bounds):
                largest=max([1.]+[abs(midpoint(x)) for x in a]+[abs(midpoint(b))])
                scales.append(F(1,2**max(0,math.ceil(math.log2(largest)))))
            rows=[[times(x,pt(s)) for x in a] for a,s in zip(rows,scales)]
            bounds=[times(b,pt(s)) for b,s in zip(bounds,scales)]
            kind=witness['kind'];assert kind in ('upper_bound','infeasible')
            pivot=witness['height_balance_pivot'];assert 0<=pivot<len(rows)
            multipliers={}
            for index,raw in witness['rational_multipliers']:
                assert index not in multipliers and index!=pivot and 0<=index<len(rows)
                q=F(raw);assert q>=0
                multipliers[index]=pt(q)
            target=c[-1] if kind=='upper_bound' else pt(0)
            numerator=minus(target,total(times(mu,rows[i][-1]) for i,mu in multipliers.items()))
            multipliers[pivot]=div(numerator,rows[pivot][-1])
            assert multipliers[pivot][0]>=0,(number,mask,'negative multiplier')
            # Its exact defining formula enforces height cancellation.
            residual=[minus(c[j] if kind=='upper_bound' else pt(0),
                            total(times(mu,rows[i][j]) for i,mu in multipliers.items()))
                      for j in range(len(c)-1)]
            constant=total(times(mu,bounds[i]) for i,mu in multipliers.items())
            bound=constant[1]+radius_upper*norm_upper(residual)
            if kind=='upper_bound':
                bound+=e[1]
            assert bound<0,(number,mask,kind,str(bound))
            used=[labels[i] for i,mu in sorted(multipliers.items()) if mu[1]>0]
            assert used==witness['verified']['used_rows']
            records.append(dict(profile=number,mask=mask,kind=kind,strict_upper_bound=str(bound),
                                height_pivot=pivot,pivot_multiplier_bounds=list(map(str,multipliers[pivot]))))
        print('independent dual profile',number,flush=True)
    counts=Counter(r['kind'] for r in records)
    assert counts==dict(infeasible=1739,upper_bound=292,zero_line_exclusion=1)
    maxima={kind:max(F(r['strict_upper_bound']) for r in records if r['kind']==kind)
            for kind in ('infeasible','upper_bound')}
    paths=[f'scripts/{TAG}.py',GEO,PRO,PRO_REVIEW,PRIMARY,
           'scripts/fresh_h20_s61_multicap_20261014_exact_review.py']
    return dict(status='PASS_INDEPENDENT_ALL_THROUGH_NULL_COXETER_ESCAPE_GATE',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                cases=len(records),counts=dict(counts),maximum_bounds={k:str(v) for k,v in maxima.items()},
                records=records,extra_center_constraint=False,
                theorem='Every norm2 wall acute to the retained135 normals, through the specified future nu and with Coxeter contacts, has strictly negative E1 pairing. All127 contact domains and16 minimum-contact branches are covered.',
                scope='No such through-nu wall cuts E1; does not exclude walls off nu, changed retained walls, or other H20 constructions. Independent physical geometry, constraint assembly and dual upper bounds; no optimizer calls or producer arithmetic imports.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
