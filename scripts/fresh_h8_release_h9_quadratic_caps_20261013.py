#!/usr/bin/env python3
"""Add genuinely mixed-angle caps over Q(phi,sqrt(d)) to a released H8 seed.

For d=2,3,2+phi the new seed pairing sqrt(d) gives label 4,6,10.
The radical coefficient of p^t G^-1 p must vanish exactly. Only accepted
caps are exact; the numerical proposal filter does not certify exhaustion.
"""
import argparse
from itertools import combinations, product
import json
from pathlib import Path
import time
import numpy as np
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,sub,neg,scale,mul,div,sign,total,dot,mv,pair,matrix,
)

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/fresh_h8_release_h9_caps_expanded_20261013.json'
ESCAPE=ROOT/'data/enumerations/fresh_h8_release_h9_volume_expanded_20261013.json'
real=lambda x:float(x[0])+(1+5**.5)/2*float(x[1])

def esign(a,b,d):
    sa,sb=sign(a),sign(b)
    if not sb:return sa
    if not sa or sa==sb:return sb
    return sa*sign(sub(mul(a,a),mul(d,mul(b,b))))

def run(d,max_support=5):
    start=time.monotonic();raw=json.loads(SOURCE.read_text());I=matrix(raw['seed_gram_inverse']);R=matrix(raw['seed_normals']);Q=matrix(raw['host_gram'])
    y=tuple(map(parse,json.loads(ESCAPE.read_text())['y_coordinates']))
    In=np.array([[real(x)for x in row]for row in I])
    aa=(O,P,scale(O,2),add(O,P),scale(P,2),Z)
    bb=(Z,Z,Z,Z,Z,O)
    vals=np.array([real(x)+real(d)**.5*real(z)for x,z in zip(aa,bb)])
    records=[];counts=[]
    for k in range(1,max_support+1):
        choices=np.array([v for v in product(range(6),repeat=k)if 5 in v],dtype=int)
        p=vals[choices];proposed=0;exact=0;cutters=0
        for support in combinations(range(10),k):
            H=In[np.ix_(support,support)];squares=np.einsum('bi,ij,bj->b',p,H,p)
            for ci in np.flatnonzero(abs(squares-2)<1e-7):
                proposed+=1;a=[Z]*10;b=[Z]*10
                for j,m in zip(support,choices[ci]):a[j]=neg(aa[m]);b[j]=neg(bb[m])
                a=tuple(a);b=tuple(b);ca=mv(I,a);cb=mv(I,b)
                norm_a=add(dot(a,ca),mul(d,dot(b,cb)));norm_b=scale(dot(a,cb),2)
                if norm_a!=scale(O,2) or norm_b!=Z:continue
                ra=tuple(total(mul(c,w[j])for c,w in zip(ca,R))for j in range(10))
                rb=tuple(total(mul(c,w[j])for c,w in zip(cb,R))for j in range(10))
                assert add(pair(ra,Q,ra),mul(d,pair(rb,Q,rb)))==scale(O,2) and pair(ra,Q,rb)==Z
                assert tuple(pair(ra,Q,w)for w in R)==a and tuple(pair(rb,Q,w)for w in R)==b
                for x,z in zip(a,b):
                    assert esign(x,z,d)<=0
                    va=add(mul(x,x),mul(d,mul(z,z)));vb=scale(mul(x,z),2)
                    assert vb==Z and (va in (Z,O,scale(O,2),scale(O,3),add(O,P),add(scale(O,2),P)) or sign(sub(va,scale(O,4)))>=0)
                pairing_a=neg(dot(ca,y));pairing_b=neg(dot(cb,y));cuts=esign(pairing_a,pairing_b,d)>0
                exact+=1;cutters+=int(cuts)
                encode=lambda v:[enc(x)for x in v]
                records.append({'normal_rational':encode(ra),'normal_radical':encode(rb),
                    'seed_pairings_rational':encode(a),'seed_pairings_radical':encode(b),
                    'seed_coefficients_rational':encode(ca),'seed_coefficients_radical':encode(cb),
                    'old_escape_pairing_rational':enc(pairing_a),'old_escape_pairing_radical':enc(pairing_b),
                    'cuts_old_common_escape':cuts,'nonzero_pairings':k})
        row={'support_size':k,'mixed_patterns':len(choices)*len(list(combinations(range(10),k))),
             'numeric_norm_proposals':proposed,'exact_caps':exact,'old_escape_cutters':cutters}
        counts.append(row);print(row,flush=True)
    out={'source':str(SOURCE.relative_to(ROOT)),'old_escape_source':str(ESCAPE.relative_to(ROOT)),
         'extension_radicand':enc(d),'coordinate_convention':'rational + sqrt(d) * radical, each coefficient in Q(phi)',
         'maximum_pairing_support':max_support,'counts':counts,'caps':records,
         'old_escape_cutters':sum(c['cuts_old_common_escape']for c in records),
         'scope':'Exact accepted mixed-angle caps; bounded numerical proposal scan, not an exhaustion theorem.',
         'seconds':round(time.monotonic()-start,3)}
    name='2'if d==scale(O,2)else'3'if d==scale(O,3)else'2phi'
    path=ROOT/f'data/enumerations/fresh_h8_release_h9_quadratic_caps_{name}_20261013.json'
    path.write_text(json.dumps(out,indent=2)+'\n');print({'field':name,'caps':len(records),'cutters':out['old_escape_cutters'],'seconds':out['seconds']},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--radicand',choices=('2','3','2phi'),required=True);p.add_argument('--max-support',type=int,default=5);a=p.parse_args();run({'2':scale(O,2),'3':scale(O,3),'2phi':add(scale(O,2),P)}[a.radicand],a.max_support)
