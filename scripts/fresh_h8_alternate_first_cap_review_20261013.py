#!/usr/bin/env python3
"""Independent exact twelve-wall replay of the first quartic coupled cap."""
from fractions import Fraction as F
import json
from pathlib import Path
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,sub,neg,scale,mul,div,sign,total,dot,mv,pair,matrix,inverse,rref)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fresh_h8_alternate_first_cap_independent_audit_20261013.json'

def ea(x,y):return add(x[0],y[0]),add(x[1],y[1])
def en(x):return neg(x[0]),neg(x[1])
def es(x,y):return ea(x,en(y))
def em(x,y,d):return add(mul(x[0],y[0]),mul(d,mul(x[1],y[1]))),add(mul(x[0],y[1]),mul(x[1],y[0]))
def ed(x,y,d):
    den=sub(mul(y[0],y[0]),mul(d,mul(y[1],y[1])))
    z=em(x,(y[0],neg(y[1])),d);return div(z[0],den),div(z[1],den)
def sgn(x,d):
    a,b=x;sa,sb=sign(a),sign(b)
    if not sb:return sa
    if not sa or sa==sb:return sb
    return sa*sign(sub(mul(a,a),mul(d,mul(b,b))))
def etotal(xs):
    result=(Z,Z)
    for x in xs:result=ea(result,x)
    return result
def epair(x,g,y,d):return etotal(em(em(x[i],(g[i][j],Z),d),y[j],d) for i in range(len(x)) for j in range(len(y)) if g[i][j]!=Z)
def pe(x):return parse(x['rational']),parse(x['radical'])
def ve(x):return tuple(map(parse,x))
def ee(x):return {'rational':enc(x[0]),'radical':enc(x[1])}

def run():
    base=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_3_6_20261013.json').read_text())
    polar=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_3_6_polar_20261013.json').read_text())
    raw=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_polar_coupled_20261013.json').read_text());cap=raw['caps'][0]
    Q=matrix(base['host_gram']);R=matrix(base['seed_normals']);G=matrix(base['seed_gram']);I=inverse(G)
    assert G==tuple(tuple(pair(r,Q,s) for s in R) for r in R) and matrix(base['seed_gram_inverse'])==I
    assert len(rref(R)[1])==10
    # Q is the original diagonal time direction plus positive E8 and scalar.
    assert Q[0][0]==scale(P,-2) and all(Q[0][i]==Q[i][0]==Z for i in range(1,10))
    positive=[list(row[1:]) for row in Q[1:]]
    for i in range(9):
        t=positive[i][i];assert sign(t)>0
        for j in range(i+1,9):
            for k in range(i+1,9):positive[j][k]=sub(positive[j][k],div(mul(positive[j][i],positive[i][k]),t))
    d=parse(cap['discriminant']);assert d==(F(49),F(-24)) and sign(d)>0
    normd=d[0]*d[0]+d[0]*d[1]-d[1]*d[1];assert normd==649
    x=pe(cap['free_x']);A,B,C=ve(cap['norm_polynomial_descending'])
    assert (A,B,C)==(scale(P,F(1,2)),add(scale(O,-7),scale(P,4)),add(scale(O,-16),scale(P,8)))
    assert es(ea(ea(em((A,Z),em(x,x,d),d),em((B,Z),x,d)),(C,Z)),(Z,Z))==(Z,Z)
    assert sgn(es(x,(scale(O,2),Z)),d)>0
    pa,pb=ve(cap['pairings_rational']),ve(cap['pairings_radical']);ca,cb=mv(I,pa),mv(I,pb)
    assert ca==ve(cap['coefficients_rational']) and cb==ve(cap['coefficients_radical'])
    assert tuple(zip(pa,pb))==tuple(en(x) if i==0 else (neg(add(scale(O,2),P)),Z) if i==8 else (Z,Z) for i in range(10))
    pro=polar['polar_caps'][0];pc=ve(pro['seed_coefficients']);pq=parse(pro['norm'])
    assert len(polar['polar_caps'])==1 and pro['seed_index']==5 and pro['original_wall']==5
    assert pc==tuple(neg(t) for t in I[5]) and pq==I[5][5]==(F(6,5),F(-2,5))
    # Work in the original independent wall basis. The Gram rank and inertia
    # of all twelve normals follow from the leading nonsingular Lorentz block.
    normals=[tuple((O if i==j else Z,Z) for i in range(10)) for j in range(10)]
    normals.extend((tuple((t,Z) for t in pc),tuple(zip(ca,cb))))
    norms=[epair(r,G,r,d) for r in normals]
    assert norms[-1]==(scale(O,2),Z) and norms[-2]==(pq,Z)
    paircp=epair(normals[-1],G,normals[-2],d);assert paircp==(neg(O),Z)
    permitted={Z:2,O:3,scale(O,2):4,scale(O,3):6,add(O,P):5,add(scale(O,2),P):10}
    pair_rows=[]
    for i,r in enumerate(normals):
        assert sgn(norms[i],d)>0
        for j,s in enumerate(normals[:i]):
            p=epair(r,G,s,d);assert sgn(p,d)<=0
            ratio=ed(em((scale(O,4),Z),em(p,p,d),d),em(norms[i],norms[j],d),d)
            if sgn(es(ratio,(scale(O,4),Z)),d)>=0:label='parallel' if ratio==(scale(O,4),Z) else 'dotted'
            else:
                assert ratio[1]==Z and ratio[0] in permitted;label=permitted[ratio[0]]
            pair_rows.append({'i':i,'j':j,'pairing':ee(p),'four_cos_squared':ee(ratio),'label':label})
    assert next(r['label'] for r in pair_rows if (r['i'],r['j'])==(11,10))==10
    y=ve(base['strict_center_y']);h=tuple((neg(t),Z) for t in mv(I,y))
    hnorm=epair(h,G,h,d);assert sgn(hnorm,d)<0
    products=[epair(h,G,r,d) for r in normals];assert all(sgn(v,d)<0 for v in products)
    facets=[]
    for i,r in enumerate(normals):
        foot=tuple(es(v,em(ed(products[i],norms[i],d),w,d)) for v,w in zip(h,r))
        fqs=epair(foot,G,foot,d);assert sgn(fqs,d)<0 and epair(foot,G,r,d)==(Z,Z)
        assert sgn(epair(foot,G,h,d),d)<0
        assert all(sgn(epair(foot,G,s,d),d)<0 for j,s in enumerate(normals) if j!=i)
        facets.append({'wall':i,'interior_point_coefficients':[ee(t) for t in foot],'norm':ee(fqs)})
    target=ve(raw['rounds'][0]['target_escape_y']);cut=(dot(ca,target),dot(cb,target))
    assert cut==pe(cap['target_escape_pairing']) and sgn(cut,d)<0
    out={'status':'PASS','cap_source':'data/enumerations/fresh_h8_alternate_polar_coupled_20261013.json',
        'cap_index':0,'coefficient_field':'Q(phi,sqrt(49-24phi))','extension_degree':4,
        'discriminant':enc(d),'discriminant_field_norm':str(normd),'free_x':ee(x),
        'old_nonzero_pairing_indices':[0,8],'old_nonzero_original_wall_indices':[0,9],
        'polar_seed_index':5,'wall_count':12,'gram_rank':10,'gram_inertia':[9,1,2],
        'cap_norm':ee(norms[-1]),'polar_norm':enc(pq),'cap_polar_pairing':ee(paircp),
        'cap_polar_label':10,'common_interior_norm':ee(hnorm),'pair_checks':pair_rows,
        'strict_facet_witnesses':facets,'target_cut_dual_pairing':ee(cut),
        'scope':'Exactly realized twelve-facet Coxeter polyhedron in H9 with a common strict timelike interior. No finite-volume or compactness claim.'}
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k in ('status','wall_count','gram_rank','gram_inertia','cap_polar_label','scope')}))
    return out

if __name__=='__main__':run()
