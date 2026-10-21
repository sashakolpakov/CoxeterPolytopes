#!/usr/bin/env python3
"""Free dotted-distance caps for the released ten-wall H9 seed.

Finite seed pairings are 0,-1,-phi; one other pairing is -x, x>=2.
The cap norm equation determines x quadratically over Q(phi). Exact checks
certify accepted caps, mutual conflicts and strict rational escapes. The
numerical candidate proposal stage has bounded support and no exhaustion
claim. Eligible old caps form a SUPERSET of all Coxeter-compatible ones.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations,product
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,sub,neg,scale,mul,div,sign,total,dot,mv,pair,matrix,
)

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/enumerations/fresh_h8_release_h9_caps_expanded_20261013.json'
ESCAPE=ROOT/'data/enumerations/fresh_h8_release_h9_volume_expanded_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h8_release_h9_free_caps_20261013.json'
real=lambda x:float(x[0])+(1+5**.5)/2*float(x[1])
encvec=lambda v:[enc(x)for x in v]
EZ=(Z,Z);EO=(O,Z)

def eadd(x,y):return add(x[0],y[0]),add(x[1],y[1])
def eneg(x):return neg(x[0]),neg(x[1])
def esub(x,y):return eadd(x,eneg(y))
def escale(x,k):return scale(x[0],k),scale(x[1],k)
def emul(x,y,d):return add(mul(x[0],y[0]),mul(d,mul(x[1],y[1]))),add(mul(x[0],y[1]),mul(x[1],y[0]))
def esign(x,d):
    a,b=x;sa,sb=sign(a),sign(b)
    if not sb:return sa
    if not sa or sa==sb:return sb
    return sa*sign(sub(mul(a,a),mul(d,mul(b,b))))
def ereal(x,d):return real(x[0])+real(d)**.5*real(x[1])
def eencode(x):return {'rational':enc(x[0]),'radical':enc(x[1])}
def polysign(coeffs,x,d):
    p=EZ
    for c in coeffs:p=eadd(emul(p,x,d),(scale(O,c),Z))
    return esign(p,d)

def potentially_compatible(p,d):
    if esign(p,d)>0:return False
    square=emul(p,p,d)
    if esign(esub(square,(scale(O,4),Z)),d)>=0:return True
    # If q=4cos²(pi/m) belongs to a field of degree dividing four, its
    # cyclotomic degree phi(m)/2 divides four. These are all possibilities.
    # Polynomial conjugates are intentionally retained, making this an
    # overinclusive gate rather than a false rejection of a valid angle.
    for q in (0,1,2,3):
        if esign(esub(square,(scale(O,q),Z)),d)==0:return True
    t=esub(square,(scale(O,2),Z))
    polynomials=((1,1,-1),(1,-1,-1),(1,0,-2),(1,0,-3),
                 (1,-1,-4,4,1),(1,0,-4,0,2),(1,0,-5,0,5),
                 (1,0,-4,0,1),(1,1,-4,-4,1))
    return any(polysign(cs,t,d)==0 for cs in polynomials)

def bsign(x,d,e):
    a=(x[0],x[1]);b=(x[2],x[3]);sa,sb=esign(a,d),esign(b,d)
    if not sb:return sa
    if not sa or sa==sb:return sb
    difference=esub(emul(a,a,d),(mul(e,emul(b,b,d)[0]),mul(e,emul(b,b,d)[1])))
    return sa*esign(difference,d)
def bmul(x,y,d,e):
    out=[Z]*4
    for i,a in enumerate(x):
        for j,b in enumerate(y):
            v=mul(a,b)
            if i&j&1:v=mul(v,d)
            if i&j&2:v=mul(v,e)
            out[i^j]=add(out[i^j],v)
    return tuple(out)

def proposals(I,y,max_support=4):
    In=np.array([[real(x)for x in row]for row in I]);yn=np.array([real(x)for x in y]);Iy=In@yn
    accepted=[];count=0
    for free in range(10):
        rest=[i for i in range(10)if i!=free]
        for k in range(max_support+1):
            for ids in combinations(rest,k):
                for labs in product((O,P),repeat=k):
                    fixed=[Z]*10
                    for j,m in zip(ids,labs):fixed[j]=m
                    fixed=tuple(fixed);pn=np.array([real(x)for x in fixed])
                    a=I[free][free];b=scale(dot(I[free],fixed),2);c=sub(pair(fixed,I,fixed),scale(O,2))
                    da=real(a);db=real(b);dc=real(c);discn=db*db-4*da*dc
                    if discn<=0:continue
                    for branch in (1,-1):
                        xn=(-db+branch*discn**.5)/(2*da)
                        if xn<2-1e-8:continue
                        count+=1;pn[free]=xn
                        if pn@Iy<=1e-8:pn[free]=0;continue
                        disc=sub(mul(b,b),scale(mul(a,c),4));assert sign(disc)>0
                        x=(div(neg(b),scale(a,2)),div(scale(O,branch),scale(a,2)))
                        assert esign(esub(x,(scale(O,2),Z)),disc)>0
                        equation=eadd(eadd(emul((a,Z),emul(x,x,disc),disc),emul((b,Z),x,disc)),(c,Z))
                        assert esign(equation,disc)==0
                        pa=tuple(neg(add(v,x[0]if j==free else Z))for j,v in enumerate(fixed))
                        pb=tuple(neg(x[1])if j==free else Z for j in range(10))
                        ca,cb=mv(I,pa),mv(I,pb)
                        norm=(add(dot(pa,ca),mul(disc,dot(pb,cb))),scale(dot(pa,cb),2))
                        assert esign(esub(norm,(scale(O,2),Z)),disc)==0
                        cut=(neg(dot(ca,y)),neg(dot(cb,y)));assert esign(cut,disc)>0
                        accepted.append({'free':free,'fixed':fixed,'polynomial':(a,b,c),'disc':disc,'x':x,
                                         'pa':pa,'pb':pb,'ca':ca,'cb':cb,'old_escape_pairing':cut})
                        pn[free]=0
    return accepted,count

def load_extra():
    out=[]
    for name in ('2','2phi'):
        path=ROOT/f'data/enumerations/fresh_h8_release_h9_quadratic_caps_{name}_20261013.json'
        d=json.loads(path.read_text());rad=parse(d['extension_radicand'])
        for cap in d['caps']:
            out.append((rad,tuple(map(parse,cap['seed_coefficients_rational'])),tuple(map(parse,cap['seed_coefficients_radical']))))
    return out

def strict_escape(I,old,cap,extra,center):
    d=cap['disc'];a,b=cap['ca'],cap['cb']
    qn=np.array([[real(x)for x in row]for row in I]);cn=np.array([real(x)for x in center])
    normals=[[real(x)for x in row]for row in old]
    normals.append([real(x)+real(d)**.5*real(z)for x,z in zip(a,b)])
    normals.extend([real(x)+real(rad)**.5*real(z)for x,z in zip(aa,bb)]for rad,aa,bb in extra)
    A=np.array(normals);A=A/np.linalg.norm(A,axis=1)[:,None]
    rng=np.random.default_rng(20261013)
    for trial in range(100):
        prop=linprog(rng.normal(size=10),A_ub=-A,b_ub=np.zeros(len(A)),A_eq=np.ones((1,10)),b_eq=[1],bounds=(0,None),method='highs')
        assert prop.success,prop.message
        oldq=prop.x@qn@prop.x
        if oldq<=1e-8:continue
        epsilon=.01
        while (prop.x+epsilon*cn)@qn@(prop.x+epsilon*cn)<oldq/2:epsilon/=2
        target=prop.x+epsilon*cn
        for den in (10**8,10**10,10**12,10**14):
            y=tuple(scale(O,F(round(float(t)*den),den))for t in target)
            q=pair(y,I,y)
            if sign(q)<=0 or any(sign(x)<=0 for x in y):continue
            if any(sign(dot(row,y))<=0 for row in old):continue
            value=(dot(a,y),dot(b,y))
            if esign(value,d)<=0:continue
            if any(esign((dot(aa,y),dot(bb,y)),rad)<=0 for rad,aa,bb in extra):continue
            assert sign(pair(y,I,center))<0
            return {'status':'EXACT_STRICT_SPACELIKE_ESCAPE','y':encvec(y),'norm':enc(q),
                    'center_pairing':enc(pair(y,I,center)),'free_cap_constraint_value':eencode(value),
                    'extra_cap_constraint_values':[eencode((dot(aa,y),dot(bb,y)))for rad,aa,bb in extra],
                    'all_seed_and_selected_cap_constraints_strict':True}
    return {'status':'UNDECIDED_NO_CERTIFIED_ESCAPE'}

def build():
    raw=json.loads(BASE.read_text());esc=json.loads(ESCAPE.read_text());I=matrix(raw['seed_gram_inverse']);R=matrix(raw['seed_normals']);Q=matrix(raw['host_gram'])
    old=[tuple(map(parse,c['seed_coefficients']))for c in raw['caps']];oldp=[tuple(map(parse,c['seed_pairings']))for c in raw['caps']]
    y=tuple(map(parse,esc['y_coordinates']));center=tuple(map(parse,esc['strict_center_y']));extra=load_extra()
    assert sign(pair(y,I,y))>0 and all(sign(dot(a,y))>=0 for a in old)
    assert sign(pair(center,I,center))<0 and all(sign(x)>0 for x in center)
    assert all(sign(dot(a,center))>0 for a in old)
    assert all(esign((dot(a,y),dot(b,y)),d)>=0 and esign((dot(a,center),dot(b,center)),d)>0 for d,a,b in extra)
    caps,count=proposals(I,y);cases=[]
    for index,cap in enumerate(caps):
        d=cap['disc'];assert esign((dot(cap['ca'],center),dot(cap['cb'],center)),d)>0
        eligible=[j for j,p in enumerate(oldp)if potentially_compatible((dot(cap['ca'],p),dot(cap['cb'],p)),d)]
        gate=strict_escape(I,[old[j]for j in eligible],cap,extra,center)
        ra=tuple(total(mul(c,r[j])for c,r in zip(cap['ca'],R))for j in range(10))
        rb=tuple(total(mul(c,r[j])for c,r in zip(cap['cb'],R))for j in range(10))
        assert esign((sub(add(pair(ra,Q,ra),mul(d,pair(rb,Q,rb))),scale(O,2)),scale(pair(ra,Q,rb),2)),d)==0
        case={'index':index,'free_pairing_seed_index':cap['free'],'free_pairing_original_wall':raw['retained_original_walls'][cap['free']],
              'fixed_pairing_magnitudes':encvec(cap['fixed']),'norm_equation_coefficients_descending':encvec(cap['polynomial']),
              'discriminant':enc(d),'free_magnitude':eencode(cap['x']),'free_magnitude_approx':ereal(cap['x'],d),
              'normal_rational':encvec(ra),'normal_radical':encvec(rb),
              'seed_coefficients_rational':encvec(cap['ca']),'seed_coefficients_radical':encvec(cap['cb']),
              'old_escape_pairing':eencode(cap['old_escape_pairing']),
              'potentially_compatible_old_cap_indices':eligible,'volume_gate':gate}
        cases.append(case);print({'cap':index,'eligible_old_caps':len(eligible),'gate':gate['status']},flush=True)
    conflicts=[]
    for i,j in combinations(range(len(caps)),2):
        a,b=caps[i],caps[j];d,e=a['disc'],b['disc']
        p=(dot(a['pa'],b['ca']),dot(a['pb'],b['ca']),dot(a['pa'],b['cb']),dot(a['pb'],b['cb']))
        if bsign(p,d,e)>0:reason='positive_outward_pairing'
        else:
            assert bsign(p,d,e)<0 and bsign((add(p[0],P),*p[1:]),d,e)>0
            sq=bmul(p,p,d,e)
            assert bsign((sub(sq[0],O),*sq[1:]),d,e)!=0
            assert bsign((sub(sq[0],scale(O,2)),*sq[1:]),d,e)!=0
            reason='finite_nonright_angle_below_pi5_cosine_but_neither_label3_nor4'
        conflicts.append({'pair':[i,j],'reason':reason,'biquadratic_pairing_coefficients':encvec(p)})
    assert len(conflicts)==len(caps)*(len(caps)-1)//2
    return {'schema':'released-H8-free-dotted-H9-caps-v1','base_source':str(BASE.relative_to(ROOT)),
            'old_escape_source':str(ESCAPE.relative_to(ROOT)),'maximum_fixed_pairing_support':4,
            'numeric_admissible_free_distance_proposals':count,'exact_old_escape_cutters':len(caps),
            'cases':cases,'all_free_cap_pairs_conflict':True,'free_cap_conflicts':conflicts,
            'all_eligible_supersystems_have_exact_escapes':all(c['volume_gate']['status']=='EXACT_STRICT_SPACELIKE_ESCAPE'for c in cases),
            'eligible_gate_note':'All possible finite Coxeter angle degrees dividing four are retained through their minimal polynomials, including conjugate values; thus eligibility is an overinclusive exact gate.',
            'scope':'No finite-volume Coxeter system using the ten seed walls and any subset of the archived434 base caps, two quadratic caps and16 free caps. The finite numerical proposal scan is not an exhaustion theorem for other caps or distances.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=build()
    if a.verify:assert out==json.loads(OUTPUT.read_text())
    else:OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
    print({'free_caps':out['exact_old_escape_cutters'],'pair_conflicts':len(out['free_cap_conflicts']),'all_escape_gates_pass':out['all_eligible_supersystems_have_exact_escapes']},flush=True)

if __name__=='__main__':main()
