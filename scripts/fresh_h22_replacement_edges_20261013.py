#!/usr/bin/env python3
"""Exact complete22-direction star of the813 replacement's finite vertex."""
import argparse,json
from fractions import Fraction as F
import fresh_h22_core_replacement_20261013 as P
import fresh_h22_multiquadratic_20261013_review as R
from fresh_h22_contact_bridge_20261013_verify import inverse,combine

TAG='fresh_h22_replacement_edges_20261013';OUT='results/'+TAG+'.json'
SOURCE='results/fresh_h22_core_replacement_20261013.json'

def setup():
    src=R.read(SOURCE)
    for p,s in src['sources'].items():assert R.sha(p)==s,p
    Q,W,C,h,M,old=P.setup();ids=src['retained_base_indices']+['new_plus','new_minus']
    u=list(map(F,src['old_projection_numerator_u']));U=F(src['u_square'])
    plus=R.vm(R.ratvec(u+[0]),R.sqrtq(1/U));plus[-1]=R.sqrtq(F(1,2))
    walls={i:R.ratvec(old[i])for i in src['retained_base_indices']}
    walls.update(new_plus=plus,new_minus=R.flip(plus));cov={i:R.mv(M,w)for i,w in walls.items()}
    point=R.ratvec(list(map(F,src['finite_vertex']['old_coordinates']))+[0])
    hh=R.ratvec(h+[0]);inc=src['twenty_retained_incident']+['new_plus','new_minus']
    return src,M,walls,cov,hh,point,inc,U

def run(verify=False):
    src,M,W,cov,h,v,inc,U=setup();hc=R.mv(M,h)
    assert R.pair(v,R.mv(M,v))==R.C(-U)
    A=[]
    for i in inc:
        rr=[]
        for j in inc:
            a=R.pair(W[i],cov[j]);assert set(a)<={1};rr.append(a.get(1,F(0)))
        A.append(rr)
    R.pd(A);AI=inverse(A)
    rows=[];outside=[i for i in W if i not in inc]
    for new in ('new_plus','new_minus'):
        for i in outside:
            a=R.pair(W[new],cov[i]);q=R.pair(W[i],cov[i])
            assert R.sign(a)<0 and R.sign(R.add(R.mul(a,a),R.scale(q,-2)))>0
    for k,released in enumerate(inc):
        edge=[i for i in inc if i!=released]
        dual=combine(AI[k],[W[i]for i in inc]);dn=AI[k][k]
        assert R.pair(dual,R.mv(M,dual))==R.C(dn)
        assert R.pair(v,R.mv(M,dual))=={}
        assert all(R.pair(dual,cov[i])==R.C(int(i==released))for i in inc)
        alpha=R.sqrtq(U/dn)
        null=R.vadd(v,R.vm(dual,R.scale(alpha,-1)))
        assert R.pair(null,R.mv(M,null))=={}
        assert R.sign(R.pair(null,hc))<0
        pp={i:R.pair(null,cov[i])for i in W}
        assert all(pp[i]=={}for i in edge)
        assert all(R.sign(pp[i])<0 for i in W if i not in edge)
        for power in range(4,61,4):
            eps=F(1,2**power)
            x=R.vadd(v,R.vm(dual,R.scale(alpha,-1-eps)))
            xp={i:R.pair(x,cov[i])for i in W}
            if all(R.sign(a)<0 for i,a in xp.items()if i not in edge):break
        else:raise AssertionError(('no beyond-null witness',released))
        q=R.pair(x,R.mv(M,x));assert q==R.C(U*(2*eps+eps*eps))
        assert R.sign(R.pair(x,hc))<0 and all(xp[i]=={}for i in edge)
        retained_new=[i for i in edge if isinstance(i,str)];assert retained_new
        rows.append({'released':released,'edge':edge,'retained_new_caps':retained_new,
                     'dual_square':str(dn),'dual':list(map(R.enc,dual)),
                     'null_parameter':R.enc(alpha),'null_point':list(map(R.enc,null)),
                     'null_center_pairing':R.enc(R.pair(null,hc)),
                     'null_zero_walls':edge,'strict_nonedge_null_pairings':len(W)-len(edge),
                     'beyond_null_factor':str(1+eps),'spacelike_point':list(map(R.enc,x)),
                     'spacelike_square':R.enc(q),'spacelike_center_pairing':R.enc(R.pair(x,hc)),
                     'status':'CERTIFIED_OPEN_EDGE_NO_EXISTING_FINITE_OR_IDEAL_ENDPOINT'})
    assert len(rows)==22
    out={'status':'PASS_EXACT_COMPLETE22_EDGE_STAR_ALL_OPEN',
         'sources':{p:R.sha(p)for p in(SOURCE,'scripts/'+TAG+'.py','scripts/fresh_h22_core_replacement_20261013.py','scripts/fresh_h22_multiquadratic_20261013_review.py','scripts/fresh_h22_contact_bridge_20261013_verify.py')},
         'chamber_facets':813,'known_actual_start_vertex':1,'incident_walls':inc,
         'vertex_square':str(-U),'directions':22,'finite_other_endpoints':0,'ideal_other_endpoints':0,'open_directions':22,
         'rows':rows,'scope':'Complete outgoing edge star of this one actual finite vertex, not the complete chamber vertex graph. Each null endpoint has exactly21 incident walls and strict signs against all others; explicit feasible spacelike continuation certifies an unbounded direction, not a bounded cusp. Every edge retains a newcap and all nonincident oldwalls are strictlydotted to thatcap, giving a structural proof that none can terminate theedge even at infinity. No newcap or finite-volume example is claimed.'}
    if verify:assert out==R.read(OUT)
    else:(R.P/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'22open',flush=True);return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');run(ap.parse_args().verify)
