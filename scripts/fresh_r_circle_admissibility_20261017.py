#!/usr/bin/env python3
"""Exact finite Coxeter-label bound for the rightS6+8+63 R circle."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import numpy as np
import fresh_old8_triangle_double_flip_20261017 as SC
import fresh_old30_connector_native_gate_20261017 as E
Q,P,U,I,G,ROOT=E.Q,E.P,E.U,E.I,E.G,E.ROOT
TAG='fresh_r_circle_admissibility_20261017'

def run():
    pp=ROOT/'results'/f'{SC.TAG}.json';parent=json.loads(pp.read_text())
    _,labels,ns,h,pro,_=SC.H.setup();r=np.array(parent['prefixes'][0]['normals'][1])
    angle=np.pi/15;k=(r+np.cos(angle)*ns[44])/np.sin(angle)
    base=Q.Engine();g=E.Gram(base,ns);g.add('K',[8,44,63],0,[2,2,2],k)
    assert g.pair(44,10).l==g.pair(44,10).h==0
    assert g.pair(8,63).h<-2
    assert g.pair(44,8).l==g.pair(44,8).h==g.pair(44,63).l==g.pair(44,63).h==0
    assert g.pair('K',10).h<0 and g.pair('K',53).h<0 and g.pair(44,53).h<0
    co=Q.theta(46)/2;si=Q.sqrt_closed(1-co*co)
    value=U.rnd(-co*g.pair(44,53)+si*g.pair('K',53));assert value.l>0
    finite_rejections=[]
    for m in range(16,46):
        cm=Q.theta(m)/2;sm=Q.sqrt_closed(1-cm*cm)
        vm=U.rnd(-cm*g.pair(44,30)+sm*g.pair('K',30));proof=Q.exclusion(vm)
        assert proof is not None and proof['kind']=='strict_consecutive_angle_gap'
        finite_rejections.append(dict(label=m,wall=30,rejection=proof))
    # The other orientation has positive product with old10 for every
    # alpha in(0,pi/2]. On this orientation f(alpha)=Acos(alpha)-Bsin(alpha)
    # has strictly negative derivative; f(pi/46)>0 excludes every m>=46.
    rows=[]
    for m in range(2,46):
        v=-np.cos(np.pi/m)*ns[44]+np.sin(np.pi/m)*k
        error=max(abs(float(v@SC.M@v)-2),abs(float(v@SC.M@ns[8])),abs(float(v@SC.M@ns[63])),abs(float(v@SC.M@ns[44])+SC.B.theta(m)))
        assert error<1e-7
        vals=SC.T.native_values(v,0,labels,ns,pro);bad=SC.T.gates(vals,[8,44,63])
        center=float(v@SC.M@h)
        rows.append(dict(label=m,normal=v.tolist(),native_products=[[i,z]for i,z in vals.items()],bad_native=bad,center_product=center,native_and_center_pass=not bad and center< -1e-8))
    paths=[Path(__file__),Path(SC.__file__),pp,Path(E.__file__),Path(Q.__file__),Path(G.__file__),Path(G.E.X.__file__),G.CONTROL,G.PARENT,SC.H.PARENT,SC.H.BOUNDARY]
    return dict(status='PASS_EXACT_R_CIRCLE_LABEL_BOUND_WITH_NUMERICAL_NATIVE_FILTER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      exact_circle=dict(K_recipe=g.recipes[0],old8_63=g.pair(8,63).enc(),old44_10=g.pair(44,10).enc(),K_10=g.pair('K',10).enc(),old44_53=g.pair(44,53).enc(),K_53=g.pair('K',53).enc(),cos_pi46=co.enc(),sin_pi46=si.enc(),R46_53=value.enc(),finite_rejections=finite_rejections,necessary_labels=list(range(2,16))),
      numerical_R_candidates=rows,numerical_native_and_center_labels=[r['label']for r in rows if r['native_and_center_pass']],
      theorem='For norm2 R rightS6+8+63 and a coherent acute Coxeter contact with retained44, R=-cos(pi/m)44 plus or minus sin(pi/m)K for integer m>=2. The minus sheet is strictly obtuse to10. On the plus sheet, old53 product=Acos(alpha)-Bsin(alpha), A,B>0, is strictly decreasing on[0,pi/2]; its exact value atpi/46 is positive, excluding allm>=46 and the endpointR=-44. Each remainingm16..45 is excluded by an exact old30 Coxeter gap. Thus every such R has m in2..15. Dotted contact with44 is impossible in the positive circle plane.',
      scope='Exact necessary finite label bound for the specified zero-profile R circle, using inherited physical outward rational intervals. The subsequent all-old native/center filter is numerical discovery only; it is not an exact classification of the retained labels. No new chamber or volume result is asserted.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['numerical_native_and_center_labels'],flush=True)
