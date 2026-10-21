#!/usr/bin/env python3
"""Exact low-Z8 complement to the frozen continuous high-Z8 gate.

The union is complete for the stated fixed prefix and RZ<=5 polar-Y duty,
not for arbitrary high RZ or other prefix geometry.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_old8_triangle_double_flip_20261017 as SC
import fresh_old30_connector_native_gate_20261017 as E
Q,P,U,I,G,ROOT=E.Q,E.P,E.U,E.I,E.G,E.ROOT
TAG='fresh_joint_rz_low_z8_gate_20261017'
HIGH='fresh_joint_rz_high_z8_gate_20261017'

def run():
    pp=ROOT/'results'/f'{SC.TAG}.json';cp=ROOT/'results/fresh_r_circle_admissibility_20261017.json'
    hp=ROOT/'results'/f'{HIGH}.json'
    parent=json.loads(pp.read_text());circle=json.loads(cp.read_text());high=json.loads(hp.read_text())
    assert high['status']=='PASS_ALL_42_HIGH_Z8_FAMILIES_EXCLUDED' and not high['unresolved']
    for path,h in high['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    assert circle['exact_circle']['necessary_labels']==list(range(2,16))
    a,_,c,x,_,_=map(np.array,parent['prefixes'][0]['normals'])
    rdata={r['label']:r for r in circle['numerical_R_candidates']}
    # The polar-frame signature depends only on these literal contacts,
    # not on the R-circle coordinate.  Certify this small table once.
    domain=[];dc=Counter();negative=[]
    for pi,m8,mr,mx in product(range(7),range(2,6),range(2,6),range(2,6)):
        if SC.P.corner3(0,0,pi,2,m8,mx)is None:continue
        if SC.P.corner3(0,0,pi,2,m8,2)is None or SC.P.corner3(0,0,pi,2,mr,2)is None:continue
        prof=G.profile(pi);q=U.rnd(2-G.white(prof,prof));t8,tr=Q.theta(m8),Q.theta(mr)
        rem=U.rnd(q-(t8*t8+tr*tr)/2)
        if rem.l>0:kind='positive_definite_polar_frame'
        elif rem.h<0:
            kind='lorentzian_polar_frame';negative.append((pi,m8,mr,mx))
        else:
            assert (pi,m8,mr,mx)in[(0,4,4,2),(0,4,4,3)]
            # Literal affine C2 Gram: q=2, theta4^2=2; positive kernel
            # (1,1,sqrt2). A strict common center prevents this positive
            # null combination from vanishing, so its orthogonal line
            # is null and cannot contain norm2 Y.
            assert q.l==q.h==2 and rem.l<=0<=rem.h
            kind='affine_polar_frame_no_spacelike_orthogonal_Y'
        dc[kind]+=1;domain.append(dict(profile=pi,labels=[m8,mr,mx],projected_Z_square=q.enc(),polar_schur=rem.enc(),kind=kind))
    assert tuple(dc[k]for k in['positive_definite_polar_frame','affine_polar_frame_no_spacelike_orthogonal_Y','lorentzian_polar_frame'])==(96,2,28)
    _,_,ns,_,_,_=SC.H.setup();base=Q.Engine();records=[];counts=Counter()
    for m in circle['exact_circle']['necessary_labels']:
        g=E.Gram(base,ns)
        g.add('C',[7,30,10],0,[2,2,9],c)
        g.add('A',[7,8,'C'],0,[2,2,3],a)
        g.add('R',[8,44,63],0,[2,m,2],np.array(rdata[m]['normal']))
        assert g.pair('R',10).h<0
        assert g.pair(8,'R').l==g.pair(8,'R').h==0
        g.add('X',[8,'A','C'],0,[2,3,2],x)
        assert g.pair(8,'X').l==g.pair(8,'X').h==0
        assert all(F(r['remaining_square'][0])>0 for r in g.recipes)
        cases=[]
        for pi,m8,mr,mx in negative:
            rows,info=g.norm_frame([8,'R','X'],pi,[-Q.theta(m8),-Q.theta(mr),-Q.theta(mx)])
            counts['frame_recipes']+=1
            if not rows:
                cases.append(dict(profile=pi,labels=[m8,mr,mx],reason='negative_remaining_square',frame_certificate=info));continue
            branches=[];prof=G.profile(pi)
            for sg,row in rows:
                bad=[]
                for wall in E.TARGETS:
                    value=U.rnd(row[wall]+G.white(prof,g.cp[wall]));rejection=Q.exclusion(value)
                    if rejection:bad.append(dict(wall=wall,**rejection))
                assert bad,('unresolved_low_Z_native_branch',m,pi,m8,mr,mx,sg)
                counts['native_rejected_branches']+=1
                branches.append(dict(sign=sg,rejection=bad[0]))
            cases.append(dict(profile=pi,labels=[m8,mr,mx],reason='both_native_norm_branches_excluded',frame_certificate=info,branches=branches))
        records.append(dict(R44=m,exact_R10=g.pair('R',10).enc(),exact_prefix_recipes=g.recipes,cases=cases))
        print('LOW_Z8',m,dict(counts),flush=True)
    assert counts['frame_recipes']==392 and counts['native_rejected_branches']==784
    paths=[Path(__file__),Path(SC.__file__),Path(SC.P.__file__),Path(E.__file__),Path(Q.__file__),Path(P.__file__),Path(G.__file__),Path(G.E.X.__file__),
        pp,cp,ROOT/'scripts/fresh_r_circle_admissibility_20261017.py',G.CONTROL,G.PARENT,SC.H.PARENT,SC.H.BOUNDARY]
    paths +=[ROOT/f'{folder}/{HIGH}.{ext}'for folder,ext in[('scripts','py'),('results','json'),('results','md')]]
    return dict(status='PASS_EXACT_LOW_Z8_AND_COMPLETE_FIXED_PREFIX_LOW_RZ_UNION',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},counts=dict(counts),profile_domain_counts=dict(dc),profile_domain=domain,records=records,
        exact_affine_limit=dict(profile=0,Z8=4,ZR=4,ZX=[2,3],projected_gram=[['2','0','-sqrt2'],['0','2','-sqrt2'],['-sqrt2','-sqrt2','2']],positive_kernel=['1','1','sqrt2'],common_center_required=True),
        union_domain='Fixed C10=9,AC3,X8/XA/XC2/3/2, profile0; zero-profile R rightS6+8+63 withR44 in2..15; RZ<=5; actual finite8/X/Z and Y/8/Z,Y/R/Z corners; norm2 zero-profileY right8/R/Z. All seven Z A2 profiles and every finite Z8/ZX label are covered, using this finite low table and the42 frozen high-Z8 families. HighZX is impossible by finite-corner isolation plus the positive-definite8/R/Z frame.',
        completeness='The low table has96 positive polar frames,2 affine polar frames,28 Lorentzian frames independent ofR. Positive frames have timelike orthogonal complement. Affine frames have null complement under a common strict center. Every Lorentzian frame is solved in both exact norm branches, with native A2 offsets restored, and a strict old-native rejection. All14 permittedR sheets are explicitly checked byR10<0. HighZ8 is the frozen42-family continuous certificate. No candidate label cutoff is used in the union.',
        scope='Exact necessary-condition obstruction for this fixed-prefix polar-Y incidence withRZ<=5. It does not exclude highRZ>=6, changedA/C/X, nonzero-profileR, nonpolarY or different disk incidences. No parent numerical chamber is promoted; no new finite-volume or improved-control result is claimed. Shared physical interval and finite-profile implementations are inherited.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],flush=True)
