#!/usr/bin/env python3
"""A norm-two cap from a coupledD4 polar repair with a free dotted anchor."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from fresh_h20_coupled_completion_20261013 import setup as setup98,read,sha,ROOT,F10,num
from fresh_h20_coupled_completion_20261013_exact import OUT as CUSP_SOURCE
from fresh_h20_gap_next_20261013_exact import root_checks,cusp_signs,base_point,strictify
from fresh_h20_edge_family_completion_20261013_cusp_repairs import base_root,sumvec
from fresh_dimension_bridge_20261013_verify import E,Z,O,vec,parse,dot,mv
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot,same
from fresh_h20_targeted_escape_20261013_volume import check
SOURCE='data/enumerations/fresh_h20_gap_next_20261013_exact.json'
SCOUT='data/enumerations/fresh_h20_coupled_completion_20261013_dotted_volume_scout.json'
OUT='data/enumerations/fresh_h20_coupled_completion_20261013_dotted_exact.json'
def setup():
 Q,h,hc,roots=setup98()
 for c in read(CUSP_SOURCE)['new_caps']:roots.append((F10,tuple(map(F10.parse,c['normal']))))
 assert len(roots)==104
 return Q,h,hc,roots

def run():
 Q,h,hc,roots=setup();p=tuple(E(F(15,4))*x for x in vec(read(SOURCE)['final_hyperideal_point']['vector']))
 leaves=[tuple(x/E(F(3))for x in base_root(*roots[i]))for i in(26,28,42)];C=tuple(E(Z.a,F(1,3))*x for x in base_root(*roots[49]));S=sumvec([(O,x)for x in leaves])
 bs=leaves+[C];assert [[dot(a,mv(Q,b))for b in bs]for a in bs]==[[E(F(2))if i==j else -O if i==3 or j==3 else Z for j in range(4)]for i in range(4)]
 assert dot(p,mv(Q,p))==E(F(97,16))and all(dot(p,mv(Q,l))==-O for l in leaves)
 assert dot(p,mv(Q,C))==E(F(-17,4),F(11,12))
 u=sumvec([(E(F(3)),p),(O,S)]);v=sumvec([(O,C),(-O,p)])
 A=dot(v,mv(Q,v));BB=E(F(2))*dot(u,mv(Q,v));CC=dot(u,mv(Q,u))-E(F(2));D=BB*BB-E(F(4))*A*CC
 assert A.sign()>0 and D.sign()>0
 f=Extension(B,D);b=(-BB/(E(F(2))*A),-O/(E(F(2))*A))
 assert f.sign(f.sub(b,f.rational(F(119,100))))>0 and f.sign(f.sub(f.rational(F(6,5)),b))>0
 r=tuple(f.add(f.lift(x),f.mul(b,f.lift(y)))for x,y in zip(u,v,strict=True));q,center,pairs=root_checks(f,r,Q,hc,roots);assert same(f,q,f.rational(2))
 target=read(SOURCE)['strict_escape'];old=check(vec(target['vector']),Q,hc,roots);assert old['norm']==target['norm'];cut=fdot(f,r,tuple(map(f.lift,mv(Q,vec(target['vector'])))));assert f.sign(cut)>0
 assert all(pairs[i]['kind']=='3'for i in(26,28,42,7,11,13,17,18));assert pairs[23]['kind']=='dotted'and pairs[49]['kind']=='dotted'
 cusp=cusp_signs(f,r,Q);roots.append((f,r))
 ids=read(SCOUT)['full_bank_escape']['active_indices'];vpoint,pv=base_point(Q,hc,roots,ids);escape=strictify(vpoint,Q,h,hc,roots)
 return {'status':'EXACT_105_WALL_H20_FREE_DOTTED_ANCHOR_CUTTER_WITH_STRICT_ESCAPE','source':CUSP_SOURCE,'source_sha256':sha(CUSP_SOURCE),'polar_source':SOURCE,'polar_sha256':sha(SOURCE),'active_set_source':SCOUT,'active_set_sha256':sha(SCOUT),'retained_walls':104,'wall_count':105,'formula':{'polar':'(15/4)*old_hyperideal_point','leaves':[26,28,42],'leaf_scale':'1/3','center_wall':49,'center_scale':['0','1/3'],'u':'3*p+sum(leaves)','v':'center-p','root':'u+b*v','b_quadratic':[A.enc(),BB.enc(),CC.enc()],'b_rational_isolation':['119/100','6/5'],'b':f.encode(b)},'new_cap':{'wall':104,'field':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':D.enc()},'normal':list(map(f.encode,r)),'norm':f.encode(q),'strict_center_pairing':f.encode(center),'old_pair_checks':pairs,'pair_counts':dict(sorted(Counter(t['kind']for t in pairs).items())),'inherited_cusp_strict_pairings':cusp},'strict_target':dict(old,new_pairing=f.encode(cut)),'final_hyperideal_point':pv,'strict_escape':escape,'scope':'Exact norm-two target cutter retains104 walls. It repairs three D4 leaf angles to pi/3 while the wall23 and central-wall49 anchor relations are strictly dotted. All104 new pair gates,107 inherited cusp signs,strict target cut,and final105 escape are exact. Infinite volume; no complete search or vertex census.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print(d['new_cap']['pair_counts']);print('escape',d['strict_escape']['norm'])
