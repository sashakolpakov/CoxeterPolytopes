#!/usr/bin/env python3
"""Exact all-real/all-label exclusion for 88 specified right-edge families.

Sparse squarefree radicals reconstruct frames; integer outward intervals
certify strict inequalities. No floating geometry or scout module is imported.
The numerical continuous file only chooses short proof witnesses.
"""
import argparse,json,math,sys
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from collections import Counter
sys.path.insert(0,str(Path(__file__).resolve().parent))
import fresh_h22_multiquadratic_20261013_review as R
import fresh_h22_incidence_prism_20261013_verify as I
TAG='fresh_h22_contact_bridge_20261013';OUT='results/'+TAG+'_verify.json';SCOUT='results/'+TAG+'_continuous.json';SOURCE='results/fresh_h22_multiquadratic_20261013.json'
def neg(a):return R.scale(a,-1)
def sub(a,b):return R.add(a,neg(b))
def vsub(a,b):return[sub(x,y)for x,y in zip(a,b,strict=True)]
def div(a,b):
 assert b[0]>0 or b[1]<0
 return I.mul(a,(I.S*I.S//b[1],-((-I.S*I.S)//b[0])))
def insub(a,b):return I.add(a,I.neg(b))
@lru_cache(None)
def sqrtiv(n):return I.sqrt(I.rat(n))
def iv(a):return I.isum(I.scale(sqrtiv(n),c)for n,c in a.items())
def enc(a):return list(map(str,a))
def inverse(A):
 n=len(A);B=[list(map(F,row))+[F(i==j)for j in range(n)]for i,row in enumerate(A)]
 for j in range(n):
  assert B[j][j]>0
  q=B[j][j];B[j]=[x/q for x in B[j]]
  for i in range(n):
   if i!=j:
    c=B[i][j];B[i]=[x-c*y for x,y in zip(B[i],B[j],strict=True)]
 return [row[n:]for row in B]
def combine(co,rows):return[R.psum(R.scale(r[j],c)for r,c in zip(rows,co,strict=True))for j in range(len(rows[0]))]
def setup():
 d=R.read(SOURCE)
 for p,s in d['sources'].items():assert R.sha(p)==s,p
 raw=R.read(I.RAW)['coxeter_input'];Q=[[-int(x)for x in row]for row in raw['basis_gram']];M=[row+[0]for row in Q]+[[0]*22+[2]];W=[[int(x)for x in row]+[0]for row in raw['wall_coordinates_as_columns_transposed']];core=R.read(I.CORE);points=core['incidence_parts'][0];vp=core['old_projections'][1];parents=[]
 for rec in R.read(I.BANK)['candidates']:
  if rec['twice_height_square']!=4:continue
  bits=[rec['mask']>>i&1 for i in range(21)];twice=[rec['A']*vp[j]-sum(bits[i]*W[points[i]][j]for i in range(21))for j in range(22)];assert twice==rec['twice_old_projection']and all(x%2==0 for x in twice);parents.append([x//2 for x in twice])
 assert len(parents)==301;old=W+[r+[1]for r in parents]+[r+[-1]for r in parents];h=core['center']+[0];recipes=R.read(I.REC)['recipes'];caps={};norms={};T=R.sqrtq(3)
 for k in d['selected_original_types']:
  rec=recipes[k];u=list(map(F,rec['u']));v=list(map(F,rec['v']));q=F(rec['norm']);qu=I.pair(M,u,u);qv=I.pair(M,v,v);rad=q/qv-3*q*q*qu/(4*qv);r=[R.add(R.scale(T,-q*x/2),R.scale(R.sqrtq(rad),rec['radical_branch']*y))for x,y in zip(u,v,strict=True)];caps[k]=r if R.sign(r[-1])>0 else R.flip(r);norms[k]=q
 coupled=R.read('results/fresh_h22_pi7_coupled_caps_20261013.json');a,b=[list(map(F,recipes[17][s]))for s in('u','v')];Z=[R.add(R.scale(R.sqrtq(31),x/3),R.scale(T,-3*y/2))for x,y in zip(a,b,strict=True)];Z=R.vs(Z,d['normalized_frame_orientation']);du,dv=[list(map(F,row))for row in coupled['rational_duals']]
 caps[23]=R.vs(R.vadd(R.vs(Z,209),R.ratvec([-241*x for x in du])),F(1,120));caps[24]=R.vm(R.vadd(R.vs(Z,241),R.ratvec([-418*x for x in dv])),R.sqrtq(F(1,29281)));norms[23]=norms[24]=F(2)
 walls=list(map(R.ratvec,old));nn=[I.pair(M,r,r)for r in old]
 for k in d['selected_original_types']+[23,24]:walls.extend((caps[k],R.flip(caps[k])));nn.extend((norms[k],norms[k]))
 assert len(walls)==len(nn)==860;cov=[R.mv(M,r)for r in walls]
 for r,c,q in zip(walls,cov,nn,strict=True):assert R.pair(r,c)==R.C(q)
 vertices=[R.vadd(R.vs(Z,A),R.ratvec([B*x+C*y for x,y in zip(du,dv,strict=True)]))for A,B,C in [v['frame_coefficients']for v in d['compact_H2_vertices']]]
 return d,M,walls,cov,nn,h,vertices
def cutoff(a,b,lower):
 # a*T-b<=0, with a<0 for a lower bound, a>0 for upper.
 ai,bi=iv(a),iv(b);assert ai[1]<0 if lower else ai[0]>0
 return div(bi,ai)
def finite_upper(ai,bi,q):
 # Since a<0, a nonpositive anchor pairing with normalized cosine<=1
 # lies between roots of a²*T²+(-2ab-4q)*T+b²<=0.
 assert ai[1]<0
 A=I.mul(ai,ai);B=insub(I.scale(I.mul(ai,bi),-2),I.rat(4*q));C=I.mul(bi,bi);D=insub(I.mul(B,B),I.scale(I.mul(A,C),4));assert D[0]>0
 return div(I.add(I.neg(B),I.sqrt(D)),I.scale(A,2))
def cosine(ai,bi,q,T):
 root=I.sqrt(T)
 return div(I.add(I.neg(I.mul(ai,root)),div(bi,root)),I.sqrt(I.rat(4*q)))
def run(verify=False):
 source,M,walls,cov,nn,h,points=setup();scout=R.read(SCOUT);assert len(scout['rows'])==88;by={(r['vertex'],r['released']):r for r in scout['rows']};assert len(by)==88;counts=Counter();rows=[]
 for vi,src in enumerate(source['compact_H2_vertices']):
  inc=src['incident'];assert len(inc)==22;A=[]
  for i in inc:
   row=[]
   for j in inc:
    p=R.pair(walls[i],cov[j]);assert set(p)<= {1};row.append(p.get(1,F(0)))
   A.append(row)
  R.pd(A);AI=inverse(A);assert all(x>=0 for row in AI for x in row)
  p=points[vi];assert R.pair(p,R.mv(M,p))==R.C(src['norm'])and src['norm']<0
  v=R.vm(p,R.sqrtq(F(-1,src['norm'])));assert R.pair(v,R.mv(M,v))==R.C(-1)
  for k,released in enumerate(inc):
   info=by[(vi,released)];edge=[i for i in inc if i!=released];assert info['edge']==edge
   dual=combine(AI[k],[walls[j]for j in inc]);e=R.vm(dual,R.sqrtq(1/AI[k][k]));ec=R.mv(M,e);assert R.pair(e,ec)==R.C(1)and R.pair(v,ec)=={}
   assert all(R.pair(e,cov[j])=={}and R.pair(v,cov[j])=={}for j in edge)
   assert R.sign(R.pair(e,cov[released]))>0
   minus=vsub(v,e);plus=R.vadd(v,e);cache={}
   def coefficients(tag):
    if tag not in cache:
     if tag=='upper_height':aa,bb=neg(minus[-1]),neg(plus[-1])
     elif tag=='center':aa=R.pair(minus,R.ratvec(I.mv(M,h)));bb=R.pair(plus,R.ratvec(I.mv(M,h)))
     else:assert isinstance(tag,int)and 0<=tag<860;aa=R.pair(minus,cov[tag]);bb=R.pair(plus,cov[tag])
     cache[tag]=(aa,bb)
    return cache[tag]
   lower=info['active_lower'];upper=info['active_upper'];assert lower!='center'and upper!='center'
   li=I.rat(1)if lower=='outgoing'else cutoff(*coefficients(lower),True);L=max(I.S,li[0]);rec={'vertex':vi,'released':released,'edge':edge,'inverse_nonnegative':True,'lower_source':lower,'lower_cutoff':enc(li),'upper_source':upper,'anchors':[]}
   if upper is not None:
    ui=cutoff(*coefficients(upper),False);assert L>ui[1];rec['upper_cutoff']=enc(ui);rec['conclusion']='EMPTY_ACUTE_UPPER_CAP_INTERVAL';counts['empty_edges']+=1;rows.append(rec);print('empty',vi,released,flush=True);continue
   counts['nonempty_edges']+=1;proofs={q['anchor']:q for q in info['anchor_intervals']}
   for anchor in range(210,511):
    assert anchor not in edge;aa,bb=coefficients(anchor);ai,bi=iv(aa),iv(bb);U=finite_upper(ai,bi,nn[anchor]);qrec={'anchor':anchor,'finite_upper_bound':enc(U)};counts['anchor_tests']+=1
    if U[1]<L:qrec['conclusion']='FINITE_ANCHOR_INTERVAL_EMPTY';counts['empty_anchor_intervals']+=1
    else:
     witness=proofs[anchor]['forbidden_range_witness'];bad=witness['wall'];m=witness['gap_labels'][0];assert witness['gap_labels']==[m,m+1]and m in(2,3);ca,cb=coefficients(bad);ci,di=iv(ca),iv(cb);Tinterval=(L,U[1]);assert L>0
     # Derivative sign of c(T) is sign(-a*T-b). Certify monotonicity
     # on the whole retained interval, then enclose the two endpoint values.
     derivative=insub(I.neg(I.mul(ci,Tinterval)),di);assert derivative[0]>0 or derivative[1]<0
     cL=cosine(ci,di,nn[bad],(L,L));cU=cosine(ci,di,nn[bad],(U[1],U[1]));cr=(min(cL[0],cU[0]),max(cL[1],cU[1]));gaplo=I.rat(0)if m==2 else I.rat(F(1,2));gaphi=I.rat(F(1,2))if m==2 else I.sqrt(I.rat(F(1,2)));assert cr[0]>gaplo[1]and cr[1]<gaphi[0]
     qrec.update({'conclusion':'FORBIDDEN_COXETER_COSINE_ON_WHOLE_INTERVAL','interval':[str(L),str(U[1])],'witness_wall':bad,'gap_labels':[m,m+1],'cosine_bounds':enc(cr),'derivative_bounds':enc(derivative)});counts['forbidden_angle_intervals']+=1
    rec['anchors'].append(qrec)
   rows.append(rec);print('edge',vi,released,'anchors301',flush=True)
 assert counts=={'empty_edges':10,'nonempty_edges':78,'anchor_tests':23478,'empty_anchor_intervals':23003,'forbidden_angle_intervals':475}
 out={'status':'PASS_EXACT_ALL_REAL_RIGHT21_EDGE_FINITE_PARENT_BRIDGE_OBSTRUCTION','sources':{p:R.sha(p)for p in(SOURCE,SCOUT,I.RAW,I.CORE,I.BANK,I.REC,'results/fresh_h22_pi7_coupled_caps_20261013.json','scripts/'+TAG+'_verify.py','scripts/fresh_h22_multiquadratic_20261013_review.py','scripts/fresh_h22_incidence_prism_20261013_verify.py')},'interval_denominator':str(I.S),'counts':dict(counts),'edge_families':88,'anchors_per_feasible_edge':301,'rows':rows,'normalization':'r=(t(v-e)-t^-1(v+e))/sqrt2; T=t²>1. v²=-1,e²=1,v.e=0. Nonnegative inverse of the old positive rank22 Gram excludes t<=1 for any strict compatible completion; opposite e orientation has positive pairing with the released old wall.','scope':'No positive-height norm-two cap orthogonal to any21-wall edge support at the four specified actual finite vertices can retain the860 walls acutely and meet any of the301 positive-height rational parent caps at a finite Coxeter angle (even allowing parallel contact). All real parameters and all finite angle labels are covered in these88 families. No claim covers nonright edge profiles, other vertices, different anchor families, wall replacements or general H22 completion.'}
 if verify:assert out==R.read(OUT)
 else:(R.P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],dict(counts),flush=True);return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
