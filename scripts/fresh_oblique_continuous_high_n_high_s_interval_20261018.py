#!/usr/bin/env python3
"""All finite high-N parameters joined to the nine symbolic high-S families.

Directed physical intervals cover both N families on two fixed parents.
This is a necessary obstruction, not a new wall/chamber certificate.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_oblique_nonright_high_s_interval_20261018 as Z
P,Q,C,ROOT=Z.P,Z.Q,Z.C,Z.ROOT
TAG='fresh_oblique_continuous_high_n_high_s_interval_20261018'
V,add,sub,mul,div,neg,sqrt,dot=Z.V,Z.add,Z.sub,Z.mul,Z.div,Z.neg,Z.sqrt,Z.dot
def ng(g,n,a):
 frame=[10,67,69];inv,det=P.inverse([[g[i][j]for j in frame]for i in frame]);cv=[g[i][64]for i in frame]
 sig=sub(g[64][64],dot(cv,inv,cv));rhs=[V(-a),V(0),neg(n)];rho=sub(V(2),dot(rhs,inv,rhs));assert rho[0]>0 and sig[0]>0
 fac=sqrt(div(rho,sig));contacts={}
 positive8=None
 for j in g:
  cj=[g[i][j]for i in frame];lin=dot(rhs,inv,cj);res=mul(fac,sub(g[64][j],dot(cv,inv,cj)))
  contacts[j]=rhs[frame.index(j)]if j in frame else add(lin,res)
  if j==8:positive8=sub(lin,res);assert positive8[0]>0
 out={i:row.copy()for i,row in g.items()}
 for j in out:out[j][73]=contacts[j]
 out[73]={**contacts,73:V(2)}
 del out[70]
 for row in out.values():del row[70]
 return out,dict(frame=frame,inverse=inv,determinant=det,witness=64,witness_residual=sig,
  N_domain=n,N10_magnitude=a,rhs=rhs,norm_remainder=rho,opposite_N_sheet_old8=positive8,new_products=contacts)

def highrm(g):
 frame=[52,68,71];inv,det=P.inverse([[g[i][j]for j in frame]for i in frame]);rhs=[V(-1),V(0),V(0)];cv=[g[i][8]for i in frame]
 alpha=neg(dot(rhs,inv,rhs));sig=sub(g[8][8],dot(cv,inv,cv));assert alpha[0]>0 and sig[0]>0
 aa={j:dot(rhs,inv,[g[i][j]for i in frame])for j in g}
 bb={j:div(sub(g[8][j],dot(cv,inv,[g[i][j]for i in frame])),sqrt(sig))for j in g}
 z=(div(P.theta(6),sqrt(add(V(2),mul(V(3),alpha))))[0],div(V(2),sqrt(add(V(2),mul(V(4),alpha))))[1])
 vv=sub(V(1),mul(alpha,mul(z,z)));vx=sub(mul(z,aa[69]),bb[69]);vn=sub(mul(z,aa[73]),bb[73])
 mat=[[g[69][69],g[69][73],vx],[g[73][69],g[73][73],vn],[vx,vn,vv]]
 det3=P.determinant(mat)
 minor=sub(mul(g[69][69],vv),mul(vx,vx))
 H=sub(mul(sub(V(3),mul(V(2),g[71][71])),neg(det3)),mul(V(2),mul(mul(g[71][73],g[71][73]),minor)))
 # X.M and S.M are defining right duties; do not compare a rounded
 # inverse reconstruction of the latter with the singleton interval zero.
 assert g[69][71]==V(0) and g[71][71]==V(2) and rhs[2]==V(0)
 assert det3[1]<0,('highrm_det',det3)
 assert H[0]>0,('highrm_H',H)
 return dict(S_frame=frame,S_inverse=inv,S_determinant=det,alpha=alpha,residual8_square=sig,
  transformed_z_domain=z,V_square=vv,VX=vx,VN=vn,full_Gram=mat,determinant=det3,
  M_covector=[V(0),g[71][73],V(0)],middle_cofactor=minor,strict_obstruction=H)

def check(g,wc,a,n,inh):
 gn,nproof=ng(g,n,a);wn={**wc,73:profiles[0]};del wn[70]
 cube=C.complete_cube(gn,wn,wi,profiles,inh)
 # Remaining eight high-S boxes: fixed S frame independent of N.
 frame=[52,68,71];inv,det=P.inverse([[gn[i][j]for j in frame]for i in frame]);aa=[V(-1),V(0),V(0)];cc=[V(0),V(0),V(-1)];cv=[gn[i][8]for i in frame]
 sig=sub(gn[8][8],dot(cv,inv,cv))
 basis=dict(AA=dot(aa,inv,aa),AC=dot(aa,inv,cc),CC=dot(cc,inv,cc),A8=dot(aa,inv,cv),C8=dot(cc,inv,cv),residual8_square=sig)
 coeffs=[]
 for r in cube['recipes']:
  c=[r['negative_R_sheet_products'][i]for i in frame]
  coeffs.append(dict(index=r['index'],AR=dot(aa,inv,c),CR=dot(cc,inv,c),BR=div(sub(r['negative_R_sheet_products'][8],dot(cv,inv,c)),sqrt(sig))))
 rect=[Z.rectangle(basis,coeffs,f)for f in Z.FAMILIES]
 polar=highrm(gn)
 return dict(N=nproof,R_cube=cube,S_basis=basis,nonright_families=rect,high_RM=polar)


def run(saved=None):
 global profiles,wi
 prior=json.loads((ROOT/f'results/{P.TAG}.json').read_text())
 old=json.loads((ROOT/f'results/{P.H.TAG}.json').read_text())
 previous=json.loads((ROOT/f'results/{Z.TAG}.json').read_text())
 for d in(prior,old,previous):
  for path,h in d['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
 profiles=[[V(q)for q in row]for row in prior['profiles']]
 wi=[[V(q)for q in row]for row in prior['primary_projection']['white_inverse']]
 g0={int(a):{int(b):V(q)for b,q in row.items()}for a,row in prior['primary_projection']['projected_Gram'].items()}
 wc0={int(a):[V(q)for q in row]for a,row in prior['primary_projection']['old_white_covectors'].items()}
 assert len(prior['parents'])==len(old['records'])==2
 rows=[];attempts=0
 for parent,source in zip(prior['parents'],old['records']):
  assert parent['index']==source['index'] and parent['instruction_key']==source['instruction_key']
  g={i:r.copy()for i,r in g0.items()};wc=wc0.copy()
  for step in parent['steps']:Q.restore(g,wc,profiles,step)
  assert len(parent['steps'])==5
  inherited=next(r['continuous']for r in prior['records']if r['parent']==parent['index'])
  for a in(0,1):
   pending=[(P.theta(6)[0],F(2),0)]
   while pending:
    lo,hi,depth=pending.pop();attempts+=1
    try:out=check(g,wc,a,(lo,hi),inherited)
    except AssertionError:
     assert depth<8,('unresolved interval',parent['index'],a,lo,hi)
     mid=(lo+hi)/2;pending.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
    else:
     rows.append(dict(index=len(rows),parent=parent['index'],N10_magnitude=a,depth=depth,**out))
     print('CONTINUOUS_N_BOX',parent['index'],a,depth,flush=True)
 for p in range(2):
  for a in(0,1):
   parts=sorted(tuple(r['N']['N_domain'])for r in rows if r['parent']==p and r['N10_magnitude']==a)
   assert parts[0][0]==P.theta(6)[0] and parts[-1][1]==2
   assert all(x[1]==y[0]for x,y in zip(parts,parts[1:])) and all(x<y for x,y in parts)
 recipes=[r for b in rows for r in b['R_cube']['recipes']]
 families=[f for b in rows for f in b['nonright_families']]
 joins=[j for f in families for j in f['negative_S_R_joins']]
 assert len(rows)==66 and attempts==128 and len(recipes)==66*28 and len(joins)==66*8*28
 debit=old['exact_deletion_debit']
 assert debit['lost_finite_count']==19 and debit['compact_faces_by_dimension']=={'2':12,'3':1}
 paths={Path(__file__),Path(Z.__file__),Path(C.__file__),Path(Q.__file__),Path(P.__file__),Path(P.I.__file__)}
 paths.update(ROOT/path for d in(prior,old,previous)for path in d['sources'])
 paths.update(ROOT/f'{folder}/{stem}.{ext}'for stem in(Z.TAG,P.TAG,P.H.TAG)for folder,ext in(('scripts','py'),('results','json'),('results','md')))
 result=P.I.enc(dict(status='PASS_CONTINUOUS_HIGH_N_ALL_NINE_HIGH_S_FAMILIES_INTERVAL_OBSTRUCTION',
  sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
  interval_scale=P.I.PREC,N_parameter_domain=[P.theta(6)[0],F(2)],N_families=[[2,2,'m>=6'],[3,2,'m>=6']],
  nonright_high_S_families=Z.FAMILIES,parents=prior['parents'],records=rows,
  exact_deletion_debit=debit,deleted_walls=old['deleted_walls'],topology=old['topology'],
  summary=dict(parent_classes=2,parent_origins=4,continuous_N_families=4,attempted_boxes=attempts,covering_boxes=len(rows),
   low_R_recipe_boxes=len(recipes),low_R_sheet_boxes=2*len(recipes),SM_right_continuous_joins=len(recipes),
   nonright_SM_continuous_joins=len(joins),high_RM_polar_boxes=len(rows),
   minimum_opposite_N_old8=min(b['N']['opposite_N_sheet_old8'][0]for b in rows),
   minimum_positive_R_old8=min(r['positive_R_sheet_old8'][0]for r in recipes),
   maximum_SM_right_SR=max(r['continuous']['full_product_range'][1]for r in recipes),
   maximum_nonright_SM_SR=max(j['full_SR_range'][1]for j in joins),
   maximum_high_RM_determinant=max(b['high_RM']['determinant'][1]for b in rows),
   minimum_high_RM_polynomial=min(b['high_RM']['strict_obstruction'][0]for b in rows),
   unresolved=0,accepted=0,new_cones=0),
  scope='All finite high-N labels in both frozen N families, all nine symbolic high-S families and all native RM labels on two exact fixed parent chains. Physical rational intervals, no numerical R-bank rejection. Low-S joins for arbitrary N and other parents/incidences remain open. No accepted new wall or volume claim.'))
 if saved:assert result==saved
 return result


if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
 assert args.verify or not path.exists()
 out=run(json.loads(path.read_text())if args.verify else None)
 if not args.verify:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['summary'],flush=True)
