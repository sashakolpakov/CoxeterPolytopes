#!/usr/bin/env python3
"""All-real coupling bounds for the 51 isolated dihedral endpoint tails.

The scout supplies only names of two supporting inequalities. Each dual
identity is reconstructed by an exact two-by-two solve with outward
intervals. No optimizer, numerical coordinates, norm bound or angle-label
cutoff occurs in the proof.
"""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_component_bridge_20261013_verify as N
V=N.V;K=N.K;BASE=N.PREFIX;OUT=BASE+'_dihedral_verify.json'

def run(verify=False):
 d=V.read(BASE+'_caps.json');cert=N.check(d);assert cert==V.read(BASE+'_verify.json')
 fc=V.read(BASE+'_continuous_verify.json')
 for path,hsh in fc['sources'].items():assert V.sha(path)==hsh
 assert fc['status']=='PASS_EXACT_FINITE_PROFILE_ALL_HEIGHT_FRONTIER_OBSTRUCTION'
 proposals=V.read(BASE+'_dihedral_scout.json')['rows']
 G,I,h,g,old,names,center,pairs,_=N.old_state()
 caps=[N.make(s,b,g,old,names,center,pairs)[0]for s,b in zip(N.SPECS,d['branches'],strict=True)]
 for spec,new in zip(N.SPECS,N.NEW,strict=True):
  ff=dict(zip(spec['edge'],spec['labels'],strict=True));ff[spec['anchor']]=spec['anchor_label']
  for n in names:pairs['|'.join(sorted((new,n))) ]={'kind':'finite','label':ff[n]}if n in ff else{'kind':'dotted'}
 pairs['|'.join(sorted(N.NEW))]={'kind':'dotted'};old+=caps;names+=N.NEW
 target=tuple(map(K.parse,d['escape']));yi=list(map(V.golden,target));w=list(map(V.golden,[K.neg(x)for x in K.mv(I,target)]))
 edges=[r['edge']for r in cert['new_vertex_frontier']if len(r['end_walls'])==1]
 expected=[];rows=[]
 for ei,ss in enumerate(edges):
  M=N.C.gram(ss,pairs);Mi=K.inverse(M);assert N.T.P.T.psd(M)=='positive_definite'
  isolated=[j for j in range(8)if all(M[j][k]==K.Z for k in range(8)if j!=k)]
  expected.extend((ei,ss[j])for j in isolated)
  A=[old[names.index(n)]for n in ss];mi=[[V.golden(x)for x in row]for row in Mi]
  def proj(v):
   z=N.combine(V.mv(mi,[V.pair(v,g,r)for r in A]),A)
   return[V.sub(a,b)for a,b in zip(v,z,strict=True)]
  tv=proj(center);qt=V.pair(tv,g,tv);assert qt[1]<0;tv=[V.div(x,V.sqrt(V.neg(qt)))for x in tv]
  ev=proj(w);lam=V.pair(ev,g,tv);ev=[V.add(a,V.mul(lam,b))for a,b in zip(ev,tv,strict=True)];qe=V.pair(ev,g,ev);assert qe[0]>0;ev=[V.div(x,V.sqrt(qe))for x in ev]
  # Orthogonal normalized timelike/spacelike vectors span E^perp exactly.
  target_linear=[V.dot(tv,yi),V.dot(ev,yi)]
  for j in isolated:
   prop=[r for r in proposals if(r['edge'],r['wall'])==(ei,ss[j])];assert len(prop)==1
   support=[n for n,_ in prop[0]['dual_support']];assert len(support)==len(set(support))==2
   assert all(n in names and n not in ss for n in support)
   base=[V.div(V.neg(x),V.rat(2))for x in A[j]]
   rr=[old[names.index(n)]for n in support]
   mat=[[V.pair(tv,g,r),V.pair(ev,g,r)]for r in rr]
   a,b=mat[0];c,e=mat[1];x,y=target_linear
   det=V.sub(V.mul(a,e),V.mul(b,c));assert det[0]>0 or det[1]<0
   # mat^T lambda = -target_linear; this exact inverse formula proves
   # cancellation of both free coordinates, not merely intervals about0.
   l1=V.div(V.sub(V.mul(c,y),V.mul(e,x)),det)
   l2=V.div(V.sub(V.mul(b,x),V.mul(a,y)),det)
   assert l1[0]>0 and l2[0]>0
   constants=[V.pair(base,g,r)for r in rr];constant=V.add(V.dot(base,yi),V.add(V.mul(l1,constants[0]),V.mul(l2,constants[1])))
   assert constant[0]>0,(ei,ss[j],support,constant)
   rows.append({'edge':ei,'isolated_wall':ss[j],'support':support,'dual_weights':[V.enc(l1),V.enc(l2)],
                'nonzero_determinant':V.enc(det),'uniform_target_pairing_lower_bound':V.enc(constant)})
 assert len(expected)==len(set(expected))==51
 assert set(expected)=={(r['edge'],r['wall'])for r in fc['uncovered_I2_tails']}
 assert set(expected)=={(r['edge'],r['wall'])for r in proposals}
 out={'status':'PASS_EXACT_ALL_REAL_ISOLATED_DIHEDRAL_FRONTIER_OBSTRUCTION',
      'sources':{p:V.sha(p)for p in(BASE+'_caps.json',BASE+'_verify.json',BASE+'_continuous_verify.json',BASE+'_dihedral_scout.json')},
      'interval_scale':str(V.SCALE),'tail_cases':51,'certificates':rows,
      'identity':'For f_j=(r,old_j)<=0 and r=-a/2+x*v+y*w, target_pairing = positive_constant - lambda1*f1 - lambda2*f2, with both lambdas positive.',
      'proof':'Scale any candidate with its sole nonzero edge pairing -mu, mu>0, by1/mu. The isolated edge normal a has norm2 and is orthogonal to the other7 normals, hence r=-a/2+z with z in their two-dimensional Lorentz complement. The exact two-wall dual identity bounds its target pairing strictly positively for ALL real z; no norm or angle restriction is needed.',
      'combined_scope':'Together with the independently replayed finite-profile theorem, no acute-compatible new wall retaining all104 walls can both cut this exact escape and supply a finite endpoint on any of these16 open spherical8 edges, for ANY integer Coxeter angle labels. Other edges/components, wall replacement, caps not creating such an endpoint, and multiple noncutting caps are not excluded. No global compactH9 nonexistence claim.'}
 if verify:assert out==V.read(OUT)
 else:(V.ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],len(rows),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();run(args.verify)
