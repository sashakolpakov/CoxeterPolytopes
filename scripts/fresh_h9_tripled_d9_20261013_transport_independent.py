#!/usr/bin/env python3
"""Separate Fraction replay of the tripled-D9 chamber transport.

No producer/K-class implementation is imported. The old independent field
engine is reused with the NEW form <-3(T-1)>+I9, and the new matrix is derived
again from the two new oriented corners, never from the old-host matrix.
"""
import argparse,json
import fresh_h9_d9_symmetry_transfer_20261013 as D
a=D.a;ROOT=D.ROOT;PREFIX='fresh_h9_tripled_d9_20261013'
D.G=tuple(tuple(a.neg(a.scale(a.ALPHA,3))if i==j==0 else a.O if i==j else a.Z for j in range(10))for i in range(10))
OUT=ROOT/'results'/f'{PREFIX}_transport_independent.json'

def primitive(r):
 assert D.lattice(r)
 for p in(a.T,a.scale(a.O,3)):assert not D.lattice(tuple(a.divide(x,p)for x in r))
 q=D.pair(r,r);assert all(a.sign(q,e)>0 for e in range(4))
 assert all(a.integral(a.divide(a.scale(D.pair(r,b),2),q))for b in D.BASIS)
D.primitive_reflective=primitive

def run(verify=False):
 paths=[f'results/{PREFIX}_{s}.json'for s in('edge1','step2','step3','edge3','symmetry','independent','chain_independent')]
 sources=D.checked_sources(paths+['scripts/fresh_h9_glued_host_20261013_independent.py','scripts/fresh_h9_d9_symmetry_transfer_20261013.py'])
 first,second,third,new,cert,ia,ca=[json.loads((ROOT/p).read_text())for p in paths]
 assert ia['status']=='PASS_INDEPENDENT_TRIPLED_D9_COMPLETE_INITIAL_EDGES'
 assert ca['status']=='PASS_INDEPENDENT_TRIPLED_D9_CHAIN_EDGES'
 assert (ia['total_profiles'],ia['total_interval_points'],ca['profiles'],ca['interval_points'])==(68,1121,68,40)
 assert ia['norm_classes']==26
 assert cert['status']=='PASS_EXACT_TRIPLED_D9_INFINITE_CHAMBER_TRANSPORT'
 assert a.dec(cert['coefficient'])==a.scale(a.ALPHA,3)
 origin=(a.O,*(a.Z for _ in range(9)));frames=[D.INITIAL];points=[origin]
 for d in(first,second,third):
  frames.append(D.matrix(d['endpoint_incident']));points.append(a.vec(d['endpoint']))
 newroots=D.matrix(new['endpoint_incident']);newpoint=a.vec(new['endpoint'])
 steps=[]
 for i,d in enumerate((first,second,third)):
  steps.append(D.path_step(frames[i],points[i],frames[i+1],points[i+1],a.vec(d['chosen_root']),f'chain {i+1}'))
 steps.append(D.path_step(D.INITIAL,origin,newroots,newpoint,a.vec(new['chosen_root']),'other initial edge'))
 # All four edge endpoints lie in the same retained finite prefix. The
 # mirror-free open segments are the explicitly hash-bound first-edge proof.
 walls=tuple(dict.fromkeys((*D.matrix(third['all_walls']),*D.matrix(new['all_walls']))))
 for p in(*points,newpoint):assert all(a.sign(D.pair(p,r))<=0 for r in walls)
 perm=tuple(cert['permutation']);assert perm==(8,1,2,3,4,5,6,7,0)
 oldroots=frames[-1];oldpoint=points[-1];scale=a.dec(cert['positive_point_scale']);assert a.sign(scale)>0
 assert D.pair(oldpoint,oldpoint)==a.mul(a.sq(scale),D.pair(newpoint,newpoint))
 C,pivots=D.simple_system(oldroots,oldpoint);CN,pivotsn=D.simple_system(newroots,newpoint)
 assert all(C[i][j]==CN[perm[i]][perm[j]]for i in range(9)for j in range(9))
 F=D.mm(D.transpose((*(newroots[i]for i in perm),D.vmul(scale,newpoint))),a.minverse(D.transpose((*oldroots,oldpoint))))
 FI=a.minverse(F)
 assert F==D.matrix(cert['matrix'])and FI==D.matrix(cert['inverse'])
 assert D.mm(F,FI)==D.mm(FI,F)==D.IDENT and D.mm(D.mm(D.transpose(F),D.G),F)==D.G
 assert all(a.mv(F,r)==newroots[perm[i]]for i,r in enumerate(oldroots))
 assert a.mv(F,oldpoint)==D.vmul(scale,newpoint)
 assert a.sign(a.mv(F,origin)[0])>0
 module=D.mm(D.mm(D.BI,F),D.B);inverse_module=D.mm(D.mm(D.BI,FI),D.B)
 assert module==D.matrix(cert['module_matrix'])and inverse_module==D.matrix(cert['module_inverse'])
 assert all(a.integral(x)for M in(module,inverse_module)for row in M for x in row)
 assert all(D.lattice(a.mv(M,b))for M in(F,FI)for b in D.BASIS)
 trace=D.tr(F);assert trace==a.dec(cert['trace'])and a.sign(a.sub(trace,a.scale(a.O,10)))>0
 assert D.mirrors(oldroots)==D.mirrors(newroots)==45
 assert cert['finite_mirror_counts']==[45,45]
 forged=[list(row)for row in F];forged[0][0]=a.add(forged[0][0],a.O)
 assert D.mm(D.mm(D.transpose(forged),D.G),forged)!=D.G
 out={'status':'PASS_INDEPENDENT_FRACTION_TRIPLED_D9_GLOBAL_TRANSPORT','sources':sources,
      'host':'<-3(T-1)>+D9, T^4-4T^2+2=0, T=2cos(pi/8)',
      'arithmetic':'Separate Fraction polynomial arithmetic; no producer field-class imports.',
      'required_complete_edges':4,'path':steps,'all_path_prefix_walls':len(walls),
      'total_profiles':136,'total_interval_points':1161,'safe_norm_classes':26,
      'corner_mirrors':[45,45],'corner_pivots':[list(map(a.enc,pivots)),list(map(a.enc,pivotsn))],
      'matrix':D.encmat(F),'inverse':D.encmat(FI),'module_matrix':D.encmat(module),'module_inverse':D.encmat(inverse_module),
      'positive_corner_scale':a.enc(scale),'trace':a.enc(trace),'trace_exceeds10':True,'forged_matrix_rejected':True,
      'scope':'The separately certified complete first edges connect the two complete oriented corners in the same full lattice reflection chamber. The independently reconstructed forward/inverse D9 module isometry transports that sector and is future preserving. It therefore stabilizes the full chamber. Trace>10 gives infinite order, excluding a finite-volume chamber. This is specific to this new tripled host, not a commensurability transfer.'}
 if verify:assert out==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['trace'],flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a0=p.parse_args();run(a0.verify)
