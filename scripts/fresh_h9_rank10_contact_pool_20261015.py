#!/usr/bin/env python3
"""Exact bounded contact pool tested against the28-wall chamber."""
import argparse,json
import fresh_h9_quartet_polar_second_20261014 as B
A,D,R,C,ROOT=B.A,B.D,B.P,B.C,B.ROOT
TAG='fresh_h9_rank10_contact_pool_20261015'
BASE='results/fresh_h9_replacement_escape_extension_20261015.json'
RECIPES='data/enumerations/fresh_h9_rank10_contact_recipes_20261015.json'
def run():
 d=B.read(BASE);A.RAD=A.A.Q(*d['field']['new_positive_radical_square'])
 Q,N,H,E=D.dm(d['form']),D.dm(d['normals']),D.dv(d['center']),D.dv(d['strict_escape']['point'])
 cov=[A.mv(Q,n)for n in N];norms=[A.pair(n,n,Q)for n in N]
 phi=A.L(A.A.Q(A.F(1,2),A.F(1,2)));vals={2:A.Z,3:-A.ONE,5:-phi};out=[]
 for r in json.loads((ROOT/RECIPES).read_text()):
  cs=dict(r['contacts']);basis=[]
  for i in cs:
   assert norms[i]==2
   if C._row_rank([cov[j]for j in basis+[i]],10)>len(basis):basis.append(i)
   if len(basis)==10:break
  assert len(basis)==10
  n=A.mv(A.inv([cov[j]for j in basis]),[vals[cs[j]]for j in basis]);q=A.pair(n,n,Q)
  assert q==2 and all(A.pair(n,N[i],Q)==vals[m]for i,m in cs.items())
  old=next((i for i,v in enumerate(N)if n==v),None)
  relations=[('existing',None,None)if i==old else ('forbidden',None,None)if A.pair(n,v,Q).sign()>0 else R.relation(A.pair(n,v,Q),q,norms[i])for i,v in enumerate(N)]
  bad=[i for i,x in enumerate(relations)if x[0]=='forbidden']
  hp,ep=A.pair(n,H,Q),A.pair(n,E,Q);assert hp.sign()<0
  out.append(dict(id=r['id'],contacts=r['contacts'],basis=basis,normal=n,norm=q,
                  existing_wall=old,forbidden_indices=bad,relations=relations,
                  center_pairing=hp,escape_pairing=ep,
                  compatible_cutter=old is None and not bad and ep.sign()>0))
  print('recipe',r['id'],'existing',old,'bad',bad,'escape_sign',ep.sign(),flush=True)
 assert len(out)==33 and not any(r['compatible_cutter']for r in out)
 paths=[BASE,RECIPES,f'scripts/{TAG}.py','scripts/fresh_h9_quartet_polar_second_20261014.py','scripts/fresh_h9_quartet_polar_20261014.py']
 return A.enc(dict(status='PASS_BOUNDED_33_EXACT_CONTACT_RECIPES',sources={p:R.sha(p)for p in paths},
  rows=out,recipe_count=len(out),existing_wall_count=sum(r['existing_wall']is not None for r in out),
  compatible_new_count=sum(r['existing_wall']is None and not r['forbidden_indices']for r in out),
  compatible_escape_cutter_count=0,
  scope='Only33explicit rank-ten finite-contact recipes drawn from bounded earlier scouts. Each normal is reconstructed exactly from linear contacts, then all28 Coxeter relations and the current strict-escape sign are tested. This is not a complete wall census or an all-label nonexistence result.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['compatible_new_count'],d['compatible_escape_cutter_count'])
