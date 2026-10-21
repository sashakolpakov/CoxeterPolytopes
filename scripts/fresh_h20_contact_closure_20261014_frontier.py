#!/usr/bin/env python3
"""Exact local20-edge stars of the16 specified pi/3 bridge vertices."""
import argparse,json
from collections import Counter
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
import fresh_h20_contact_closure_20261014 as p
K=p.K;Z=K(0);ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_contact_closure_20261014_frontier';SOURCE='results/fresh_h20_contact_closure_20261014.json'
def sha(f):return sha256((ROOT/f).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/SOURCE).read_text())
 for f,h in d['sources'].items():assert sha(f)==h
 Q=[[K.from_json(x)for x in r]for r in d['form']];N=[list(map(K.from_json,r))for r in d['normals']];cov=[p.b.matvec(Q,r)for r in N];G=[[sum((a*b for a,b in zip(x,y)),Z)for y in cov]for x in N];norms=[G[i][i]for i in range(76)]
 dotted={(i,j)for i in range(76)for j in range(76)if G[i][j]<0 and G[i][j]*G[i][j]>norms[i]*norms[j]}
 counts=Counter();rows=[];ninth=0
 for vi,vertex in enumerate(d['finite_vertices']):
  inc=set(vertex['incident'])
  for released in sorted(inc):
   S=sorted(inc-{released});B=[[G[i][j]for j in S]for i in S];candidates=[a for a in range(76)if a not in inc and not any((a,i)in dotted for i in S)];hits=[]
   if candidates:
    I=p.b.inverse(B)
    for a in candidates:
     q=[G[a][i]for i in S];rho=G[a][a]-sum((x*y for x,y in zip(q,p.b.matvec(I,q))),Z);ninth+=1
     if rho>=0:hits.append([a,rho])
   if hits:
    assert len(hits)==1 and hits[0][1]>0
    dest=[j for j,w in enumerate(d['finite_vertices'])if set(w['incident'])==set(S+[hits[0][0]])];assert len(dest)==1 and dest[0]!=vi
    kind='finite';counts['finite_germs']+=1
   else:kind='open';dest=[];counts['open_germs']+=1
   rows.append(dict(vertex=vi,released=released,support=S,kind=kind,other_vertex=dest,nondotted_candidates=candidates,psd_extensions=hits))
 assert counts=={'finite_germs':16,'open_germs':304}
 return p.enc(dict(status='PASS_16_SPECIFIED_CORNER_STARS_8_FINITE_EDGES_304_OPEN_GERMS',sources={f:sha(f)for f in[SOURCE,f'scripts/{TAG}.py','scripts/fresh_h20_contact_closure_20261014.py']},counts=dict(counts),exact_schur_tests=ninth,edge_germs=rows,scope='Complete local stars at the16 specified actual vertices only. Each open ray reaches the ideal boundary with no other wall active there, so it has no second finite or bounded-link ideal endpoint. No all-chamber vertex/component census.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();f=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(f.read_text())
 else:f.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['exact_schur_tests'])
