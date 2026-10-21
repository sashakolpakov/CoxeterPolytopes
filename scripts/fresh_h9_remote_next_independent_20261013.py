#!/usr/bin/env python3
"""Source audit and100-digit replay, reusing the remote intervalengine."""
import argparse,copy,json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts import fresh_h9_remote_next_20261013_verify as N
V=N.V
OUT='results/fresh_h9_remote_next_independent_20261013.json'

def strip(x,path=()):
 if isinstance(x,list):
  if len(x)==2 and all(isinstance(z,str)for z in x):
   try:a,b=map(int,x)
   except ValueError:pass
   else:
    assert a<=b;return 'INTERVAL',{path:(a,b)}
  rows=[];parts={}
  for i,y in enumerate(x):v,p=strip(y,path+(i,));rows.append(v);parts.update(p)
  return rows,parts
 if isinstance(x,dict):
  rows={};parts={}
  for k,y in x.items():v,p=strip(y,path+(k,));rows[k]=v;parts.update(p)
  return rows,parts
 return x,{}

def run():
 data=V.read(N.PREFIX+'_caps.json');old=V.read(N.PREFIX+'_verify.json')
 V.SCALE=10**100;V.PHI=V.div(V.add(V.rat(1),V.sqrt(V.rat(5))),V.rat(2));V.radical.cache_clear()
 new=N.check(data);new['artifact_sha256']=V.sha(N.PREFIX+'_caps.json')
 om,oi=strip(old);nm,ni=strip(new);assert om==nm and oi.keys()==ni.keys()
 for p,(a,b)in ni.items():lo,hi=oi[p];assert lo*10**20<=a<=b<=hi*10**20,p
 return {'status':'PASS_REMOTE86_H4_SOURCE_REVIEW_AND_100_DIGIT_REPLAY',
  'source':N.PREFIX+'_caps.json','source_sha256':V.sha(N.PREFIX+'_caps.json'),
  'checked_enclosures':len(ni),'new_pair_counts':new['new_pair_counts'],
  'actual_finite_vertex':new['actual_vertex']['type'],'closed_old_edge':new['closed_old_edge'],
  'strict_escape_norm':new['remaining_escape']['norm'],
  'scope':'Sourceproofreview andhigherprecisionreplay ofsameoutwardintervalimplementation,notaseparatearithmeticimplementation. Norm/frame/threeanchoridentities,all76strictdottedrelations,actualH4endpoint,oldendpointpreservation,andstrict86escapePASS. Noexhaustiveorfinitevolumeclaim.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==V.read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out))
