#!/usr/bin/env python3
"""Exact identification of joint proposals against the old1084 mirrors and peers."""
import argparse,json
from collections import defaultdict
from scripts.fresh_h9_joint_completion_20261013 import setup,ROOT,read
from scripts.fresh_h8_threefree_verify_20261013 import (
 node,sq,qsign,qs,bsign,Z,O,P,add,sub,neg,mul,div,scale)
OUTPUT=ROOT/'data/enumerations/fresh_h9_joint_completion_20261013_mirrors.json'
def run(verify=False):
 G,I,center,old,*_=setup();data=read('data/enumerations/fresh_h9_joint_completion_20261013.json');new=[node(c,G)for c in data['caps']];groups=defaultdict(list);squares=[];ds=[];origins=[];tests=0
 def key(n):
  vals=[tuple(mul(div(scale(O,2),n['q']),t)for t in sq(a,b,n['d']))for a,b in zip(n['pa'],n['pb'])];dotted=tuple(j for j,v in enumerate(vals)if qsign(qs(v,(scale(O,4),Z)),n['d'])>0)
  if len(dotted)!=3:return None,vals
  fixed=[]
  for j,v in enumerate(vals):
   if j in dotted:continue
   k=next((k for k,z in enumerate((Z,O,add(O,P)))if qsign(qs(v,(z,Z)),n['d'])==0),None)
   if k is None:return None,vals
   fixed.append(k)
  return(dotted,tuple(fixed)),vals
 for i,n in enumerate(old):
  k,v=key(n)
  if k is None:continue
  groups[k].append(len(squares));squares.append(v);ds.append(n['d']);origins.append({'dictionary':'old','index':i})
 aliases=[];novel=[]
 for i,n in enumerate(new):
  k,v=key(n);assert k is not None;alias=None
  for j in groups[k]:
   tests+=1
   if all(bsign((sub(a[0],b[0]),a[1],neg(b[1]),Z),n['d'],ds[j])==0 for a,b in zip(v,squares[j])):alias=origins[j];break
  if alias is not None:aliases.append({'candidate':i,**alias});continue
  novel.append(i);groups[k].append(len(squares));squares.append(v);ds.append(n['d']);origins.append({'dictionary':'joint','index':i})
 out={'status':'EXACT_JOINT_MIRROR_IDENTIFICATION','source':'data/enumerations/fresh_h9_joint_completion_20261013.json','old_dictionary_size':len(old),'proposal_entries':len(new),'distinct_new_mirrors':len(novel),'new_representative_indices':novel,'aliases':aliases,'exact_same_pattern_pair_checks':tests,'scope':'All normals are square-normalized through pairing squares; nonpositive seed pairings and independent seed normals make equality of all normalized squares equivalent to equality of oriented mirrors.'}
 if verify:assert read(str(OUTPUT.relative_to(ROOT)))==out
 else:OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
 print({k:v for k,v in out.items()if k not in('new_representative_indices','aliases')},flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
