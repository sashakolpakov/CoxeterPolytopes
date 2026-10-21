#!/usr/bin/env python3
"""All-pair exact mirror equality audit, filtered only by exact dotted support."""
import argparse,json
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from scripts.fresh_h8_alternate_core_angles_independent_20261013 import parse,qm,bsign,sub,neg,Z
ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_h8_threefree_combined_caps_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h8_threefree_mirror_aliases_20261013.json'
def run(verify=False):
 caps=json.loads((ROOT/SOURCE).read_text())['caps'];groups=defaultdict(list);sq=[];ds=[]
 for i,c in enumerate(caps):
  r=parse(c['discriminant']);ds.append(r);a=list(map(parse,c['pairings_rational']));b=list(map(parse,c['pairings_radical']));sq.append([qm((x,y),(x,y),r)for x,y in zip(a,b)]);groups[tuple(sorted(c['free_seed_indices']))].append(i)
 alias=[];tests=0
 for g in groups.values():
  for i,j in combinations(g,2):
   tests+=1
   if all(bsign((sub(a[0],b[0]),a[1],neg(b[1]),Z),ds[i],ds[j])==0 for a,b in zip(sq[i],sq[j])):alias.append([i,j])
 leaders=list(range(len(caps)))
 for i,j in alias:
  a,b=leaders[i],leaders[j];leaders=[a if t==b else t for t in leaders]
 out={'source':SOURCE,'equal_dotted_support_pair_checks':tests,'duplicate_mirror_pairs':alias,'distinct_new_mirrors':len(set(leaders)),'scope':'All norms equal2 and all seed pairings nonpositive. Equality of their ten pairing squares therefore exactly characterizes equality of oriented mirrors. Different dotted supports cannot coincide.'}
 if verify:assert json.loads(OUTPUT.read_text())==out
 else:OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
