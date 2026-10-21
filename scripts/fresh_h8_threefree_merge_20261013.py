#!/usr/bin/env python3
"""Merge adaptive accepted cap proposals, preserving exact witness provenance."""
import json
from copy import deepcopy
from scripts.fresh_h8_alternate_caps_20261013 import *
from scripts.fresh_h8_release_h9_free_caps_20261013 import esign,eencode

def run():
 files=['data/enumerations/fresh_h8_threefree_caps_20261013.json','data/enumerations/fresh_h8_threefree_adaptive_caps_20261013.json'];data=[json.loads((ROOT/f).read_text())for f in files];proof=json.loads((ROOT/data[0]['witness_source']).read_text());rows=[{'source':f'old-core-{i}','certificate':{'y':r['escape_y'],'norm':r['escape_norm'],'center_pairing':r['escape_center_pairing']}}for i,r in enumerate(proof['five_complete_escape_clauses'])];rows+=json.loads((ROOT/data[1]['witness_source']).read_text())['clauses'];wsource='data/enumerations/fresh_h8_threefree_combined_targets_20261013.json';(ROOT/wsource).write_text(json.dumps({'clauses':rows},indent=2)+'\n');ys=[tuple(map(parse,r['certificate']['y']))for r in rows];caps=[];seen=set()
 for source,raw in zip(files,data):
  for j,c0 in enumerate(raw['caps']):
   c=deepcopy(c0);key=json.dumps([c['discriminant'],c['coefficients_rational'],c['coefficients_radical']],sort_keys=True)
   if key in seen:continue
   seen.add(key);a=tuple(map(parse,c['coefficients_rational']));b=tuple(map(parse,c['coefficients_radical']));d=parse(c['discriminant']);vals=[(dot(a,y),dot(b,y))for y in ys];c['cut_core_witness_indices']=[i for i,v in enumerate(vals)if esign(v,d)<0];c['core_constraint_values']=[eencode(v)for v in vals];c['original_proposal_source']=source;c['original_proposal_index']=j;caps.append(c)
 out=deepcopy(data[0]);out['caps']=caps;out['witness_source']=wsource;out['target_core_witnesses']=list(range(len(ys)));out['merged_sources']=files;out['proposal_counts']=[d['proposal_counts']for d in data];(ROOT/'data/enumerations/fresh_h8_threefree_combined_caps_20261013.json').write_text(json.dumps(out,indent=2)+'\n');print({'caps':len(caps),'additional':len(caps)-len(data[0]['caps'])})
if __name__=='__main__':run()
