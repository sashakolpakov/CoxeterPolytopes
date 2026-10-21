#!/usr/bin/env python3
"""Reject the already-cut old escape and a forged source binding."""
import copy,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_prism_completion_20261013_verify as T

def run():
    d=T.V.read(T.PREFIX+'_caps.json');rows=[]
    for name in('forged_source_hash','old_target_as_remaining_escape'):
        bad=copy.deepcopy(d)
        if name=='forged_source_hash':bad['source_sha256']='0'*64
        else:bad['escape']=bad['target']
        try:T.check(bad)
        except AssertionError:rows.append({'mutation':name,'rejected':True})
        else:raise AssertionError('mutation accepted: '+name)
    out={'status':'PASS_REJECTION_TESTS','source_sha256':T.V.sha(T.PREFIX+'_caps.json'),'tests':rows}
    (T.V.ROOT/'results/fresh_h9_prism_completion_20261013_regressions.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'])

if __name__=='__main__':run()
