#!/usr/bin/env python3
"""Reject the wrong quadratic branch and a previously cut escape."""
import copy,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_face_bridge_20261013_verify as T

def run():
 d=T.V.read(T.PREFIX+'_caps.json');tests=[]
 for name,mutate in [('opposite cap branch',lambda x:x['caps'][0].update(branch=1)),('old escape cut by faceGap',lambda x:x.update(escape=x['target']))]:
  bad=copy.deepcopy(d);mutate(bad)
  try:T.check(bad)
  except AssertionError:tests.append({'case':name,'rejected':True})
  else:raise AssertionError('forgery accepted: '+name)
 out={'status':'PASS_FACE_BRIDGE_REGRESSIONS','source':T.PREFIX+'_caps.json','source_sha256':T.V.sha(T.PREFIX+'_caps.json'),'tests':tests}
 (T.V.ROOT/'results/fresh_h9_face_bridge_20261013_regressions.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'])

if __name__=='__main__':run()
