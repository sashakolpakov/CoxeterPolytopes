#!/usr/bin/env python3
"""Archive the bounded numerical test of the already saved opposite escape."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts.fresh_h9_face_bridge_20261013_scout import setup,PREFIX,V,K

def run():
 Q,_,_,_,_,_,_=setup();bank=V.read(PREFIX+'_scout.json');d=V.read(PREFIX+'_caps.json')
 R=np.array([r['normal']for r in bank['candidates']]);y=np.array([sum(V.golden(K.parse(x)))/(2*V.SCALE)for x in d['escape']]);chosen=R[[0,12219]]
 pair=R@Q@chosen.T;cut=-R@y;cutters=cut>1e-7;rows={}
 for parallel in(False,True):
  mask=cutters.copy()
  for j in range(2):
   finite=(abs(pair[:,j])<1e-6)|np.logical_or.reduce([abs(pair[:,j]+2*np.cos(np.pi/m))<1e-6 for m in(3,4,5,6,8,10,12)])
   ok=(pair[:,j]<-2-1e-6)|finite
   if parallel:ok|=abs(pair[:,j]+2)<1e-6
   mask&=ok
  rows['with_parallel'if parallel else'without_parallel']=[int(i)for i in np.flatnonzero(mask)]
 out={'status':'NUMERICAL_BOUNDED_OPPOSITE_ESCAPE_GATE','source_sha256':{p:V.sha(p)for p in(PREFIX+'_scout.json',PREFIX+'_caps.json')},'bank_entries':len(R),'old_selected_indices':[0,12219],'target_cut_proposals':int(cutters.sum()),'compatible_cutters':rows,'scope':'Floating test of the fixed saved endpoint bank only. No exact all-cap obstruction or wider-family claim.'}
 (V.ROOT/(PREFIX+'_opposite_scout.json')).write_text(json.dumps(out,indent=2)+'\n');print(out)

if __name__=='__main__':run()
