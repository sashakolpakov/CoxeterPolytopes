#!/usr/bin/env python3
"""Reject invalid volume certificates for the new finite-face construction."""
import copy,json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_compact_next_20261013_verify as T
V=T.V;K=T.K

def run():
    d=V.read(T.PREFIX+'_caps.json');T.check(d);_,_,h,*_=T.reconstruct(d)
    corruptions=[]
    a=copy.deepcopy(d);a['escape']=list(map(K.enc,h));corruptions.append(('timelike replacement escape',a))
    a=copy.deepcopy(d);a['escape']=d['target'];corruptions.append(('previous escape cut by compact0 claimed surviving',a))
    out={'status':'PASS','source':T.PREFIX+'_caps.json','source_sha256':V.sha(T.PREFIX+'_caps.json'),'rejected':[]}
    for label,bad in corruptions:
        try:T.check(bad)
        except AssertionError:out['rejected'].append(label)
        else:raise AssertionError(label+' accepted')
    (V.ROOT/'results/fresh_h9_compact_next_20261013_regressions.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))

if __name__=='__main__':run()
