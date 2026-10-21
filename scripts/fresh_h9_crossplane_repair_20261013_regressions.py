#!/usr/bin/env python3
"""Meaningful false-certificate rejection checks for the two exact checkpoints."""
import copy
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_crossplane_repair_20261013_pair_verify as P
V=P.V;K=P.K

def run():
    d=V.read(P.PREFIX+'_caps.json');P.check(d)
    G,I,h,*_=P.reconstruct(d)
    corruptions=[]
    a=copy.deepcopy(d);a['escape']=list(map(K.enc,h));corruptions.append(('timelike point substituted for escape',a))
    a=copy.deepcopy(d);a['escape']=d['targets'][0];corruptions.append(('old escape cut by both new walls substituted as surviving',a))
    a=copy.deepcopy(d);a['branch']=-d['branch'];corruptions.append(('wrong norm-equation branch',a))
    out={'status':'PASS','source':P.PREFIX+'_caps.json','source_sha256':V.sha(P.PREFIX+'_caps.json'),'rejected':[]}
    for name,bad in corruptions:
        try:P.check(bad)
        except AssertionError:out['rejected'].append(name)
        else:raise AssertionError('Accepted false certificate: '+name)
    path='results/fresh_h9_crossplane_repair_20261013_regressions.json'
    (V.ROOT/path).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))

if __name__=='__main__':run()
