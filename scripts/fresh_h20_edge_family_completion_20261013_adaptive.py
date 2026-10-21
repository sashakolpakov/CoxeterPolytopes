#!/usr/bin/env python3
"""Bounded adaptive mixed-bank completion scout; every epoch is frozen."""
import argparse,json
from fresh_h20_edge_family_completion_20261013_bank import ROOT,PREFIX,run as bank
from fresh_h20_edge_family_completion_20261013_scout import run as scout

def run(epochs):
    extras=['data/enumerations/'+PREFIX+'_scout0.json'];log=[]
    for i in range(1,epochs+1):
        bp='data/enumerations/'+PREFIX+f'_bank{i}.json';sp='data/enumerations/'+PREFIX+f'_scout{i}.json'
        b=bank(4,extras,bp);g=b['whole_bank_geometry']
        item={'epoch':i,'bank':bp,'entries':len(b['caps']),'whole_bank_status':g['status'],'whole_bank_norm':g.get('norm'),'selected_walls':[x['walls']for x in b['selections']]};log.append(item)
        (ROOT/('data/enumerations/'+PREFIX+'_adaptive.json')).write_text(json.dumps({'status':'NUMERICAL_BOUNDED_ADAPTIVE_MIXED_BANK','epoch_limit':epochs,'epochs':log},indent=2)+'\n')
        if g['status']!='NUMERICAL_SPACELIKE_ESCAPE':break
        try:s=scout(bp,sp)
        except AssertionError as e:item['stopped']=str(e);break
        item['new_scout']=sp;item['new_entries']=len(s['caps'])
        if not s['caps']:break
        extras.append(sp)
    (ROOT/('data/enumerations/'+PREFIX+'_adaptive.json')).write_text(json.dumps({'status':'NUMERICAL_BOUNDED_ADAPTIVE_MIXED_BANK','epoch_limit':epochs,'epochs':log},indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--epochs',type=int,default=8);a=p.parse_args();run(a.epochs)
