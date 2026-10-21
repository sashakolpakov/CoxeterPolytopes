#!/usr/bin/env python3
"""Finite mixed-bank growth with an18-short/17-short direction fallback."""
import argparse,json
from fresh_h20_edge_family_completion_20261013_bank import ROOT,PREFIX,read,run as bank
from fresh_h20_edge_family_completion_20261013_scout import run as scout18
from fresh_h20_edge_family_completion_20261013_three_angles import run as scout17

def run(epochs):
    current='data/enumerations/'+PREFIX+'_bank5.json';extras=[p for p in read(current)['sources']if not p.endswith('_vertex_bank.json')and not p.endswith('_edge_directions.json')]
    log=[];out='data/enumerations/'+PREFIX+'_face_walk.json'
    for j in range(6,6+epochs):
        step={'target_bank':current};log.append(step)
        for k,scout in[(18,scout18),(17,scout17)]:
            sp='data/enumerations/'+PREFIX+f'_pass{j}_{k}.json'
            try:s=scout(current,sp)
            except AssertionError as e:step[str(k)+'_scope_limit']=str(e);continue
            step[str(k)+'_scout']=sp;step[str(k)+'_entries']=len(s['caps'])
            if s['caps']:extras.append(sp);break
        else:step['stopped']='No targeted cutter in either scoped family';break
        bp='data/enumerations/'+PREFIX+f'_bank{j}.json';b=bank(4,extras,bp);g=b['whole_bank_geometry']
        step.update({'bank':bp,'entries':len(b['caps']),'whole_bank_geometry':g,'selected_walls':[x['walls']for x in b['selections']]});current=bp
        (ROOT/out).write_text(json.dumps({'status':'NUMERICAL_BOUNDED_OLD_FACE_DIRECTION_WALK','epoch_limit':epochs,'steps':log},indent=2)+'\n')
        if g['status']!='NUMERICAL_SPACELIKE_ESCAPE':break
        if len(b['caps'])>5000:step['stopped']='Explicit5000-entry bank-size cap';break
    (ROOT/out).write_text(json.dumps({'status':'NUMERICAL_BOUNDED_OLD_FACE_DIRECTION_WALK','epoch_limit':epochs,'steps':log},indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--epochs',type=int,default=12);a=p.parse_args();run(a.epochs)
