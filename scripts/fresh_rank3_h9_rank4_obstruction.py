#!/usr/bin/env python3
"""Bound-independent obstruction to completing the two-weight slice pair.

Fraction-only certificate producer/replay; no project arithmetic imports.
"""
import argparse
import json
from pathlib import Path
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,sub,scale,mul,div,conj,sign,integral,total,dot,mv,tr,mm,pair,matrix,eye,inverse,
)
from scripts.fresh_multibasepoint_h9_verify import ROOT,read
from scripts.fresh_multibasepoint_h9_targeted_verify import escape

OUTPUT=ROOT/'results/fresh_rank3_h9_rank4_obstruction_20261013.json'
def projected(w):return mul(scale(P,-2),w[0]),w[5],w[8],w[9]

def verify(payload):
    base=read(payload['base_source']);completion=read(payload['completion_source']);dictionary=read(completion['dictionary'])
    gram=matrix(base['gram']);old=matrix(base['cases'][0]['all_walls']);U=matrix(base['standard_block_swap_involution'])
    assert len(old)==382 and mm(U,U)==eye(10) and mm(tr(U),mm(gram,U))==gram
    assert all(integral(x) for row in U for x in row) and set(mv(U,r) for r in old)==set(old)
    row=completion['cases'][payload['completion_case']];assert payload['completion_case']==1 and row['selection_method']=='decreasing_height'
    roots=matrix(dictionary['roots']);selected=tuple(roots[i] for i in row['selected_root_indices'])
    v=escape(row['hyperideal_witness'],(*old,*selected),gram)
    assert v==tuple(map(parse,payload['escape'])) and pair(v,gram,v)==parse(payload['escape_norm'])
    assert pair(v,gram,v)==parse({'a':'901480/358801','b':'-2767458/1794005'})
    inside=tuple(map(parse,dictionary['strict_interior']))
    assert sign(pair(inside,gram,inside))<0 and sign(inside[0])>0
    assert all(sign(pair(inside,gram,r))<0 for r in (*old,*roots))
    for record,target in zip(payload['identities'],(v,mv(U,v)),strict=True):
        ids=record['old_wall_indices'];coeff=tuple(map(parse,record['nonnegative_coefficients']))
        assert len(ids)==len(coeff)==4 and all(0<=i<382 for i in ids)
        assert all(sign(c)>0 for c in coeff)
        expected=projected(target);assert tuple(map(parse,record['target_functional']))==expected
        assert tuple(total(mul(c,projected(old[i])[k]) for i,c in zip(ids,coeff)) for k in range(4))==expected
    # Independent finite-dictionary corroboration; the functional identities
    # prove more and have no parameter-box or reflective-norm restriction.
    assert all(sign(pair(v,gram,r))<=0 for r in (*old,*roots))
    return {'status':'EXACT_BOUND_INDEPENDENT_TWO_SLICE_OBSTRUCTION_PASS',
        'base_walls':382,'completion_case':1,'escape_norm':enc(pair(v,gram,v)),
        'strictly_positive_farkas_coefficients':8,'all_1052_dictionary_pairings_nonpositive':True,
        'scope':'This spacelike vector survives every original382-compatible real wall in S or U(S), where S=span(time,omega4,omega7,scalar). It excludes a finite-volume full-dimensional cone formed only by adding such walls; not ambient reflectivity.'}

def produce():
    base_path='results/fresh_rank3_h9_cones_20261013.json';comp_path='results/fresh_rank4_h9_completion_20261013.json'
    base=read(base_path);comp=read(comp_path);old=matrix(base['cases'][0]['all_walls']);U=matrix(base['standard_block_swap_involution'])
    v=tuple(map(parse,comp['cases'][1]['hyperideal_witness']['vector']));records=[]
    for label,target,ids in (('escape',v,(4,28,91,277)),('swapped_escape',mv(U,v),(29,186,226,252))):
        cols=tuple(projected(old[i]) for i in ids);coefs=mv(inverse(tr(cols)),projected(target))
        records.append({'target':label,'old_wall_indices':list(ids),'nonnegative_coefficients':list(map(enc,coefs)),
            'target_functional':list(map(enc,projected(target)))})
    payload={'base_source':base_path,'base_case':0,'completion_source':comp_path,'completion_case':1,
        'fundamental_weight_indices':[4,7],'escape':list(map(enc,v)),
        'escape_norm':comp['cases'][1]['hyperideal_witness']['norm'],
        'restricted_pairing_coordinates':'(-2phi*time, E8_coordinate5, E8_coordinate8, last_scalar)',
        'identities':records,'scope':'Nonnegative restricted-pairing identities against original382 walls, for S and its integral block-swap image. No arithmetic, height or coefficient bound in the resulting obstruction.'}
    payload['audit']=verify(payload);OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');return payload['audit']

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify',action='store_true');a=p.parse_args()
    print(json.dumps(verify(json.loads(OUTPUT.read_text())) if a.verify else produce(),sort_keys=True))
