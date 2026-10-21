#!/usr/bin/env python3
"""Finite perturbation scout around the remaining H5 escape.

Only exact isotropic-direction inequalities select PARI norm-equation inputs.
These are proposal directions, never asserted reflective roots.
"""
import json
from fresh_cubic9_positive_20261013 import *

def gp(x):
    return f'Mod(({x.a})+({x.b})*y+({x.c})*y^2,f)'

def run():
    source=ROOT/'results/fresh_h9_joint_arithmetic_20261013_pell2.json'
    data=json.loads(source.read_text());G=Host(T,0).G;sp=tuple(row[1:]for row in G[1:])
    walls=[tuple(map(K.from_json,r))for r in data['forced_h5_face']['all_walls']]
    target=tuple(map(K.from_json,data['forced_h5_face']['geometry']['hyperideal_witness']['vector']))
    d=tuple(30*v/(T*T-1)for v in target[1:]);assert all(t.is_integral for t in d)
    directions=[(*col,Z)for col in zip(*(tuple(5*t for t in row)for row in AINV))]
    directions.append((Z,Z,Z,Z,K(5)))
    def feasible(e):
        delta=pair(e,sp,e);b=pair(e,sp,target[1:])
        if b.sign()<=0 or (b*b-T*delta).sign()<=0:return None
        for w in walls:
            p=pair(e,sp,w[1:])
            if w[0]==Z:
                if p.sign()>0:return None
            elif p.sign()>0 and (T*delta*w[0]*w[0]-p*p).sign()<=0:return None
        return delta
    rows=[]
    for j in range(1,5):
        for k,u in enumerate(directions):
            for sign in (-1,1):
                e=tuple(a+sign*D**j*b for a,b in zip(d,u));delta=feasible(e)
                if delta is None:continue
                assert all(t.is_integral for t in e)
                rows.append({'unit_power':j,'direction':k,'sign':sign,'spatial_vector':enc(e),
                    'spatial_norm':delta.to_json(),'absolute_norm':str(delta.norm())})
    rows.sort(key=lambda row:abs(F(row['absolute_norm'])))
    out={'status':'EXACT_FINITE_DIRECTION_SCREEN','source':str(source.relative_to(ROOT)),
        'tested':40,'passed':len(rows),'target':enc(target),'base_spatial_vector':enc(d),'directions':rows,
        'scope':'Forty explicitly prescribed unit-scaled spatial perturbations. Exact open inequalities certify that a sufficiently long positive-norm root on a passed direction would cut the target and pair acutely with all old walls. No root existence or angle conclusion follows from the direction screen.'}
    output=ROOT/'results/fresh_h9_joint_arithmetic_20261013_pell_scout.json'
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    lines=['f=y^3-3*y-1;nf=nfinit(f);t=Mod(y,f);']
    for i,row in enumerate(rows):
        delta=K.from_json(row['spatial_norm'])
        lines.append(f'delta={gp(delta)};print({i}," ",vector(4,j,nfhilbert(nf,t*delta,-[1,2,5,10][j]*t)));')
    path=ROOT/'scripts/fresh_h9_joint_arithmetic_20261013_pell_scout.gp'
    path.write_text('\n'.join(lines)+'\n');print(json.dumps({'directions':len(rows)}))
    return out
if __name__=='__main__':run()
