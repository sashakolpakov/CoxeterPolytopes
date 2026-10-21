#!/usr/bin/env python3
"""Parameterize, then exactly recognize, the two D9 Pell transports."""
from fresh_h9_glued_host_20261013_arithmetic import *
from pathlib import Path
import argparse,hashlib,json

ROOT=Path(__file__).resolve().parents[1]
PREFIX='fresh_h9_tripled_d9_20261013_mechanism'
OUT=ROOT/'results'/f'{PREFIX}.json'
IDENT=tuple(tuple(O if i==j else Z for j in range(10))for i in range(10))
def mm(A,B):return tuple(tuple(dot(r,c)for c in zip(*B))for r in A)


def check():
    rows=[];sources={}
    for c in(1,3):
        if c==1:
            fp='fresh_h9_changed_spatial_20261013';firstp=f'{fp}.json'
            targetp='fresh_h9_changed_spatial_next_20261013.json';symp='fresh_h9_d9_symmetry_transfer_20261013.json';key='transport_matrix'
        else:
            fp='fresh_h9_tripled_d9_20261013';firstp=f'{fp}_edge1.json';targetp=f'{fp}_edge3.json';symp=f'{fp}_symmetry.json';key='matrix'
        pp=[ROOT/'results'/p for p in(firstp,targetp,symp)]
        first,target,sym=[json.loads(p.read_text())for p in pp]
        for p in pp:sources[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
        a=c*ALPHA;zz,xx=map(K.from_json,first['relative_unit']);r=tuple(map(K.from_json,target['chosen_root']))
        ww=sum(r[1:4],Z)/2;yy=r[0]/2
        assert zz*zz-a*xx*xx==ww*ww-3*a*yy*yy==1
        assert zz>1 and ww>1 and xx>0 and yy>0
        assert all(x.is_integral for x in(zz,xx,ww,yy))and((ww+1)/3).is_integral
        A=IDENT
        for j in range(1,4):
            B=[list(row)for row in IDENT];B[0][0]=B[j][j]=zz;B[0][j]=-xx;B[j][0]=-a*xx
            A=mm(tuple(map(tuple,B)),A)
        C=[list(row)for row in IDENT];C[0][0]=ww
        for i in range(1,4):
            C[0][i]=yy;C[i][0]=a*yy
            for j in range(1,4):C[i][j]=(ww+1)/3-(1 if i==j else 0)
        P=[list(row)for row in IDENT];P[1]=list(IDENT[3]);P[3]=list(IDENT[1])
        F=mm(mm(tuple(map(tuple,C)),tuple(map(tuple,P))),A)
        assert F==tuple(tuple(map(K.from_json,row))for row in sym[key])
        qs=(a*(3*yy*zz**3-ww*xx*(1+zz+zz**2)),
            ww*(zz+a*xx**2*(1+zz))-3*a*xx*yy*zz**2,
            ww*(zz+a*xx**2)-3*a*xx*yy*zz,
            ww*zz-3*a*xx*yy)
        assert qs==tuple(sum((F[i][j]for i in range(1,4)),Z)for j in range(4))
        residues=(qs[0]/2,*((x-1)/2 for x in qs[1:]));assert all(x.is_integral for x in residues)
        trace=(ww*a*xx**2*zz+2*ww*a*xx**2+3*ww*zz**3+3*ww*zz
               -2*a*xx**2*zz+2*a*xx**2-6*a*xx*yy*(zz**2+zz+1)+18)/3
        assert trace==sum((F[i][i]for i in range(10)),Z)>10
        rows.append({'c':c,'axis_ZX':enc((zz,xx)),'triple_WY':enc((ww,yy)),
            'exact_factorization':'F=C*P*B3*B2*B1','D9_parity_coefficients':enc(qs),
            'parity_quotients':enc(residues),'trace':trace.to_json()})
    return {'status':'PASS_EXACT_COMMON_D9_TRANSPORT_MECHANISM','sources':sources,'examples':rows,
        'conditional_scope':'This is a symbolic sufficient mechanism once the four complete first-edge/corner hypotheses hold. It does not prove nonreflectivity for every coefficient c. In a new coefficient, all primitive norm classes, first-contact profiles, parity congruences and trace must be checked again.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();out=check()
    if args.verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
