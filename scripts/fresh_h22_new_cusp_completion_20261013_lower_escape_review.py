#!/usr/bin/env python3
"""Separate Fraction check of the lowered orbit's target cut and escape norm.

The full 5,397,548-inequality escape check is the primary verifier's replay;
this independently reconstructs the cutting orbit member and rational point.
"""
import argparse,hashlib,json
from fractions import Fraction as F
from math import gcd
from functools import reduce
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[1];PRE='fresh_h22_new_cusp_completion_20261013'
SOURCE='results/'+PRE+'_lower_escape.json';OLD='results/'+PRE+'_escape.json';OUT='results/'+PRE+'_lower_escape_review.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json';FRAMES='data/enumerations/'+PRE+'_frames.npz';HIGHEST='data/enumerations/'+PRE+'_highest_roots.npz'
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True)),F(0))
def mv(M,x):return[dot(row,x)for row in M]
def run():
 d=read(SOURCE);old=read(OLD);raw=read(RAW)['coxeter_input'];Q=[[-F(x)for x in r]for r in raw['basis_gram']];W=[[F(x)for x in r]for r in raw['wall_coordinates_as_columns_transposed']]
 assert d['previous_escape_source']==OLD and d['previous_escape_sha256']==sha(OLD)
 assert d['source_sha256']==sha(d['source'])
 arrays=np.load(P/HIGHEST);frames=np.load(P/FRAMES);cut=d['old_target_cut'];j=cut['entry'];ix=cut['cusp_index'];den=cut['old_normal_denominator'];v=cut['twenty_shift'];assert cut['family']=='highest_roots'and int(arrays['cusp_indices'][j])==ix
 root=[F(int(a))+F(v,120)*int(b)for a,b in zip(arrays['old_bases'][j],frames['six_cusp_old'][ix],strict=True)]
 assert root==[F(x,den)for x in cut['old_normal_numerator']]
 q=dot(root,mv(Q,root))+11*F(v,10)**2;assert q==2
 xold=[F(x,old['denominator'])for x in old['numerator'][:22]];zold=F(old['numerator'][22],old['denominator'])
 A=dot(root,mv(Q,xold));B=F(v,10)*zold;assert B>0 and(A>=0 or 11*B*B>A*A)
 assert A==F(int(cut['pairing_numerator_rational']),cut['pairing_denominator'])and B==F(int(cut['pairing_numerator_sqrt11']),cut['pairing_denominator'])
 x=[F(x,d['denominator'])for x in d['numerator'][:22]];z=F(d['numerator'][22],d['denominator']);xc=mv(Q,x);square=dot(x,xc)+z*z;assert square==F(d['square'])and square>0
 h=[int(sum(r[i]for r in W))for i in range(22)];g=reduce(gcd,h);h=[F(a//g)for a in h];assert dot(h,mv(Q,h))==-14
 hp=dot(h,xc);assert hp==F(d['center_pairing'])and hp<0
 oldmaximum=max(dot(w,xc)for w in W);assert oldmaximum<0
 return {'status':'PASS_SEPARATE_FRACTION_LOWER_TARGET_CUT_AND_ESCAPE_POINT','sources':{p:sha(p)for p in(SOURCE,OLD,RAW,FRAMES,HIGHEST)},'cut_wall_norm':str(q),'previous_escape_pairing_A_plus_B_sqrt11':[str(A),str(B)],'positive_sign':True,'new_escape_square':str(square),'new_escape_center_pairing':str(hp),'new_escape_old210_maximum':str(oldmaximum),'scope':'Separate exact reconstruction of the cutting orbit entry, positive target pairing, rational escape norm/time orientation and original210 strict inequalities. All remaining inequalities are covered by the separately rerun primary lower-escape certificate.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(P/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
