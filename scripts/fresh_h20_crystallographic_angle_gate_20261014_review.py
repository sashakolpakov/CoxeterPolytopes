#!/usr/bin/env python3
"""Separate Fraction replay of the six new uniform section-angle tables.

This complements source/proof review and a fresh primary-checker replay.
It does not independently repeat the previous H3 irreducibility census.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_crystallographic_angle_gate_20261014'
SOURCE=ROOT/f'results/{TAG}.json'
OUT=ROOT/f'results/{TAG}_review.json'

def linear_positive(a,b,lo,hi):
    assert a+b*lo>=0 and a+b*hi>=0 and (a or b)
    return True

def ksign(a,b=F(0)):
    a,b=F(a),F(b)
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if (a>0)==(b>0):return 1 if a>0 else -1
    d=a*a-5*b*b
    assert d
    return ((d>0)-(d<0))*(1 if a>0 else -1)

def run():
    d=json.loads(SOURCE.read_text())
    for path,h in d['sources'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
    tables=[]
    for table in d['uniform_pair_tables']:
        lo,hi=map(F,table['q_interval']);counts={}
        for row in table['rows']:
            (da,ka),(db,kb)=row['types'];g=row['old_pairing'];v=row['verdict']
            linear_positive(F(da),F(-ka*ka),lo,hi)
            linear_positive(F(db),F(-kb*kb),lo,hi)
            if v=='RIGHT_TO_UNCHANGED_SOURCE_ROOT':assert g==ka*kb==0 and (ka==0)!=(kb==0)
            elif v=='STRICTLY_DIVERGENT':
                linear_positive(F(-g),F(ka*kb),lo,hi)
                # q² cancels identically in p²-na*nb.
                linear_positive(F(g*g-da*db),F(da*kb*kb+db*ka*ka-2*g*ka*kb),lo,hi)
            else:
                assert v=='SPECIAL_FINITE_RATIO' and g==-1 and (ka==0)!=(kb==0)
            counts[v]=counts.get(v,0)+1
        tables.append({'name':table['name'],'row_counts':counts})
    assert len(tables)==6
    # Independent rational endpoint gates for f(C)=(2C-1)/(6C-4).
    f=lambda c:(2*c-1)/(6*c-4)
    assert f(F(13,16))==F(5,7) and f(F(4,5))==F(3,4)
    assert f(F(16,19))==F(13,20)
    assert ksign(F(5,7)-F(3,8),-F(1,8))>0
    assert ksign(F(3,8)-F(13,20),F(1,8))>0
    assert ksign(F(3,8)-F(5,8),F(1,8))>0
    assert F(2)>F(26,19)**2  # C8>16/19
    assert F(36,25)>F(5,4)  # sqrt5/2<6/5
    # The largest cos²(pi/7) root lies in (4/5,13/16): its cubic
    # is increasing on [3/4,1], and changes sign at these endpoints.
    P=lambda c:64*c**3-80*c*c+24*c-1
    assert P(F(4,5))<0<P(F(13,16)) and 192*F(3,4)**2-160*F(3,4)+24>0
    # The m=5 long-A2 determinant is linear in q=(1+sqrt5)/8.
    for row in d['long_A2_m5_table']:
        (da,ka),(db,kb)=row['types'];g=row['old_pairing']
        assert ksign(F(da)-F(ka*ka,8),-F(ka*ka,8))>0
        assert ksign(F(db)-F(kb*kb,8),-F(kb*kb,8))>0
        if row['verdict']=='RIGHT_TO_UNCHANGED_SOURCE_ROOT':assert g==ka*kb==0
        else:
            A=da*kb*kb+db*ka*ka-2*g*ka*kb
            assert ksign(F(g*g-da*db)+F(A,8),F(A,8))>0
    assert ksign(F(1),-F(1,2))<0  # 2-qBshort at m5
    assert ksign(F(2),-F(1))<0    # 4-4qBlong at m5
    paths=[SOURCE,ROOT/f'scripts/{TAG}.py',Path(__file__).resolve()]
    return {'status':'PASS_SEPARATE_FRACTION_TABLE_REPLAY_AND_SOURCE_PROOF_REVIEW',
            'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
            'uniform_tables':tables,'label5_long_A2_rows':len(d['long_A2_m5_table']),
            'scope':'Separate Fraction/Qsqrt5 arithmetic for the six new interval tables and gap constants. Primary script also freshly replayed. Mathematical review covers arbitrary two-wall contacts and transfer of short-A2 proof to an assumed actual ridge; preceding H3/short-A2 theorem certificates remain dependencies.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    if a.verify:assert d==json.loads(OUT.read_text());print('PASS separate Fraction angle-table review')
    else:OUT.write_text(json.dumps(d,indent=2)+'\n');print('WROTE',OUT)
if __name__=='__main__':main()
