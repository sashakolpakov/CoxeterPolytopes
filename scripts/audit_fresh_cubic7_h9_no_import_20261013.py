#!/usr/bin/env python3
"""Independent Fraction-polynomial audit of fresh cubic H9 escape files.

No repository modules, search code, field arithmetic, or numeric libraries
are imported. This verifies saved roots, level records, and geometry, not
search exhaustion.
"""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path

ZERO=(F(0),F(0),F(0)); ONE=(F(1),F(0),F(0)); THETA=(F(0),F(1),F(0))
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def neg(x):return tuple(-a for a in x)
def sub(x,y):return add(x,neg(y))
def scale(x,n):return tuple(n*a for a in x)
def mul(x,y):
    p=[F(0)]*5
    for i,a in enumerate(x):
        for j,b in enumerate(y):p[i+j]+=a*b
    # Independent generic reduction by theta^3 + theta^2 -2 theta -1.
    for d in (4,3):
        p[d-3]+=p[d];p[d-2]+=2*p[d];p[d-1]-=p[d]
    return tuple(p[:3])
@lru_cache(None)
def inverse(x):
    cols=[mul(x,b) for b in (ONE,THETA,mul(THETA,THETA))]
    matrix=[[cols[j][i] for j in range(3)]+[F(i==0)] for i in range(3)]
    for k in range(3):
        pivot=next(i for i in range(k,3) if matrix[i][k])
        matrix[k],matrix[pivot]=matrix[pivot],matrix[k]
        t=matrix[k][k];matrix[k]=[z/t for z in matrix[k]]
        for i in range(3):
            if i!=k:
                t=matrix[i][k];matrix[i]=[a-t*b for a,b in zip(matrix[i],matrix[k])]
    return tuple(matrix[i][3] for i in range(3))
def divide(x,y):return mul(x,inverse(y))
def decode(x):return tuple(F(x[k]) for k in ('a','b','c'))
def encoded(x):return dict(zip(('a','b','c'),map(str,x)))
def polynomial(t):return t**3+t*t-2*t-1
def interval_mul(x,y):
    p=[a*b for a in x for b in y];return min(p),max(p)
@lru_cache(None)
def sign(x,place=0):
    if x==ZERO:return 0
    low,high=((F(1),F(2)),(F(-1),F(0)),(F(-2),F(-1)))[place]
    low_polynomial_positive=polynomial(low)>0
    while True:
        p=interval_mul((x[2],x[2]),(low,high))
        p=(p[0]+x[1],p[1]+x[1]);p=interval_mul(p,(low,high))
        p=(p[0]+x[0],p[1]+x[0])
        if p[0]>0:return 1
        if p[1]<0:return -1
        mid=(low+high)/2
        if (polynomial(mid)>0)==low_polynomial_positive:low=mid
        else:high=mid

EDGES=((1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(3,8))
def pairing(x,y,a):
    out=neg(mul(a,mul(x[0],y[0])))
    out=add(out,mul(x[9],y[9]))
    for i in range(1,9):out=add(out,scale(mul(x[i],y[i]),2))
    for i,j in EDGES:out=sub(out,add(mul(x[i],y[j]),mul(x[j],y[i])))
    return out
def gram_vector(v,a):
    out=[neg(mul(a,v[0])),*(scale(v[i],2) for i in range(1,9)),v[9]]
    for i,j in EDGES:out[i]=sub(out[i],v[j]);out[j]=sub(out[j],v[i])
    return out

def verify(p):
    assert p['search_complete'] is True
    assert p['completed_control_levels']==p['total_control_levels']==len(p['levels'])
    delta=sub(scale(ONE,2),THETA)
    expected_a=(THETA,scale(THETA,2),mul(THETA,delta))[p['seed_index']]
    a=decode(p['timelike_coefficient']);assert a==expected_a
    assert [sign(a,i) for i in range(3)]==[1,-1,-1]
    norms=[ONE,scale(ONE,2)]+([delta,scale(delta,2)] if p['seed_index']==2 else [])
    assert [decode(x) for x in p['norm_classes']]==norms
    walls=[tuple(decode(x) for x in row) for row in p['all_walls']]
    assert len(walls)==p['wall_count']==9+sum(row['new_walls'] for row in p['levels'])
    assert all(len(row)==10 for row in walls)
    bound=scale(ONE,F(p['height_bound']))
    assert sign(bound)>0
    previous_height=ZERO;cursor=9;seen_controls=set()
    for index,row in enumerate(p['levels']):
        assert row['index']==index
        control=decode(row['control']);q=decode(row['norm']);height=decode(row['height'])
        assert all(c.denominator==1 for c in control) and sign(control)>0
        assert q in norms
        assert (control,q) not in seen_controls
        seen_controls.add((control,q))
        assert height==divide(mul(a,mul(control,control)),q)
        assert sign(height)>0 and sign(sub(height,previous_height))>=0
        assert sign(sub(bound,height))>=0
        previous_height=height
        count=row['new_walls']
        assert type(count) is int and count>=0
        assert type(row['candidates']) is int and row['candidates']>=count
        assert type(row['nodes']) is int and row['nodes']>=0
        end=cursor+count
        assert end<=len(walls) and row['walls']==end
        for w in walls[cursor:end]:
            assert w[0]==control
            assert pairing(w,w,a)==q
            assert divide(mul(a,mul(w[0],w[0])),pairing(w,w,a))==height
        cursor=end
    assert cursor==len(walls)
    for i,w in enumerate(walls):
        assert all(c.denominator==1 for x in w for c in x)
        q=pairing(w,w,a);assert q in norms
        assert all(sign(q,j)>0 for j in range(3))
        assert all(c.denominator==1 for x in gram_vector(w,a) for c in divide(scale(x,2),q))
        if i<9:assert w==tuple(neg(ONE) if j==i+1 else ZERO for j in range(10))
        else:
            assert sign(w[0])>0
            height=divide(mul(a,mul(w[0],w[0])),q)
            assert sign(sub(bound,height))>=0
        assert all(sign(pairing(w,z,a))<=0 for z in walls[:i])
    inside=tuple(decode(x) for x in p['strict_interior'])
    assert len(inside)==10 and sign(inside[0])>0
    assert sign(pairing(inside,inside,a))<0
    assert all(sign(pairing(inside,w,a))<0 for w in walls)
    geometry=p['geometry'];assert geometry['verdict']=='EXACT_SPACELIKE_ESCAPE'
    escape=tuple(decode(x) for x in geometry['vector'])
    assert len(escape)==10 and sign(escape[0])>0
    q=pairing(escape,escape,a)
    assert q==decode(geometry['norm']) and sign(q)>0
    products=[pairing(escape,w,a) for w in walls]
    assert products==[decode(x) for x in geometry['wall_pairings']]
    assert all(sign(x)<=0 for x in products)
    return {'seed':p['seed_index'],'height_bound':p['height_bound'],'walls':len(walls),'verified_control_levels':len(p['levels']),'verified_escape_norm':encoded(q),'verdict':'PASS_INDEPENDENT_EXACT_ESCAPE_ROOTS_AND_LEVEL_RECORDS'}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('files',type=Path,nargs='+');parser.add_argument('--output',type=Path)
    args=parser.parse_args();results=[]
    for path in args.files:
        p=json.loads(path.read_text());row=verify(p)
        forged=copy.deepcopy(p);forged['geometry']['norm']=encoded(neg(decode(forged['geometry']['norm'])))
        try:verify(forged)
        except AssertionError:pass
        else:raise AssertionError('forged norm was accepted')
        forged=copy.deepcopy(p);forged['geometry']['vector']=copy.deepcopy(forged['strict_interior'])
        forged['geometry']['norm']=encoded(pairing(tuple(map(decode,forged['strict_interior'])),tuple(map(decode,forged['strict_interior'])),decode(forged['timelike_coefficient'])))
        try:verify(forged)
        except AssertionError:pass
        else:raise AssertionError('timelike fake escape was accepted')
        forged=copy.deepcopy(p)
        forged['levels'][0]['height']=encoded(add(decode(forged['levels'][0]['height']),ONE))
        try:verify(forged)
        except AssertionError:pass
        else:raise AssertionError('forged level height was accepted')
        order_test=None
        if len(p['levels'])>=2:
            forged=copy.deepcopy(p)
            # Preserve each packet and its pairing data, so this corruption
            # isolates chronology rather than failing only a wall count.
            packets=[];cursor=9
            for level in forged['levels']:
                end=cursor+level['new_walls']
                packets.append((forged['all_walls'][cursor:end],forged['geometry']['wall_pairings'][cursor:end]))
                cursor=end
            permutation=[1,0,*range(2,len(packets))]
            forged['levels']=[forged['levels'][i] for i in permutation]
            forged['all_walls']=forged['all_walls'][:9]+[w for i in permutation for w in packets[i][0]]
            forged['geometry']['wall_pairings']=forged['geometry']['wall_pairings'][:9]+[x for i in permutation for x in packets[i][1]]
            cursor=9
            for i,level in enumerate(forged['levels']):
                level['index']=i;cursor+=level['new_walls'];level['walls']=cursor
            try:verify(forged)
            except AssertionError:order_test=True
            else:raise AssertionError('nonchronological levels were accepted')
        packet_test=None
        if sum(r['new_walls']>0 for r in p['levels'])>=2:
            forged=copy.deepcopy(p)
            forged['all_walls'][9],forged['all_walls'][-1]=forged['all_walls'][-1],forged['all_walls'][9]
            products=forged['geometry']['wall_pairings']
            products[9],products[-1]=products[-1],products[9]
            try:verify(forged)
            except AssertionError:packet_test=True
            else:raise AssertionError('wall assigned to the wrong level was accepted')
        row.update({'source':str(path),'forged_norm_rejected':True,'timelike_fake_escape_rejected':True,'forged_level_height_rejected':True,'nonchronological_levels_rejected':order_test,'wrong_level_wall_rejected':packet_test});results.append(row)
    result={'status':'PASS','scope':'Independent saved-root, exact chronological level-record, strict-interior, and spacelike-escape validation. Not independent Vinberg enumeration exhaustion.','negative_control_null_meaning':'Not applicable: fewer than two levels or fewer than two nonempty wall packets.','cases':results}
    if args.output:args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
