#!/usr/bin/env python3
"""Small independent rational frame check accompanying source/proof review."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_contact_closure_20261014'

def rref(A):
    A=[r[:]for r in A];p=[];k=0
    for j in range(len(A[0])):
        z=next((i for i in range(k,len(A))if A[i][j]),None)
        if z is None:continue
        A[k],A[z]=A[z],A[k];t=A[k][j];A[k]=[x/t for x in A[k]]
        for i in range(len(A)):
            if i!=k:
                t=A[i][j];A[i]=[a-t*b for a,b in zip(A[i],A[k])]
        p.append(j);k+=1
    return A,p

def run():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    prior=json.loads((ROOT/'results/fresh_h20_coupled_lift_20261014.json').read_text())
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert all(F(x['b'])==0 for row in d['form']for x in row)
    Q=[[F(x['a'])for x in row]for row in d['form']]
    S=prior['base_recipes'][1]['support'];N=[]
    for i in S:
        r=d['normals'][i];assert all(F(x['b'])==0 for x in r);N.append([F(x['a'])for x in r])
    mv=lambda x:[sum(a*b for a,b in zip(row,x))for row in Q]
    ip=lambda x,y:sum(a*b for a,b in zip(x,mv(y)))
    A,p=rref([mv(x)for x in N]);assert len(p)==18
    free=[i for i in range(21)if i not in p];assert free==[11,19,20]
    X=[]
    for j in free:
        x=[F(0)]*21;x[j]=1
        for row,k in zip(A,p):x[k]=-row[j]
        X.append(x)
    G=[[ip(x,y)for y in X]for x in X]
    expected=[[F(x)for x in row]for row in [['143/2','-99/2','165'],['-99/2','71/2','-115'],['165','-115','381']]]
    assert G==expected
    def cap(c):return[sum(x[i]*t for x,t in zip(X,c))for i in range(21)]
    r=cap(list(map(F,['-19/2','13/2','6'])));s=cap(list(map(F,['-191','99','110'])))
    assert ip(r,r)==2 and ip(s,s)==418 and ip(r,s)==0
    a=[F(x['a'])for x in d['normals'][1]];assert all(F(x['b'])==0 for x in d['normals'][1])
    assert ip(r,a)==-1 and ip(s,a)==-37
    for v,w in zip((r,s),d['base_caps'][:2]):assert all(F(x['b'])==0 for x in w)and v==[F(x['a'])for x in w]
    assert {(z['template'],tuple(z['signs']))for z in d['positive_dual_support_obstructions']}=={(i,(1,s))for i in range(4)for s in (1,-1)}
    assert len({tuple(sorted(z['support']))for z in d['finite_edges']})==8
    assert len({tuple(sorted(z['incident']))for z in d['finite_vertices']})==16
    paths=[f'scripts/{TAG}.py',f'results/{TAG}.json','results/fresh_h20_coupled_lift_20261014.json',f'scripts/{TAG}_review.py']
    return dict(status='PASS_SCOPED_SOURCE_PROOF_REVIEW_AND_SEPARATE_RATIONAL_FRAME_CHECK',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                independent_frame_rank=18,complement_dimension=3,norms=['2','418'],mutual_pair='0',anchor_pairings=['-1','-37'],
                distinct_saved_vertices=16,distinct_saved_edges=8,dual_support_records=8,
                scope='Separate Fraction nullspace/frame and cap identities, plus source/proof review and fresh complete primary saved replay. No second full1080-pair arithmetic census or global frontier census.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();dest=ROOT/f'results/{TAG}_review.json'
    if a.verify:assert d==json.loads(dest.read_text())
    else:dest.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
