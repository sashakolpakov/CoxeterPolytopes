#!/usr/bin/env python3
"""Fixed-uniform-self-edge common quadratic-space obstruction; no vertex exclusion."""
import argparse,json,math
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_vertex_common_space_20261016'
EDGE='results/fresh_sqrt5_universal_self_edges_20261016.json'
FRONTS='results/fresh_sqrt5_h3_center_selected_diagnostic_20261016.json'
def k(a=0,b=0):return F(a),F(b)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def sub(x,y):return add(x,neg(y))
def mul(x,y):
    a,b=x;c,d=y;return a*c+b*d,a*d+b*c+b*d
def norm(x):a,b=x;return a*a+a*b-b*b
def div(x,y):
    a,b=y;n=norm(y);assert n;return mul(x,((a+b)/n,-b/n))
def enc(x):return list(map(str,x))
def dec(x):return tuple(map(F,x))
def square_root_rational(x):
    if x<0:return None
    a=math.isqrt(x.numerator);b=math.isqrt(x.denominator)
    return F(a,b)if a*a==x.numerator and b*b==x.denominator else None
def square_test(x):
    nn=norm(x);s=square_root_rational(nn)
    if s is None:return dict(square=False,reason='nonsquare rational norm',norm=str(nn))
    A=x[0]+x[1]/2;B=x[1]/2;attempts=[]
    for t in sorted(set([s,-s])):
        xx=(A+t)/2;u=square_root_rational(xx);attempts.append(str(xx))
        if u:
            v=B/(2*u);answer=(u-v,2*v)
            if mul(answer,answer)==x:return dict(square=True,root=enc(answer),norm=str(nn))
        elif u==0 and B==0:
            v=square_root_rational(A/5)
            if v is not None:
                answer=(-v,2*v);assert mul(answer,answer)==x;return dict(square=True,root=enc(answer),norm=str(nn))
    return dict(square=False,reason='rational norm square but no rational quadratic coordinates',norm=str(nn),candidate_x_squared=attempts)
def inverse(M):
    n=len(M);A=[row[:]+[k(i==j)for j in range(n)]for i,row in enumerate(M)]
    for col in range(n):
        p=next(i for i in range(col,n)if A[i][col]!=k());A[col],A[p]=A[p],A[col];d=A[col][col];A[col]=[div(x,d)for x in A[col]]
        for i in range(n):
            if i!=col:
                d=A[i][col];A[i]=[sub(a,mul(d,b))for a,b in zip(A[i],A[col])]
    assert all(A[i][j]==k(i==j)for i in range(n)for j in range(n));return[row[n:]for row in A]
def run():
    edge=json.loads((ROOT/EDGE).read_text());front=json.loads((ROOT/FRONTS).read_text());lam=k(13,21)
    assert norm(lam)==1 and 2*lam[0]+lam[1]==47
    components=[];index={}
    for rec in edge['component_records']:
        A=[[dec(x)for x in row]for row in rec['integral_gram']];J=inverse(A);ps=[];ports=[]
        for i,w in enumerate(rec['root_witnesses']):
            assert w['root']==i and J[i][i]==dec(w['finite_inverse_diagonal'])and A[i][i]==dec(w['norm'])
            d=A[i][i];tau=mul(d,lam);p=sub(mul(k(2),tau),mul(mul(tau,tau),J[i][i]));ps.append(p)
            ports.append(dict(root=i,norm=enc(d),inverse_diagonal=enc(J[i][i]),edge_Schur=enc(p),edge_Schur_norm=str(norm(p))))
        pairs=[]
        for i,j in combinations(range(len(A)),2):
            ratio=div(ps[i],ps[j]);test=square_test(ratio)
            pairs.append(dict(roots=[i,j],ratio=enc(ratio),**test))
        bad=next((p for p in pairs if not p['square']),None)
        item=dict(type=rec['type'],rank=len(A),ports=ports,pair_squareclass_checks=pairs,within_component_obstruction=bad)
        components.append(item);index[item['type']]=item
    d5=index['D5'];i=next(p for p in d5['ports']if p['inverse_diagonal']==['1','0']);j=next(p for p in d5['ports']if p['inverse_diagonal']==['5/4','0'])
    assert F(i['edge_Schur_norm'])==-720 and F(j['edge_Schur_norm'])==-899
    assert square_root_rational(F(720,899))is None
    records=[]
    for idx,rec in enumerate(front['maximizing_witnesses']):
        forest=rec['forest'];assert sum(index[name]['rank']*count for name,count in forest.items())==21
        bad=next((index[name]['within_component_obstruction']for name,count in sorted(forest.items())if count and index[name]['within_component_obstruction']),None)
        badname=next((name for name,count in sorted(forest.items())if count and index[name]['within_component_obstruction']),None)
        records.append(dict(front_index=idx,forest=forest,marked_component=rec['marked_component'],facet_root=rec['facet_root'],fixed_uniform_edge_obstruction_component=badname,witness=bad))
    assert len(records)==306
    files=[f'scripts/{TAG}.py',EDGE,'results/fresh_sqrt5_universal_self_edges_20261016_review.json',FRONTS]
    return dict(status='PASS_EXACT_FIXED_UNIFORM_SELF_EDGE_COMMON_SPACE_DIAGNOSTIC',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files},lambda_value=enc(lam),component_records=components,maximizing_forest_records=records,counts=dict(Counter('FIXED_UNIFORM_EDGE_CONFLICT'if r['witness']else'NO_WITHIN_COMPONENT_CONFLICT_DETECTED'for r in records)),d5_witness=dict(roots=[i['root'],j['root']],inverse_diagonals=['1','5/4'],Schur_norms=['-720','-899'],ratio_norm='720/899'),scope='Only the fixed lambda=phi^8 duplicate-root witnesses are tested. A within-component mismatch persists under common component rescaling and obstructs their simultaneous realization in one Qphi quadratic space, irrespective of new-to-new pairings. It does not exclude the finite vertex forest, variable edge parameters, different endpoint diagrams, or compact dimension21.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
