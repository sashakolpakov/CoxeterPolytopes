#!/usr/bin/env python3
"""High-precision geometric diagnostics for structured changed-core faces.

Numerical exploration only. No frozen predecessor is modified.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_coupled_roles_20261014'
INPUTS=['results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json',
        'results/fresh_h20_e7_bridge_face_20261014_exact.json',
        'results/fresh_h20_bridge_edge_closure_20261014.json',
        'results/fresh_h20_polar_pair_20261014.json']
S75=[1,2,3,4,5,6,7,8,9,10,12,13,14,15,17,31,75,87]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load():
    old,bridge,edge,polars=[json.loads((ROOT/p).read_text())for p in INPUTS]
    for d in[old,bridge,edge,polars]:
        for p,h in d['sources'].items():assert sha(p)==h,p
    mp.mp.dps=100
    def f(x):
        q=F(x);return mp.mpf(q.numerator)/q.denominator
    rad=[mp.sqrt(x)for x in old['field']['base_square_root_basis']]
    def v0(x):return sum(f(q)*r for q,r in zip(x,rad))
    D=mp.sqrt(v0(old['field']['discriminant']))
    def v1(x):return v0(x[0])+D*v0(x[1])
    K=mp.sqrt(v1(bridge['field']['new_positive_sqrt_radical']))
    def v2(x):return v1(x[0])+K*v1(x[1])
    Delta=mp.sqrt(v2(edge['field']['new_positive_sqrt_radical']))
    def v3(x):return v2(x[0])+Delta*v2(x[1])
    N=[mp.matrix([v1(x)for x in v])for v in old['normals']]
    N+=[mp.matrix([v2(x)for x in bridge['new_normal']])]
    N +=[mp.matrix([v3(x)for x in v])for v in edge['new_normals']]
    N +=[mp.matrix(list(map(f,v)))for v in polars['new_normals']]
    norms=list(map(f,old['norms']))+[mp.mpf(2)]+[v3(edge['new_normal_norm'])]*2+list(map(f,polars['new_norms']))
    Q=mp.matrix([[f(x)for x in row]for row in old['form']])
    def dot(v,w):return(v.T*Q*w)[0]
    for n,q in zip(N,norms):assert abs(dot(n,n)-q)<mp.mpf('1e-70')*max(1,q)
    N=[n*mp.sqrt(2/q)for n,q in zip(N,norms)]
    M=mp.matrix([[n[j]for j in range(21)]for n in N]);G=M*Q*M.T
    assert max(abs(G[i,i]-2)for i in range(100))<mp.mpf('1e-70')
    E=mp.matrix(list(map(f,polars['strict_escape']['rational_point'])))
    assert abs(dot(E,E)-f(polars['strict_escape']['norm']))<mp.mpf('1e-70')
    H=sum(N[:95],mp.zeros(21,1));assert dot(H,H)<0
    return N,Q,G,H,E,dot
def face_data(N,Q,G,H,E,dot,removed):
    S=[i for i in S75 if i not in removed];J=mp.matrix([[G[i,j]for j in S]for i in S])**-1
    def project(v):
        c=J*mp.matrix([dot(v,N[i])for i in S])
        return v-sum((a*N[i]for a,i in zip(c,S)),mp.zeros(21,1))
    C=project(H);ep=project(E);assert dot(C,C)<0
    rows=[]
    for i in range(100):
        v=project(N[i]);q=dot(v,v);ch=dot(v,C)
        typ='zero'if i in S else'spacelike'if q>mp.mpf('1e-60')else'timelike'if q< -mp.mpf('1e-60')else'near_null'
        if i in S:assert abs(q)<mp.mpf('1e-60')
        if i not in S and i not in removed:assert ch<0
        rows.append(dict(wall=i,removed=i in removed,norm=mp.nstr(q,35),center_pairing=mp.nstr(ch,35),type=typ))
    return dict(removed=removed,support=S,face_dimension=20-len(S),center_norm=mp.nstr(dot(C,C),35),projected_escape_norm=mp.nstr(dot(ep,ep),35),projected_escape_center_pairing=mp.nstr(dot(ep,C),35),retained_spacelike=[r['wall']for r in rows if not r['removed']and r['type']=='spacelike'],retained_near_null=[r['wall']for r in rows if not r['removed']and r['type']=='near_null'],rows=rows)
def run():
    data=load();faces=[face_data(*data,removed)for removed in[[8],[15],[8,15]]]
    paths=INPUTS+[f'scripts/{TAG}.py']
    return dict(status='NUMERICAL_CHANGED_CORE_FACE_DIAGNOSTIC',sources={p:sha(p)for p in paths},precision_decimal_digits=100,faces=faces,scope='Only three specified lower-rank support faces after deleting8,15,orboth. Projected normal signs are high-precision numerical diagnostics, not interval certificates. No cap acceptance, complete-face claim, or global volume conclusion.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    for f in d['faces']:print(f['removed'],'H'+str(f['face_dimension']),'spacelike',f['retained_spacelike'],'null',f['retained_near_null'],'E²',f['projected_escape_norm'],flush=True)
