#!/usr/bin/env python3
"""All-real endpoint test on the changed fifteen-wall seed's exact frontier."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
import fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate as B

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_h8_disjoint_walk_20261014_universal_gate'
BASE='results/fresh_h9_h8_disjoint_cap_20261014.json'
FRONT='results/fresh_h9_h8_disjoint_walk_20261014.json'
I,ei,qi=B.I,B.ei,B.qi
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def coefficients():
    d,f=[json.loads((ROOT/p).read_text())for p in(BASE,FRONT)]
    for obj in(d,f):
        for p,h in obj['sources'].items():assert sha(p)==h,p
    assert f['open_count']==62 and f['finite_profile_count']==1076 and f['symbolic_dihedral_tail_count']==64
    Q,N=A.dm(d['form']),A.dm(d['normals']);E=A.dv(d['strict_escape']['point']);H=A.dv(d['center'])
    assert A.pair(E,E,Q).sign()>0 and A.pair(E,H,Q).sign()<0
    norms=[A.pair(n,n,Q)for n in N];ep=[A.pair(E,n,Q)for n in N]
    assert all(x.sign()<0 for x in ep)
    scales=[(2/ei(n)).sqrt()for n in norms]
    G=[[ei(A.pair(n,m,Q))*scales[i]*scales[j]for j,m in enumerate(N)]for i,n in enumerate(N)]
    gE=[ei(x)*s for x,s in zip(ep,scales)];edges=[]
    for germ in f['queue']:
        row=germ['profile_data'];support=row['support'];outside=[i for i in range(15)if i not in support]
        M=[[A.Q(*x)for x in rr]for rr in row['normalized_Gram']]
        J=[[A.Q(*x)for x in rr]for rr in row['normalized_inverse']]
        assert A.inertia([[A.E(x)for x in rr]for rr in M])==[8,0]
        for i,ii in enumerate(support):
            for j,jj in enumerate(support):
                assert sum((M[i][k]*J[k][j]for k in range(8)),A.Q())==A.Q(i==j)
                raw=A.pair(N[ii],N[jj],Q);want=A.E(M[i][j])
                assert raw.sign()==want.sign() and 4*raw*raw==want*want*norms[ii]*norms[jj]
        X,Y=A.dv(germ['start_point']),A.dv(germ['direction'])
        xx,yy=-A.pair(X,X,Q),A.pair(Y,Y,Q)
        assert xx.sign()>0 and yy.sign()>0 and not A.pair(X,Y,Q)
        assert A.pair(X,H,Q).sign()<0 and all(A.pair(X,n,Q).sign()<=0 for n in N)
        assert all(not A.pair(Z,N[i],Q)for Z in(X,Y)for i in support)
        release=germ['released_wall']
        assert not A.pair(X,N[release],Q) and A.pair(Y,N[release],Q).sign()<0
        dx,dy=ei(xx).sqrt(),ei(yy).sqrt()
        avec=[ei(A.pair(X,N[j],Q))*scales[j]/dx for j in outside]+[ei(A.pair(X,E,Q))/dx]
        cvec=[ei(A.pair(Y,N[j],Q))*scales[j]/dy for j in outside]+[ei(A.pair(Y,E,Q))/dy]
        assert avec[-1].hi<0<cvec[-1].lo and (avec[-1]+cvec[-1]).hi<0
        assert all(a.hi<=0 and(a+c).hi<0 for a,c in zip(avec,cvec))
        cols=[[G[i][j]for i in support]for j in outside]+[[gE[i]for i in support]]
        profiles=[]
        for p in row['finite_profiles']:
            cc=[I(0)]*8;active=[];cost=A.Q()
            for wall,m in p['contacts']:
                j=support.index(wall);active.append(j)
                cc[j]={3:I(1),4:B.SQ2,5:B.PHI}[m]
                cost+=J[j][j]*{3:A.Q(1),4:A.Q(2),5:A.Q(F(3,2),F(1,2))}[m]
            assert all(not J[i][j]for i in active for j in active if i!=j)
            assert A.Q(*p['residual_norm'])==2-cost
            weights=[sum((qi(J[i][j])*cc[j]for j in range(8)),I(0))for i in range(8)]
            b=[-B.isdot(weights,col)for col in cols];rho=qi(A.Q(*p['residual_norm']))
            assert rho.lo>0
            profiles.append(dict(kind='finite',contacts=p['contacts'],b=b,rho=rho))
        for tail in row['dihedral_tail_profiles']:
            j=support.index(tail['wall'])
            assert J[j][j]==A.Q(F(1,2)) and all(not J[i][j]for i in range(8)if i!=j)
            weights=[qi(J[i][j])for i in range(8)]
            profiles.append(dict(kind='tail',wall=tail['wall'],b=[-B.isdot(weights,col)for col in cols]))
        edges.append(dict(edge_id=germ['id'],support=support,outside=outside,a=avec,c=cvec,profiles=profiles))
    return edges
def run(limit=None):
    edges=coefficients();records=[];total=Counter()
    for edge in edges:
        rows=[]
        for p in edge['profiles']:
            r=B.certify(edge,p);rows.append(r);total.update(r['counts']);total[p['kind']+'_profiles']+=1
        records.append(dict(edge_id=edge['edge_id'],profiles=rows))
        print(edge['edge_id'],len(rows),'boxes',sum(r['counts']['boxes']for r in rows),flush=True)
        if limit is not None and len(records)>=limit:break
    if limit is None:assert len(records)==62 and total['finite_profiles']==1076 and total['tail_profiles']==64
    paths=[BASE,FRONT,f'scripts/{TAG}.py','scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',
           'scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.py']
    return dict(status='PASS_CHANGED_H8_ALL_ENDPOINT_PROFILES_AND_CONTINUOUS_TAILS_MISS_ESCAPE',
        sources={p:sha(p)for p in paths},dyadic_precision_bits=B.BITS,counts=dict(total),edges=records,
        parameterization='r=p+hU+sqrt(h²+rho)V,h>=0. Tails retain the continuous closed interval sqrt3<=k<=2, with p=k*p0,rho=2-k²/2.',
        allowed_pairing_relaxation='Labels2..15 individually; every label>=16, parallel and dotted pairings retained in (-infinity,-2cos(pi/16)].',
        cosine_enclosures=dict(pi=B.PI.enc(),isolated=[x.enc()for x in B.ALLOWED_POINTS],tail_threshold=B.HIGH_THRESHOLD.enc()),
        maximum_depth=max(r['maximum_depth']for e in records for r in e['profiles']),
        scope='Every finite endpoint cap on the complete62 saved finite-corner germs, retaining all15 walls, strictly misses the literal seed escape. Exact outward dyadic arithmetic covers all real endpoint parameters and every integer Coxeter label. Complete finite-profile classification is inherited. This does not exclude intermediate caps, later supports, ideal endpoints, other wall directions or changed chambers. The inherited interval algorithm is reused, not independently reimplemented.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--limit',type=int);a=ap.parse_args()
    result=run(a.limit);out=ROOT/f'results/{TAG}.json'
    if a.limit is not None:print('PARTIAL',result['counts'])
    elif a.verify:assert result==json.loads(out.read_text());print('SAVED_REPLAY_PASS',result['counts'])
    else:out.write_text(json.dumps(result,indent=2)+'\n');print(result['status'],result['counts'])
