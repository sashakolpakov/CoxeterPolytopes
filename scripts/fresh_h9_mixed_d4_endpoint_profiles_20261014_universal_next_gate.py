#!/usr/bin/env python3
"""Extension of the frozen endpoint gate to all71 directions of the14-wall pair.

Uses the already reviewed dyadic interval engine; reconstructs every new physical
frame and pairing coefficient from the bound fourteen-wall certificate.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate as G

A,I=G.A,G.I
qi,ei,jq,isdot=G.qi,G.ei,G.jq,G.isdot
SQ2,PHI=G.SQ2,G.PHI
ESCAPE=G.ESCAPE
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_next_gate'
BASE='results/fresh_h9_mixed_d4_endpoint_pair_20261014.json'
PROFILES='results/fresh_h9_mixed_d4_endpoint_pair_20261014_next_profiles.json'
ENGINE='scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.py'
FIRST='results/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.json'

def coefficients():
    data=json.loads((ROOT/BASE).read_text())
    pp=json.loads((ROOT/PROFILES).read_text())
    for obj in (data,pp):
        assert obj['sources']
        for path,digest in obj['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest()==digest, path
    assert len(pp['profiles'])==pp['open_edge_count']==71
    assert sum(len(r['finite_profiles']) for r in pp['profiles'].values())==pp['finite_profile_count']==2221
    assert sum(len(r['dihedral_tail_profiles']) for r in pp['profiles'].values())==pp['unbounded_tail_count']==92
    M,N=A.dm(data['form']),A.dm(data['normals'])
    assert len(M)==10 and len(N)==14
    E=[A.E(F(s)) for s in ESCAPE]
    H=A.dv(data['center'])
    en=A.pair(E,E,M);eh=A.pair(E,H,M)
    ep=[A.pair(E,n,M) for n in N]
    assert en.sign()>0 and eh.sign()<0 and all(x.sign()<0 for x in ep)
    norm=[A.pair(n,n,M) for n in N]
    scales=[(2/ei(q)).sqrt() for q in norm]
    cross=[[ei(A.pair(a,b,M))*scales[i]*scales[j]
            for j,b in enumerate(N)] for i,a in enumerate(N)]
    crossE=[ei(x)*s for x,s in zip(ep,scales)]
    edges=[]
    for es,row in pp['profiles'].items():
        support=row['support']; outside=[i for i in range(len(N)) if i not in support]
        B=[[jq(v) for v in rr] for rr in row['normalized_Gram']]
        J=[[jq(v) for v in rr] for rr in row['normalized_inverse']]
        assert all(sum((B[i][k]*J[k][j] for k in range(8)),A.Q())==A.Q(i==j)
                   for i in range(8) for j in range(8))
        # Verify the saved normalized Gram by signs and squared exact values.
        for i,ii in enumerate(support):
            for j,jj in enumerate(support):
                raw=A.pair(N[ii],N[jj],M); want=A.E(B[i][j])
                assert raw.sign()==want.sign()
                assert 4*raw*raw==want*want*norm[ii]*norm[jj]
        assert A.inertia([[A.E(v) for v in rr] for rr in B])==[8,0]
        X=A.dv(data['complete_ray_census']['rays'][row['start_ray']]['vector'])
        other=A.dv(data['complete_ray_census']['rays'][row['other_cone_ray']]['vector'])
        assert A.pair(other,other,M).sign()>0
        xn=-A.pair(X,X,M); xy=A.pair(X,other,M)
        Y=[xn*y+xy*x for x,y in zip(X,other)]
        xx=-A.pair(X,X,M); yy=A.pair(Y,Y,M)
        assert xx.sign()>0 and yy.sign()>0 and not A.pair(X,Y,M)
        assert A.pair(X,H,M).sign()<0
        assert all(A.pair(X,n,M).sign()<=0 for n in N)
        assert all(not A.pair(Z,N[i],M) for i in support for Z in (X,Y))
        assert A.pair(X,N[row['released_wall']],M)==0
        assert A.pair(Y,N[row['released_wall']],M).sign()<0
        dx,dy=ei(xx).sqrt(),ei(yy).sqrt()
        avec=[ei(A.pair(X,N[j],M))*scales[j]/dx for j in outside]
        cvec=[ei(A.pair(Y,N[j],M))*scales[j]/dy for j in outside]
        avec.append(ei(A.pair(X,E,M))/dx)
        cvec.append(ei(A.pair(Y,E,M))/dy)
        assert avec[-1].hi<0<cvec[-1].lo and (avec[-1]+cvec[-1]).hi<0
        # Every scalar pairing strictly decreases with h>=0.
        for a,c in zip(avec,cvec):
            assert a.hi<=0 and (a+c).hi<0
        cols=[[cross[i][j] for i in support] for j in outside]
        cols.append([crossE[i] for i in support])
        profiledata=[]
        for p in row['finite_profiles']:
            cc=[I(0)]*8
            active=[];cost=A.Q()
            for wall,m in p['contacts']:
                assert m in (3,4,5)
                j=support.index(wall);active.append(j)
                cc[j]={3:I(1),4:SQ2,5:PHI}[m]
                cost+=J[j][j]*{3:A.Q(1),4:A.Q(2),5:A.Q(F(3,2),F(1,2))}[m]
            assert all(not J[i][j] for i in active for j in active if i!=j)
            assert jq(p['residual_norm'])==2-cost
            weights=[sum((qi(J[i][j])*cc[j] for j in range(8)),I(0)) for i in range(8)]
            b=[-isdot(weights,col) for col in cols]
            rho=qi(jq(p['residual_norm']));assert rho.lo>0
            profiledata.append(dict(kind='finite',contacts=p['contacts'],b=b,rho=rho))
        for tail in row['dihedral_tail_profiles']:
            j=support.index(tail['wall'])
            assert J[j][j]==A.Q(F(1,2)) and all(not J[i][j] for i in range(8) if i!=j)
            weights=[qi(J[i][j]) for i in range(8)]
            b=[-isdot(weights,col) for col in cols]
            profiledata.append(dict(kind='tail',wall=tail['wall'],b=b))
        edges.append(dict(edge_id=int(es),support=support,outside=outside,a=avec,c=cvec,profiles=profiledata,from_new_finite_vertex=row['from_new_finite_vertex']))
    return edges,dict(coordinates=ESCAPE,norm=A.enc(en),future_center_pairing=A.enc(eh),old_pairings=A.enc(ep))

def run():
    first=json.loads((ROOT/FIRST).read_text())
    for path,digest in first['sources'].items():
        assert sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    edges,escape=coefficients();records=[];total=Counter();priority=Counter()
    for edge in edges:
        rows=[]
        for p in edge['profiles']:
            r=G.certify(edge,p);rows.append(r);total.update(r['counts']);total[p['kind']+'_profiles']+=1
            if edge['from_new_finite_vertex']:priority[p['kind']+'_profiles']+=1
        records.append(dict(edge_id=edge['edge_id'],from_new_finite_vertex=edge['from_new_finite_vertex'],profiles=rows))
        print(edge['edge_id'],len(rows),'boxes',sum(r['counts']['boxes'] for r in rows),flush=True)
    assert len(records)==71 and total['finite_profiles']==2221 and total['tail_profiles']==92
    assert sum(e['from_new_finite_vertex'] for e in edges)==21
    assert priority==dict(finite_profiles=459,tail_profiles=21)
    sources=[BASE,PROFILES,ENGINE,FIRST,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',f'scripts/{TAG}.py']
    return dict(status='PASS_ALL_71_ENDPOINT_FAMILIES_MISS_SAME_RATIONAL_GAP',
        sources={s:sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
        dyadic_precision_bits=G.BITS,escape=escape,counts=dict(total),new_vertex_counts=dict(priority),edges=records,
        max_depth=max(r['maximum_depth'] for e in records for r in e['profiles']),
        terminal_boxes=sum(v for k,v in total.items() if k=='escape' or k.startswith('wall:')),
        parameterization=first['parameterization'],allowed_pairing_relaxation=first['allowed_pairing_relaxation'],
        cosine_enclosures=first['cosine_enclosures'],cosine_error_proof=first['cosine_error_proof'],
        scope='All71 open directions of the14-wall chamber, including all21 at its3new vertices. Every finite-endpoint cap preserving those14 Coxeter walls strictly misses E. No arbitrary-cap, affine-endpoint or later-support exclusion.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result==json.loads(path.read_text());print('SAVED_REPLAY_PASS',result['counts'])
    else:
        path.write_text(json.dumps(result,indent=2)+'\n');print(result['status'],result['counts'])
