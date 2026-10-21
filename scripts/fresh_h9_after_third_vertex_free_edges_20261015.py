#!/usr/bin/env python3
"""Exact rank-eight faces with no finite endpoint in the28-wall chamber."""
import argparse,json
from collections import Counter
from itertools import combinations
import fresh_h9_quartet_polar_second_20261014 as B
A,D,R,ROOT=B.A,B.D,B.P,B.ROOT
TAG='fresh_h9_after_third_vertex_free_edges_20261015'
BASE='results/fresh_h9_replacement_escape_extension_20261015.json'

def run():
    d=B.read(BASE);A.RAD=A.A.Q(*d['field']['new_positive_radical_square'])
    Q,N,H=D.dm(d['form']),D.dm(d['normals']),D.dv(d['center']);G=D.dm(d['gram'])
    funcs=[A.mv(Q,n)for n in N]
    rays=[dict(active=a,normal=B.ray(tuple(a),funcs)[0])for a in d['complete_census']['hyperideal_rays']]
    V=[r['normal']for r in rays]
    labels={}
    for i,j in combinations(range(len(N)),2):
        typ,m,_=R.relation(G[i][j],G[i][i],G[j][j])
        if typ in('finite','right'):labels[i,j]=m
    counts=Counter();out=[]
    for i,j in combinations(range(len(rays)),2):
        common=sorted(set(rays[i]['active'])&set(rays[j]['active']))
        if len(common)<8:continue
        counts['shared_rank8_candidates']+=1
        supports=[S for S in combinations(common,8)if all(tuple(sorted(p))in labels for p in combinations(S,2))]
        if not supports:counts['not_spherical_clique']+=1;continue
        pd=[S for S in supports if A.pd([[G[a][b]for b in S]for a in S])]
        if not pd:counts['not_positive_definite']+=1;continue
        assert len(common)==8 and len(pd)==1;S=pd[0]
        u,v=V[i],V[j];q=A.pair(u,u,Q);r=A.pair(v,v,Q);g=A.pair(u,v,Q)
        assert q.sign()>0 and r.sign()>0 and (q*r-g*g).sign()<0
        if g.sign()>=0:counts['misses_hyperbolic_space']+=1;continue
        p=[x-g/r*y for x,y in zip(u,v)];norm=A.pair(p,p,Q)
        assert norm.sign()<0 and A.pair(p,H,Q).sign()<0
        vals=[A.pair(p,n,Q)for n in N]
        assert [k for k,z in enumerate(vals)if not z]==list(S)
        assert all(z.sign()<=0 for z in vals)
        direction=[y-A.pair(v,p,Q)/norm*x for x,y in zip(p,v)]
        assert not A.pair(p,direction,Q)and A.pair(direction,direction,Q).sign()>0
        out.append(dict(id=len(out),support=list(S),hyperideal_endpoints=[i,j],
                        first_ray=u,second_ray=v,start_point=p,direction=direction,
                        norm_polynomial=[norm,A.Z,A.pair(direction,direction,Q)],
                        start_pairings=vals,source_pairing=g))
        counts['complete_geodesic_faces']+=1
        print('line',len(out),'support',S,'rays',i,j,flush=True)
    paths=[BASE,f'scripts/{TAG}.py','scripts/fresh_h9_quartet_polar_second_20261014.py',
           'scripts/fresh_h9_quartet_polar_20261014.py','scripts/fresh_h9_h8_disjoint_walk_20261014_exact.py',
           'scripts/fresh_h9_h8_disjoint_walk_20261014_next_frontier.py']
    return A.enc(dict(status='EXACT_VERTEX_FREE_GEODESIC_FACE_CENSUS',sources={p:R.sha(p)for p in paths},
                     counts=dict(counts),faces=out,
                     scope='Every pair of exterior extreme rays sharing at least8 activewalls is checked for a rank8 positive support. Future timelike interior witnesses certify complete geodesic one-faces having no finite endpoints. Cone completeness inherited from independently verified28-wall census. This is not a cap or volume claim.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'])
