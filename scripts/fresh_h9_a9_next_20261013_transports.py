#!/usr/bin/env python3
"""Complete oriented diagram-frame transports among nine A4+A4+A1 corners."""
import argparse,hashlib,json
from itertools import product
from fresh_h9_a9_next_20261013_profiles import ROOT,K,Z,O,G,inverse,mv,dot,pair,menc
from fresh_h9_quartic_mixed_20261013_endpoint import square_root

SOURCE=ROOT/'results/fresh_h9_a9_next_20261013_geometry.json'
FRONTIER=ROOT/'results/fresh_h9_a9_next_20261013_frontier.json'
OUT=ROOT/'results/fresh_h9_a9_next_20261013_transports.json'

def paths(roots):
    adj={i:{j for j in range(9)if i!=j and pair(roots[i],G,roots[j])!=Z}for i in range(9)}
    assert all(pair(roots[i],G,roots[j])==-O for i in adj for j in adj[i])
    single=[i for i in adj if not adj[i]];assert len(single)==1
    todo=set(adj)-set(single);ps=[]
    while todo:
        start=min(i for i in todo if len(adj[i])==1);path=[start]
        while True:
            nxt=adj[path[-1]]-set(path)
            if not nxt:break
            assert len(nxt)==1;path.append(nxt.pop())
        assert len(path)==4;todo-=set(path);ps.append(path)
    assert len(ps)==2
    return ps,single[0]

def permutations(E,F):
    ep,ei=paths(E);fp,fi=paths(F)
    for swap,r0,r1 in product((False,True),repeat=3):
        targets=fp[::-1]if swap else fp;targets=[targets[0][::-1]if r0 else targets[0],targets[1][::-1]if r1 else targets[1]]
        iso={ei:fi}
        for a,b in zip(ep,targets):iso.update(zip(a,b))
        yield [iso[i]for i in range(9)]

def run(verify=False):
    data,fd=[json.loads(p.read_text())for p in(SOURCE,FRONTIER)]
    walls=[tuple(map(K.from_json,r))for r in data['all_walls']];frames=[]
    for row in fd['vertices'][1:]:frames.append(([walls[i]for i in row['incident_indices']],tuple(map(K.from_json,row['point'])),row['name']))
    entries=[];accepted=[];I=tuple(tuple(O if i==j else Z for j in range(10))for i in range(10))
    J=tuple(tuple(O if (i==0 and j==0)or(i>0 and j==10-i)else Z for j in range(10))for i in range(10))
    for ai,(E,p,name)in enumerate(frames):
        eg=tuple(tuple(pair(r,G,s)for s in E)for r in E);ei=inverse(eg)
        duals=tuple(mv(ei,tuple(mv(G,r)[j]for r in E))for j in range(10));pcov=tuple(z/pair(p,G,p)for z in mv(G,p))
        for bi,(F,q,label)in enumerate(frames[:ai+1]):
            rho=square_root(pair(p,G,p)/pair(q,G,q));assert rho is not None and rho>0
            for perm in permutations(E,F):
                f=[F[i]for i in perm]
                assert tuple(tuple(pair(r,G,s)for s in f)for r in f)==eg
                H=tuple(tuple(sum((f[k][i]*duals[j][k]for k in range(9)),Z)+rho*q[i]*pcov[j]for j in range(10))for i in range(10))
                row={'source_frame':ai,'target_frame':bi,'permutation':perm,'rho':rho.to_json()}
                bad=next(((i,j)for i in range(10)for j in range(10)if not H[i][j].is_integral),None)
                if bad is not None:
                    i,j=bad;row.update(status='NOT_INTEGRAL',nonintegral_entry={'row':i,'column':j,'value':H[i][j].to_json()})
                else:
                    hi=inverse(H);assert all(z.is_integral for rr in hi for z in rr)
                    assert all(mv(H,r)==s for r,s in zip(E,f))and mv(H,p)==tuple(rho*z for z in q)
                    assert tuple(tuple(pair(tuple(H[k][i]for k in range(10)),G,tuple(H[k][j]for k in range(10)))for j in range(10))for i in range(10))==G
                    kind='IDENTITY'if H==I else'KNOWN_REVERSAL'if H==J else'NEW_INTEGRAL_MAP'
                    row.update(status='INTEGRAL',kind=kind,matrix=menc(H),inverse=menc(hi),trace=sum((H[i][i]for i in range(10)),Z).to_json());accepted.append(row)
                entries.append(row)
        print('checked frame',ai,name,flush=True)
    assert len(entries)==360
    out={'status':'COMPLETE_ORIENTED_A4_A4_A1_FRAME_TRANSPORT_CENSUS','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,FRONTIER)},
      'frames':[name for E,p,name in frames],'frame_pairs_including_self':45,'diagram_maps_per_pair':8,'transports':entries,'integral_count':len(accepted),
      'new_integral_count':sum(r['kind']=='NEW_INTEGRAL_MAP'for r in accepted),
      'scope':'Every oriented equal-norm simple-frame transport between these nine known corners, including each corner automorphism, is tested for lattice integrality in both directions. This finite list is not the full lattice isometry group and cannot establish reflectivity or nonreflectivity by failure to find a new map.'}
    if verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'integral',len(accepted),'new',out['new_integral_count'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
