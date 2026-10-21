#!/usr/bin/env python3
"""Remove M's star and test the changed L-X-R finite duty exactly.

The H4 duty and high N-X bond force R perpendicular to L,X,N.
An original retained wall52 excludes both polar signs at every native
high-N label. No S bank, cone reconstruction or volume claim is made.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_oblique_low_s_high_r_polar_interval_20261018 as P

ROOT=P.ROOT; A=P.A
TAG='fresh_oblique_m_star_retriangulation_interval_20261018'
V,add,sub,mul,div,neg,sqrt,dot=P.V,P.add,P.sub,P.mul,P.div,P.neg,P.sqrt,P.dot
GAPS={0:(19,20),1:(16,17)}


def topology(old):
    oldtri=old['triangles']
    removed=[[68,69,71],[69,70,71],[70,71,72],[68,71,72]]
    keys={tuple(sorted(x))for x in removed}
    assert keys<={tuple(sorted(x))for x in oldtri}
    tri=[list(x)for x in oldtri if tuple(sorted(x))not in keys]+[[68,69,70],[68,70,72]]
    vertices=sorted({i for f in tri for i in f})
    edge=Counter(tuple(sorted(e))for f in tri for e in combinations(f,2))
    boundary=[8,10,73,70,72,52]
    be={tuple(sorted((a,b)))for a,b in zip(boundary,boundary[1:]+boundary[:1])}
    assert set(e for e,c in edge.items()if c==1)==be and set(edge.values())=={1,2}
    links={}
    for v in vertices:
        le=[sorted(set(f)-{v})for f in tri if v in f];adj={i:[]for e in le for i in e}
        for x,y in le:adj[x].append(y);adj[y].append(x)
        seen=set();todo=[next(iter(adj))]
        while todo:
            x=todo.pop()
            if x in seen:continue
            seen.add(x);todo+=adj[x]
        assert seen==set(adj)
        deg=sorted(map(len,adj.values()))
        assert deg==([1,1]+[2]*(len(deg)-2)if v in boundary else[2]*len(deg))
        links[str(v)]=le
    assert(len(vertices),len(edge),len(tri))==(9,18,10)
    assert 71 not in vertices and (68,70)in edge
    return dict(triangles=tri,removed_triangles=removed,boundary=boundary,vertices=vertices,
                V=9,E=18,F=10,links=links,euler=1)


def parents(prior):
    base={int(i):{int(j):V(z)for j,z in row.items()}for i,row in prior['refined_primary']['projected_Gram'].items()}
    out=[]
    for row in prior['records']:
        g={i:r.copy()for i,r in base.items()}
        steps=row['refined_parent_steps'][:3]
        assert [s['wall']for s in steps]==[67,68,69]
        for s in steps:P.restore(g,s)
        assert 70 not in g and 71 not in g
        assert steps[0]['profile']==steps[2]['profile']==0 and steps[1]['profile']in(1,4)
        assert steps[2]['frame']==[8,67,68]and steps[2]['labels']==[2,2,5]
        lp=[V(z)for z in prior['profiles'][steps[1]['profile']]]
        assert lp[:4]==[V(0)]*4 and sorted(lp[4:])==[V(-1),V(0)]
        assert g[68][68][0]<=F(4,3)<=g[68][68][1]and g[69][69]==V(2)
        h4minor=sub(mul(g[68][68],g[69][69]),mul(g[68][69],g[69][68]))
        assert h4minor[0]>0
        out.append((row,g,dict(L_profile=steps[1]['profile'],L_white_covector=lp,
            H4_path_labels=[3,3,5],projected_LX_minor=h4minor,
            forced_R_profile=0,forced_R_frame=[69,73,68],forced_R_products=[0,0,0])))
    assert len(out)==2
    return out


def check(g,n,a,side):
    gn,np=P.n_extension(g,n,a)
    frame=[69,73,68];gram=[[gn[i][j]for j in frame]for i in frame]
    inv,det=A.P.inverse(gram)
    assert gram[0][0][0]>0 and det[1]<0
    c=[gn[i][52]for i in frame]
    residual=sub(gn[52][52],dot(c,inv,c))
    square=mul(V(2),residual)
    # R has profile0, so its projected product with old52 is physical.
    if side=='below_gap':assert 1<square[0]<=square[1]<2
    else:assert side=='above_gap'and 0<square[0]<=square[1]<1
    return dict(N=np,frame=frame,Gram=gram,inverse=inv,determinant=det,
        old52_covector=c,old52_projected_square=gn[52][52],old52_residual_square=residual,
        both_R_signs_old52_product_square=square,rejection=side)


def run(saved=None):
    prior=json.loads((ROOT/f'results/{A.TAG}.json').read_text())
    for p,h in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert prior['summary']['unresolved']==0
    records=[];domains=[];attempts=0
    for parent,g,duty in parents(prior):
        for a,(k,l)in GAPS.items():
            lo=A.theta(6)[0];kb=A.theta(k);lb=A.theta(l)
            assert k+1==l and kb[1]<lb[0]<2
            ranges=[('below_gap',(lo,kb[1])),('above_gap',(lb[0],F(2)))]
            domains.append(dict(parent=parent['index'],N10_magnitude=a,gap_labels=[k,l],
                gap_left_theta=kb,gap_right_theta=lb,ranges=ranges,duty=duty))
            for side,(left,right)in ranges:
                todo=[(left,right,0)]
                while todo:
                    x,y,depth=todo.pop();attempts+=1
                    try:proof=check(g,(x,y),a,side)
                    except(AssertionError,ZeroDivisionError):
                        assert depth<18,('unresolved',parent['index'],a,side,x,y)
                        mid=(x+y)/2;todo.extend([(x,mid,depth+1),(mid,y,depth+1)])
                    else:records.append(dict(index=len(records),parent=parent['index'],N10_magnitude=a,
                        side=side,domain=(x,y),depth=depth,proof=proof))
                parts=sorted(r['domain']for r in records if r['parent']==parent['index']and r['N10_magnitude']==a and r['side']==side)
                assert parts[0][0]==left and parts[-1][1]==right and all(x[1]==y[0]for x,y in zip(parts,parts[1:]))
                print('RANGE_PASS',parent['index'],a,side,len(parts),flush=True)
    paths={Path(__file__),Path(P.__file__),Path(A.__file__),Path(A.P.__file__),Path(A.P.I.__file__)}
    paths.update(ROOT/p for p in prior['sources'])
    paths.update(ROOT/f'{folder}/{A.TAG}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    sq=[r['proof']['both_R_signs_old52_product_square']for r in records]
    out=A.P.I.enc(dict(status='PASS_M_STAR_RETRIANGULATION_NATIVE52_POLAR_PREFLIGHT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        interval_scale=A.P.I.PREC,domains=domains,records=records,topology=topology(prior['topology']),
        exact_deletion_debit=prior['exact_deletion_debit'],deleted_original_walls=prior['deleted_walls'],
        removed_new_wall=71,retained_new_walls_before_R=[67,68,69,73],
        summary=dict(parents=2,parent_origins=4,N_families=2,range_domains=8,attempted_boxes=attempts,
            certified_boxes=len(records),minimum_square=min(x[0]for x in sq),maximum_square=max(x[1]for x in sq),
            minimum_unit_gap=min((x[0]-1 if r['side']=='below_gap'else 1-x[1])for r,x in zip(records,sq)),
            unresolved=0,accepted_R=0,S_banks=0,new_cones=0),
        scope='Changed incidence: delete new M71 and replace its four-triangle star by LXR/LRS. Complete native high-N labels on two fixed J/L/X parents; H4 and I2 duties force a profile0 polar R, whose physical old52 product is nonnative for both signs. No S bank, full cone or new chamber/volume claim; old exact19/12/1 deletion debit unchanged.'))
    if saved is not None:assert saved==out
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists()
    out=run(json.loads(path.read_text())if args.verify else None)
    if not args.verify:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
