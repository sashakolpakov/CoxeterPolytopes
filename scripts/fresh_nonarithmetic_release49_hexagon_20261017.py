#!/usr/bin/env python3
"""Release49 and target its entire six-sided A-region with three pentagons.

Exact old incidence/deletion and abstract disk checks; one bounded numerical
search in the unchanged exact64-wall control coordinates. No promotion.
"""
import argparse
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import numpy as np
import mpmath as mp

import fresh_nonarithmetic_release47_cycle_patch_20261017 as A

ROOT,B,M,J=A.ROOT,A.B,A.M,A.J
TAG='fresh_nonarithmetic_release49_hexagon_20261017'
PARENT=ROOT/'results/fresh_nonarithmetic_replace47_close12_20261017.json'
BOUNDARY=ROOT/'results/fresh_nonarithmetic_replace47_close12_boundary_20261017.json'
CYCLE=[7,33,30,59,60,58]
ALPHABET=list(range(2,13))+[16,24,32]
mp.mp.dps=70


def mp_interval_mid(x):return (mp.mpf(x[0])+mp.mpf(x[1]))/2


def precise_space(z):
    g=mp.matrix([[mp_interval_mid(x)for x in row]for row in z['residual_gram']])
    ns={r['wall']:mp.matrix([mp_interval_mid(x)for x in r['coordinates']])for r in z['projections']}
    h=mp.matrix([mp_interval_mid(x)for x in z['projected_common_center']]);pair=lambda x,y:(x.T*g*y)[0]
    h=h/mp.sqrt(-pair(h,h))
    gg=np.array([[float(g[i,j])for j in range(4)]for i in range(4)])
    vals,vecs=np.linalg.eigh((gg+gg.T)/2);transform=vecs@np.diag(1/np.sqrt(abs(vals)))
    return g,ns,h,pair,transform


def topology():
    old=list('abcdef');new=list('LMN')
    triangles=[list(s)for s in ['abL','bcL','cdM','deM','efN','faN','cLM','eMN','aNL','LMN']]
    edges=Counter(tuple(sorted(e))for t in triangles for e in combinations(t,2))
    boundary={tuple(sorted(e))for e in ['ab','bc','cd','de','ef','fa']}
    assert {e for e,n in edges.items()if n==1}==boundary
    assert all(n==(1 if e in boundary else 2)for e,n in edges.items())
    links={}
    for v in old+new:
        g=defaultdict(set)
        for t in triangles:
            if v in t:
                i,j=[x for x in t if x!=v];g[i].add(j);g[j].add(i)
        seen=set();todo=[next(iter(g))]
        while todo:
            x=todo.pop()
            if x not in seen:seen.add(x);todo.extend(g[x]-seen)
        assert seen==set(g)
        assert sorted(map(len,g.values()))==([1,1]+[2]*(len(g)-2)if v in old else[2]*len(g))
        links[v]={x:sorted(y)for x,y in sorted(g.items())}
    assert len(old+new)-len(edges)+len(triangles)==1
    assert all(len(links[v])==5 for v in new)
    return dict(old_cycle=CYCLE,triangles=triangles,links=links,vertices=9,edges=18,triangles_count=10,
                new_facet_degrees=[5,5,5],old_new_edges=9,new_new_edges=3,
                finite_triangle_types=dict(two_old=6,one_old=3,all_new=1),
                scope='Conditional genuine six-sided disk with ten finite vertices; no geometric existence assertion.')


def deletion_cost():
    d=json.loads(PARENT.read_text());bd=json.loads(BOUNDARY.read_text())['after64']
    vertices=sorted(tuple(v)for v in d['retained_old_finite_supports'])+sorted(tuple(v['support'])for v in d['new_finite_vertices'])
    vertices=sorted(vertices);assert len(vertices)==len(set(vertices))==198
    edges=defaultdict(list);faces=defaultdict(list)
    for i,v in enumerate(vertices):
        assert len(v)==9
        for e in combinations(v,8):edges[e].append(i)
        for size in range(8):
            for s in combinations(v,size):faces[s].append(i)
    assert all(len(v)<=2 for v in edges.values())
    lost={i for i,v in enumerate(vertices)if 49 in v};touched=[]
    for s,vs in faces.items():
        if not lost.intersection(vs)or len(vs)<10-len(s):continue
        if all(all(len(edges[tuple(j for j in vertices[i]if j!=k)])==2 for k in vertices[i]if k not in s)for i in vs):
            touched.append(dict(support=list(s),dimension=9-len(s),finite_vertices=[list(vertices[i])for i in vs]))
    incident=[e for e in bd['edges']if 49 in e['sides']]
    rim={tuple(sorted(x for x in r['active']if x!=49))for r in bd['rays']if 49 in r['active']}
    assert rim=={tuple(sorted([CYCLE[i],CYCLE[(i+1)%6]]))for i in range(6)}
    fv=[r['active']for r in bd['rays']if 49 in r['active']and r['kind']=='finite_vertex']
    exterior=[r['active']for r in bd['rays']if 49 in r['active']and r['kind']=='exterior_ray']
    assert fv==[[7,49,58]] and len(exterior)==5
    return dict(deleted_wall=49,lost_finite_count=len(lost),lost_finite_supports=[list(vertices[i])for i in sorted(lost)],
        compact_faces_touched=touched,compact_faces_touched_by_dimension=dict(Counter(r['dimension']for r in touched)),
        S6_link_cycle=CYCLE,S6_finite_corners=fv,S6_exterior_corners=exterior,
        S6_deleted_side_edges=[dict(sides=e['sides'],kind=e['kind'])for e in incident],
        S6_compact_side_faces_touched=[r['wall']for r in bd['facets']if r['compact']and any(49 in v for v in r['vertex_supports'])],
        scope='Exact old-certificate loss from the complete198 finite census. Deleting a wall alone enlarges the region; fewer ledger edges alone are not finite-volume progress.')


def solve_frame(ids,data):
    z,labels,ns,h,profiles=data;retained=[i for i in labels if i!=49]
    _,precise,hh,pair,transform=precise_space(z)
    ff=[precise[i]for i in ids]
    gm=mp.matrix([[pair(x,y)for y in ff]for x in ff]);im=gm**-1
    g=np.array([[float(gm[i,j])for j in range(3)]for i in range(3)])
    counts=Counter();out=[];seen=set()
    if abs(np.linalg.det(g))<1e-10:return [],dict(singular_frame=1)
    inv=np.linalg.inv(g);kernel=None
    for j in range(4):
        e=mp.matrix([int(i==j)for i in range(4)]);vv=im*mp.matrix([pair(e,f)for f in ff])
        kk=e-sum((vv[i]*ff[i]for i in range(3)),mp.zeros(4,1))
        if abs(pair(kk,kk))>mp.mpf('1e-25'):kernel=kk;break
    assert kernel is not None
    kqm=pair(kernel,kernel);kq=float(kqm)
    tuples=list(product(ALPHABET,repeat=3));desired=-np.array([[B.theta(m)for m in row]for row in tuples])
    for pi,pr in enumerate(profiles):
        off=np.array(pr['offsets']);q=pr['square'];cc=desired-np.array([off[labels.index(i)]for i in ids])
        good=np.ones(len(tuples),dtype=bool)
        for ab in [(0,1),(1,2)]:
            gg=g[np.ix_(ab,ab)]
            if np.linalg.eigvalsh(gg)[0]<=1e-9:good[:]=False;break
            good &=q-np.einsum('ij,jk,ik->i',cc[:,ab],np.linalg.inv(gg),cc[:,ab])>1e-8
        rem=(q-np.einsum('ij,jk,ik->i',cc,inv,cc))/kq
        ix0=np.where(good&(rem>1e-9))[0];counts['real_finite_endpoint_branches']+=2*len(ix0)
        for row in ix0:
            ms=tuples[row]
            cm=mp.matrix([-2*mp.cos(mp.pi/m)-mp.mpf(str(off[labels.index(i)]))for i,m in zip(ids,ms)])
            vv=im*cm;part=sum((vv[i]*ff[i]for i in range(3)),mp.zeros(4,1))
            qm=mp.mpf(str(q));residual=(qm-pair(part,part))/kqm
            assert residual>0
            for sign in [-1,1]:
                root=part+sign*mp.sqrt(residual)*kernel
                assert abs(pair(root,root)-qm)<mp.mpf('1e-40')
                native=[float(pair(root,precise[i])+mp.mpf(str(off[labels.index(i)])))for i in retained]
                for j,i in enumerate(ids):native[retained.index(i)]=-B.theta(ms[j])
                center=float(pair(root,hh)+mp.mpf(str(pr['center_offset'])))
                if center>=-1e-8 or max(native)>=2e-7:continue
                counts['acute_common_center_branches']+=1
                relations=[B.relation(v)for v in native]
                if any(r is None or r=='parallel_or_unresolved'for r in relations):counts['native_angle_gap']+=1;continue
                y=np.linalg.solve(transform,np.array([float(x)for x in root]));key=(pi,)+tuple(np.round(y,7))
                if key in seen:counts['duplicate']+=1;continue
                seen.add(key)
                out.append(dict(index=len(out),frame=ids,profile=pi,labels=list(ms),branch=sign,
                    projection=y.tolist(),residual_coordinates=[mp.nstr(x,65)for x in root],
                    expanded_double_norm_error=float(y@M@y-q),center_product=center,
                    finite_old_contacts=[[i,r]for i,r in zip(retained,relations)if isinstance(r,int)],
                    deleted49_product=float(pair(root,precise[49])+mp.mpf(str(off[labels.index(49)])))))
    return out,dict(counts)


def mutual(a,b,profiles):
    return A.native(np.array(a['projection']),a['profile'],np.array(b['projection']),b['profile'],profiles)


def run():
    data=A.setup();z,labels,ns,h,profiles=data
    baseline=B.ledger(np.array([ns[i]for i in z['actual_sides']]),z['actual_sides'])
    assert (baseline['finite'],baseline['exterior'],baseline['open_ends'])==(25,7,15)
    actual=[i for i in z['actual_sides']if i!=49]
    after_delete=B.ledger(np.array([ns[i]for i in actual]),actual)
    costs=deletion_cost();disk=topology()
    print('RELEASE49_CONTROL',baseline['finite'],baseline['open_ends'],'DELETE',after_delete['finite'],after_delete['exterior'],after_delete['open_ends'],'lostambient',costs['lost_finite_count'],'compactcost',costs['compact_faces_touched_by_dimension'],flush=True)
    stages=[];complete=[];improved=[];hist=Counter()
    for parity in [0,1]:
        cycle=CYCLE[parity:]+CYCLE[:parity]
        frames=[[cycle[(s+i)%6]for i in range(3)]for s in [0,2,4]]
        pools=[];pool_counts=[]
        for frame in frames:
            pool,count=solve_frame(frame,data);pools.append(pool);pool_counts.append(count)
            print('HEXAGON_POOL',parity,frame,len(pool),count,flush=True)
        pairs={};pair_counts={};relations={}
        for a,b,shared in [(0,1,cycle[2]),(1,2,cycle[4]),(0,2,cycle[0])]:
            edge_set=set();counts=Counter();relation={}
            for x in pools[a]:
                for y in pools[b]:
                    g=mutual(x,y,profiles);r=B.relation(g);counts[str(r)]+=1
                    if not isinstance(r,int):continue
                    ff=np.array([ns[shared],np.array(x['projection']),np.array(y['projection'])])
                    least=float(np.linalg.eigvalsh(ff@M@ff.T)[0])
                    if least<=1e-8:counts['nonfinite_shared_old_corner']+=1;continue
                    if x['profile']==y['profile']and np.max(abs(np.array(x['projection'])-np.array(y['projection'])))<1e-7:
                        counts['duplicate_root']+=1;continue
                    edge_set.add((x['index'],y['index']));relation[(x['index'],y['index'])]=dict(label=r,product=g,finite_corner_min_eigenvalue=least)
            pairs[a,b]=edge_set;pair_counts[f'{a},{b}']=dict(counts);relations[a,b]=relation
            print('HEXAGON_PAIRS',parity,a,b,len(edge_set),dict(counts),flush=True)
        triplets=[];trianglecounts=Counter()
        for ia,ib in sorted(pairs[0,1]):
            for ic in range(len(pools[2])):
                if (ib,ic)not in pairs[1,2]or(ia,ic)not in pairs[0,2]:continue
                trianglecounts['three_finite_join_edges']+=1
                roots=[pools[i][j]for i,j in enumerate([ia,ib,ic])]
                vec=[np.array(x['projection'])for x in roots];ff=np.array(vec)
                least=float(np.linalg.eigvalsh(ff@M@ff.T)[0])
                if least<=1e-8:trianglecounts['nonfinite_central_corner']+=1;continue
                boundary=B.ledger(np.array([ns[i]for i in actual]+vec),actual+[66,67,68])
                score=dict(indices=[ia,ib,ic],central_min_eigenvalue=least,boundary=boundary,
                    mutual=[relations[a,b][u,v]for a,b,u,v in [(0,1,ia,ib),(1,2,ib,ic),(0,2,ia,ic)]])
                triplets.append(score);complete.append(dict(parity=parity,roots=roots,**score))
                hist[boundary['finite'],boundary['exterior'],boundary['open_ends']]+=1
                if boundary['open_ends']<15:improved.append(dict(parity=parity,roots=roots,**score))
                print('HEXAGON_DISK',parity,[ia,ib,ic],boundary['finite'],boundary['exterior'],boundary['open_ends'],flush=True)
        stages.append(dict(parity=parity,cycle=cycle,frames=frames,pools=pools,pool_counts=pool_counts,
            pair_counts=pair_counts,pair_finite_edges={f'{a},{b}':sorted(map(list,es))for(a,b),es in pairs.items()},
            triple_counts=dict(trianglecounts),complete_disk_candidates=triplets))
    paths=[Path(__file__),Path(A.__file__),Path(A.C.__file__),Path(B.__file__),PARENT,BOUNDARY,A.INPUT,A.C.PARENT]
    return dict(status='BOUNDED_NUMERICAL_RELEASE49_SIX_SIDED_THREE_CAP_PATCH',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deletion_cost=costs,abstract_disk=disk,baseline=baseline,after_delete49=after_delete,
        finite_alphabet=ALPHABET,white_profiles=profiles,stages=stages,
        complete_boundary_histogram=[dict(finite=f,exterior=e,open_ends=n,count=c)for(f,e,n),c in sorted(hist.items())],
        complete_disks=complete,improvements=improved,
        scope='Exact old incidence/deletion and abstract triangulated hexagon; bounded numerical discovery only. Native root construction and products use70-digit midpoint arithmetic; preliminary filters and boundary scoring use double precision. Both alternating placements and all7 A2 profiles, both norm branches and all63 retained native gates/common center are checked at the stated finite labels. Every complete ten-corner disk candidate receives the full numerical S6 boundary. This scout does not certify high-label tails or other three-cap topologies; no numerical geometry, ambient finite census or finite-volume assertion is promoted.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'complete',len(out['complete_disks']),'improved',len(out['improvements']),flush=True)
