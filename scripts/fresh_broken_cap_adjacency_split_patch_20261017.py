#!/usr/bin/env python3
"""Two quadrilateral disks sharing the new X side: bounded numerical search."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_broken_cap_adjacency_events_20261017 as E

T,S,H,B,G,M,ROOT=E.T,E.S,E.H,E.B,E.G,E.M,E.ROOT
TAG='fresh_broken_cap_adjacency_split_patch_20261017'
X=69

def finite3(pi,frame,ms,anchor_profiles,old_pairs):
    ans=[]
    for j,k in[(0,1),(1,2)]:
        g=[[2]*5 for _ in range(5)]
        for i in range(5):g[i][i]=1
        g[0][1]=g[1][0]=3
        for col,node in[(2,frame[j]),(3,frame[k])]:
            for i,m in enumerate(G.profile_labels(anchor_profiles[node])):g[i][col]=g[col][i]=m
        for i,m in enumerate(G.profile_labels(pi)):g[i][4]=g[4][i]=m
        g[2][3]=g[3][2]=old_pairs[frozenset((frame[j],frame[k]))]
        g[2][4]=g[4][2]=ms[j];g[3][4]=g[4][3]=ms[k]
        v=T.N.A.components(tuple(v for row in g for v in row))
        if v is None:return None
        ans.append(v)
    return ans

def pair_corner(anchor,a,b,m,anchor_profiles):
    g=[[2]*5 for _ in range(5)]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    for i,z in enumerate(G.profile_labels(anchor_profiles[anchor])):g[i][2]=g[2][i]=z
    for col,cap in[(3,a),(4,b)]:
        for i,z in enumerate(G.profile_labels(cap['profile'])):g[i][col]=g[col][i]=z
        z=cap['labels'][cap['frame'].index(anchor)];g[2][col]=g[col][2]=z
    g[3][4]=g[4][3]=m
    return T.N.A.components(tuple(v for row in g for v in row))

def run():
    data=json.loads((ROOT/'results'/f'{E.TAG}.json').read_text())
    z,labels,ns,h,profiles,primary=H.setup();actual=[i for i in z['actual_sides']if i not in S.DELETE]
    oldsq={r['wall']:B.mid(r['square'])for r in z['projections']}
    anchor_profiles={7:4,8:0,30:0,44:0,63:0,X:0}
    frames={'A':[7,8,X],'C':[7,30,X],'R1':[44,8,X],'R2':[44,63,X]}
    counts=Counter();rows=[];hits=[];unresolved=[]
    for xx in data['X_bank']:
        if T.corner(0,xx['old_labels'][1],xx['old_labels'][2],0,2)is None:
            counts['X_nonfinite_old30_63_corner']+=1;continue
        xy=np.array(xx['projection']);nd=dict(ns);nd[X]=xy
        old_pairs={frozenset(k):v for k,v in[((7,8),5),((7,30),5),((44,8),2),((44,63),2),
                   ((8,X),xx['old_labels'][0]),((30,X),xx['old_labels'][1]),((63,X),xx['old_labels'][2])]}
        banks={};domain_counts={};blockers={}
        for role,frame in frames.items():
            accepted=[];dc=Counter();bc=Counter()
            for pi,ms0 in product(range(7),product(range(2,7),repeat=3)):
                ty=finite3(pi,frame,ms0,anchor_profiles,old_pairs)
                if ty is None:continue
                high=[j for j,m in enumerate(ms0)if m==6];dc['high_shapes'if high else'low_shapes']+=1
                for hh in product(T.HIGH_LABELS,repeat=len(high)):
                    ms=list(ms0)
                    for j,m in zip(high,hh):ms[j]=m
                    offsets=[H.A.native(np.zeros(4),pi,np.zeros(4),anchor_profiles[i],profiles)for i in frame]
                    rhs=[-B.theta(m)-o for m,o in zip(ms,offsets)]
                    roots,norm=T.N.norm_roots([nd[i]for i in frame],rhs,profiles[pi]['square']);dc['norm_equations']+=1
                    if norm['status'].startswith('unresolved'):
                        unresolved.append(dict(X=xx['index'],role=role,profile=pi,labels=ms,norm=norm));continue
                    for sg,cy in roots:
                        dc['real_branches']+=1
                        vals=T.native_values(cy,pi,labels,ns,profiles)
                        vals[X]=H.A.native(cy,pi,xy,0,profiles)
                        bad=T.gates(vals,frame);center=float(cy@M@h+profiles[pi]['center_offset'])
                        if center>=-1e-8:continue
                        dc['future_branches']+=1;bc[tuple(i for i,v in bad)]+=1
                        if bad:continue
                        assert abs(cy@M@cy-profiles[pi]['square'])<2e-6
                        assert max(abs(float(cy@M@nd[i])-r)for i,r in zip(frame,rhs))<2e-6
                        accepted.append(dict(index=len(accepted),role=role,frame=frame,profile=pi,labels=ms,symbolic_high_positions=high,
                                             branch=sg,projection=cy.tolist(),center_product=center))
                        dc['native_coherent_caps']+=1
            banks[role]=accepted;domain_counts[role]=dict(dc);blockers[role]=[[list(k),v]for k,v in sorted(bc.items())]
            counts.update(dc)
        joined={};relations={}
        for tag,ra,rb,anchors in[('left','A','C',[7,X]),('right','R1','R2',[44,X])]:
            jj=[];rc=Counter()
            for a,b in product(banks[ra],banks[rb]):
                value=H.A.native(np.array(a['projection']),a['profile'],np.array(b['projection']),b['profile'],profiles)
                rel=B.relation(value);rc[str(rel)]+=1
                if not isinstance(rel,int):continue
                types=[pair_corner(i,a,b,rel,anchor_profiles)for i in anchors]
                if any(t is None for t in types):continue
                jj.append(dict(indices=[a['index'],b['index']],product=value,label=rel,finite_corner_types=types))
            joined[tag]=jj;relations[tag]=dict(rc);counts[f'{tag}_coherent_pairs']+=len(jj)
        for left,right in product(joined['left'],joined['right']):
            caps=[banks['A'][left['indices'][0]],banks['C'][left['indices'][1]],banks['R1'][right['indices'][0]],banks['R2'][right['indices'][1]]]
            cross=[];valid=True
            for i,j in product(range(2),range(2,4)):
                a,b=caps[i],caps[j];v=H.A.native(np.array(a['projection']),a['profile'],np.array(b['projection']),b['profile'],profiles);rel=B.relation(v)
                cross.append(dict(roles=[a['role'],b['role']],product=v,relation=rel))
                if rel is None or rel=='parallel_or_unresolved':valid=False
            counts['joint_pair_trials']+=1
            if not valid:continue
            nd2=dict(nd);ids=[70,71,72,73]
            nd2.update({i:np.array(c['projection'])for i,c in zip(ids,caps)})
            bd=S.C.detailed_cone(nd2,actual+[X]+ids)
            inc=H.incidence_score([nd2[i]for i in actual+[X]+ids],actual+[X]+ids,
                 [oldsq[i]for i in actual]+[2]+[profiles[c['profile']]['square']for c in caps])
            hits.append(dict(X=xx['index'],left=left,right=right,caps=caps,cross_contacts=cross,boundary=bd,PD_incidence=inc))
            counts['native_coherent_whole_patches']+=1
        rows.append(dict(X=xx['index'],X_labels=xx['old_labels'],domain_counts=domain_counts,blockers=blockers,
                         coherent_cap_banks=banks,pair_relations=relations,coherent_pairs=joined))
        print('SPLIT_PATCH_X',xx['index'],'banks',{k:len(v)for k,v in banks.items()},'pairs',{k:len(v)for k,v in joined.items()},'hits',len(hits),flush=True)
    paths=[Path(__file__),Path(E.__file__),Path(T.__file__),Path(S.__file__),ROOT/'results'/f'{E.TAG}.json',H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_NUMERICAL_SHARED_X_TWO_QUADRILATERAL_DISK_SEARCH',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},deleted_walls=S.DELETE,
        boundary_cycle=S.CYCLE,exact_old_deletion_cost=data['exact_old_deletion_cost'],
        topology=dict(X_old_contacts=[8,30,63],X_old_edge_covered=[30,63],quadrilateral_disks=[[7,8,'X',30],[8,44,63,'X']],
            cap_frames=frames,central_triangles=[['A','C','X'],['R1','R2','X']],vertices=10,edges=22,triangles=13,x=5,y=6,z=2),
        counts=dict(counts,X_tested=len(rows)),records=rows,whole_patch_candidates=hits,unresolved_events=unresolved,
        scope='Numerical search in the57 X contact-event roots, retaining only finite S6+30+63+X. All finite-corner shape recipes and both norm branches are included, but isolated high labels are sampled6,10,24. All retained native contacts/common center and all cross-pair products are tested. Complete intrinsic boundary and independent PD-incidence diagnostics are required for any surviving whole patch. No exact promotion, permanent-cusp completeness or ambient finite-volume theorem is asserted.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(run()));dest=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert out==json.loads(dest.read_text())
    else:dest.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],flush=True)
