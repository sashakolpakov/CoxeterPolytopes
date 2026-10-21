#!/usr/bin/env python3
"""A directed disk-patch diagnostic and one-old-wall bridge scout.

No new chamber is certified.  The exact combinatorial disk identities are
separated from the floating-point native-contact search.
"""
import argparse,json,math
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_b_release52_patch_20261017 as P

ROOT=P.ROOT;B=P.B;M=B.METRIC;J=P.J
TAG='fresh_nonarithmetic_b_disk_bridge_20261017'
CYCLE=[7,47,12,8,65,64,30]
DUTIES=[(CYCLE[i],CYCLE[(i+1)%7])for i in range(7)]
ALPHABET=[2,3,4,5,6,8,10,12]

def disk_certificate():
    old=['o'+str(i)for i in CYCLE];new=list('LQM RNP'.replace(' ',''))
    triangles=[['o30','o7','L'],['o7','o47','L'],['o47','o12','M'],['o12','o8','M'],
        ['o8','o65','N'],['o65','o64','N'],['o64','o30','P'],
        ['o47','L','Q'],['o47','Q','M'],['o8','M','R'],['o8','R','N'],
        ['o64','N','P'],['o30','P','L'],['L','Q','P'],['M','Q','R'],['N','R','P'],['Q','R','P']]
    edge_counts=Counter(tuple(sorted(e))for t in triangles for e in combinations(t,2))
    boundary={tuple(sorted(('o'+str(i),'o'+str(j))))for i,j in DUTIES}
    assert {e for e,m in edge_counts.items()if m==1}==boundary
    assert all(m==(1 if e in boundary else 2)for e,m in edge_counts.items())
    links=[]
    for v in old+new:
        graph=defaultdict(set)
        for t in triangles:
            if v not in t:continue
            a,b=[w for w in t if w!=v];graph[a].add(b);graph[b].add(a)
        todo=[next(iter(graph))];seen=set()
        while todo:
            a=todo.pop()
            if a in seen:continue
            seen.add(a);todo.extend(graph[a]-seen)
        assert len(seen)==len(graph)
        degrees=sorted(map(len,graph.values()))
        if v in old:assert degrees==[1,1]+[2]*(len(graph)-2)
        else:assert degrees==[2]*len(graph)
        links.append(dict(vertex=v,degree_sequence=degrees))
    assert len(old)+len(new)-len(edge_counts)+len(triangles)==1
    a=sum((u in old)!=(v in old)for u,v in edge_counts)
    b=sum(u in new and v in new for u,v in edge_counts)
    xyz=Counter(sum(v in new for v in t)for t in triangles)
    assert (a,b,xyz[1],xyz[2],xyz[3])==(13,9,7,6,4)
    assert 2*(a+b)-3*sum(xyz.values())==-7
    return dict(old_boundary=old,new_vertices=new,triangles=triangles,vertex_links=links,
                Euler_characteristic=1,counts=dict(a=a,b=b,x=xyz[1],y=xyz[2],z=xyz[3]),delta_ends=-7,
                status='PASS_EXACT_ABSTRACT_TRIANGULATED_DISK_ONLY')

def duty_mask(y,q,ns):
    mask=[]
    for e,(i,j)in enumerate(DUTIES):
        frame=np.array([ns[i],ns[j]]);g=frame@M@frame.T;cc=frame@M@y
        assert np.linalg.eigvalsh(g)[0]>1e-8
        if q-cc@np.linalg.inv(g)@cc>1e-7:mask.append(e)
    return tuple(mask)

def seed_analysis(bank,ns):
    groups=defaultdict(list);masks={}
    for c in bank['candidates']:
        mask=duty_mask(np.array(c['projection']),bank['profiles'][c['profile']]['square'],ns)
        if mask:groups[mask].append(c['index']);masks[c['index']]=mask
    multi={i for i,mask in masks.items()if len(mask)>1}
    finite_multi=[p for p in bank['finite_pairs']if all(i in multi for i in p['indices'])]
    assert not finite_multi
    # The four-sided seed layout L,M,N,P requires finite L--M and M--N
    # contacts. Neither exists in this complete saved finite-pair graph.
    assert len(groups[(0,6)])==3 and len(groups[(1,2)])==136 and len(groups[(3,4)])==132
    return groups,dict(duties=[list(x)for x in DUTIES],
        endpoint_mask_counts=[[list(k),len(v)]for k,v in sorted(groups.items())],
        multi_endpoint_seed_count=len(multi),finite_pairs_between_multi_endpoint_seeds=finite_multi,
        four_cap_2221_scope='Only the 2,2,2,1 boundary-arc layout with its three multiple-endpoint seeds from the frozen finite sample is screened here. No all-label or arbitrary four-cap exclusion.',
        six_cap_layout=dict(boundary_order=['L','Q','M','R','N','P'],
            old_neighbors={'L':[30,7,47],'Q':[47],'M':[47,12,8],'R':[8],'N':[8,65,64],'P':[64,30]},
            new_edges=[['L','Q'],['Q','M'],['M','R'],['R','N'],['N','P'],['P','L'],['Q','R'],['R','P'],['P','Q']],
            new_triples=[['L','Q','P'],['M','Q','R'],['N','R','P'],['Q','R','P']],
            counts=dict(f=6,a=13,b=9,x=7,y=6,z=4,delta_ends=-7)),
        disk_identity='For a triangulated cap disk with q old boundary edges and f new dual vertices, x=q, y=a-q, z=2f+q-2-a and b=3f+q-3-a; hence 2(a+b)-3(x+y+z)=-q. This assumes a simple no-ideal disk patch, retained old vertices/edges, no extra intersections or holes. It is not an existence proof.',
        two_old_neighbor_necessity='Every cap with at most two retained old intrinsic neighbors supplies at most one boundary endpoint, so f>=q=7 is necessary to close this entire seven-duty disk. This does not exclude other topologies or partial net improvement.')

def native(y,pi,i,labels,ns,profiles):
    return y@M@ns[i]+profiles[pi]['offsets'][labels.index(i)]

def bridges(L,mc,data):
    z,labels,ns,h,profiles=data;retained=[i for i in labels if i!=52]
    ly=np.array(L['projection']);my=np.array(mc['projection'])
    lp=profiles[L['profile']];mp=profiles[mc['profile']]
    lc=np.array(lp['contacts']);mct=np.array(mp['contacts'])
    frame=np.array([ns[47],ly,my]);gram=frame@M@frame.T
    counts=Counter();out=[]
    lm=float(ly@M@my+lc@J@mct);lmrel=B.relation(lm)
    if lmrel is None or lmrel=='parallel_or_unresolved':return out,Counter(incompatible_seeds=1)
    if abs(np.linalg.det(gram))<1e-10:return out,Counter(singular_frame=1)
    inv=np.linalg.inv(gram);k=np.linalg.svd(frame@M)[2][-1];kq=k@M@k
    if abs(kq)<1e-10:return out,Counter(null_kernel=1)
    desired=-np.array([[B.theta(m)for m in ms]for ms in product(ALPHABET,repeat=3)])
    mss=list(product(ALPHABET,repeat=3));old=np.array([ns[i]for i in retained]);ix=[labels.index(i)for i in retained]
    for pi,pr in enumerate(profiles):
        c=np.array(pr['contacts']);off=np.array(pr['offsets']);q=pr['square']
        shifts=np.array([off[labels.index(47)],c@J@lc,c@J@mct])
        cc=desired-shifts;p=cc@inv@frame
        finite=np.ones(len(cc),dtype=bool)
        for inds in [(0,1),(0,2)]:
            gg=gram[np.ix_(inds,inds)]
            if np.linalg.eigvalsh(gg)[0]<=1e-8:finite[:]=False;break
            finite &= q-np.einsum('ij,jk,ik->i',cc[:,inds],np.linalg.inv(gg),cc[:,inds])>1e-8
        rem=(q-np.einsum('ij,jk,ik->i',p,M,p))/kq
        good=np.where(finite&(rem>1e-9))[0];counts['real_finite_endpoint_branches']+=2*len(good)
        for sign in[-1,1]:
            yy=p[good]+sign*np.sqrt(rem[good])[:,None]*k
            products=yy@M@old.T+off[ix]
            products[:,retained.index(47)]=desired[good,0]
            centers=yy@M@h+pr['center_offset']
            for local in np.where((centers<-1e-8)&(products.max(axis=1)<=2e-7))[0]:
                rel=[B.relation(v)for v in products[local]]
                if any(r is None or r=='parallel_or_unresolved'for r in rel):counts['native_gap']+=1;continue
                y=yy[local];ms=mss[good[local]]
                assert abs(y@M@y-q)<1e-5
                # Check both prescribed new contacts independently after solve.
                assert abs(y@M@ly+c@J@lc+B.theta(ms[1]))<1e-5
                assert abs(y@M@my+c@J@mct+B.theta(ms[2]))<1e-5
                out.append(dict(left_seed=L['index'],middle_seed=mc['index'],profile=pi,labels=list(ms),
                    branch=sign,projection=y.tolist(),finite_old_contacts=[[i,r]for i,r in zip(retained,rel)if isinstance(r,int)],
                    old_center_product=float(centers[local]),seed_mutual_relation=lmrel))
    return out,counts

def right_completion(left,ngroup,bank,data):
    """Necessary P endpoint cap for the full six-cap disk, before R coupling.

    P must meet old30/64 and L, then N.  The three first contacts determine
    at most two normals.  Only the assigned finite label lists are exhausted.
    """
    z,labels,ns,h,profiles=data;retained=[i for i in labels if i!=52]
    ly=np.array(left['projection']);lc=np.array(profiles[left['profile']]['contacts'])
    frame=np.array([ns[64],ns[30],ly]);gram=frame@M@frame.T;inv=np.linalg.inv(gram)
    k=np.linalg.svd(frame@M)[2][-1];kq=k@M@k
    longlabels=list(range(2,33))+[48,64,96,128]
    mss=list(product(longlabels,[2,3,4,5],longlabels))
    desired=-np.array([[B.theta(m)for m in ms]for ms in mss])
    old=np.array([ns[i]for i in retained]);ix=[labels.index(i)for i in retained]
    counts=Counter();out=[];joins=[]
    for pi,pr in enumerate(profiles):
        c=np.array(pr['contacts']);off=np.array(pr['offsets']);q=pr['square']
        cc=desired-np.array([off[labels.index(64)],off[labels.index(30)],c@J@lc])
        pp=cc@inv@frame;finite=np.ones(len(cc),dtype=bool)
        for inds in [(0,1),(1,2)]:
            gg=gram[np.ix_(inds,inds)];assert np.linalg.eigvalsh(gg)[0]>1e-8
            finite &= q-np.einsum('ij,jk,ik->i',cc[:,inds],np.linalg.inv(gg),cc[:,inds])>1e-8
        rem=(q-np.einsum('ij,jk,ik->i',pp,M,pp))/kq
        good=np.where(finite&(rem>1e-9))[0];counts['real_finite_endpoint_branches']+=2*len(good)
        for sign in[-1,1]:
            yy=pp[good]+sign*np.sqrt(rem[good])[:,None]*k
            vv=yy@M@old.T+off[ix]
            vv[:,retained.index(64)]=desired[good,0];vv[:,retained.index(30)]=desired[good,1]
            center=yy@M@h+pr['center_offset']
            for local in np.where((center<-1e-8)&(vv.max(axis=1)<=2e-7))[0]:
                rel=[B.relation(v)for v in vv[local]]
                if any(r is None or r=='parallel_or_unresolved'for r in rel):counts['native_gap']+=1;continue
                y=yy[local];rec=dict(index=len(out),left_seed=left['index'],profile=pi,labels=list(mss[good[local]]),
                    branch=sign,projection=y.tolist(),finite_old_contacts=[[i,r]for i,r in zip(retained,rel)if isinstance(r,int)])
                out.append(rec)
                for ni in ngroup:
                    nc=bank['candidates'][ni];ny=np.array(nc['projection']);ncontact=np.array(profiles[nc['profile']]['contacts'])
                    v=float(y@M@ny+c@J@ncontact);relation=B.relation(v);counts['N_relation_'+str(relation)]+=1
                    if not isinstance(relation,int):continue
                    projected=np.array([ns[64],ny,y]);cg=projected@M@projected.T
                    if np.linalg.eigvalsh(cg)[0]<=1e-8:counts['N_nonfinite_corner']+=1;continue
                    joins.append(dict(P=rec['index'],N_seed=ni,mutual_label=relation,product=v))
    return dict(left_seed=left['index'],frame=['old64','old30','L'],label_lists=[longlabels,[2,3,4,5],longlabels],
        counts=dict(counts),candidates=out,necessary_P_N_joins=joins,
        scope='Necessary P-side condition only, labels as explicitly listed; nonrightA2 profiles included. Failure blocks this sampled six-cap template, not all labels or other patch graphs.')

def compute():
    data=P.setup();z,labels,ns,h,profiles=data
    controlpath=ROOT/'results/fresh_nonarithmetic_replace47_close12_20261017.json'
    boundarypath=ROOT/'results/fresh_nonarithmetic_replace47_close12_boundary_20261017.json'
    costpath=ROOT/'results/fresh_nonarithmetic_b_release_costs_20261017.json'
    gatepath=ROOT/'results/fresh_nonarithmetic_b_release_incidence_20261017.json'
    control=json.loads(controlpath.read_text());boundary=json.loads(boundarypath.read_text())['after64']
    costs=json.loads(costpath.read_text());gate=json.loads(gatepath.read_text())
    assert control['total_wall_count']==64 and control['complete_finite_vertex_count']==198
    assert boundary['open_ends']==15 and costs['finite_count']==198
    assert gate['after_deletion']==dict(S6_geodesic_edges=39,S6_finite_vertices=20,nonfinite_end_incidences=18)
    assert gate['minimum_change']==-3
    bankpath=ROOT/'results'/(P.TAG+'.json');bank=json.loads(bankpath.read_text())
    groups,analysis=seed_analysis(bank,ns);counts=Counter();rows=[];seen=set()
    for li in groups[(0,6)]:
        for mi in groups[(1,2)]:
            out,cc=bridges(bank['candidates'][li],bank['candidates'][mi],data);counts.update(cc)
            for row in out:
                key=(li,mi,row['profile'])+tuple(np.round(row['projection'],7))
                if key not in seen:seen.add(key);row['index']=len(rows);rows.append(row)
        print('LEFT_SEED_COMPLETE',li,'bridges',len(rows),'counts',dict(counts),flush=True)
    actual=[i for i in z['actual_sides']if i!=52];batches=[]
    for row in rows:
        L=bank['candidates'][row['left_seed']];mc=bank['candidates'][row['middle_seed']]
        yy=[L['projection'],row['projection'],mc['projection']]
        ledger=B.ledger(np.array([ns[i]for i in actual]+yy),actual+[66,67,68])
        batches.append(dict(bridge=row['index'],boundary=ledger))
    batches.sort(key=lambda r:(r['boundary']['open_ends'],r['boundary']['exterior'],-r['boundary']['finite']))
    right=[]
    for li in sorted({r['left_seed']for r in rows}):
        r=right_completion(bank['candidates'][li],groups[(3,4)],bank,data);right.append(r)
        print('RIGHT_NECESSITY_COMPLETE',li,'P',len(r['candidates']),'PN',len(r['necessary_P_N_joins']),flush=True)
    pq=[]
    for r in right:
        histogram=Counter();joins=[];minimum=math.inf;maximum=-math.inf
        for pc in r['candidates']:
            py=np.array(pc['projection']);cp=np.array(profiles[pc['profile']]['contacts'])
            for qc in rows:
                if qc['left_seed']!=r['left_seed']:continue
                qy=np.array(qc['projection']);cq=np.array(profiles[qc['profile']]['contacts'])
                value=float(py@M@qy+cp@J@cq);rel=B.relation(value)
                histogram[str(rel)]+=1;minimum=min(minimum,value);maximum=max(maximum,value)
                if rel is not None and rel!='parallel_or_unresolved':joins.append([pc['index'],qc['index'],rel])
        pq.append(dict(left_seed=r['left_seed'],relation_counts=dict(histogram),compatible_pairs=joins,
                       numerical_product_range=[minimum,maximum]))
    paths=[Path(__file__),Path(P.__file__),bankpath,controlpath,boundarypath,costpath,gatepath]
    return dict(status='BOUNDED_NUMERICAL_ONE_OLD_NEIGHBOR_DISK_BRIDGE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_wall=52,control=dict(walls=64,finite_vertices=198,S6_open_ends=15),
        deletion_cost=dict(finite_vertices=5,compact_polygons=[8,12,64,65],higher_compact_faces=0),
        topology=analysis,abstract_disk_certificate=disk_certificate(),alphabet=ALPHABET,profiles=profiles,counts=dict(counts),bridges=rows,batches=batches,right_necessary_stage=right,
        necessary_P_Q_compatibility=pq,
        scope='Directed six-cap disk template. Its previously missing Q bridge is tested: one old anchor47, two new multiple-endpoint seedsL/M, all seven A2 profiles, sampled labels as saved, all63 retained native gates/common old strict center and both required finite endpoint triples. Then necessary P-side joins to L and N are tested before any R solve. No six-cap completion, all-label nonexistence, independent arithmetic replay, or finite-volume claim follows.')

def main():
    a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');args=a.parse_args()
    out=json.loads(json.dumps(compute()));p=ROOT/'results'/(TAG+'.json')
    if args.verify:assert json.loads(p.read_text())==out
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['bridges']),dict(out['counts']),flush=True)
    if out['batches']:print('BEST',out['batches'][0],flush=True)
if __name__=='__main__':main()
