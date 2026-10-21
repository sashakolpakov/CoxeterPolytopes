#!/usr/bin/env python3
"""A remote bridge closes old60 after three bounded local bridge failures.

Complete declared low recipe cubes, native-first. All28 frozen T caps are
intersected with every surviving remote normal by the complete old cone's
edge clipping rule. Fifteen independent direct cone enumerations corroborate
one full T parent. New normals and cones remain numerical.
"""
import argparse,json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_oblique_n_r_s_polar_cap_20261018 as B
ROOT=B.ROOT;T=B.T;C=B.C;P=B.P;E=B.E;Q=B.Q
TAG='fresh_remote60_compact_face_bridge_20261018'
WALL=75

def remote_topology():
    tri=[[58,59,60],[58,60,WALL],[59,60,WALL]];boundary=[58,59,WALL]
    ee=Counter(tuple(sorted(e))for f in tri for e in combinations(f,2))
    assert len(ee)==6 and{e for e,n in ee.items()if n==1}=={tuple(sorted(e))for e in combinations(boundary,2)}
    return dict(triangles=tri,boundary=boundary,interior_wall=60,interior_link=[58,59,WALL],V=4,E=6,F=3,euler=1,
        separate_from_main_T_disk=True,forbidden_finite_duty=[60,66,WALL])

def offset(s,pi,i,tt):
    return s.profiles[pi]['offsets'][s.index[i]]if i in s.index else s.off(pi,tt[i])

def bank(s,nd,tt,frame):
    ps=[s.ps(i,tt)for i in frame]
    bonds=[T.G.relation(s.dot(nd[frame[0]],ps[0],nd[frame[k]],ps[k]))for k in(1,2)]
    assert all(isinstance(m,int)for m in bonds)
    rows=[]
    for pi,ms in product(range(7),product(range(2,6),repeat=3)):
        duties=[T.P.corner3(ps[0],ps[k],pi,bonds[k-1],ms[0],ms[k])for k in(1,2)]
        r=dict(index=len(rows),profile=pi,frame=frame,labels=list(ms),finite_duty_types=duties,
            required_corners=[[frame[0],frame[k],WALL]for k in(1,2)],norm=None,branches=[]);rows.append(r)
        if any(x is None for x in duties):r['status']='finite_duty_rejected';continue
        rhs=[-T.B.theta(m)-offset(s,pi,i,tt)for i,m in zip(frame,ms)]
        roots,info=T.N.norm_roots([nd[i]for i in frame],rhs,s.profiles[pi]['square'])
        r['norm']=info;r['status']=info['status']
        for sign,v in roots:
            vv,bb,cp=s.native(v,pi,nd,tt)
            err=max(abs(v@T.M@v-s.profiles[pi]['square']),max(abs(v@T.M@nd[i]-q)for i,q in zip(frame,rhs)))
            duplicates=[i for i in s.cur['actual_sides']if s.ps(i,{})==pi and np.max(abs(v-s.old[i]))<2e-5]
            r['branches'].append(dict(branch=sign,normal=v.tolist(),native_products=sorted(vv.items()),
                blockers=[[i,vv[i]]for i in bb],center_product=cp,norm_residual=float(err),
                native_admissible=bool(not bb and cp<-1e-8 and err<2e-5),original_duplicates=duplicates))
    assert len(rows)==448
    return dict(frame=frame,frame_profiles=ps,frame_bonds=bonds,recipes=rows,
        summary=dict(recipes=448,status_counts=dict(Counter(r['status']for r in rows)),
            branches=sum(len(r['branches'])for r in rows),native_branches=sum(b['native_admissible']for r in rows for b in r['branches'])))

def assemble_rays(vectors,full):
    rays=[];pairs=defaultdict(list)
    for v in vectors:
        v=np.array(v);v=v/np.linalg.norm(v)
        products={i:float(w@T.M@v)/np.linalg.norm(w@T.M)for i,w in full.items()}
        assert max(products.values())<2e-8
        active=[i for i,q in products.items()if abs(q)<2e-10];assert len(active)==3
        q=float(v@T.M@v);assert abs(q)>1e-12
        r=dict(index=len(rays),active=active,coordinates=v.tolist(),square=q,kind='finite'if q<0 else'exterior');rays.append(r)
        for pair in combinations(active,2):pairs[pair].append(r['index'])
    assert len({tuple(r['active'])for r in rays})==len(rays)
    edges=[]
    for e,ii in pairs.items():
        assert len(ii)==2
        a,b=[rays[i]for i in ii];ab=float(np.array(a['coordinates'])@T.M@np.array(b['coordinates']));nf=sum(r['kind']=='finite'for r in(a,b))
        kind='compact'if nf==2 else'one_open'if nf==1 else'two_open'if ab<-1e-12 and a['square']*b['square']-ab*ab<-1e-15 else'outside'
        edges.append(dict(sides=list(e),rays=ii,kind=kind))
    return dict(rays=rays,edges=edges,warnings=[],counts=dict(Counter(r['kind']for r in rays)),
        open_ends=sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in edges))

def clip(old,v,full):
    cov=v@T.M;cov=cov/np.linalg.norm(cov)
    products=[float(cov@np.array(r['coordinates']))for r in old['rays']]
    assert min(map(abs,products))>1e-10
    keep=[i for i,q in enumerate(products)if q<0];removed=[i for i,q in enumerate(products)if q>0]
    vectors=[old['rays'][i]['coordinates']for i in keep];cross=[]
    for ei,e in enumerate(old['edges']):
        i,j=e['rays'];a,b=products[i],products[j]
        if a*b>=0:continue
        if a<0:i,j=j,i;a,b=b,a
        w=a*np.array(old['rays'][j]['coordinates'])-b*np.array(old['rays'][i]['coordinates'])
        cross.append(dict(old_edge=ei,positive_ray=i,negative_ray=j,new_ray=len(vectors)))
        vectors.append(w)
    raw=assemble_rays(vectors,full)
    return raw,dict(normalized_ray_products=products,retained_old_rays=keep,removed_old_rays=removed,crossing_edges=cross,
        complete_input_cone_inherited=True,rule='Keep negative rays and intersect every sign-changing old edge with the new wall')

def signature(raw):return sorted((tuple(sorted(r['active'])),r['kind'])for r in raw['rays'])

def run():
    read=lambda tag:json.loads((ROOT/f'results/{tag}.json').read_text())
    d=read(B.TAG);d0=read(P.TAG)
    for x in(d,d0):
        for p,h in x['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    assert P.H.H.old_frame_interface(s)==d['old_frame_interface']
    debit=T.H.deletion(s.primary,s.cur,C.DELETED);assert json.loads(json.dumps(debit))==d['exact_deletion_debit']
    parents=[];failed=[]
    for p in d0['records']:
        nd,tt,steps=E.reconstruct(s,p)
        parents.append(dict(parent=p['index'],origins=p['origins'],reconstruction=steps,
            normals={str(i):nd[i].tolist()for i in tt},profiles=tt))
        for name,frame in [('old10',[10,33,73]),('old52',[52,47,72]),('old55',[55,7,59])]:
            bb=bank(s,nd,tt,frame);assert bb['summary']['native_branches']==0
            failed.append(dict(parent=p['index'],target=name,**bb))
    base={i:v for i,v in s.old.items()if i not in C.DELETED}
    remote=bank(s,base,{},[60,58,59]);roots=[]
    for r in remote['recipes']:
        for b in r['branches']:
            if b['native_admissible']:roots.append((r,b))
    assert remote['summary']['status_counts']=={'two_real_norm_branches':15,'no_real_norm_root':2,'finite_duty_rejected':431}
    assert len(roots)==15 and all(r['profile']==0 and r['labels'][0]==2 and b['branch']==-1 for r,b in roots)
    configs=[];accepted=[]
    for ai,a in enumerate(d['accepted']):
        old=d['records'][a['configuration']];nd={**base};tt={}
        for j,ch in enumerate(old['reconstruction']):
            w=ch['wall'];v,proof=C.solve_saved(s,ch['frame'],ch['profile'],ch['labels'],ch['branch'],s.old if j<2 else nd,{}if j<2 else tt,old['normals'][str(w)])
            nd[w]=v;tt[w]=ch['profile']
        tr=old['recipes'][a['profile']];tb=next(x for x in tr['branches']if x['branch']==a['branch'])
        v,proof=C.solve_saved(s,tr['frame'],a['profile'],tr['labels'],a['branch'],nd,tt,tb['normal']);nd[74]=v;tt[74]=a['profile']
        config=dict(index=ai,source_configuration=a['configuration'],T_profile=a['profile'],T_branch=a['branch'],
            normals={str(i):nd[i].tolist()for i in tt},profiles=tt,pairs=[]);configs.append(config)
        for r,b in roots:
            v=np.array(b['normal']);vv,bb,cp=s.native(v,0,nd,tt)
            z=dict(recipe=r['index'],branch=b['branch'],native_products=sorted(vv.items()),blockers=[[i,vv[i]]for i in bb],
                center_product=cp,native_admissible=bool(not bb and cp<-1e-8),accepted_index=None);config['pairs'].append(z)
            if not z['native_admissible']:continue
            full={**nd,WALL:v};ft={**tt,WALL:0};raw,clipping=clip(a['full_projected_boundary'],v,full)
            ids=[i for i in s.cur['actual_sides']if i not in C.DELETED]+list(ft)
            gram=np.array([[full[i]@T.M@full[j]for j in ids]for i in ids])
            pd=sorted(sorted(ids[i]for i in ix)for ix in combinations(range(len(ids)),3)if np.linalg.eigvalsh(gram[np.ix_(ix,ix)]).min()>1e-8)
            finite=sorted(sorted(q['active'])for q in raw['rays']if q['kind']=='finite');assert finite==pd
            direct=None
            if ai==0:
                direct=T.S.C.detailed_cone(full,list(full));assert not direct['warnings']and signature(direct)==signature(raw)
                assert direct['open_ends']==raw['open_ends']
            oldf={tuple(sorted(q['active']))for q in a['full_projected_boundary']['rays']if q['kind']=='finite'}
            olde={tuple(sorted(q['active']))for q in a['full_projected_boundary']['rays']if q['kind']=='exterior'}
            newf={tuple(q)for q in finite};newe={tuple(sorted(q['active']))for q in raw['rays']if q['kind']=='exterior'}
            faces=[dict(wall=w,rays=[q for q in raw['rays']if w in q['active']],all_rays_finite=all(q['kind']=='finite'for q in raw['rays']if w in q['active']))for w in ids]
            closed={f['wall']for f in faces if f['rays']and f['all_rays_finite']};assert{8,11,12,60,67,68,69,70}<=closed
            assert newf-oldf=={(58,60,WALL),(59,60,WALL)}and not(oldf-newf)
            assert olde-newe=={(58,60,66),(59,60,66)}and newe-olde=={(58,66,WALL),(59,66,WALL)}
            assert raw['counts']=={'finite':27,'exterior':15}and raw['open_ends']==13
            z['accepted_index']=len(accepted)
            accepted.append(dict(configuration=ai,recipe=r['index'],branch=b['branch'],full_projected_boundary=raw,clipping=clipping,
                direct_cone=direct,PD_finite_supports=pd,whole_faces=faces,
                finite_supports_added=sorted(newf-oldf),finite_supports_removed=[],exterior_supports_removed=sorted(olde-newe),
                exterior_supports_added=sorted(newe-olde),finite_corner_types=[s.corner(f,full,ft)for f in remote_topology()['triangles']]))
        print('T_CONFIGURATION',ai,'native',sum(z['native_admissible']for z in config['pairs']),'caps',len(accepted),flush=True)
    assert len(configs)==28 and len(accepted)==420
    paths={Path(__file__)}|{ROOT/p for x in(d,d0)for p in x['sources']}
    paths.update(ROOT/f'{folder}/{tag}.{ext}'for tag in(B.TAG,P.TAG)for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    summary=dict(preflight_families=6,preflight_recipes=2688,preflight_native=0,remote_recipes=448,remote_eligible=17,
        remote_real_roots=30,remote_native=15,T_configurations=28,complete_pairs=420,accepted=420,direct_cones=15,
        finite_rays=27,exterior_rays=15,open_ends=13,new_whole_old_face=60,unresolved=0)
    return dict(status='BOUNDED_REMOTE60_WHOLE_FACE_BRIDGE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        predecessor=B.TAG,new_wall=WALL,deleted_original_walls=list(C.DELETED),exact_deletion_debit=debit,
        numerical_snapshot=d['numerical_snapshot'],old_frame_interface=d['old_frame_interface'],main_topology=d['topology'],remote_topology=remote_topology(),
        wall_count_breakdown=dict(projected_originals=52,fixed_original_whites=6,retained_new_before_U=7,all_facets_after_U=66),
        preflight_parents=parents,failed_preflights=failed,remote_bank=remote,configurations=configs,accepted=accepted,summary=summary,
        scope='Complete declared low2..5/seven-profile necessary bridge domains on the saved parents. Remote60 bank is independent of new walls; all15 native roots joined to all28 frozen T caps. Complete new numerical cones use the frozen full cone and every sign-changing old edge;15 direct cone enumerations corroborate one T configuration. Every combination has complete finite-support PD comparison and whole-face ledgers. Exact old19/12/1 debit unchanged. No high-label, exact-new-normal or finite-volume claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert a.verify or not p.exists();out=json.loads(json.dumps(run()))
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['summary'],flush=True)
