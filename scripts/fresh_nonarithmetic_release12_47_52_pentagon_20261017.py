#!/usr/bin/env python3
"""Changed five-side cavity after joint release12/47/52.

Complete finite endpoint profile shapes; bounded numerical realization test.
This is not a new exact chamber or an all-parameter geometric exclusion.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_alternative66_endcaps_20261017 as A

H,B,G,M=A.H,A.B,A.G,A.M
ROOT=H.ROOT
TAG='fresh_nonarithmetic_release12_47_52_pentagon_20261017'
DELETE=[12,47,52]
CYCLE=[7,8,65,64,30]
OLD_LABELS={frozenset((7,8)):5,frozenset((8,65)):2,
            frozenset((65,64)):2,frozenset((64,30)):2,frozenset((30,7)):5}

def finite_corner(a,b,pi,ma,mb):
    ids=[5,6,a,b,67];g=[[2]*5 for _ in ids]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    for j,k in enumerate((a,b)):
        for i,m in enumerate(G.profile_labels(4 if k==7 else 0)):
            g[i][j+2]=g[j+2][i]=m
    for i,m in enumerate(G.profile_labels(pi)):g[i][4]=g[4][i]=m
    g[2][3]=g[3][2]=OLD_LABELS[frozenset((a,b))]
    g[2][4]=g[4][2]=ma;g[3][4]=g[4][3]=mb
    return A.components(tuple(x for row in g for x in row))

def domains(frame):
    low=[];tails=[]
    for pi in range(7):
        for ms in product(range(2,7),repeat=len(frame)):
            types=[finite_corner(a,b,pi,ms[j],ms[j+1])for j,(a,b)in enumerate(zip(frame,frame[1:]))]
            if any(t is None for t in types):continue
            high=[j for j,m in enumerate(ms)if m==6]
            row=dict(profile=pi,labels=list(ms),corner_types=types,high_positions=high)
            (tails if high else low).append(row)
    return dict(low=low,tails=tails)

def full_product(y,pi,z,pj,profiles):
    return H.A.native(y,pi,z,pj,profiles)

def score(y,pi,ns,labels,h,profiles,actual,squares,assigned):
    p=profiles[pi];off=p['offsets'];native={i:float(y@M@ns[i]+off[labels.index(i)])for i in labels}
    rel={i:B.relation(x)for i,x in native.items()}
    for i,m in assigned.items():rel[i]=m
    bad=[dict(wall=i,product=native[i])for i in labels if i not in DELETE and (rel[i]is None or rel[i]=='parallel_or_unresolved')]
    center=float(y@M@h+p['center_offset'])
    out=dict(projection=y.tolist(),profile=pi,center_product=center,blockers=bad,
             finite_contacts=[[i,r]for i,r in rel.items()if isinstance(r,int)and i not in DELETE],
             norm_error=float(y@M@y-p['square']))
    if not bad and center<-1e-8:
        rr=[ns[i]for i in actual]+[y];ll=actual+[67];ss=[squares[i]for i in actual]+[p['square']]
        out['incidence']=H.incidence_score(rr,ll,ss)
        out['boundary']=B.ledger(rr,ll)
    return out

def run():
    z,labels,ns,h,profiles,primary=H.setup();actual=[i for i in z['actual_sides']if i not in DELETE]
    squares={r['wall']:B.mid(r['square'])for r in z['projections']}
    cost=H.deletion(primary,z,DELETE);removed=B.ledger([ns[i]for i in actual],actual)
    assert cost['lost_finite_count']==7 and cost['compact_faces_by_dimension']=={2:4}
    assert (removed['finite'],removed['exterior'],removed['open_ends'])==(19,7,11)
    # Every retained cycle edge is geodesic and its pair is exactly native 2 or5.
    assert all(B.relation(float(ns[a]@M@ns[b]))==OLD_LABELS[frozenset((a,b))]
               for a,b in zip(CYCLE,CYCLE[1:]+CYCLE[:1]))
    stages=[];total=Counter();improvements=[]
    for at in range(5):
        left=[CYCLE[(at+j)%5]for j in(-1,0,1)]
        right=[CYCLE[(at+j)%5]for j in(1,2,3,4)]
        dl,dr=domains(left),domains(right)
        frame=np.array([ns[i]for i in right]);gram=frame@M@frame.T;inv=np.linalg.inv(gram)
        right_rows=[];right_pool=[];counts=Counter();unresolved=[]
        for recipe in dr['low']+dr['tails']:
            pi=recipe['profile'];p=profiles[pi];off=np.array(p['offsets'])
            contact=-np.array([B.theta(m)if m!=6 else 0. for m in recipe['labels']])-off[[labels.index(i)for i in right]]
            high=recipe['high_positions'];trials=[]
            if not high:
                y=contact@inv@frame;error=float(y@M@y-p['square'])
                counts['low_norm_equations']+=1
                if abs(error)<1e-7:trials.append(([],y))
                else:counts['low_wrong_norm']+=1
                row=dict(recipe=recipe,norm_residual=error)
            elif len(high)==1:
                j=high[0];v=np.zeros(4);v[j]=-1
                poly=[float(v@inv@v),float(2*v@inv@contact),float(contact@inv@contact-p['square'])]
                roots=np.roots(poly if abs(poly[0])>1e-12 else poly[1:])
                kept=[]
                for root in roots:
                    if abs(root.imag)>1e-8 or not np.sqrt(3)-1e-8<=root.real<=2+1e-8:continue
                    c=float(root.real);cc=contact+c*v;y=cc@inv@frame;trials.append(([c],y));kept.append(c)
                counts['one_tail_norm_equations']+=1
                row=dict(recipe=recipe,cosine_norm_polynomial=poly,real_domain_roots=kept)
            else:
                # Two or more independent continuous labels are not silently cut off.
                row=dict(recipe=recipe,status='unresolved_multivariate_norm_locus',constant=contact.tolist(),inverse=inv.tolist())
                unresolved.append(row);counts['multivariate_tail_loci']+=1
            row['root_trials']=[]
            for cs,y in trials:
                assigned=dict(zip(right,recipe['labels']));tail_rel=[]
                for j,c in zip(high,cs):
                    r=B.relation(-c);tail_rel.append(r);assigned[right[j]]=r if isinstance(r,int)else 'unresolved_tail_angle'
                r=score(y,pi,ns,labels,h,profiles,actual,squares,assigned)
                r['tail_cosines']=cs;r['tail_relations']=tail_rel;row['root_trials'].append(r)
                counts['real_right_norm_roots']+=1
                if not r['blockers'] and r['center_product']<-1e-8:
                    r['all_assigned_angles_coxeter']=all(isinstance(v,int)for v in assigned.values())
                    r['recipe']=recipe;right_pool.append(r);counts['acute_right_candidates']+=1
            right_rows.append(row)
        # Low left recipes plus every symbolic high-tail family remain recorded.
        left_pool=[];left_rows=[]
        for recipe in dl['low']:
            for r in A.frame_solver(left,labels,ns,h,profiles,recipe):
                y=np.array(r['projection']);s=score(y,recipe['profile'],ns,labels,h,profiles,actual,squares,dict(zip(left,recipe['labels'])))
                s['recipe']=recipe;s['branch']=r['branch'];left_rows.append(s)
                if not s['blockers'] and s['center_product']<-1e-8:left_pool.append(s)
        pairs=[]
        for l,r in product(left_pool,right_pool):
            yl,yr=np.array(l['projection']),np.array(r['projection'])
            value=full_product(yl,l['profile'],yr,r['profile'],profiles);rel=B.relation(value)
            if not isinstance(rel,int)or not r['all_assigned_angles_coxeter']:continue
            rr=[ns[i]for i in actual]+[yl,yr];ll=actual+[67,68];ss=[squares[i]for i in actual]+[profiles[l['profile']]['square'],profiles[r['profile']]['square']]
            inc=H.incidence_score(rr,ll,ss);bound=B.ledger(rr,ll)
            rec=dict(left=l,right=r,mutual_product=value,mutual_label=rel,incidence=inc,boundary=bound)
            pairs.append(rec)
            if inc.get('open_ends',99)<10:improvements.append(dict(at=at,**rec))
        stage=dict(at=at,left_frame=left,right_frame=right,left_domain=dl,right_domain=dr,
            right_norm_trials=right_rows,left_low_trials=left_rows,accepted_left_low=left_pool,
            accepted_right=right_pool,finite_join_pairs=pairs,counts=dict(counts),unresolved_tail_loci=unresolved)
        stages.append(stage);total.update(counts)
        print('PENTAGON_SPLIT',at,left,right,'domains',len(dl['low']),len(dl['tails']),len(dr['low']),len(dr['tails']),
              'accepted',len(left_pool),len(right_pool),'pairs',len(pairs),'counts',dict(counts),flush=True)
    paths=[Path(__file__),Path(A.__file__),Path(H.__file__),Path(G.__file__),H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_NUMERICAL_TWO_CENTER_PENTAGON_REPLACEMENT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        old_control=dict(walls=63,ambient_finite=199,S6_open_ends=10),deleted=DELETE,retained_cycle=CYCLE,
        deletion_cost=cost,deletion_boundary=removed,white_profiles=profiles,stages=stages,counts=dict(total),improvements=improvements,
        topology='Five retained boundary edges; two interior vertices L,R. L covers two consecutive old edges; R covers the other three; L,R share a new edge. The two additional triangles use the two shared old endpoints. Seven new finite triangles would close this cavity if all literal geometric gates held.',
        finite_domain='All seven A2 profiles and labels2..5 plus symbolic6 for every isolated I2(m),m>=6. Complete required finite-corner shapes, but numerical norm/gate tests only. One-tail norm roots are computed without a label cutoff; multivariate tails and left-tail joins remain explicitly unresolved.',
        scope='No exact chamber promotion and no all-parameter exclusion. Deletion costs are inherited exact support incidences; all new geometry checks are numerical. Full native angles/common center and complete boundary are required before an improvement is nominated.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/'results'/(TAG+'.json')
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'improvements',len(out['improvements']),flush=True)
