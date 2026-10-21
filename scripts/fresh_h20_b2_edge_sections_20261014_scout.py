#!/usr/bin/env python3
"""Finite necessary B2-face closure census; numerical geometry proposals.

All three spherical extensions are included. This is not a parameter grid.
The output remains numerical until discrete recipes are replayed exactly.
"""
from pathlib import Path
from collections import Counter
import argparse, gzip, hashlib, json, math, os, time
import numpy as np
from fresh_h20_h4_edge_sections_20261013_scout import load, comps

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_b2_edge_sections_20261014'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

def contacts(ga,gb,d):
    ans=[]
    if gb==0:
        if d-ga*ga>0:ans.append(('ORTHOGONAL',float(ga)))
        disc=5*d-4*ga*ga
        if disc>=0:
            for branch in [-1,1]:
                x=(4*ga+branch*math.sqrt(disc))/5
                if x>ga+1e-10 and d-x*x>1e-10:
                    ans.append((f'B3_AT_MODIFIED_END_{branch}',x))
    if gb<0 and d-ga*ga==2*gb*gb:
        ans.append(('B3_AT_UNCHANGED_END',float(ga)))
    return ans

def defect(E,cc,ai,I,dual,outside,G,pp,qp,zp,ep,gamma):
    counts=Counter();neighbors=[]
    for release in [i for i in E if i not in cc]+cc:
        ri=E.index(release);dp=dual[ri]-I[ri,ai]*zp/gamma
        ids=np.array([j for j in outside if dp[j]<-1e-10],int)
        if not len(ids):counts['no_collision']+=1;continue
        j=int(ids[np.argmin(pp[ids]/dp[ids])]);step=float(pp[j]/dp[j]);xp=pp-step*dp
        qd=I[ri,ri]+I[ri,ai]**2/gamma;qx=float(qp+step*step*qd)
        active=np.flatnonzero(abs(xp)<1e-7).tolist()
        if qx>=-1e-8:
            counts['not_finite_first_contact']+=1;continue
        if len(active)!=20 or max(xp)>1e-7:
            counts['nonsimple_or_conditioned_contact']+=1;continue
        H=G[np.ix_(active,active)]-np.outer(ep[active],ep[active])
        diag=np.diag(H)
        if min(diag)<=1e-8 or np.linalg.eigvalsh(H)[0]<=1e-8:
            counts['conditioned_Gram']+=1;continue
        counts['checked_actual_finite_neighbors']+=1
        ii,jj=np.triu_indices(20,1);p=H[ii,jj];sq=p*p/(diag[ii]*diag[jj])
        bad=np.flatnonzero(p>1e-7)
        if len(bad):
            k=int(bad[np.argmax(p[bad])]);kind='OBTUSE';extra={}
        else:
            proposals=[]
            for k in np.flatnonzero((p<-1e-7)&(sq>1e-8)&(sq<1-1e-9)):
                m=math.pi/math.acos(math.sqrt(float(sq[k])));nearest=round(m)
                if abs(m-nearest)<1e-5:continue
                lower=max(2,int(math.floor(m)));upper=lower+1
                c0=math.cos(math.pi/lower)**2;c1=math.cos(math.pi/upper)**2
                margin=min(float(sq[k])-c0,c1-float(sq[k]))
                if margin>1e-6:proposals.append((margin,int(k),lower,upper))
            if not proposals:
                neighbors.append({'release':release,'new_wall':j,'active':active,'norm':qx});continue
            _,k,lower,upper=max(proposals);kind='STRICT_CONSECUTIVE_LABEL_GAP';extra={'adjacent_labels':[lower,upper]}
        return {'status':'NUMERICAL_ACTUAL_RIDGE_DEFECT','witness':{'release':release,'new_wall':j,
          'active':active,'vertex_norm':qx,'pair':[active[int(ii[k])],active[int(jj[k])]],
          'kind':kind,'squared_cosine':float(sq[k]),'projected_pairing':float(p[k]),**extra},'neighbor_counts':dict(counts)}
    return {'status':'UNRESOLVED_AFTER_ALL_INITIAL_EDGES','neighbor_counts':dict(counts),'finite_neighbors':neighbors}

def run(limit):
    Q,W,G,h=load();C=W@Q;corners=json.loads((ROOT/WALK).read_text())['corner_orbits']
    gamma=1/3;seen=set();families=[];counts=Counter();status_counts=Counter();survivors=[];sample=[]
    archive=ROOT/f'results/{TAG}_scout.jsonl.gz';start=time.monotonic();index=0;stop=False
    print('pid',os.getpid(),'limit',limit,flush=True)
    with gzip.open(archive,'wt')as stream:
      for corner in corners:
        if stop:break
        for ei,edge in enumerate(corner['edge_orbits']):
          E=edge['edge_walls']
          if tuple(E)in seen:continue
          seen.add(tuple(E));components=[cc for cc in comps(E,G)if len(cc)==2 and all(G[i,i]==2 for i in cc)and G[cc[0],cc[1]]==-1]
          if not components:continue
          I=np.linalg.inv(G[np.ix_(E,E)].astype(float));dual=I@G[E]
          left=np.array(corner['primitive_ray'],float);right=np.array(edge['endpoint_primitive_ray'],float)
          left/=-(left@Q@h);right/=-(right@Q@h);pv,pb=C@left,C@right
          plane=np.array([[left@Q@left,left@Q@right],[right@Q@left,right@Q@right]])
          assert np.linalg.det(plane)<0 and max(abs(pv[E]))<1e-8 and max(abs(pb[E]))<1e-8
          outside=[l for l in range(210)if l not in E]
          for cc in components:
            if stop:break
            for a in cc:
              if len(families)>=limit:stop=True;break
              b=next(i for i in cc if i!=a);ai=E.index(a);fid=len(families)
              fam={'family':fid,'corner':corner['id'],'edge_index':ei,'edge_walls':E,'component':cc,'anchor':a,
                'endpoint_types':[corner['type'],corners[edge['endpoint_orbit']]['type']]}
              families.append(fam);unique=[];kinds=Counter();ncases=0
              for l in outside:
                ga,gb,d=map(int,[G[l,a],G[l,b],G[l,l]])
                for contact_type,x in contacts(ga,gb,d):
                  kinds[contact_type]+=1
                  p=np.array([pv[l],pb[l]]);assert p@p>1e-12
                  target=x-(2*ga+gb)/3;z0=target*p/(p@p);y=np.array([p[1],-p[0]])
                  aa=y@plane@y;bb=2*z0@plane@y;dd=z0@plane@z0-gamma
                  if abs(aa)<1e-13:
                    if abs(bb)<1e-13:counts['degenerate_quadratic']+=1;continue
                    roots=[(0,-dd/bb)]
                  else:
                    disc=bb*bb-4*aa*dd
                    if disc< -1e-13:counts['negative_discriminant']+=1;continue
                    if abs(disc)<1e-13:counts['near_double_root']+=1
                    roots=[(s,(-bb+s*math.sqrt(max(0,disc)))/(2*aa))for s in[-1,1]]
                  for branch,root in roots:
                    zc=z0+root*y;zv,zb=plane@zc
                    if zv*zb>=-1e-12:counts['not_transverse_interior']+=1;continue
                    if any(np.max(abs(zc-old))<1e-7 for old in unique):counts['duplicate_within_family']+=1;continue
                    unique.append(zc);t=-zv/(zb-zv);pp=(1-t)*pv+t*pb;qp=np.array([1-t,t])@plane@np.array([1-t,t])
                    zp=zc[0]*pv+zc[1]*pb;ep=dual[ai]+zp
                    row={'index':index,'family':fid,'closure_wall':l,'contact_type':contact_type,
                      'normal_contact':float(x),'quadratic_branch':branch,'edge_plane_coefficients':zc.tolist(),'position':float(t)}
                    if qp>=-1e-10 or not 0<t<1 or max(pp[outside])>=-1e-8 or abs(ep[l]-x)>1e-7:
                      row['status']='NUMERICAL_BOUNDARY_OR_CONDITIONING_UNRESOLVED'
                    else:row.update(defect(E,cc,ai,I,dual,outside,G,pp,qp,zp,ep,gamma))
                    status_counts[row['status']]+=1;stream.write(json.dumps(row,separators=(',',':'))+'\n');index+=1;ncases+=1
                    if row['status']!='NUMERICAL_ACTUAL_RIDGE_DEFECT':survivors.append(row)
                    elif len(sample)<100 or ncases==1:sample.append(row)
              fam['contact_recipe_counts']=dict(kinds);fam['distinct_candidate_count']=ncases
              if fid%50==0:
                summary={'status':'RUNNING','pid':os.getpid(),'families_completed':len(families),'candidates':index,
                  'status_counts':dict(status_counts),'survivors':survivors,'elapsed_seconds':time.monotonic()-start}
                (ROOT/f'results/{TAG}_partial.json').write_text(json.dumps(summary,indent=2)+'\n')
                print('families',len(families),'candidates',index,'unresolved',len(survivors),'elapsed',round(time.monotonic()-start,1),flush=True)
    result={'status':'NUMERICAL_B2_FINITE_CLOSURE_CENSUS_COMPLETE_AT_STATED_FAMILY_LIMIT',
      'sources':{p:sha(p)for p in[RAW,WALK,f'scripts/{TAG}_scout.py','scripts/fresh_h20_h4_edge_sections_20261013_scout.py']},
      'family_limit':limit,'families_completed':len(families),'all5048_families_complete':len(families)==5048,
      'families':families,'candidate_count':index,'status_counts':dict(status_counts),'enumeration_counts':dict(counts),
      'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'survivors':survivors,'sample_defects':sample,
      'elapsed_seconds':time.monotonic()-start,
      'scope':'All orthogonal and both B3 spherical face extensions on the stated source families; numerical roots, deduplication and actual-neighbor gates only. No exact global exclusion or positive Coxeter example is claimed.'}
    (ROOT/f'results/{TAG}_scout.json').write_text(json.dumps(result,indent=2)+'\n')
    partial=ROOT/f'results/{TAG}_partial.json'
    if partial.exists():partial.unlink()
    print(result['status'],len(families),index,dict(status_counts),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--families',type=int,default=5048);args=ap.parse_args()
    assert 1<=args.families<=5048;run(args.families)
