#!/usr/bin/env python3
"""Exactify selected no-removal edge caps from the frozen split-seed scout."""
import json
import time
from itertools import combinations
import mpmath as mp
import numpy as np
from scripts.fresh_h9_split_endpoint_walk_exact_20261013_base import *
from scripts.fresh_h9_vertex_walk_20261013_caps import solve_line, finite_vertex
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K
from scripts.fresh_h9_replaced_core_20261013_cores import affine_census
from scripts.fresh_h9_replaced_core_20261013_geometry import volume_gate

def k_gram(profiles):
    result=[]
    for a in profiles:
        row=[]
        for b in profiles:
            F,p=pair_profiles(a,b)
            z=p
            while F is not K:
                assert F.base.sign(z[1])==0
                z=z[0];F=F.base
            row.append(z)
        result.append(row)
    return result

def target_proposal(scout,I,h,old):
    snap=scout['source_snapshot'];N=np.array(snap['N']);x=np.r_[1.,snap['leading_escape']['point']]
    metric=np.r_[-1.,np.ones(9)]
    y=-N[list(BASIS)]@(metric*x)
    hn=np.array([float(value(K,z))for z in h]);hn/=hn.sum()
    for epsilon in (1e-5,1e-7,1e-3):
        for den in (10**12,10**15):
            yy=tuple(scale(O,Q(round(float(z)*den),den))for z in y+epsilon*hn)
            try:return strict_escape(yy,I,h,list(old.values()))
            except AssertionError:pass
    raise AssertionError('Could not reconstruct old leading escape')

def make(proposal,index,I,h,old,target,pool):
    names=[wall_name(j)for j in proposal['edge']];anchor_name=wall_name(proposal['anchor'])
    anchor=old[anchor_name];E=anchor['field'];assert E is not K and E.base is K
    rows=[]
    for name in names+[anchor_name]:
        common,ee=embed(E,old[name]['field']);assert key(common)==key(E)
        rows.append(tuple(ee(x)for x in old[name]['c']))
    values=[Z]*8
    for j,m in proposal['finite_profile']:values[j]=K.neg({3:O,5:P}[m])
    rhs=[E.lift(x)for x in values]+[E.lift(K.neg({2:Z,3:O,5:P}[proposal['anchor_label']]))]
    free,u,v=solve_line(E,rows,rhs)
    IE=tuple(tuple(E.lift(x)for x in row)for row in I)
    A=pairing(E,v,IE,v);B=E.mul(E.rational(2),pairing(E,u,IE,v))
    C=E.sub(pairing(E,u,IE,u),E.rational(2))
    D=E.sub(E.mul(B,B),E.mul(E.rational(4),E.mul(A,C)))
    assert E.sign(A)and E.sign(D)>0
    F=Quadratic(E,D);choices=[]
    for branch in (-1,1):
        t=(E.div(E.neg(B),E.mul(E.rational(2),A)),E.div(E.rational(branch),E.mul(E.rational(2),A)))
        p=tuple(F.add(F.lift(x),F.mul(t,F.lift(y)))for x,y in zip(u,v,strict=True))
        error=sum((float(value(F,z))-proposal['seed_covector'][i])**2 for i,z in zip(BASIS,p,strict=True))
        choices.append((error,branch,t,p))
    error,branch,t,p=min(choices,key=lambda row:row[0]);assert error<1e-9
    c=matvec(F,tuple(tuple(lift_K(F,x)for x in row)for row in I),p)
    name=f'walk{index}';row={'name':name,'field':F,'c':c,'p':p,'q':scale(O,2)}
    assert F.sign(F.sub(dot(F,c,p),F.rational(2)))==0
    assert F.sign(dot(F,c,tuple(lift_K(F,x)for x in h)))>0
    cut=dot(F,c,tuple(lift_K(F,x)for x in target));assert F.sign(cut)<0
    checks={n:relation(row,r,pool)for n,r in old.items()}
    dropped=[n for n,r in checks.items()if r['kind']not in ('finite','dotted')]
    assert not dropped
    gram=k_gram([old[n]for n in names]);full=[r+[x]for r,x in zip(gram,values,strict=True)]+[values+[scale(O,2)]]
    assert inertia_K(gram)==(8,0,0)and inertia_K(full)==(9,0,0)
    rec={'proposal_index':index,'base_radicand':enc(E.d),'discriminant':E.encode(D),
         'block_names':names,'block_support':[r[0]for r in proposal['finite_profile']],
         'block_labels':[r[1]for r in proposal['finite_profile']],
         'anchor_name':anchor_name,'anchor_label':proposal['anchor_label'],
         'free_coordinate':free,'line_constant':list(map(E.encode,u)),'line_direction':list(map(E.encode,v)),
         'norm_polynomial':list(map(E.encode,(A,B,C))),'parameter':F.encode(t),'branch':branch,
         'seed_pairings':list(map(F.encode,p)),'seed_coefficients':list(map(F.encode,c)),
         'target_pairing':F.encode(cut),'dropped':dropped,'pair_checks':checks,
         'spherical_vertex_gram':[[enc(x)for x in r]for r in full]}
    rec['finite_vertex']=finite_vertex(row,rec,I,h,old|{name:row})
    return row,rec

def run():
    mp.mp.dps=90;start=time.monotonic()
    G,I,h,old,scout,banks=source();optional=read(OPTIONAL_SCOUT)
    assert optional['source_snapshot']==scout['source_snapshot']
    assert inertia_K(G)==(9,1,0)and K.sign(pairing(K,h,I,h))<0
    proof=target_proposal(scout,I,h,old);target=tuple(map(parse,proof['y']))
    pool=angle_pool(ANGLE_SOURCE);caps=[];new=[]
    result={'status':'EXACT_NO_REMOVAL_SPLIT_EDGE_CAPS','source':SCOUT,'proposal_source':OPTIONAL_SCOUT,
            'source_sha256':{p:sha(p)for p in [SCOUT,OPTIONAL_SCOUT,*banks]},
            'independent_seed_basis':list(BASIS),'gram':[[enc(x)for x in r]for r in G],
            'inverse':[[enc(x)for x in r]for r in I],'strict_center':list(map(enc,h)),
            'old_walls':list(old),'target':proof,'caps':caps,
            'scope':'Exact accepted no-removal cap equations, all old pair relations, an actual finite vertex and target cut. Geometry cases separately decide the tested chambers; no exhaustive family claim.'}
    for index in (2747,2735):
        row,rec=make(optional['candidates'][index],index,I,h,old,target,pool)
        caps.append(rec);new.append(row)
        (ROOT/(PREFIX+'.json')).write_text(json.dumps(result,indent=2)+'\n')
        print({'proposal':index,'retained_old':72,'actual_finite_vertex':True,'seconds':round(time.monotonic()-start,1)},flush=True)
    oldpairs={}
    for a,b in combinations(old.values(),2):
        rr=relation(a,b,pool);assert rr['kind']in('finite','dotted')
        oldpairs['|'.join(sorted((a['name'],b['name'])))]=rr
    newpair=relation(*new,pool)
    selections=[(0,),(1,)]
    if newpair['kind']in('finite','dotted'):selections.append((0,1))
    cases=[]
    for ids in selections:
        selected=old|{new[i]['name']:new[i]for i in ids};checks=dict(oldpairs)
        for i in ids:
            checks.update({'|'.join(sorted((new[i]['name'],n))):r for n,r in caps[i]['pair_checks'].items()})
        if len(ids)==2:checks['|'.join(sorted(r['name']for r in new))]=newpair
        census=affine_census(list(selected),checks)
        vertices=[finite_vertex(new[i],caps[i],I,h,selected)for i in ids]
        geometry=volume_gate(I,h,list(selected.values()))
        cases.append({'new_caps':list(ids),'walls':list(selected),'facet_count':len(selected),
                      'pair_checks':checks,'affine_census':census,'actual_finite_vertices':vertices,
                      'geometry':geometry,'status':geometry['status']})
        (ROOT/(PREFIX+'_geometry.json')).write_text(json.dumps({'source':PREFIX+'.json','new_pair_relation':newpair,'cases':cases},indent=2)+'\n')
        print({'caps':ids,'facets':len(selected),'affine':len(census['affine_circuits']),'status':geometry['status'],
               'maximum_norm':geometry['maximum_numeric_vertex_norm'],'seconds':round(time.monotonic()-start,1)},flush=True)

if __name__=='__main__':run()
