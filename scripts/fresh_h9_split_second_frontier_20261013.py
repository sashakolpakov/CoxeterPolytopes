#!/usr/bin/env python3
"""Exact simultaneous no-removal second-frontier caps."""
import json
import time
from itertools import combinations
import mpmath as mp
import numpy as np
from scripts.fresh_h9_split_second_frontier_20261013_base import *
from scripts.fresh_h9_vertex_walk_20261013_caps import solve_line,finite_vertex
from scripts.fresh_h9_split_endpoint_walk_exact_20261013 import k_gram
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K
from scripts.fresh_h9_replaced_core_20261013_cores import affine_census
from scripts.fresh_h9_replaced_core_20261013_geometry import volume_gate
from scripts.fresh_h9_split_edge_frontier_20261013 import frontier

CHOSEN=(33,45,4,136,5)

def make(proposal,index,G,I,h,old,pool):
    E=old['walk2720']['field'];block=proposal['edge_names'];anchor=proposal['anchor_name']
    rows=[]
    for n in block+[anchor]:
        common,er=embed(E,old[n]['field']);assert key(common)==key(E)
        rows.append(tuple(er(x)for x in old[n]['c']))
    values=[Z]*8
    for j,m in proposal['finite_profile']:values[j]=K.neg({3:O,5:P}[m])
    rhs=[lift_K(E,x)for x in values]+[lift_K(E,K.neg({2:Z,3:O,5:P}[proposal['anchor_label']]))]
    free,u,v=solve_line(E,rows,rhs)
    IE=tuple(tuple(lift_K(E,x)for x in row)for row in I)
    A=pairing(E,v,IE,v);B=E.mul(E.rational(2),pairing(E,u,IE,v));C=E.sub(pairing(E,u,IE,u),E.rational(2))
    D=E.sub(E.mul(B,B),E.mul(E.rational(4),E.mul(A,C)));assert E.sign(A)and E.sign(D)>0
    F=Quadratic(E,D);choices=[]
    for branch in(-1,1):
        t=(E.div(E.neg(B),E.mul(E.rational(2),A)),E.div(E.rational(branch),E.mul(E.rational(2),A)))
        p=tuple(F.add(F.lift(x),F.mul(t,F.lift(y)))for x,y in zip(u,v,strict=True))
        error=sum((float(value(F,x))-y)**2 for x,y in zip(p,proposal['basis_covector'],strict=True))
        choices.append((error,branch,t,p))
    error,branch,t,p=min(choices,key=lambda row:row[0]);assert error<1e-9
    c=matvec(F,tuple(tuple(lift_K(F,x)for x in r)for r in I),p)
    name=f'second{index}';row={'name':name,'field':F,'c':c,'p':p,'q':scale(O,2)}
    assert F.sign(F.sub(dot(F,c,p),F.rational(2)))==0
    assert F.sign(dot(F,c,tuple(lift_K(F,x)for x in h)))>0
    for rr,want in zip(rows,rhs,strict=True):assert F.sign(F.sub(dot(F,p,tuple(F.lift(x)for x in rr)),F.lift(want)))==0
    checks={n:relation(row,r,pool)for n,r in old.items()}
    assert all(r['kind']in('finite','dotted')for r in checks.values()),checks
    gram=k_gram([old[n]for n in block]);full=[r+[x]for r,x in zip(gram,values,strict=True)]+[values+[scale(O,2)]]
    assert inertia_K(full)==(9,0,0)
    record={'name':name,'proposal_index':index,'field':field_record(F),'base_field':field_record(E),
            'block_names':block,'block_support':[r[0]for r in proposal['finite_profile']],
            'block_labels':[r[1]for r in proposal['finite_profile']],
            'anchor_name':anchor,'anchor_label':proposal['anchor_label'],
            'free_coordinate':free,'line_constant':list(map(E.encode,u)),'line_direction':list(map(E.encode,v)),
            'norm_polynomial':list(map(E.encode,(A,B,C))),'parameter':F.encode(t),'branch':branch,
            'covector':list(map(F.encode,p)),'coefficients':list(map(F.encode,c)),
            'pair_checks':checks,'spherical_vertex_gram':[[enc(x)for x in r]for r in full],
            'numerical_closed_old_edge_indices':proposal['open_edge_cut_indices']}
    record['finite_vertex']=finite_vertex(row,record,I,h,old|{name:row})
    return row,record

def numerical_frontier(selected,snapshot,I,records):
    N=np.array(snapshot['N']);NB=N[list(BASIS)];names=list(selected)
    W=np.array([[float(value(r['field'],x))for x in r['c']]for r in selected.values()])@NB
    # Discover all vertices in the original six-vertex component from that component alone.
    result=frontier(W,names,snapshot['finite_vertices'])
    return result

def run():
    mp.mp.dps=90;start=time.monotonic();G,I,h,old,_,oldcase,_=current();scout=read(PREFIX+'_scout.json');snapshot=read(PREFIX+'_source.json')
    assert sha(PREFIX+'_source.json')==scout['source_sha256']
    pool=angle_pool(ANGLE_SOURCE);records=[];new=[]
    output={'status':'EXACT_SECOND_FRONTIER_CAPS','source':PREFIX+'_source.json','scout':PREFIX+'_scout.json',
            'source_sha256':{p:sha(p)for p in(PREFIX+'_source.json',PREFIX+'_scout.json',INITIAL+'_family.json',INITIAL+'_family_geometry.json')},
            'old_walls':list(old),'caps':records,'scope':'Exact no-removal second-frontier roots with norm, angles, center and an actual finite vertex. Joint geometry is archived separately.'}
    for index in CHOSEN:
        row,record=make(scout['candidates'][index],index,G,I,h,old,pool);records.append(record);new.append(row)
        (ROOT/(PREFIX+'_caps.json')).write_text(json.dumps(output,indent=2)+'\n')
        print({'cap':index,'actual_finite_vertex':True,'seconds':round(time.monotonic()-start,1)},flush=True)
    pairs=dict(oldcase['pair_checks'])
    for record in records:pairs.update({'|'.join(sorted((record['name'],n))):r for n,r in record['pair_checks'].items()})
    for a,b in combinations(new,2):
        rr=relation(a,b,pool);assert rr['kind']in('finite','dotted'),(a['name'],b['name'],rr)
        pairs['|'.join(sorted((a['name'],b['name'])))]=rr
    selected=old|{r['name']:r for r in new};assert len(pairs)==len(selected)*(len(selected)-1)//2
    census=affine_census(list(selected),pairs)
    vertices=[finite_vertex(r,rec,I,h,selected)for r,rec in zip(new,records,strict=True)]
    geometry=volume_gate(I,h,list(selected.values()))
    front=numerical_frontier(selected,snapshot,I,records)
    result={'source':PREFIX+'_caps.json','source_sha256':sha(PREFIX+'_caps.json'),'walls':list(selected),
            'facet_count':len(selected),'pair_checks':pairs,'affine_census':census,
            'actual_finite_vertices':vertices,'geometry':geometry,'frontier':front,
            'scope':'Exact pair/affine and specified actual-vertex checks, with exact strict escapes. Existing-wall frontier classification is numerical until its endpoint certificates are replayed.'}
    (ROOT/(PREFIX+'_geometry.json')).write_text(json.dumps(result,indent=2)+'\n')
    print({'facets':len(selected),'affines':len(census['affine_circuits']),'geometry':geometry['status'],
           'frontier':front['counts'],'finite_vertices_reached':front['finite_vertices_reached'],'seconds':round(time.monotonic()-start,1)},flush=True)

if __name__=='__main__':run()
