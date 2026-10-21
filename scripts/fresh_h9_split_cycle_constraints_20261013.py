#!/usr/bin/env python3
"""Exact simultaneous closure of two edges of the split-seed H9 chamber.

All arithmetic and branch selection are exact. Candidate2009 suggested the
small pairing system, but its decimal coordinates are not used to solve it.
The script does not enumerate every possible cap or assert finite volume.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

from scripts.fresh_h9_split_endpoint_walk_exact_20261013_base import (
    ROOT,K,Quadratic,Z,O,P,source,read,enc,parse,scale,add,sub,mul,div,
    dot,matvec,pairing,lift_K,embed,lift_existing,key,relation,angle_pool,
    ANGLE_SOURCE,strict_escape)
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import inertia_K
from scripts.fresh_h4_blocks_independent_verify import inverse

PREFIX='data/enumerations/fresh_h9_split_cycle_constraints_20261013'
OUT=ROOT/(PREFIX+'.json')
SCOUT='data/enumerations/fresh_h9_split_endpoint_walk_optional_anchor_20261013_scout.json'
BASE_SCOUT='data/enumerations/fresh_h9_split_endpoint_walk_20261013_scout.json'
OLD_CAPS='data/enumerations/fresh_h9_split_endpoint_walk_exact_20261013.json'
COMMON=[f'seed{i}'for i in(0,1,2,4,5,6,7)]


def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def affine_line(F,rows,rhs):
    assert len(rows)==len(rhs)==9 and all(len(row)==10 for row in rows)
    a=[list(row)+[x]for row,x in zip(rows,rhs,strict=True)];pivots=[]
    for j in range(10):
        i=next((i for i in range(len(pivots),9)if F.sign(a[i][j])),None)
        if i is None:continue
        k=len(pivots);a[k],a[i]=a[i],a[k]
        pivot=a[k][j];a[k]=[F.div(x,pivot)for x in a[k]]
        for i in range(9):
            if i==k:continue
            c=a[i][j];a[i]=[F.sub(x,F.mul(c,y))for x,y in zip(a[i],a[k],strict=True)]
        pivots.append(j)
    assert len(pivots)==9
    free=next(i for i in range(10)if i not in pivots)
    u=[F.zero]*10;v=[F.zero]*10;v[free]=F.one
    for i,j in enumerate(pivots):u[j]=a[i][-1];v[j]=F.neg(a[i][free])
    for row,want in zip(rows,rhs,strict=True):
        assert F.sign(F.sub(dot(F,row,u),want))==0 and F.sign(dot(F,row,v))==0
    return free,tuple(u),tuple(v)


def pair_with_old(F,c,old):
    E,er=embed(F,old['field'])
    return E,dot(E,tuple(lift_existing(E,F,x)for x in c),tuple(er(x)for x in old['p']))


def exact_vertex(F,incident,old,new,I,h,embed_root):
    profiles=old|{new['name']:new};roots=[]
    for name in incident:
        row=profiles[name]
        if row is new:c,p=row['c'],row['p']
        else:c,p=embed_root(row['c'],row['field']),embed_root(row['p'],row['field'])
        roots.append((c,p))
    gram=[]
    for c,p in roots:
        row=[]
        for cc,pp in roots:
            val=dot(F,c,pp)
            # All specified vertex Gram entries are 0,2,−1,−phi.
            match=[z for z in(Z,scale(O,2),K.neg(O),K.neg(P))if F.sign(F.sub(val,lift_K(F,z)))==0]
            assert len(match)==1;row.append(match[0])
        gram.append(row)
    assert inertia_K(gram)==(9,0,0)
    gh=tuple(dot(F,c,tuple(lift_K(F,x)for x in h))for c,p in roots)
    weights=matvec(F,tuple(tuple(lift_K(F,x)for x in row)for row in inverse(gram)),gh)
    y=tuple(lift_K(F,x)for x in h)
    for(c,p),a in zip(roots,weights,strict=True):y=tuple(F.sub(x,F.mul(a,z))for x,z in zip(y,p,strict=True))
    If=tuple(tuple(lift_K(F,x)for x in row)for row in I)
    square=pairing(F,y,If,y);center=pairing(F,y,If,tuple(lift_K(F,x)for x in h))
    assert F.sign(square)<0 and F.sign(center)<0
    signs={}
    for name,row in profiles.items():
        if row is new:sf=F;val=dot(F,row['c'],y)
        else:
            sf,er=embed(F,row['field'])
            val=dot(sf,tuple(er(x)for x in row['c']),tuple(lift_existing(sf,F,x)for x in y))
        signs[name]=sf.sign(val)
    assert all(s==(0 if n in incident else 1)for n,s in signs.items())
    return {'incident':incident,'y':list(map(F.encode,y)),'gram':[[enc(x)for x in row]for row in gram],
            'square':F.encode(square),'center_pairing':F.encode(center),'all_wall_signs':signs}


def build():
    G,I,h,old,scout,banks=source();pool=angle_pool(ANGLE_SOURCE)
    assert inertia_K(G)==(9,1,0)
    assert all(matvec(K,I,row)==tuple(O if i==j else Z for j in range(10))for i,row in enumerate(G))
    c17,c22=old['wall17'],old['wall22']
    E17=c17['field'];E22=c22['field']
    assert E17.base is E22.base is K
    assert E17.d==add(scale(O,11),scale(P,Q(-23,16)))
    assert E22.d==add(scale(O,-4),scale(P,Q(169,16)))
    E,ee22=embed(E17,E22)
    def embed_E(values,field):
        if field is K:return tuple(lift_K(E,x)for x in values)
        if key(field)==key(E17):return tuple(E.lift(x)for x in values)
        if key(field)==key(E22):return tuple(ee22(x)for x in values)
        raise AssertionError('Unexpected local endpoint field')
    equations=COMMON+['wall17','wall22']
    rows=[embed_E(old[n]['c'],old[n]['field'])for n in equations]
    rhs=[E.zero]*8+[E.rational(-1)]
    free,u,v=affine_line(E,rows,rhs)
    IE=tuple(tuple(lift_K(E,x)for x in row)for row in I)
    A=pairing(E,v,IE,v);B=E.mul(E.rational(2),pairing(E,u,IE,v))
    C=E.sub(pairing(E,u,IE,u),E.rational(2))
    D=E.sub(E.mul(B,B),E.mul(E.rational(4),E.mul(A,C)))
    assert E.sign(A)and E.sign(D)>0
    F=Quadratic(E,D)
    IF=tuple(tuple(lift_K(F,x)for x in row)for row in I)
    roots=[];rejected=[]
    for branch in(-1,1):
        t=(E.div(E.neg(B),E.mul(E.rational(2),A)),E.div(E.rational(branch),E.mul(E.rational(2),A)))
        p=tuple(F.add(F.lift(x),F.mul(t,F.lift(z)))for x,z in zip(u,v,strict=True));c=matvec(F,IF,p)
        assert F.sign(F.sub(dot(F,c,p),F.rational(2)))==0
        cp=dot(F,c,tuple(lift_K(F,x)for x in h))
        row={'name':'cycle2009','field':F,'c':c,'p':p,'q':scale(O,2)}
        seeds={n:relation(row,old[n],pool)for n in old if n.startswith('seed')}
        if F.sign(cp)<=0 or any(x['kind']not in('finite','dotted')for x in seeds.values()):
            rejected.append({'branch':branch,'center_sign':F.sign(cp),'seed_relations':seeds});continue
        checks={n:relation(row,r,pool)for n,r in old.items()}
        assert all(x['kind']in('finite','dotted')for x in checks.values())
        roots.append((row,branch,t,cp,checks))
    assert len(roots)==1
    row,branch,t,cp,checks=roots[0]
    for name,want in zip(equations,rhs,strict=True):
        coeff=tuple(F.lift(x)for x in embed_E(old[name]['c'],old[name]['field']))
        assert F.sign(F.sub(dot(F,coeff,row['p']),F.lift(want)))==0
    def embed_F(vals,field):return tuple(F.lift(x)for x in embed_E(vals,field))
    orders=[COMMON+['seed3','wall22'],COMMON+['seed3','wall17'],
            COMMON+['wall17','cycle2009'],COMMON+['wall22','cycle2009']]
    vertices=[]
    for i,order in enumerate(orders):
        vertices.append(exact_vertex(F,order,old,row,I,h,embed_F))
        print('Exact cycle vertex',i,'with nine independent incident walls',flush=True)
    # Four consecutive pairs share exactly eight facets. Their segment is
    # future timelike, and all remaining inequalities are strict by convexity.
    edges=[]
    for i,j in((0,1),(1,2),(2,3),(3,0)):
        common=sorted(set(orders[i])&set(orders[j]));assert len(common)==8
        yi=tuple(map(F.parse,vertices[i]['y']));yj=tuple(map(F.parse,vertices[j]['y']))
        mid=tuple(F.add(a,b)for a,b in zip(yi,yj,strict=True))
        assert F.sign(pairing(F,mid,IF,mid))<0
        # Distinct endpoints cannot be proportional: their exact zero sets differ.
        assert orders[i]!=orders[j]
        edges.append({'vertices':[i,j],'common_walls':common,'midpoint_square':F.encode(pairing(F,mid,IF,mid))})
    # Verify the numerical proposal's identity by rational intervals, rather
    # than using it to choose a branch or to decide any geometric sign.
    optional_scout=read(SCOUT)
    assert optional_scout['source_snapshot']==scout['source_snapshot']
    proposed=optional_scout['candidates'][2009]
    for i,number in enumerate(proposed['seed_covector']):
        ff,value=pair_with_old(F,row['c'],old[f'seed{i}'])
        target=Q(str(number));tol=Q(1,10**7)
        assert ff.sign(ff.sub(value,ff.rational(target-tol)))>0
        assert ff.sign(ff.sub(value,ff.rational(target+tol)))<0
    previous=read(OLD_CAPS)['target'];y=tuple(map(parse,previous['y']))
    escape=strict_escape(y,I,h,list(old.values())+[row])
    local={n:checks[n]for n in COMMON+['seed3','wall17','wall22']}
    return {'status':'EXACT_TWO_EDGE_CLOSURE_AND_FINITE_FOUR_CYCLE',
        'source_sha256':{p:sha(p)for p in sorted(set(banks)|{SCOUT,BASE_SCOUT,OLD_CAPS,ANGLE_SOURCE})},
        'source_old_walls':list(old),'independent_seed_basis':[0,1,2,3,4,5,6,7,8,12],
        'seed_gram':[[enc(x)for x in r]for r in G],'seed_gram_inverse':[[enc(x)for x in r]for r in I],
        'center_covector':list(map(enc,h)),
        'local_equation_walls':equations,'local_rhs':list(map(E.encode,rhs)),
        'base_radicands':[enc(E17.d),enc(E22.d)],'norm_discriminant':E.encode(D),
        'field_degree_upper_bound':16,'line_free_coordinate':free,
        'line_constant':list(map(E.encode,u)),'line_direction':list(map(E.encode,v)),
        'norm_polynomial':list(map(E.encode,(A,B,C))),'branch':branch,'parameter':F.encode(t),
        'normal_coefficients':list(map(F.encode,row['c'])),'normal_covector':list(map(F.encode,row['p'])),
        'norm':F.encode(dot(F,row['c'],row['p'])),'strict_center_pairing':F.encode(cp),
        'pair_relations':checks,'local_relations':local,'rejected_branch':rejected,
        'finite_vertices':vertices,'four_cycle_edges':edges,
        'old_open_edge_indices_repaired':[3,11],
        'old_escape_survives':escape,'full_wall_count':73,
        'scope':'One exact cap preserving all72 old walls and closing two old outgoing edges simultaneously. Four actual finite vertices and all four connecting edges form a cycle in the common seven-wall H2 face. The resulting73-wall chamber still has the saved strict spacelike escape. No complete cap-family enumeration, new finite-volume polytope, or unrestricted seed obstruction.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();r=build()
    if a.verify:assert r==read(PREFIX+'.json')
    else:OUT.write_text(json.dumps(r,indent=2)+'\n')
    print({k:r[k]for k in('status','field_degree_upper_bound','branch','full_wall_count')},flush=True)
