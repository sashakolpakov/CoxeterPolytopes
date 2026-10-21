#!/usr/bin/env python3
"""Fraction-only verification of targeted H9 cuts and local obstructions."""
import argparse
import json
from math import isqrt
from pathlib import Path
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,sub,scale,mul,div,sign,integral,dot,mv,tr,mm,pair,
    matrix,eye,rref,e8_host,
)
from scripts.fresh_multibasepoint_h9_verify import ROOT,read,permutation,valid_coxeter_value

def norm(x):return x[0]*x[0]+x[0]*x[1]-x[1]*x[1]
def nearest(x):return (2*x.numerator+x.denominator)//(2*x.denominator)
def gcd(a,b):
    while b!=Z:
        q=div(a,b);q=(nearest(q[0]),nearest(q[1]));r=sub(a,mul(q,b))
        assert r==Z or abs(norm(r))<abs(norm(b));a,b=b,r
    return a

def escape(witness,walls,gram):
    v=tuple(map(parse,witness['vector']));q=pair(v,gram,v)
    assert q==parse(witness['norm']) and sign(q)>0 and sign(v[0])>0
    normals=tuple(mv(gram,r) for r in walls);values=tuple(dot(v,n) for n in normals)
    assert values==tuple(map(parse,witness['wall_pairings'])) and all(sign(x)<=0 for x in values)
    active=witness['active_facets'];assert len(set(active))==len(active)==9
    assert all(values[i]==Z for i in active) and len(rref([normals[i] for i in active])[1])==9
    return v

def plane(data,walls,gram,witness,inside):
    v=escape(witness,walls,gram);alpha=parse(data['alpha']);z=tuple(map(parse,data['saturated_spatial_direction']))
    beta=parse(data['spatial_norm']);c=parse(data['escape_spatial_multiplier'])
    assert all(map(integral,z)) and sign(c)>0 and tuple(div(x,c) for x in z)==v[1:]
    assert beta==pair((Z,*z),gram,(Z,*z))
    content=Z
    for x in z:content=gcd(content,x)
    assert abs(norm(content))==1
    assert data['positive_t_bound']>0
    obstruction_count=0
    for row in data['norm_searches']:
        q=parse(row['norm']);assert integral(q)
        for o in row['local_obstructions']:
            p=o['prime'];f=o['phi_residue'];residue=o['required_square_residue']
            assert p>2 and all(p%d for d in range(2,isqrt(p)+1))
            assert (f*f-f-1)%p==0 and (int(beta[0])+int(beta[1])*f)%p==0
            al=(int(alpha[0])+int(alpha[1])*f)%p;qq=(int(q[0])+int(q[1])*f)%p
            assert al and residue==(-qq*pow(al,-1,p))%p
            assert residue and pow(residue,(p-1)//2,p)==p-1
            obstruction_count+=1
        if row['status']=='NO_ROOT_IN_THIS_PLANE_BY_LOCAL_OBSTRUCTION':assert row['local_obstructions']
    roots=[];normals=tuple(mv(gram,r) for r in walls)
    for cut in data['cuts']:
        r=tuple(map(parse,cut['root']));x=parse(cut['x']);t=parse(cut['t']);q=parse(cut['norm'])
        assert r==(x,*(mul(t,y) for y in z)) and all(map(integral,r))
        assert sign(t)>0 and sign(sub(scale(O,data['positive_t_bound']),t))>=0
        assert sign(sub(x,mul(c,t)))>=0
        assert sub(mul(beta,mul(t,t)),mul(alpha,mul(x,x)))==q==pair(r,gram,r)
        assert sign(q)>0 and all(integral(div(scale(y,2),q)) for y in mv(gram,r))
        values=tuple(dot(r,n) for n in normals)
        assert values==tuple(map(parse,cut['wall_pairings'])) and all(sign(y)<=0 for y in values)
        assert pair(r,gram,v)==parse(cut['escape_pairing']) and sign(pair(r,gram,v))>0
        assert sign(pair(r,gram,inside))<0
        assert div(mul(alpha,mul(x,x)),q)==parse(cut['old_basepoint_height'])
        for w,n,value in zip(walls,normals,values):
            assert valid_coxeter_value(div(scale(mul(value,value),4),mul(q,dot(w,n))))
        roots.append(r)
    return roots,obstruction_count

def verify(path):
    d=json.loads(Path(path).read_text());base=read(d['source']);gram=matrix(d['gram']);alpha=parse(d['alpha'])
    assert alpha==scale(P,2)
    p=permutation([0,1,8,2,3,4,5,6,7,9]);assert gram==mm(p,mm(e8_host(alpha,True),tr(p)))
    walls=list(matrix(base['all_walls']));inside=tuple(map(parse,d['strict_interior']))
    assert inside==tuple(map(parse,base['strict_interior'])) and sign(pair(inside,gram,inside))<0
    u=matrix(d['standard_block_swap_involution'])
    assert u==matrix(base['standard_block_swap_involution']) and mm(u,u)==eye(10)
    assert mm(tr(u),mm(gram,u))==gram and all(integral(x) for row in u for x in row)
    added=[];orbits=0;obstructions=0
    def append_orbit(r,origin):
        nonlocal orbits
        orbit=[]
        for candidate in (r,mv(u,r)):
            if candidate in walls or candidate in orbit:continue
            assert all(sign(pair(candidate,gram,w))<=0 for w in (*walls,*orbit));orbit.append(candidate)
        assert len(orbit)==2
        for candidate in orbit:
            walls.append(candidate);added.append((candidate,origin))
        orbits+=1
    initial=read('results/fresh_multibasepoint_h9_2phi_escape_plane_20261013.json')
    roots,count=plane(initial,walls,gram,base['hyperideal_witness'],inside);obstructions+=count
    assert len(roots)==1
    append_orbit(roots[0],'initial_escape_plane')
    for index,record in enumerate(d['rounds']):
        assert record['round']==index
        cut_found=False
        for trial in record['trials']:
            if trial.get('status')=='SKIPPED_LARGE_SLICE_DISCRIMINANT':continue
            assert not cut_found
            roots,count=plane(trial,walls,gram,trial['targeted_escape'],inside);obstructions+=count
            if roots:
                assert len(roots)==1;append_orbit(roots[0],f'round_{index}');cut_found=True
        assert cut_found and record['added_walls']==2
    expected=matrix(d['all_walls']);assert tuple(walls)==expected and len(walls)==d['wall_count']==372
    assert [(tuple(map(parse,z['root'])),z['origin']) for z in d['added_walls']]==added
    assert set(mv(u,r) for r in walls)==set(walls)
    normals=tuple(mv(gram,r) for r in walls);norms=tuple(dot(r,n) for r,n in zip(walls,normals))
    new_finite_angles={}
    for i,(r,n,q) in enumerate(zip(walls,normals,norms)):
        assert q in (O,scale(O,2)) and all(map(integral,r))
        assert all(integral(div(scale(y,2),q)) for y in n) and sign(dot(inside,n))<0
        for j in range(i):
            value=dot(r,normals[j]);assert sign(value)<=0
            ratio=div(scale(mul(value,value),4),mul(q,norms[j]));assert valid_coxeter_value(ratio)
            if i>=len(base['all_walls']) and sign(sub(ratio,scale(O,4)))<0:
                label={Z:2,O:3,scale(O,2):4,scale(O,3):6,add(O,P):5,add(scale(O,2),P):10}[ratio]
                new_finite_angles[label]=new_finite_angles.get(label,0)+1
    v=escape(d['hyperideal_witness'],walls,gram)
    assert d['status']=='EXACT_ESCAPE_REMAINS_AFTER_BOUNDED_ADAPTIVE_CUTS'
    # Also replay every local obstruction from the separate sqrt5 plane.
    other=read('results/fresh_multibasepoint_h9_sqrt5_escape_plane_20261013.json')
    other_base=read(other['source'])
    _,other_obstructions=plane(other,matrix(other_base['all_walls']),matrix(other_base['gram']),
        other_base['hyperideal_witness'],tuple(map(parse,other_base['strict_interior'])))
    result={'artifact':str(path),'walls':len(walls),'cut_orbits':orbits,'new_walls':len(added),
        'new_finite_angle_counts':new_finite_angles,'verified_local_obstructions':obstructions+other_obstructions,
        'escape_norm':enc(pair(v,gram,v)),'status':'EXACT_COXETER_CUTS_AND_FINAL_ESCAPE_PASS'}
    print(json.dumps(result,sort_keys=True));return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('artifact',type=Path);p.add_argument('--output',type=Path)
    a=p.parse_args();result=verify(a.artifact)
    if a.output:a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
