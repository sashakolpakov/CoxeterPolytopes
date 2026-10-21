#!/usr/bin/env python3
"""Independent K-arithmetic replay of the D9 infinite chamber transport.

The transport producer uses a separate Fraction implementation. This checker
does not import that producer or its matrix construction.
"""
import argparse,hashlib,json
from pathlib import Path
from fresh_h9_glued_host_20261013_arithmetic import *

ROOT=Path(__file__).resolve().parents[1]
PREFIX='fresh_h9_changed_spatial_next_20261013_transport_independent'
SOURCE=ROOT/'results/fresh_h9_d9_symmetry_transfer_20261013.json'
OUT=ROOT/'results'/f'{PREFIX}.json'


def check():
    d=json.loads(SOURCE.read_text())
    for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    A=tuple(tuple(map(K.from_json,row))for row in d['transport_matrix'])
    AI=tuple(tuple(map(K.from_json,row))for row in d['inverse_matrix'])
    assert inverse(A)==AI
    G=tuple(tuple(-ALPHA if i==j==0 else O if i==j else Z for j in range(10))for i in range(10))
    basis=[(O,*(Z for _ in range(9)))]
    for i in range(1,9):
        v=[Z]*10;v[i]=O;v[i+1]=-O;basis.append(tuple(v))
    v=[Z]*10;v[8]=v[9]=O;basis.append(tuple(v))
    def lattice(v):return all(x.is_integral for x in v)and(sum(v[1:],Z)/2).is_integral
    assert len(basis)==10
    modules=[]
    B=tuple(zip(*basis));BI=inverse(B)
    for M in(A,AI):
        images=[mv(M,v)for v in basis]
        assert all(lattice(v)for v in images)
        assert all(pair(images[i],G,images[j])==pair(basis[i],G,basis[j])for i in range(10)for j in range(10))
        module=tuple(zip(*(mv(BI,v)for v in images)))
        assert all(x.is_integral for row in module for x in row)
        modules.append(module)
    assert modules[0]==tuple(tuple(map(K.from_json,row))for row in d['module_matrix'])
    assert modules[1]==tuple(tuple(map(K.from_json,row))for row in d['module_inverse'])
    sp=ROOT/'results/fresh_h9_changed_spatial_20261013_step3.json'
    tp=ROOT/'results/fresh_h9_changed_spatial_next_20261013.json'
    s=json.loads(sp.read_text());t=json.loads(tp.read_text())
    old=tuple(tuple(map(K.from_json,r))for r in s['endpoint_incident'])
    new=tuple(tuple(map(K.from_json,r))for r in t['endpoint_incident'])
    perm=d['permutation'];assert perm==[8,1,2,3,4,5,6,7,0]
    assert all(mv(A,r)==new[perm[i]]for i,r in enumerate(old))
    scale=K.from_json(d['positive_point_scale']);assert scale>0
    p=tuple(map(K.from_json,s['endpoint']));q=tuple(map(K.from_json,t['endpoint']))
    assert pair(p,G,p)<0 and pair(q,G,q)<0
    assert mv(A,p)==tuple(scale*x for x in q)
    assert mv(A,(O,*(Z for _ in range(9))))[0]>0
    tr=sum((A[i][i]for i in range(10)),Z)
    assert tr==K(-7378736,-3993344,12596300,6817068)==K.from_json(d['trace'])and tr>10
    # Bind both independent all-profile certificates. Their source hashes
    # fix the four complete first steps needed to connect these sectors.
    audits=[ROOT/'results/fresh_h9_changed_spatial_20261013_independent.json',
            ROOT/'results/fresh_h9_changed_spatial_next_20261013_independent.json']
    older=json.loads(audits[0].read_text());newer=json.loads(audits[1].read_text())
    assert older['status']=='PASS_INDEPENDENT_D9_COMPLETE_EDGE_CHAIN'
    assert newer['status']=='PASS_INDEPENDENT_COMPLETE_D9_A2_B6_FIRST_EDGE'
    for pth,h in older['sources'].items():assert hashlib.sha256((ROOT/pth).read_bytes()).hexdigest()==h
    assert newer['source_sha256']==hashlib.sha256(tp.read_bytes()).hexdigest()
    assert newer['endpoint_type']=='B3+B6'and newer['endpoint_mirrors']==45
    assert all(str(p.relative_to(ROOT))in older['sources']for p in[
        ROOT/'results/fresh_h9_changed_spatial_20261013.json',
        ROOT/'results/fresh_h9_changed_spatial_20261013_step2.json',sp])
    first=json.loads((ROOT/'results/fresh_h9_changed_spatial_20261013.json').read_text())
    second=json.loads((ROOT/'results/fresh_h9_changed_spatial_20261013_step2.json').read_text())
    parse=lambda rows:tuple(tuple(map(K.from_json,row))for row in rows)
    initial=parse(first['initial_walls']);origin=(O,*(Z for _ in range(9)))
    f1=initial[1:]+(tuple(map(K.from_json,first['first_mirror'])),)
    frames=[initial,f1,parse(second['endpoint_incident']),old]
    points=[origin,tuple(map(K.from_json,first['first_endpoint'])),
            tuple(map(K.from_json,second['endpoint'])),p]
    paths=[(frames[i],points[i],frames[i+1],points[i+1])for i in range(3)]+[(initial,origin,new,q)]
    walls=parse(t['all_walls'])
    for i,(start,sv,end,ev)in enumerate(paths):
        shared=set(start)&set(end);arrive=set(end)-set(start);release=set(start)-set(end)
        assert len(shared)==8 and len(arrive)==len(release)==1
        assert pair(sv,G,next(iter(arrive)))<0 and pair(ev,G,next(iter(release)))<0
        assert pair(sv,G,ev)<0
        assert all(pair(sv,G,w)==0 for w in start)and all(pair(ev,G,w)==0 for w in end)
        assert all(pair(v,G,w)<=0 for v in(sv,ev)for w in walls)
        row=d['oriented_path'][i]
        assert tuple(map(K.from_json,row['start']))==sv and tuple(map(K.from_json,row['end']))==ev
        assert set(parse(row['edge_normals']))==shared
    return {'status':'PASS_INDEPENDENT_D9_FULL_CHAMBER_TRANSPORT',
        'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,sp,tp,*audits)},
        'matrix_size':10,'lattice_basis_vectors_per_direction':10,'oriented_simple_normals':9,
        'complete_oriented_path_steps_reconnected':4,
        'positive_point_scale':scale.to_json(),'trace':tr.to_json(),'trace_exceeds_dimension':True,
        'proof_scope':'The D9-integral isometry and inverse preserve the form and future cone and map all nine oriented simple roots between complete B3+B6 corners connected by certified first-edge paths. Hence the same full reflection chamber is preserved. Trace>10 proves infinite order, so discreteness and interior-ball translates imply infinite chamber volume. The certified PARI dependencies are inherited from the first-edge certificates.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();out=check()
    if args.verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
