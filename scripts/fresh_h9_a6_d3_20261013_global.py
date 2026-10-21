#!/usr/bin/env python3
"""Full-chamber infinite stabilizer for the prime-seven quartic A6+D3 host."""
from fresh_h9_a6_d3_20261013 import *
TAG='fresh_h9_a6_d3_20261013_global';OUT=ROOT/'results'/f'{TAG}.json'
def mm(A,B):return tuple(tuple(sum((A[i][k]*B[k][j]for k in range(len(B))),Z)for j in range(len(B[0])))for i in range(len(A)))
def closure(frame):
    C=tuple(tuple(pair(r,G,s)for s in frame)for r in frame);n=len(frame);seen={tuple(O if i==j else Z for i in range(n))for j in range(n)};todo=list(seen)
    while todo:
        v=todo.pop()
        for j in range(n):
            c=2*sum((v[k]*C[k][j]for k in range(n)),Z)/C[j][j];w=list(v);w[j]-=c;w=tuple(w)
            if w not in seen:seen.add(w);todo.append(w)
        assert len(seen)<=60
    assert len(seen)==60 and all(tuple(-x for x in v)in seen for v in seen)
    return len(seen)//2

def run(verify=False):
    prefix='fresh_h9_a6_d3_20261013';paths=[ROOT/'results'/f'{prefix}{s}.json'for s in('','_second','_third','_a5','_independent')];docs=[json.loads(p.read_text())for p in paths];first,second,third,end,audit=docs
    assert audit['status']=='PASS_INDEPENDENT_A6_D3_FOUR_COMPLETE_FIRST_EDGES'
    for d in docs:
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    for row in audit['edges']:assert hashlib.sha256((ROOT/row['source']).read_bytes()).hexdigest()==row['source_sha256']
    assert [(r['profiles'],r['complete_interval_points'])for r in audit['edges']]==[(15,6),(31,13),(65,58),(71,2)]
    walls=[tuple(map(K.from_json,r))for r in end['combined_walls']];origin=STD[0];point=tuple(map(K.from_json,end['endpoint']));target=[walls[i]for i in[1,2,3,4,5,12,6,7,8]]
    assert end['endpoint_mirrors']==30 and len(end['primitive_reflective_survivors'])==1
    assert set(target)==set(tuple(map(K.from_json,r))for r in end['endpoint_incident'])
    scale=square_root(-BETA/pair(point,G,point));assert scale is not None and scale>0
    F=mm(tuple(zip(*(target+[tuple(scale*x for x in point)]))),inverse(tuple(zip(*(INITIAL+[origin])))));FI=inverse(F)
    I=tuple(tuple(O if i==j else Z for j in range(10))for i in range(10));B=tuple(zip(*BASIS));module=mm(BI,mm(F,B));modulei=mm(BI,mm(FI,B))
    assert mm(F,FI)==I and mm(FI,F)==I and mm(tuple(zip(*F)),mm(G,F))==G
    assert all(x.is_integral for row in module for x in row)and all(x.is_integral for row in modulei for x in row)
    assert all(mv(F,r)==s for r,s in zip(INITIAL,target))and mv(F,origin)==tuple(scale*x for x in point)
    assert point[0]>0 and pair(point,G,point)<0
    for r in (*INITIAL,*target):assert reason(r,pair(r,G,r))is None
    assert all(pair(p,G,walls[i])==0 for p in(origin,point)for i in range(1,9))
    assert pair(origin,G,walls[0])==pair(point,G,walls[12])==0
    assert pair(origin,G,walls[12])<0 and pair(point,G,walls[0])<0
    assert all(pair(p,G,r)<=0 for p in(origin,point)for r in walls)
    tr=sum((F[i][i]for i in range(10)),Z);assert tr>10
    # D3 is fixed pointwise; the transport is supported on time plus A6.
    assert all(mv(F,STD[i])==STD[i]for i in(7,8,9))
    out={'status':'PASS_FULL_PRIME_SEVEN_A6_D3_NONREFLECTIVITY','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(*paths,Path(__file__).resolve(),ROOT/'scripts/fresh_h9_glued_host_20261013_arithmetic.py',Path(square_root.__code__.co_filename).resolve())},
         'host':first['host'],'coefficient':BETA.to_json(),'ambient_gram':menc(G),'lattice_basis':menc(BASIS),'transport_matrix':menc(F),'inverse_matrix':menc(FI),'lattice_matrix':menc(module),'lattice_inverse':menc(modulei),'trace':tr.to_json(),'trace_minus_dimension':(tr-10).to_json(),
         'source_incident':menc(INITIAL),'target_incident':menc(target),'source_point':enc(origin),'target_point':enc(point),'positive_point_scale':scale.to_json(),'positive_root_scales':enc((O,)*9),'full_corner_types':['A6+B3','A6+B3'],'full_corner_mirror_counts':[closure(INITIAL),closure(target)],
         'oriented_certified_path':{'first_edge_source':str(paths[3].relative_to(ROOT)),'first_edge_audit':str(paths[4].relative_to(ROOT)),'retained_indices':list(range(1,9)),'released_index':0,'new_index':12,'source_pairing_new':pair(origin,G,walls[12]).to_json(),'target_pairing_released':pair(point,G,walls[0]).to_json(),'all_thirteen_walls_nonpositive':True},
         'complete_edges':4,'integral_profiles':sum(r['profiles']for r in audit['edges']),'interval_points':sum(r['complete_interval_points']for r in audit['edges']),'finite_prefix_facets':13,'finite_prefix_strict_escape':end['strict_escape'],'finite_prefix_escape_norm':end['escape_norm'],
         'theorem':'The full root reflection chamber of <-(T-1)(2+T)>+A6+D3 over Z[T], T^4-4T^2+2=0, has infinite volume. The initial zero chamber is complete A6+B3. The complete first-contact census on its A5+B3 edge gives the full finite endpoint sector. This exact future-preserving lattice isometry maps every oriented initial simple normal to the endpoint simple normals. The two sectors belong to the same full chamber, so the isometry preserves that chamber. Its trace exceeds ten; it has infinite order. Discreteness then forces infinite volume.',
         'scope':'A full-host nonreflectivity theorem for this determinant-squareclass-distinct lattice. The PARI-certified absolute norm/unit interface and the finite reflection classification are explicit shared dependencies. No conclusion about all compact H9 polytopes or arbitrary changed lattices is asserted.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'trace>10',tr>10,'edges',out['complete_edges'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');aa=p.parse_args();run(aa.verify)
