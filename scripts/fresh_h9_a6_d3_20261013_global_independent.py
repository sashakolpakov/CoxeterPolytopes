#!/usr/bin/env python3
"""Independent Fraction matrix and chamber-path replay; no producer K arithmetic."""
from fresh_h9_a6_d3_20261013_independent import *
BASE=TAG;TAG=BASE+'_global_independent';OUT=ROOT/'results'/f'{TAG}.json'
def mm(A,B):return tuple(tuple(old.summ(a.mul(A[i][k],B[k][j])for k in range(len(B)))for j in range(len(B[0])))for i in range(len(A)))
def closure(frame):
    C=tuple(tuple(P(r,s)for s in frame)for r in frame);seen={tuple(int(i==j)for i in range(9))for j in range(9)};todo=list(seen)
    while todo:
        v=todo.pop()
        for j in range(8,-1,-1):
            c=a.divide(a.scale(old.summ(a.scale(C[k][j],v[k])for k in range(9)),2),C[j][j]);assert c[1:]==(0,0,0)and c[0].denominator==1
            w=list(v);w[j]-=int(c[0]);w=tuple(w)
            if w not in seen:seen.add(w);todo.append(w)
        assert len(seen)<=60
    assert len(seen)==60 and all(tuple(-x for x in v)in seen for v in seen)
    return 30

def run(verify=False):
    paths=[ROOT/'results'/f'{BASE}{s}.json'for s in('_global','_a5','_independent')];glob,end,audit=[json.loads(p.read_text())for p in paths]
    assert audit['status']=='PASS_INDEPENDENT_A6_D3_FOUR_COMPLETE_FIRST_EDGES'
    for d in(glob,end,audit):
        for path,h in d['sources'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
    for row in audit['edges']:assert hashlib.sha256((ROOT/row['source']).read_bytes()).hexdigest()==row['source_sha256']
    assert [(r['profiles'],r['complete_interval_points'])for r in audit['edges']]==[(15,6),(31,13),(65,58),(71,2)]
    walls=old.matrix(end['combined_walls']);origin=STD[0];point=a.vec(end['endpoint']);target=tuple(walls[i]for i in[1,2,3,4,5,12,6,7,8])
    assert set(target)==set(old.matrix(end['endpoint_incident']))and len(end['primitive_reflective_survivors'])==1
    scales=[s for s in old.quartic_roots(a.divide(a.neg(BETA),P(point,point)))if a.sign(s)>0];assert len(scales)==1;sc=scales[0]
    src=tuple(zip(*(*INITIAL,origin)));dst=tuple(zip(*(*target,tuple(a.mul(sc,x)for x in point))));F=mm(dst,a.minverse(src));FI=a.minverse(F)
    I=tuple(tuple(a.O if i==j else a.Z for j in range(10))for i in range(10));B=tuple(zip(*BASIS));module=mm(BI,mm(F,B));modulei=mm(BI,mm(FI,B))
    assert F==old.matrix(glob['transport_matrix'])and FI==old.matrix(glob['inverse_matrix'])
    assert mm(F,FI)==mm(FI,F)==I and mm(tuple(zip(*F)),mm(G,F))==G
    assert all(a.integral(x)for row in(*module,*modulei)for x in row)
    assert module==old.matrix(glob['lattice_matrix'])and modulei==old.matrix(glob['lattice_inverse'])
    assert all(a.mv(F,r)==s for r,s in zip(INITIAL,target))and a.mv(F,origin)==tuple(a.mul(sc,x)for x in point)
    assert sc==a.dec(glob['positive_point_scale'])and all(s==a.O for s in a.vec(glob['positive_root_scales']))
    assert a.sign(point[0])>0 and a.sign(P(point,point))<0
    assert all(reason(r,P(r,r))is None for r in(*INITIAL,*target))
    assert all(a.mv(F,STD[i])==STD[i]for i in(7,8,9))
    assert all(P(p,walls[i])==a.Z for p in(origin,point)for i in range(1,9))
    assert P(origin,walls[0])==P(point,walls[12])==a.Z
    assert a.sign(P(origin,walls[12]))<0 and a.sign(P(point,walls[0]))<0
    assert all(a.sign(P(p,r))<=0 for p in(origin,point)for r in walls)
    tr=old.summ(F[i][i]for i in range(10));assert tr==a.dec(glob['trace'])and a.sign(a.sub(tr,a.scale(a.O,10)))>0
    counts=[closure(INITIAL),closure(target)];assert counts==glob['full_corner_mirror_counts']
    center,escape=map(a.vec,(end['strict_interior'],end['strict_escape']))
    assert a.sign(P(center,center))<0 and a.sign(P(escape,escape))>0 and a.sign(P(center,escape))<0
    assert all(a.sign(P(v,r))<0 for v in(center,escape)for r in walls)
    out={'status':'PASS_INDEPENDENT_FULL_A6_D3_NONREFLECTIVITY_TRANSPORT','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(*paths,Path(__file__).resolve(),ROOT/'scripts'/f'{BASE}_independent.py',Path(old.__file__).resolve(),Path(a.__file__).resolve())},'trace':a.enc(tr),'trace_gt_ten':True,'actual_lattice_and_inverse_integral':True,'positive_sector_transport':True,'full_corner_mirrors':counts,'complete_first_edge_profiles':71,'complete_first_edge_points':2,'first_contact_primitive_mirrors':1,'all_certified_edges':4,'all_integral_profiles':182,'all_interval_points':79,'finite_prefix_facets':13,'strict_escape_norm':a.enc(P(escape,escape)),
         'scope':'A separately reconstructed Fraction matrix verifies the actual lattice and inverse, form, trace, positive corner-sector map, full mirror closures, and certified same-chamber path. Its full first-edge premise comes from the separate exhaustive Fraction profile/orbit replay; certified PARI norm/unit completeness remains explicit. This establishes full-host infinite volume, not merely a prefix escape.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'trace>10',out['trace_gt_ten'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');aa=p.parse_args();run(aa.verify)
