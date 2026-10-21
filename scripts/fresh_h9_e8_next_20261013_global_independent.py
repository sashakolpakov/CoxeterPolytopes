#!/usr/bin/env python3
"""Separate Fraction reconstruction of the E8+1 infinite chamber stabilizer."""
from fresh_h9_e8_first_20261013_independent import *
TAG='fresh_h9_e8_next_20261013_global';OUT=ROOT/'results'/f'{TAG}_independent.json'

def mm(A,B):return tuple(tuple(old.summ(a.mul(A[i][k],B[k][j])for k in range(len(B)))for j in range(len(B[0])))for i in range(len(A)))
def transport(src,target):return mm(tuple(zip(*target)),a.minverse(tuple(zip(*src))))
def closure(frame):
    n=len(frame);fg=[[P(r,s)for s in frame]for r in frame]
    seen={tuple(a.O if i==j else a.Z for i in range(n))for j in range(n)};todo=list(seen)
    while todo:
        v=todo.pop()
        for j in range(n-1,-1,-1):
            c=a.divide(a.scale(old.summ(a.mul(v[k],fg[k][j])for k in range(n)),2),fg[j][j]);w=list(v);w[j]=a.sub(w[j],c);w=tuple(w)
            if w not in seen:seen.add(w);todo.append(w)
        assert len(seen)<=242
    assert all(tuple(map(a.neg,v))in seen for v in seen)
    return len(seen)//2

def run(verify=False):
    source=ROOT/'results'/f'{TAG}.json';g=json.loads(source.read_text())
    for p,h in g['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    names=[PREFIX,PREFIX+'_independent',PREFIX+'_e7',PREFIX+'_e7_independent','fresh_h9_e8_next_20261013','fresh_h9_e8_next_20261013_independent','fresh_h9_e8_next_20261013_d6','fresh_h9_e8_next_20261013_d6_independent','fresh_h9_e8_next_20261013_d7','fresh_h9_e8_next_20261013_d7_independent']
    paths=[ROOT/'results'/f'{n}.json'for n in names];axis,ai,e7,ei,scalar,si,d6,di,d7,ji=[json.loads(p.read_text())for p in paths]
    for obj in(axis,ai,e7,ei,scalar,si,d6,di,d7,ji):
        assert obj['status'].startswith('PASS')
        for p,h in obj['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    assert(si['profiles'],si['complete_interval_points'],di['profiles'],di['complete_interval_points'],ji['profiles'],ji['complete_interval_points'])==(11,2,24,17,14,18)
    basis=old.matrix(axis['lattice_basis']);B=tuple(zip(*basis));BI=a.minverse(B);walls=old.matrix(d7['combined_walls'])
    origin=(a.O,*(a.Z for _ in range(9)));initial=list(walls[:9]);dp=a.vec(d6['endpoint']);jp=a.vec(d7['endpoint']);ep=a.vec(e7['endpoint'])
    df=[walls[i]for i in[2,3,4,5,6,7,8,10,12]];jf=[walls[i]for i in[1,2,3,4,5,6,7,8,13]]
    assert tuple(df)==old.matrix(d6['endpoint_incident'])and tuple(jf)==old.matrix(d7['endpoint_incident'])
    hperm=[0,1,2,3,5,4,6,8,7];jperm=[8,0,1,2,3,4,6,5,7]
    scales=[s for s in old.quartic_roots(a.divide(a.neg(BETA),P(jp,jp)))if a.sign(s)>0];assert len(scales)==1;scale=scales[0]
    H=transport(df+[dp],[df[i]for i in hperm]+[dp]);J=transport(initial+[origin],[jf[i]for i in jperm]+[tuple(a.mul(scale,x)for x in jp)])
    I=tuple(tuple(a.O if i==j else a.Z for j in range(10))for i in range(10))
    G=tuple(tuple(a.neg(BETA)if i==j==0 else a.O if i==j else a.Z for j in range(10))for i in range(10))
    for rec,M,src,target,point,topoint,perm in zip(g['generators'],[H,J],[df,initial],[[df[i]for i in hperm],[jf[i]for i in jperm]],[dp,origin],[dp,tuple(a.mul(scale,x)for x in jp)],[hperm,jperm]):
        MI=a.minverse(M);lm=mm(BI,mm(M,B));li=mm(BI,mm(MI,B))
        assert all(a.integral(x)for row in lm for x in row)and all(a.integral(x)for row in li for x in row)
        assert mm(M,M)==I and mm(tuple(zip(*M)),mm(G,M))==G
        assert all(a.mv(M,r)==s for r,s in zip(src,target))and a.mv(M,point)==topoint
        assert a.sign(point[0])>0 and a.sign(topoint[0])>0 and a.sign(P(point,point))<0
        assert old.matrix(rec['matrix'])==M and old.matrix(rec['inverse'])==MI
        assert old.matrix(rec['lattice_matrix'])==lm and old.matrix(rec['lattice_inverse'])==li
        assert rec['permutation']==perm and a.vec(rec['positive_root_scales'])==(a.O,)*9
    F=mm(H,J);FI=mm(J,H);tr=old.summ(F[i][i]for i in range(10))
    assert F==old.matrix(g['product'])and FI==old.matrix(g['product_inverse'])and mm(F,FI)==I
    assert tr==a.dec(g['trace'])and a.sign(a.sub(tr,a.scale(a.O,10)))>0
    counts=[closure(initial),closure(df),closure(jf)];assert counts==g['complete_corner_mirror_counts']==[121,33,121]
    def integral(r):return all(map(a.integral,a.mv(BI,r)))
    for r in set(initial+df+jf):
        q=P(r,r);assert a.sign(q)>0 and integral(r)
        assert not any(integral(tuple(a.divide(x,p)for x in r))for p in(a.T,a.scale(a.O,3)))
        assert all(a.integral(a.divide(a.scale(P(r,b),2),q))for b in basis)
    path_specs=[('INITIAL_E7',origin,ep,[0,2,3,4,5,6,7,8],1,10),('E7_TO_D6',ep,dp,[2,3,4,5,6,7,8,10],0,12),('INITIAL_D7',origin,jp,list(range(1,9)),0,13)]
    for rec,(label,start,end,retained,released,new)in zip(g['oriented_same_chamber_paths'],path_specs):
        assert all(a.sign(P(p,p))<0 and a.sign(p[0])>0 for p in(start,end))
        assert all(P(p,walls[i])==a.Z for p in(start,end)for i in retained)
        assert P(start,walls[released])==P(end,walls[new])==a.Z
        assert a.sign(P(start,walls[new]))<0 and a.sign(P(end,walls[released]))<0
        assert all(a.sign(P(p,r))<=0 for p in(start,end)for r in walls)
        assert rec['name']==label and rec['retained_indices']==retained and rec['released_index']==released and rec['new_index']==new
        assert a.vec(rec['start'])==start and a.vec(rec['end'])==end
    center=a.vec(d7['strict_interior']);escape=a.vec(d7['strict_escape'])
    assert a.sign(P(center,center))<0 and a.sign(P(escape,escape))>0 and a.sign(P(center,escape))<0
    assert all(a.sign(P(v,r))<0 for v in(center,escape)for r in walls)
    result={'status':'PASS_INDEPENDENT_FULL_E8_PLUS_ONE_NONREFLECTIVITY','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(source,*paths,Path(__file__).resolve(),Path(a.__file__).resolve(),Path(old.__file__).resolve())},
        'independently_reconstructed_generators':2,'generator_involutions':True,'actual_E8_module_integral_both_directions':True,'positively_oriented_simple_frames':True,
        'full_corner_mirror_counts':counts,'certified_path_segments':3,'product_trace':a.enc(tr),'trace_exceeds_dimension':True,'finite_prefix_facets':14,
        'new_profiles':49,'new_complete_interval_points':37,'all_certified_edges':5,'all_profiles':71,'all_complete_interval_points':46,
        'scope':'No producer K arithmetic is imported. Separate Fraction arithmetic reconstructs both matrices from the source and target corner frames, checks form, inverses, actual E8 module, oriented paths, complete mirror closures, trace and finite-prefix geometry. Independent complete edge replays are source-bound; certified PARI completeness remains the shared stated dependency.'}
    if verify:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'corners',counts,'trace>10',True,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');aa=p.parse_args();run(aa.verify)
