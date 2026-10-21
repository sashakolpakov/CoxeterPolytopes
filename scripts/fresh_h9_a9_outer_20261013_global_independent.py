#!/usr/bin/env python3
"""Independent Fraction matrix/path certificate of full A9 nonreflectivity."""
import argparse,hashlib,json
from fresh_h9_new_host_next_20261013_a9_independent import a,ALPHA,G,P
from fresh_h9_quartic_first_endpoint_20261013_independent import closure

ROOT=a.ROOT
TRANSPORT=ROOT/'results/fresh_h9_a9_outer_20261013_axis_transports.json'
GEOMETRY=ROOT/'results/fresh_h9_a9_outer_20261013_axis_geometry.json'
FRONTIER=ROOT/'results/fresh_h9_a9_outer_20261013_frontier.json'
FIRST=ROOT/'results/fresh_h9_new_host_next_20261013_a9_first_independent.json'
NEXT=ROOT/'results/fresh_h9_a9_next_20261013_independent.json'
AXIS=ROOT/'results/fresh_h9_a9_outer_20261013_axis_independent.json'
OUT=ROOT/'results/fresh_h9_a9_outer_20261013_global_independent.json'

def determinant(matrix):
    m=[list(r)for r in matrix];d=a.O
    for i in range(len(m)):
        pivot=next(j for j in range(i,len(m))if m[j][i]!=a.Z)
        if pivot!=i:m[i],m[pivot]=m[pivot],m[i];d=a.neg(d)
        z=m[i][i];d=a.mul(d,z)
        for j in range(i+1,len(m)):
            c=a.divide(m[j][i],z)
            for k in range(i,len(m)):m[j][k]=a.sub(m[j][k],a.mul(c,m[i][k]))
    return d

def run(verify=False,replay_path=False):
    if replay_path:
        import fresh_h9_new_host_next_20261013_a9_first_independent as first_replay
        import fresh_h9_a9_next_20261013_independent as next_replay
        import fresh_h9_a9_outer_20261013_axis_independent as axis_replay
        first_replay.run(True);next_replay.run(True);axis_replay.run(True)
    data,geom,fd,first,nxt,axis=[json.loads(p.read_text())for p in(TRANSPORT,GEOMETRY,FRONTIER,FIRST,NEXT,AXIS)]
    source_files=(TRANSPORT,GEOMETRY,FRONTIER,FIRST,NEXT,AXIS);checked=set()
    def ancestry(path):
        if path in checked:return
        checked.add(path);record=json.loads(path.read_text())
        for name,digest in record.get('sources',{}).items():
            pp=ROOT/name;assert hashlib.sha256(pp.read_bytes()).hexdigest()==digest
            if pp.suffix=='.json':ancestry(pp)
    for path in source_files:ancestry(path)
    assert first['status']=='PASS_COMPLETE_A9_FIRST_ENDPOINT_WITH_CERTIFIED_PARI_INTERFACE'
    assert first['mirror_counts_initial_edge_endpoint']==[45,20,21]and first['primitive_reflective_mirrors_in_interval']==1
    assert nxt['status']=='PASS_COMPLETE_FOUR_A9_FIRST_ENDPOINTS'and axis['status']=='PASS_COMPLETE_OUTER_A9_AXIS_FIRST_ENDPOINT'
    for record in(nxt['cases'][0],axis['cases'][0]):
        assert record['source_corner_mirrors']==record['full_endpoint_mirrors']==21
        assert record['primitive_mirrors']==record['nearest_mirrors']==1 and record['exact_orientation_and_sector_join']
    assert a.sign(ALPHA)>0 and all(a.sign(ALPHA,e)<0 for e in(1,2,3))
    walls=[a.vec(r)for r in geom['all_walls']];assert len(walls)==35
    frames=[]
    for row in fd['vertices'][1:]:frames.append(([walls[i]for i in row['incident_indices']],a.vec(row['point']),row['name'],row['incident_indices']))
    for row in geom['actual_finite_vertices']:frames.append(([walls[i]for i in row['certificate']['indices']],a.vec(row['point']),f"axis{row['case']}{'J'if row['reversed']else''}",row['certificate']['indices']))
    assert len(frames)==25
    selected=next(r for r in data['transports']if r.get('kind')=='NEW_INTEGRAL_MAP'and r['source_frame']==23 and r['target_frame']==2)
    H=tuple(a.vec(r)for r in selected['matrix']);hi=a.minverse(H);assert hi==tuple(a.vec(r)for r in selected['inverse'])
    assert all(a.integral(x)for row in H+hi for x in row)
    columns=tuple(zip(*H));assert tuple(tuple(P(r,s)for s in columns)for r in columns)==G
    det=determinant(H);assert det in(a.O,a.neg(a.O))
    E,p,name,eids=frames[23];F,q,label,fids=frames[2];assert name=='axis0'and label=='case0_reversed'
    permutation=selected['permutation'];assert sorted(permutation)==list(range(9))
    assert all(a.mv(H,r)==F[j]for r,j in zip(E,permutation,strict=True))
    rho=a.dec(selected['rho']);assert a.sign(rho)>0 and a.mv(H,p)==tuple(a.mul(rho,x)for x in q)
    assert p[0]==q[0]==a.O and a.sign(P(p,p))<0 and a.sign(P(q,q))<0
    for frame in(E,F):
        assert all(P(r,r)==a.scale(a.O,2)for r in frame)
        a.positive_definite(tuple(tuple(P(r,s)for s in frame)for r in frame))
    trace=a.Z
    for i in range(10):trace=a.add(trace,H[i][i])
    assert trace==a.dec(selected['trace'])and a.sign(a.sub(trace,a.scale(a.O,10)))>0
    assert all(a.sign(a.sub(trace,a.scale(a.O,10)),e)<=0 and a.sign(a.add(trace,a.scale(a.O,10)),e)>=0 for e in(1,2,3))
    class Host:pair=staticmethod(P)
    host=Host();p0=(a.O,)+(a.Z,)*9
    initial=list(range(9));central=fd['vertices'][1];outer=fd['vertices'][2];reverse_outer=fd['vertices'][3]
    pc=a.vec(central['point']);po=a.vec(outer['point']);pr=a.vec(reverse_outer['point'])
    assert first['nearest_root']==geom['all_walls'][10]and nxt['cases'][0]['root']==geom['all_walls'][11]
    assert axis['cases'][0]['root']==geom['all_walls'][33]
    assert p==a.vec(geom['actual_finite_vertices'][0]['point'])and q==pr
    assert len(closure(walls[:9],walls[:9],host))==45
    corners=[('central',pc,central['incident_indices']),('outer',po,outer['incident_indices']),('axis',p,eids),('reversed_outer',pr,reverse_outer['incident_indices'])]
    for nm,point,ids in corners:
        rs=[walls[i]for i in ids];assert len(ids)==9
        assert len(closure(rs,rs,host))==21 and all(P(point,r)==a.Z for r in rs)
        assert point[0]==a.O and a.sign(P(point,point))<0 and all(a.sign(P(point,r))<=0 for r in walls)
    paths=[]
    def join(start,end,source_ids,release,new,source_name,target_name):
        edge=[i for i in source_ids if i!=release];assert len(edge)==8
        assert all(P(start,walls[i])==P(end,walls[i])==a.Z for i in edge)
        assert a.sign(P(start,walls[new]))<0 and P(end,walls[new])==a.Z
        assert P(start,walls[release])==a.Z and a.sign(P(end,walls[release]))<0
        assert start[0]==end[0]==a.O and a.sign(P(start,start))<0 and a.sign(P(end,end))<0
        a.positive_definite(tuple(tuple(P(walls[i],walls[j])for j in edge)for i in edge))
        paths.append({'from':source_name,'to':target_name,'released_wall':release,'first_new_wall':new,'edge_walls':edge,'strict_sector_orientation':True})
    join(p0,pc,initial,4,10,'initial_A9','central')
    join(pc,po,central['incident_indices'],0,11,'central','outer')
    join(po,p,outer['incident_indices'],10,33,'outer','axis')
    join(pc,pr,central['incident_indices'],8,12,'central','reversed_outer')
    # The reversed branch is obtained by an exact integral symmetry fixing
    # the initial complete A9 chamber germ and the central cap.
    J=tuple(tuple(a.O if(i==0 and j==0)or(i>0 and j==10-i)else a.Z for j in range(10))for i in range(10))
    assert a.mv(J,p0)==p0 and a.mv(J,pc)==pc and a.mv(J,po)==pr
    assert [a.mv(J,r)for r in walls[:9]]==list(reversed(walls[:9]))and a.mv(J,walls[10])==walls[10]and a.mv(J,walls[11])==walls[12]
    jcols=tuple(zip(*J));assert tuple(tuple(P(r,s)for s in jcols)for r in jcols)==G
    out={'status':'PASS_FULL_A9_HOST_INFINITE_CHAMBER_STABILIZER','sources':{str(pp.relative_to(ROOT)):hashlib.sha256(pp.read_bytes()).hexdigest()for pp in source_files},
      'host':'Q(T), T^4-4T²+2=0, T=2cos(pi/8); <- (T-1)> + A9',
      'matrix':[list(map(a.enc,row))for row in H],'inverse':[list(map(a.enc,row))for row in hi],
      'determinant':a.enc(det),'integral_both_directions':True,'form_preserving':True,'future_preserving':True,
      'source_corner':{'name':name,'point':list(map(a.enc,p)),'simple_wall_indices':eids,'complete_mirrors':21},
      'target_corner':{'name':label,'point':list(map(a.enc,q)),'simple_wall_indices':fids,'complete_mirrors':21},
      'positive_scale':a.enc(rho),'oriented_simple_root_permutation':permutation,'certified_same_chamber_paths':paths,
      'trace':a.enc(trace),'trace_strictly_greater_than_dimension':10,'other_embedding_traces_between_minus10_and10':True,
      'source_ancestry_json_files_checked':len(checked),
      'argument':'Each path edge is the first mirror of the full root arrangement, certified over the complete norm/profile/orbit interface. The complete initial A9 germ and the complete oriented endpoint stabilizers therefore belong to one full W(L) chamber. The integral future form isometry permutes all nine outward simple normals of the source germ to the target germ, so local uniqueness makes it preserve that entire chamber. Trace>10 proves infinite order; the arithmetic lattice isometry group is discrete, hence an infinite chamber stabilizer forces infinite hyperbolic volume. The full host is nonreflective.',
      'CAS_dependency':'First-edge completeness uses the explicitly archived certified PARI full-unit and absolute-norm interfaces. This matrix, path, closure, orientation, and trace replay uses separate Fraction arithmetic and no PARI, NumPy, SymPy, or producer transport arithmetic.'}
    if verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'trace>10','path joins',len(paths),'complete mirrors',21,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--replay-path',action='store_true');args=p.parse_args();run(args.verify,args.replay_path)
