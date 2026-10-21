#!/usr/bin/env python3
"""Fraction-only infinite full-chamber stabilizer, including complete connecting edges."""
from fresh_h9_quartic_first_endpoint_20261013_independent import *
from fresh_h9_quartic_edge_walk_20261013_independent import mm,ident

SYMMETRIES='results/fresh_h9_quartic_next_frontier_20261013_symmetries.json'
OUT=ROOT/'results/fresh_h9_quartic_next_frontier_20261013_global_independent.json'


def transpose(a):return tuple(zip(*a))


def lift_matrix(F):
    inds=(0,1,2,3,4,9);out=[list(row)for row in ident(10)]
    for i,ii in enumerate(inds):
        for j,jj in enumerate(inds):out[ii][jj]=F[i][j]
    return tuple(map(tuple,out))


def global_replay():
    started=time.monotonic()
    # Re-enumerate every root meeting the two complete connecting intervals.
    step_audit=[replay(1),replay(2)]
    paths=['results/fresh_h9_quartic_first_endpoint_20261013'+s+'.json'for s in('','_second')]
    steps=[json.loads((ROOT/p).read_text())for p in paths]
    assert steps[0]['nearest_point']==steps[1]['start_corner']
    rfirst=[*matrix(steps[0]['edge_roots']),vec(steps[0]['nearest_root'])]
    base=json.loads((ROOT/steps[1]['base_source']).read_text())['cases'][0]
    rsecond=[*matrix(steps[1]['edge_roots']),vec(base['forced_H5_face']['all_walls'][steps[1]['start_corner_wall']])]
    assert set(rfirst)==set(rsecond)
    assert all(row['mirrors_meeting_open_segment']==0 for row in step_audit)
    sector_checks=[];host=Host(1)
    for step in steps:
        previous=json.loads((ROOT/step['base_source']).read_text())
        if step['base_source_case']is not None:previous=previous['cases'][step['base_source_case']]
        start_wall=vec(previous['forced_H5_face']['all_walls'][step['start_corner_wall']])
        next_wall=vec(step['nearest_root']);k0=vec(step['start_corner']);k1=vec(step['nearest_point'])
        start_on_next=host.pair(k0,next_wall);next_on_start=host.pair(k1,start_wall)
        assert sign(start_on_next)<0 and sign(next_on_start)<0
        assert all(host.pair(k0,e)==Z and host.pair(k1,e)==Z for e in matrix(step['edge_roots']))
        sector_checks.append({'edge_indices':step['edge_indices'],'same_oriented_edge_normals':True,
                              'start_point_on_next_wall_negative_side':enc(start_on_next),
                              'next_point_on_released_wall_negative_side':enc(next_on_start),
                              'open_edge_has_no_other_mirror':True})
    data=json.loads((ROOT/SYMMETRIES).read_text())
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in data['sources_sha256'].items())
    host=Host(1);host9=Host(2);autos=[]
    for index in(2,4):
        saved=data['distinct_integral_automorphisms'][index];corner=data['corners'][saved['corner_index']]
        assert (corner['source'],corner['kind'])==((paths[0],'start')if index==2 else(paths[1],'nearest'))
        step=steps[0 if index==2 else 1]
        if index==2:
            base=json.loads((ROOT/step['base_source']).read_text())
            root=vec(base['forced_H5_face']['all_walls'][step['start_corner_wall']]);point=vec(step['start_corner'])
        else:root=vec(step['nearest_root']);point=vec(step['nearest_point'])
        roots=(*matrix(step['edge_roots']),root)
        assert roots==matrix(corner['roots'])and point==vec(corner['point'])
        gram=tuple(tuple(host.pair(r,s)for s in roots)for r in roots);pd(gram)
        assert sign(host.pair(point,point))<0 and sign(point[0])>0
        assert all(host.pair(point,r)==Z for r in roots)
        permutation=saved['permutation'];assert sorted(permutation)==list(range(5))
        F=matrix(saved['matrix']);Fi=matrix(saved['inverse'])
        assert all(integral(z)for row in F+Fi for z in row)
        assert mm(F,Fi)==ident(6)and mm(F,F)==ident(6)
        assert mm(mm(transpose(F),host.g),F)==host.g
        assert mv(F,point)==point and all(mv(F,r)==roots[permutation[i]]for i,r in enumerate(roots))
        frame=transpose((*roots,point));assert mm(F,frame)==transpose((*(roots[i]for i in permutation),point))
        full=lift_matrix(F);fulli=lift_matrix(Fi)
        assert mm(full,fulli)==ident(10)and mm(mm(transpose(full),host9.g),full)==host9.g
        assert all(integral(z)for row in full+fulli for z in row)
        lift=lambda r:(r[0],*r[1:5],*(Z,)*4,r[5])
        assert mv(full,lift(point))==lift(point)
        other=[tuple(neg(O)if i==j else Z for i in range(10))for j in range(5,9)]
        assert all(mv(full,r)==r for r in other)
        stabilizer=closure([*map(lift,roots),*other],[*map(lift,roots),*other],host9)
        expected=step['H9_start_stabilizer_mirrors']if index==2 else step['H9_end_stabilizer_mirrors']
        assert len(stabilizer)==expected
        autos.append({'source_corner':corner['source'],'kind':corner['kind'],'permutation':permutation,
                      'H5_matrix':saved['matrix'],'H9_matrix':[[enc(z)for z in row]for row in full],
                      'full_corner_reflection_count':expected,'integral_inverse':True,'future_preserving':True,
                      'fixes_corner_and_permutes_full_simple_system':True,'involution':True})
    # The factors are local chamber symmetries at corners connected by the two
    # completely enumerated mirror-free edges. Thus they preserve one full W(L)
    # chamber, not merely the finite saved wall prefix.
    A,B=(matrix(row['H5_matrix'])for row in autos);product5=mm(B,A)
    trace=Z
    for i in range(6):trace=add(trace,product5[i][i])
    expected=(Q(-14),Q(-8),Q(26),Q(14));assert trace==expected and sign(sub(trace,scale(O,6)))>0
    product9=lift_matrix(product5);trace9=add(trace,scale(O,4))
    assert sign(sub(trace9,scale(O,10)))>0
    assert mm(mm(transpose(product9),host9.g),product9)==host9.g
    assert all(integral(z)for row in product9+a.minverse(product9)for z in row)
    forged=[list(row)for row in product9];forged[0][0]=add(forged[0][0],O);forged=tuple(map(tuple,forged))
    assert mm(mm(transpose(forged),host9.g),forged)!=host9.g
    return {'status':'PASS_INDEPENDENT_INFINITE_FULL_CHAMBER_STABILIZER',
            'sources_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in [*paths,SYMMETRIES]},
            'complete_connecting_edge_audits':[{k:v for k,v in row.items()if k!='seconds'}for row in step_audit],
            'oriented_sector_path':sector_checks,'middle_corner_oriented_simple_roots_match':True,
            'automorphisms':autos,'H5_product':[[enc(z)for z in row]for row in product5],
            'H9_product':[[enc(z)for z in row]for row in product9],'H5_trace':enc(trace),'H9_trace':enc(trace9),
            'trace_above_dimension':True,'same_full_chamber_via_certified_edges':True,'forged_matrix_rejected':True,
            'seconds':round(time.monotonic()-started,3),
            'scope':'A discrete integral isometry group contains an infinite-order stabilizer of its full reflection chamber: two future-preserving integral involutions fix certified corners and permute their complete local simple systems; two complete nearest-edge censuses connect these sectors. Therefore the full H5 and H9 reflection chambers have infinite volume. No finite-volume reflection-subgroup polytope in these lattices can exist.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=global_replay()
    if args.verify:
        saved=json.loads(OUT.read_text());saved.pop('seconds');result.pop('seconds');assert saved==result
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print({k:result[k]for k in('status','H5_trace','H9_trace','trace_above_dimension','same_full_chamber_via_certified_edges')},flush=True)
