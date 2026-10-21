#!/usr/bin/env python3
"""Exact finite-anchor mixed-Lorentz-face search and changed fourteen-wall cones.

No positive-definite inequality is applied to a Lorentzian two-plane. The
finite search covers spacelike projections with two independent finite-angle
contacts, not every possible projected normal. Full changed-cone checks prove
that the four surviving roots give infinite-volume, not compact, intersections.
"""
import argparse,json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations,product
import fresh_raw23_four_face_contact_preflight_20261018 as F
P,H,Q=F.P,F.H,F.Q
ROOT=H.ROOT
TAG='fresh_raw23_lorentz_two_face_native_roots_20261018'
PINS={
 'scripts/fresh_raw23_four_face_contact_preflight_20261018.py':'6d5421bbd340c77b9a052360d47c653b5af3bc0a66b4a33d3d360d915915fb02',
 'scripts/fresh_raw23_union_polar_completion_preflight_20261018.py':'801198f128ee0e8d3ff9587fbd65d656edaedc0d413e4b4665c094c290010ee8',
 'results/fresh_raw23_union_polar_completion_preflight_20261018.json':'eee8c3836a16c6342b2882c3b65d84999bb7fbc7363413fbc59626c739ff8361',
 'results/fresh_raw23_union_polar_completion_preflight_20261018.md':'42bd1f4baeecd9c8188a1590adf1e689c5f9582ae5df378a75d94dc5116da973',
 'scripts/check_raw23_same_family_saved_20261018.py':'5b1ab9873b2e069f18966a7052856d8a9df8a8c1b0eef4794522a812a54616a7',
}
def add(x,y):return[a+b for a,b in zip(x,y)]
def scale(x,q):return[q*a for a in x]
def products(rows,x):return[sum((a*b for a,b in zip(row,x)),Q())for row in rows]
def integral(q):return(2*q.a).denominator==1 and(q.a*q.a-5*q.b*q.b).denominator==1

def inventory(o):
    rays=o['initial_cone']['rays'];xyz=[F.vec(r['coordinates'])for r in rays]
    g=H.mat(o['metric']);rows=H.mat(o['initial_cone']['normal_rows']);out=[]
    assert len(rays)==50 and[r['index']for r in rays]==list(range(50))
    for i,j in combinations(range(50),2):
        common=sorted(set(rays[i]['zero_facets'])&set(rays[j]['zero_facets']))
        if len(common)<8 or F.rank([rows[k]for k in common])!=8:continue
        full=[k for k,r in enumerate(rays)if set(common)<=set(r['zero_facets'])]
        assert full==[i,j]
        gram=[[P.inner(xyz[a],g,xyz[b])for b in(i,j)]for a in(i,j)]
        contacts=[[-x for x in F.vec(rays[k]['products'])]for k in(i,j)]
        out.append(dict(rays=[i,j],active=common,signs=[rays[k]['norm_sign']for k in(i,j)],
            gram=gram,inertia=H.inertia(gram),contacts=contacts,
            nonzero=[k for k in range(12)if contacts[0][k]!=0 or contacts[1][k]!=0]))
    cc=Counter(tuple(r['signs'])for r in out)
    assert cc=={(1,1):112,(-1,1):64,(-1,-1):49}
    assert all(r['inertia']==([2,0,0]if r['signs']==[1,1]else[1,1,0])for r in out)
    # All causal generators lie on the same future sheet as the strict center.
    center=F.vec(o['initial_cone']['strict_timelike_center'])
    assert all(P.inner(x,g,center).sign()<0 for x,r in zip(xyz,rays)if r['norm_sign']<0)
    return out

def anchors(o,faces):
    fields={};sc=F.vec(o['original_squared_scales']);counter=Counter();records=[];survivors=[]
    for face in faces:
        if face['signs']!=[-1,1]:continue
        h,p=face['gram'],face['contacts'];assert len(face['nonzero'])==4
        pairs=[]
        for ia,ib in combinations(face['nonzero'],2):
            det=p[0][ia]*p[1][ib]-p[0][ib]*p[1][ia]
            if det==0:continue
            pairs.append([ia,ib]);decisions=[]
            for fi,fj in product(range(5),repeat=2):
                fa,fb=P.FINITE_SQUARES[fi+1],P.FINITE_SQUARES[fj+1]
                qa,qb=sc[ia]*fa,sc[ib]*fb;key=(qa.a,qa.b,qb.a,qb.b)
                if key not in fields:fields[key]=H.Field([qa,qb])
                k=fields[key];ca,cb=k.sqrt(qa),k.sqrt(qb)
                a=scale(add(scale(ca,p[1][ib]),scale(cb,-p[1][ia])),Q(1)/det)
                b=scale(add(scale(cb,p[0][ia]),scale(ca,-p[0][ib])),Q(1)/det)
                status=None;counter['tuples']+=1
                if k.sign(a)<=0 or k.sign(b)<=0:status='outside_relative_interior'
                else:
                    norm=add(add(scale(k.mul(a,a),h[0][0]),scale(k.mul(a,b),2*h[0][1])),scale(k.mul(b,b),h[1][1]))
                    if k.sign(norm)<=0:status='causal_projection'
                    elif k.sign(add(norm,k.scalar(-2)))>0:status='norm_above2'
                    else:
                        contacts=[add(scale(a,p[0][j]),scale(b,p[1][j]))for j in range(12)]
                        squares=[scale(k.mul(v,v),Q(1)/s)for v,s in zip(contacts,sc)]
                        if any(any(v!=0 for v in z[1:])for z in squares):status='nonbase_two_cycle'
                        elif any(P.square_gate(z[0])for z in squares):status='native_failure'
                        else:
                            # No conjugate-projection/norm gate is used to exclude candidates.
                            assert not k.basis and norm==[Q(2)]
                            status='survivor'
                            survivors.append(dict(rays=face['rays'],anchors=[ia,ib],finite_square_indices=[fi,fj],
                                coefficients=[a[0],b[0]],norm=norm[0],physical_squares=[z[0]for z in squares]))
                counter[status]+=1;decisions.append(status)
            records.append(dict(rays=face['rays'],anchors=[ia,ib],decisions=decisions))
        assert len(pairs)==5
    assert counter['tuples']==8000
    assert len(records)==320 and len(survivors)==(0 if o['kinds']==[1,6]else 2)
    return dict(counts=dict(counter),tests=records,survivors=survivors)

def cycles(gram,scales):
    n=len(gram);adj={i:[j for j in range(n)if i!=j and gram[i][j]!=0]for i in range(n)};out=[]
    for start in range(n):
        def visit(path,value):
            i=path[-1]
            for j in adj[i]:
                if j==start and len(path)>=3 and path[1]<path[-1]:
                    z=value*gram[i][start]/scales[i];assert integral(z)
                    out.append(dict(path=path,value=z))
                elif j>start and j not in path:visit(path+[j],value*gram[i][j]/scales[i])
        visit([start],Q(1))
    return out

def roots(o,tests):
    g=H.mat(o['metric']);rows=H.mat(o['initial_cone']['normal_rows'])
    rays=o['initial_cone']['rays'];sc=F.vec(o['original_squared_scales']);out=[]
    for s in tests['survivors']:
        a,b=s['coefficients'];u,v=[F.vec(rays[i]['coordinates'])for i in s['rays']]
        x=add(scale(u,a),scale(v,b));assert P.inner(x,g,x)==2
        raw=products(rows,x);assert all(v.sign()<=0 for v in raw)
        assert [v*v/z for v,z in zip(raw,sc)]==s['physical_squares']
        finite=[P.inner(x,g,F.vec(r['coordinates']))for r in rays if r['norm_sign']<0]
        assert len(finite)==18 and all(v.sign()<0 for v in finite)
        cuts=[r['index']for r in rays if r['norm_sign']>0 and P.inner(x,g,F.vec(r['coordinates'])).sign()>0]
        assert cuts==[r['index']for r in rays if r['norm_sign']>0]
        pole=rays[o['initial_polar_tests']['passing_rays'][0]];pn=Q(*pole['norm'])
        z=P.inner(x,g,F.vec(pole['coordinates']));assert z.sign()>0 and 2*z*z/pn==1
        assert sum(0<v.sign()and(v-4).sign()<0 for v in s['physical_squares'])==2
        out.append(dict(source=s,coordinates=x,raw_contacts=raw,finite_vertex_products=finite,
            cut_original_spacelike_rays=cuts,old_pole_product_sign=1,old_pole_physical_product=Q(1)))
    assert len({tuple((v.a,v.b)for v in r['coordinates'])for r in out})==len(out)
    return out

def joint(o,rs):
    if not rs:return None
    g=H.mat(o['metric']);oldrows=H.mat(o['initial_cone']['normal_rows']);xs=[r['coordinates']for r in rs]
    rows=oldrows+[[sum((g[i][j]*x[j]for j in range(10)),Q())for i in range(10)]for x in xs]
    co=P.cone(rows,g);assert co['ray_counts']=={'-1':34,'0':0,'1':72}
    oldgram=H.mat(o['original_scaled_gram']);sc=F.vec(o['original_squared_scales'])+[Q(1),Q(1)]
    gg=[r+[products([row],x)[0]for x in xs]for r,row in zip(oldgram,oldrows)]
    for i,x in enumerate(xs):gg.append([gg[j][12+i]for j in range(12)]+[P.inner(x,g,y)for y in xs])
    assert gg[12][13]==Q(-1,-1)
    assert all(gg[i][j]==gg[j][i]for i,j in product(range(14),repeat=2))
    assert all(gg[i][i]==2*sc[i]for i in range(14))
    pair=[]
    for i,j in combinations(range(14),2):
        assert gg[i][j].sign()<=0;z=gg[i][j]*gg[i][j]/sc[i]/sc[j]
        assert P.square_gate(z)is None;pair.append(dict(walls=[i,j],physical_square=z))
    cy=cycles(gg,sc);assert len(cy)==93
    sig=H.inertia(gg);conj=H.inertia([[Q(x.a,-x.b)for x in r]for r in gg])
    assert sig==[9,1,4]and conj==[10,0,4]
    # Exact relative-interior timelike point on EACH wall, including the new ones.
    facets=[]
    for i in range(14):
        allr=[r for r in co['rays']if i in r['zero_facets']]
        finite=[r for r in allr if r['norm_sign']<0]
        x=[sum((r['coordinates'][k]for r in finite),Q())for k in range(10)]
        z=[sum((r['coordinates'][k]for r in allr),Q())for k in range(10)]
        assert F.rank([r['coordinates']for r in allr])==9
        n=0
        while True:
            eps=Q(Fraction(1,2**n));w=add(x,scale(z,eps));wn=P.inner(w,g,w)
            if wn.sign()<0:break
            n+=1;assert n<100
        wp=products(rows,w)
        assert wp[i]==0 and all(v.sign()<0 for j,v in enumerate(wp)if j!=i)
        facets.append(dict(wall=i,ray_indices=[r['index']for r in allr],finite_ray_indices=[r['index']for r in finite],
            span_rank=9,epsilon=eps,relative_interior_point=w,products=wp,norm=wn))
    # A strict spacelike point and strict timelike center certify an open ideal patch.
    r=next(r for r in co['rays']if r['norm_sign']>0);center=co['strict_timelike_center'];n=0
    while True:
        eps=Q(Fraction(1,2**n));y=add(r['coordinates'],scale(center,eps));yn=P.inner(y,g,y)
        if yn.sign()>0:break
        n+=1;assert n<100
    yp=products(rows,y);assert all(v.sign()<0 for v in yp)
    delta=add(y,scale(center,-1))
    poly=[P.inner(center,g,center),2*P.inner(center,g,delta),P.inner(delta,g,delta)]
    assert poly[0].sign()<0 and sum(poly,Q()).sign()>0
    escape=dict(ray_index=r['index'],ray_coordinates=r['coordinates'],ray_norm=r['norm'],ray_products=r['products'],
        epsilon=eps,strict_spacelike_point=y,strict_products=yp,strict_norm=yn,
        center_norm=poly[0],segment_null_polynomial=poly)
    return dict(cone=co,gram=gg,squared_scales=sc,physical_inertia=sig,conjugate_inertia=conj,
        native_pair_checks=pair,simple_cycles=cy,facet_centers=facets,escape=escape,
        infinite_hyperbolic_volume=True,all_fourteen_walls_are_facets=True)

def run():
    for p,h in PINS.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    old=H.read(P.TAG);H.bind(old);records=[];total=Counter()
    for o in old['cases']:
        faces=inventory(o);tests=anchors(o,faces);rs=roots(o,tests)
        print('ANCHORS',o['kinds'],tests['counts'],flush=True)
        jo=joint(o,rs)
        if jo:print('JOINT',o['kinds'],jo['cone']['ray_counts'],'infinite_volume',True,flush=True)
        records.append(dict(kinds=o['kinds'],relative_index=o['relative_index'],two_faces=faces,
            finite_anchor_domain=tests,new_roots=rs,joint=jo));total.update(tests['counts'])
    assert total==dict(tuples=24000,outside_relative_interior=7160,causal_projection=2210,
        norm_above2=11156,nonbase_two_cycle=2834,native_failure=636,survivor=4)
    paths=set(PINS)|{f'scripts/{TAG}.py','results/fresh_raw23_spacelike_pole_classification_20261018.md'}
    return P.enc(dict(status='PASS_RAW23_MIXED_TWO_FACE_NATIVE_ROOTS_AND_INFINITE_JOINT_CONES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},cases=records,
        two_face_count=675,mixed_faces=192,timelike_pair_faces=147,finite_anchor_tuples=24000,
        counts=dict(total),new_norm2_roots=4,joint_cones=2,joint_active_sets=4004,
        joint_spacelike_rays=144,joint_finite_vertices=68,joint_simple_cycles=186,
        candidate_ground_field='Q(sqrt(5))',arbitrary_passive_lift_allowed_in_search=True,
        surviving_passive_norm_zero=True,all_mixed_face_normals_classified=False,
        arbitrary_completion_excluded=False,compact_or_finite_volume_example_found=False,
        exported_columns=[],old_matrix_changed=False,optimizer_invocations=0,
        scope='Complete original two-face inventory; complete spacelike mixed-face candidate domain with two independent nonright finite Coxeter contacts over Q(sqrt(5)). Four exact pure roots survive. Each same-kind pair gives a genuine fourteen-facet infinite-volume hyperbolic intersection, with exact strict ideal escape. No claim for the untested mixed-face contact domains, arbitrary completion, capacity or model improvement.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert a.verify or not path.exists();d=run()
    if a.verify:assert d==json.loads(path.read_text())
    else:path.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
