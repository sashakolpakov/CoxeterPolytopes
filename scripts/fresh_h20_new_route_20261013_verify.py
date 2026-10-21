#!/usr/bin/env python3
"""Independent Fraction replay of cusp-adapted H20 section proposals.

Checks the affine cusp decomposition, full constant Coxeter link, every real
angle-equation root, exact interiors or cusp-cone obstructions, and actual
forbidden ridges. A modular filter followed by rational rank tests determines
whether each rational normal plane is contained in any four-wall span.
"""
from collections import Counter,defaultdict
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations
from math import gcd,lcm
import json
import fresh_h20_fourwall_verify_20261013 as C

A=C.A;V=C.V;U=C.U;B=C.B;ROOT=C.ROOT;Field=C.Field
SOURCE=ROOT/'data/enumerations/fresh_h20_new_route_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_new_route_20261013_verify.json'


def groups(M,indices):
    left=set(indices);out=[]
    while left:
        first=min(left);left.remove(first);todo=[first];seen={first}
        while todo:
            a=todo.pop();neighbors={b for b in left if M[a][b]}
            left-=neighbors;seen|=neighbors;todo.extend(neighbors)
        out.append(sorted(seen))
    return out


def kernel(M):
    rows,pivots=C.rref(M);n=len(M);free=[i for i in range(n) if i not in pivots]
    assert len(free)==1;v=[F(0)]*n;v[free[0]]=1
    for i,p in enumerate(pivots):v[p]=-rows[i][free[0]]
    if all(x<0 for x in v):v=[-x for x in v]
    assert all(x>0 for x in v)
    denominator=lcm(*(x.denominator for x in v));v=[int(x*denominator) for x in v];d=gcd(*v)
    return [x//d for x in v]


def modrref(rows,p):
    rows=[[int(x.numerator)*pow(int(x.denominator),-1,p)%p for x in map(F,row)] for row in rows]
    rank=0;pivots=[]
    for col in range(len(rows[0])):
        pivot=next((i for i in range(rank,len(rows)) if rows[i][col]),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank];d=pow(rows[rank][col],-1,p);rows[rank]=[x*d%p for x in rows[rank]]
        for i in range(len(rows)):
            if i!=rank:
                d=rows[i][col];rows[i]=[(a-d*b)%p for a,b in zip(rows[i],rows[rank])]
        pivots.append(col);rank+=1
        if rank==len(rows):break
    return rows,pivots


def plane_four_support(E,W):
    """Complete <=4-wall containment test; mod-p ranks only discard impossibles."""
    p=next(p for p in (1009,1013,1019,1021) if all(F(x).denominator%p for r in E for x in r)
           and len(modrref(E,p)[1])==2)
    er,pivots=modrref(E,p);raygroups=defaultdict(list);zero=[]
    for i,w in enumerate(W):
        row=[int(x)%p for x in w]
        for basis,pivot in zip(er,pivots):
            d=row[pivot];row=[(a-d*b)%p for a,b in zip(row,basis)]
        first=next((x for x in row if x),None)
        if first is None:zero.append(i)
        else:
            d=pow(first,-1,p);raygroups[tuple(x*d%p for x in row)].append(i)
    rays=list(raygroups);planes={}
    for i,j in combinations(range(len(rays)),2):
        a=list(rays[i]);b=list(rays[j]);pa=next(k for k,x in enumerate(a) if x);pb=next(k for k,x in enumerate(b) if x)
        if pa>pb:a,b=b,a;pa,pb=pb,pa
        d=b[pa];b=[(x-d*y)%p for x,y in zip(b,a)];pb=next(k for k,x in enumerate(b) if x)
        d=pow(b[pb],-1,p);b=[x*d%p for x in b];d=a[pb];a=[(x-d*y)%p for x,y in zip(a,b)]
        planes.setdefault((tuple(a),tuple(b)),set()).update((i,j))
    candidates=set()
    buckets=[zero]+[list(js)+zero for js in raygroups.values()]
    buckets += [zero+[i for r in rs for i in raygroups[rays[r]]] for rs in planes.values()]
    for bucket in buckets:
        for ids in combinations(sorted(bucket),4):candidates.add(ids)
    rejections=[]
    for ids in sorted(candidates):
        rows=[W[i] for i in ids];rank=len(C.rref(rows)[1])
        enlarged=len(C.rref(rows+E)[1])
        if enlarged==rank:
            return {'spanning_walls_one_based':[i+1 for i in ids],'prime':p,'modular_candidate_four_sets':len(candidates)}
        if len(candidates)<=20:rejections.append({'walls_one_based':[i+1 for i in ids],'wall_rank':rank,'rank_after_adding_normal_plane':enlarged})
    return {'spanning_walls_one_based':None,'prime':p,'modular_candidate_four_sets':len(candidates),
            'small_exact_rejection_list':rejections}


def section(row,Q,W,G,h):
    u=list(map(F,row['cusp_vector']));e0=list(map(F,row['base_normal']))
    assert len(u)==len(e0)==22 and A.dot(u,u,Q)==0 and A.dot(u,h,Q)<0
    up=[A.dot(w,u,Q) for w in W];assert max(up)<=0
    incident=[i for i,x in enumerate(up) if x==0]
    ambient=groups(G,incident);ids=[i-1 for i in row['ambient_component_walls_one_based']]
    assert ids in ambient
    M=[[G[i][j] for j in ids] for i in ids];marks=kernel(M)
    a,b=[ids.index(i-1) for i in row['bridge_walls_one_based']];assert M[a][b]<0
    cut=[r[:] for r in M];cut[a][b]=cut[b][a]=0;assert len(groups(cut,range(len(ids))))==2
    expected=[F(0)]*len(ids);expected[a]=marks[b];expected[b]=-marks[a]
    assert [A.dot(W[i],e0,Q) for i in ids]==expected
    assert len(C.rref([W[i] for i in ids]+[e0])[1])==len(C.rref([W[i] for i in ids])[1])
    q=A.dot(e0,e0,Q);assert q==F(row['normal_square']) and q>0
    assert A.dot(e0,h,Q)==A.dot(e0,u,Q)==0
    e=[A.trim((e0[d],u[d])) for d in range(22)]
    q=(q,);vals=[A.trim((A.dot(w,e0,Q),A.dot(w,u,Q))) for w in W]
    he=A.trim((0,A.dot(h,u,Q)))
    def gr(a,b):return A.sub(A.scale(q,G[a][b]),A.mul(vals[a],vals[b]))
    a,b=[i-1 for i in row['target_pair_one_based']];N=A.mul(gr(a,b),gr(a,b));den=A.mul(gr(a,a),gr(b,b))
    equations={2:N,3:A.sub(A.scale(N,4),den),4:A.sub(A.scale(N,2),den),
               5:A.add(A.sub(A.scale(A.mul(N,N),16),A.scale(A.mul(N,den),12)),A.mul(den,den)),
               6:A.sub(A.scale(N,4),A.scale(den,3))}
    assert row['target_labels']==[2,3,4,5,6]
    assert {int(m):A.readp(p) for m,p in row['target_equations'].items()}==equations
    return (e,q,vals,gr,N,den,equations),he,u,e0,incident


def cusp_link(row,Q,W,G,h,sd,u,e0,incident):
    e,q,vals,gr,*_=sd;q=q[0]
    H=[[F(G[i][j])-A.dot(W[i],e0,Q)*A.dot(W[j],e0,Q)/q for j in incident] for i in incident]
    assert min(H[i][i] for i in range(len(incident)))>=0
    live=[i for i in range(len(incident)) if H[i][i]>0];active=[incident[i] for i in live]
    data=row['central_cusp_link'];assert data['active_incident_walls_one_based']==[i+1 for i in active]
    assert data['inactive_incident_walls_one_based']==[i+1 for i in incident if i not in active]
    assert not data['inactive_incident_walls_one_based'] # True for this saved batch.
    H=[[H[i][j] for j in live] for i in live];assert len(C.rref(H)[1])==data['rank']==19
    for a,b in combinations(range(len(active)),2):
        assert H[a][b]<=0;ratio=H[a][b]**2/(H[a][a]*H[b][b])
        assert ratio in (0,F(1,4),F(1,2),F(3,4),1)
    components=groups(H,range(len(active)));assert len(components)==len(data['components'])
    multipliers=[]
    for js,proof in zip(components,data['components']):
        ids=[active[j] for j in js];M=[[H[i][j] for j in js] for i in js];marks=kernel(M)
        assert proof['walls_one_based']==[i+1 for i in ids]
        assert [[F(x) for x in r] for r in proof['gram']]==M
        assert proof['positive_kernel_marks']==marks and proof['rank']==len(ids)-1
        # Nonpositive off-diagonals and a positive kernel turn this Gram by
        # diagonal congruence into a connected graph Laplacian: it is PSD.
        null=[V.sum_polys(A.scale(A.sub(A.scale((q,),W[i][d]),A.mul(vals[i],e[d])),mark)
                          for mark,i in zip(marks,ids)) for d in range(22)]
        d=next(d for d,x in enumerate(u) if x);mu=A.scale(null[d],1/u[d])
        assert all(null[d]==A.scale(mu,u[d]) for d in range(22))
        assert mu[0]==q*F(proof['null_sum_cusp_multiple']) and mu[0]>0
        multipliers.append(mu)
    return multipliers


def check(c,row,Q,W,G,h,sd,he,u,multipliers):
    e,q,vals,gr,N,den,equations=sd;field=Field(A.readp(c['minimal_polynomial']),*c['isolating_interval'])
    status=c['status']
    if status=='CENTRAL_PARAMETER':assert field.sign(A.T)==0;return False
    assert field.sign(A.T)!=0
    a,b=[i-1 for i in row['target_pair_one_based']]
    if status=='TARGET_NOT_ACUTE_ELLIPTIC_PAIR':
        assert min(field.sign(gr(a,a)),field.sign(gr(b,b)))<=0 or field.sign(gr(a,b))>0;return False
    assert min(field.sign(gr(a,a)),field.sign(gr(b,b)))>0 and field.sign(gr(a,b))<=0
    labels=[m for m,p in equations.items() if field.rem(p)==A.Z and (m!=5 or field.sign(A.sub(A.scale(N,2),den))>0)]
    if status=='CONJUGATE_ANGLE_BRANCH':assert not labels;return False
    assert labels and labels==c['target_labels']
    if status=='NO_INTERIOR_BY_CUSP_CONE_IDENTITY':
        proof=c['cusp_cone_obstruction'];mu=field.rem(multipliers[proof['component_index']])
        assert mu==A.readp(proof['cusp_multiple']) and field.sign(mu)<=0;return False
    assert all(field.sign(mu)>0 for mu in multipliers)
    assert [field.rem(mu) for mu in multipliers]==[A.readp(p) for p in c['cusp_component_null_multiples']]
    Cpoint=[A.readp(p) for p in c['strict_interior']];dual,cp,square=U.field_products(Cpoint,Q,W,field)
    assert field.sign(square)<0 and all(field.sign(p)<0 for p in cp)
    assert field.sign(V.sum_polys(A.scale(dual[d],h[d]) for d in range(22)))<0
    assert field.rem(V.sum_polys(A.mul(dual[d],e[d]) for d in range(22)))==A.Z
    if c['interior_method']=='projected_original_center':
        assert Cpoint==[field.rem(A.sub(A.scale(q,h[d]),A.mul(he,e[d]))) for d in range(22)]
    else:assert c['interior_method']=='rational_LP_point_projected_exactly'
    if c['repaired_target_ridge'] is not None:
        assert c['repaired_target_ridge']['ridge_walls_one_based']==row['target_pair_one_based']
        B.actual(c['repaired_target_ridge'],Q,W,G,h,sd,field)
    assert status=='NON_COXETER_ACTUAL_RIDGE'
    proof=c['actual_forbidden_ridge'];aa,ab,bb=B.actual(proof,Q,W,G,h,sd,field)
    N=field.rem(A.mul(ab,ab));den=field.rem(A.mul(aa,bb))
    if proof['angle_rejection']=='obtuse':assert field.sign(ab)>0
    else:
        assert field.sign(ab)<=0 and field.sign(N)>0 and field.sign(A.sub(A.scale(den,3),A.scale(N,4)))>0
        for p in (A.sub(A.scale(N,4),den),A.sub(A.scale(N,2),den),
                  A.add(A.sub(A.scale(A.mul(N,N),16),A.scale(A.mul(N,den),12)),A.mul(den,den))):assert field.sign(p)!=0
    if row['cusp_source_index']==1 and row['bridge_walls_one_based']==[5,6] and c['repaired_target_ridge'] is not None:
        assert field.sign(A.T)>0 and c['target_labels']==[3]
        assert field.rem(A.sub(A.mul(N,(F(16),F(30))),A.scale(den,5)))==A.Z
        assert field.sign(A.sub(A.scale(N,4),den))>0 and field.sign(A.sub(A.scale(N,2),den))<0
    return True


def main():
    Q,W,G,h=A.source();data=json.loads(SOURCE.read_text());assert data['complete']
    assert len(data['pencils'])==data['pencil_count']==8
    total=accepted=repaired=minimum_five=0;statuses=Counter();degrees=Counter();supports=[];cases=[]
    for index,row in enumerate(data['pencils']):
        sd,he,u,e0,incident=section(row,Q,W,G,h);multipliers=cusp_link(row,Q,W,G,h,sd,u,e0,incident)
        total+=C.root_census(row,sd)
        support=plane_four_support([e0,u],W);supports.append({'cusp_index':row['cusp_source_index'],'bridge':row['bridge_walls_one_based'],**support})
        for c in row['candidates']:
            statuses[c['status']]+=1;yes=check(c,row,Q,W,G,h,sd,he,u,multipliers);accepted+=yes
            if yes:
                repaired+=c['repaired_target_ridge'] is not None;degrees[len(c['minimal_polynomial'])-1]+=1
                if len(c['minimal_polynomial'])>2 and support['spanning_walls_one_based'] is None:
                    assert C.no_rational_root(A.readp(c['minimal_polynomial']));minimum_five+=1
                cases.append((c,row,sd,he,u,multipliers))
        print('verified cusp pencil',index+1,'four-wall plane support',support,flush=True)
    assert dict(statuses)==data['statuses'] and cases
    c,row,sd,he,u,multipliers=cases[0];forged=deepcopy(c)
    forged['actual_forbidden_ridge']['point_coordinate_polynomials'][0][-1]=str(F(forged['actual_forbidden_ridge']['point_coordinate_polynomials'][0][-1])+1)
    try:check(forged,row,Q,W,G,h,sd,he,u,multipliers)
    except AssertionError:mutation=True
    else:raise AssertionError('forged actual ridge accepted')
    forged=deepcopy(row);forged['central_cusp_link']['components'][0]['positive_kernel_marks'][0]+=1
    try:cusp_link(forged,Q,W,G,h,sd,u,list(map(F,row['base_normal'])),[i for i,w in enumerate(W) if A.dot(w,u,Q)==0])
    except AssertionError:link_mutation=True
    else:raise AssertionError('forged cusp kernel accepted')
    out={'status':'PASS','standard_library_only':True,'pencils':len(data['pencils']),
         'complete_real_angle_equation_roots':total,'statuses':dict(statuses),
         'certified_finite_volume_sections_with_full_Coxeter_cusp':accepted,
         'actual_target_angle_repairs':repaired,'actual_forbidden_ridge_rejections':accepted,
         'accepted_parameter_polynomial_degrees':dict(degrees),'normal_plane_support_tests':supports,
         'irrational_normals_proved_to_need_at_least_five_original_walls':minimum_five,
         'forged_ridge_rejected':mutation,'forged_cusp_kernel_rejected':link_mutation,
         'scope':'Literal eight cusp pencils and labels2..6; no completeness claim over cusp types, bridge choices, all coefficients or labels. Modular support pruning is followed by exact rational rank checks.'}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)


if __name__=='__main__':main()
