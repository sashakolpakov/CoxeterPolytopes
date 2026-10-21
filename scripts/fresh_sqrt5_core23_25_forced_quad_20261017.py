#!/usr/bin/env python3
"""Exact neighbouring-face incidence forced by two literal prism profiles.

No optimization or mass repair.  This checks the actual common simplex star,
the forced baseline-two quadrilateral, and a new incidence inequality.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_sqrt5_multispectator_cap_capacity_20261017 as M

V,K=M.V,M.K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_core23_25_forced_quad_20261017'
MODEL='fresh_sqrt5_a3_d5_cap_joint_ports_20261017'
FORESTS='fresh_sqrt5_rank13_combined_disjoint_ports_20261017'
FULL='fresh_sqrt5_cross_coface_cap_full_model_20261017'
PHI=K(F(1,2),F(1,2))


def sub(G,inds):return [[G[i][j]for j in inds]for i in inds]
def out(G):return [[x.out()for x in row]for row in G]
def dot(v,G,w):return sum((v[i]*G[i][j]*w[j]for i in range(len(v))for j in range(len(w))),K())
def conjugate(G):return [[K(x.a,-x.b)for x in row]for row in G]


def diagonal_pivots(G):
    """Exact fixed-order congruence; a zero pivot is allowed only on a zero row."""
    H=[[K(x)for x in row]for row in G];piv=[]
    for k in range(len(H)):
        z=H[k][k];piv.append(z)
        if not z:
            assert all(not H[k][i]for i in range(k,len(H)))
            continue
        for i in range(k+1,len(H)):
            for j in range(k+1,len(H)):H[i][j]-=H[i][k]*H[k][j]/z
    return piv


def inertia(G):
    p=diagonal_pivots(G)
    return [sum(x.sign()==s for x in p)for s in (1,-1,0)]


def identify_corner(G,inds):
    """Recognize exactly the elementary A/H components used in this witness."""
    unseen=set(inds);components=[]
    while unseen:
        todo=[min(unseen)];part=set()
        while todo:
            i=todo.pop()
            if i in part:continue
            part.add(i);todo.extend(j for j in unseen if j!=i and G[i][j])
        unseen-=part;part=sorted(part)
        assert all(G[i][i]==2 for i in part)
        edges=[(i,j,3 if G[i][j]==-1 else 5 if G[i][j]==-PHI else None)
               for i,j in combinations(part,2)if G[i][j]]
        assert all(m is not None for i,j,m in edges)
        if len(part)==1:name='A1'
        elif len(part)==2:
            assert len(edges)==1 and edges[0][2]==3;name='A2'
        elif len(part)==3:
            assert sorted(e[2]for e in edges)==[3,5];name='H3'
        else:
            assert len(part)==4 and len(edges)==3 and all(e[2]==3 for e in edges)
            degrees=Counter(i for a,b,m in edges for i in (a,b))
            assert sorted(degrees.values())==[1,1,2,2];name='A4'
        components.append(dict(name=name,roots=part,edges=edges))
    return sorted(c['name']for c in components),components


def run():
    D={key:M.read(tag)for key,tag in M.TAGS.items()if key!='model'}
    S,I,A,Q,R=[D[k]for k in ('slices','incidence','actual','quad','raw')]
    model=M.read(MODEL);alpha={i:F(a)for i,a in model['face_multipliers']}
    mu={i:F(a)for i,a in model['vertex_multipliers']}
    full=M.read(FULL);pmap=dict(full['prepared']['provenance_to_face'])
    typ=I['base_types'][55301];assert typ['index']==55301 and typ['core']==1
    assert typ['spectator_forest']==['A1','A1','A4','H3'] and typ['old_face']==580
    profiles={p['base_type']:p for p in S['coface_slice_demands']}
    assert profiles[55301]['simplex_face']==580
    core=I['cores'][1];cg=V.core_graph(core)
    Gc,scales=V.physical_core(cg);assert scales is None
    J,det=V.inverse(Gc)
    assert len(cg)==5 and core['spherical_deletions']==list(range(5))
    assert V.automorphisms(cg)==[list(range(5))]
    assert inertia(Gc)==[4,1,0] and V.pd(conjugate(Gc))
    assert all(V.pd(sub(Gc,[j for j in range(5)if j!=i]))for i in range(5))
    templates={t['canonical_core']:t for t in I['prism_templates']}
    records=[]
    specs=[(1513,483,26410,23,3,list(range(5)),K(F(7,2),F(1,2))),
           (1617,535,26880,25,2,[0,1,2,4],K(6,2))]
    for bi,pr,fi,raw_index,contact,vertices,tau in specs:
        B=S['actual_bases'][bi];assert B['index']==bi and B['base_type']==55301
        assert B['provenance']==pr and B['face']==fi==pmap[pr]
        assert B['coface_location']==['finite','A1',0] and B['canonical_vertices']==vertices
        ci=A['provenance_masses'][pr]['canonical_prism_core'];C=Q['canonical_cores'][ci]
        assert C['source_core']==raw_index
        raw=R['records'][raw_index];g=V.raw_graph(raw)
        T=next(t for t in templates[ci]['bases']if t['deleted_normal']==B['deleted'])
        assert T['base_core']==1 and T['finite_extra']==['A1']
        b,q=B['deleted'],B['coface_root'];mp=dict(T['core_isomorphism'])
        assert T['opposite_base_facet']==q and set(range(raw['nodes']))==set(mp)|{b,q}
        assert all(g[x][y]==cg[mp[x]][mp[y]]for x,y in combinations(mp,2))
        assert all(g[q][x]==V.RIGHT for x in mp)
        assert all(g[b][x]==(V.edge(3)if mp[x]==contact else V.RIGHT)for x in mp)
        assert V.squared(g[b][q])==tau and (tau-4).sign()>0
        pair_vertices=set()
        for left,right in C['quad_partitions']:
            if b in left:side,other=left,right
            elif b in right:side,other=right,left
            else:continue
            assert set(side)=={b,q};pair_vertices.update(other)
        assert sorted(mp[x]for x in pair_vertices)==vertices
        features=[f for f in model['actual_quad_features']if f['face']==fi]
        assert len(features)==1 and features[0]['provenance']==pr
        assert features[0]['spectator_forest']==['A1','A4','H3']
        assert features[0]['prism_dimension']==len(vertices)==B['prism_dimension']
        # This selected Q has exactly one base belonging to S580, not two.
        assert len([z for z in S['actual_bases']if z['provenance']==pr and z['base_type']==55301])==1
        records.append(dict(base_index=bi,provenance=pr,face=fi,raw_source=raw_index,
            canonical_prism_core=ci,deleted=b,coface=q,core_isomorphism=T['core_isomorphism'],
            contact_root=contact,canonical_base_vertices=vertices,tau=tau.out(),mass=str(alpha[fi])))
    ta,tb=specs[0][-1],specs[1][-1]
    assert J[3][3]==2-ta/2==K(F(1,4),F(-1,4))
    assert J[2][2]==2-tb/2==K(-1,-1) and J[3][2]==-PHI
    assert J[3][3].sign()<0 and (J[3][3]+1).sign()>0
    assert M.signrad(-PHI,K(F(1,2)),ta*tb)>0
    assert (J[2][2]+2).sign()<0
    # S-star: core0..4, isolated qA=5,qB=6, A4=7..10,H3=11..13.
    G=[[K(2)if i==j else K()for j in range(16)]for i in range(16)]
    for i in range(5):
        for j in range(5):G[i][j]=Gc[i][j]
    for i,j in [(7,8),(8,9),(9,10),(11,12)]:G[i][j]=G[j][i]=K(-1)
    G[12][13]=G[13][12]=-PHI
    G[5][5]=2*ta;G[6][6]=2*tb
    G[14][3]=G[3][14]=K(-1);G[15][2]=G[2][15]=K(-1)
    G[14][5]=G[5][14]=-ta;G[15][6]=G[6][15]=-tb
    G[14][15]=G[15][14]=-PHI
    basis=list(range(14));IG,detG=V.inverse(sub(G,basis))
    assert inertia(sub(G,basis))==[13,1,0] and V.pd(conjugate(sub(G,basis)))
    contacts=[[G[i][r]for i in basis]for r in (14,15)]
    projection=[[dot(x,IG,y)for y in contacts]for x in contacts]
    assert projection==[[K(2),-PHI],[-PHI,K(2)]]
    assert inertia(G)==[13,1,2] and inertia(conjugate(G))==[14,0,2]
    assert all(G[i][j].sign()<=0 for i,j in combinations(range(16),2))
    residuals=[]
    for i in range(5):
        W=[j for j in range(5)if j!=i]+list(range(7,14))
        IW,_=V.inverse(sub(G,W));cv=[[G[j][r]for j in W]for r in (14,15)]
        H=[[G[r][s]-dot(cv[u],IW,cv[v])for v,s in enumerate((14,15))]for u,r in enumerate((14,15))]
        determinant=H[0][0]*H[1][1]-H[0][1]*H[1][0]
        pd=H[0][0].sign()>0 and determinant.sign()>0
        assert pd==(i==2)
        residuals.append(dict(deleted_simplex_vertex=i,schur=out(H),determinant=determinant.out(),positive_definite=pd))
    assert residuals[2]['determinant']==K(F(3,2),F(-1,2)).out()
    W=[0,1,3,4]+list(range(7,14));black=[5,6,14,15]
    # Cyclic order qA,qB,A,B.  All corners are independently checked PD.
    corners=[]
    forest_source=M.read(FORESTS)['vertices'];flookup={tuple(v['forest']):v['index']for v in forest_source}
    for x,y in [(5,6),(6,14),(14,15),(15,5)]:
        selected=W+[x,y];assert len(selected)==13 and V.pd(sub(G,selected))
        # Convert the two scaled isolated q roots back to norm2 for recognition.
        normal=[[K(z)for z in row]for row in G]
        normal[5][5]=normal[6][6]=K(2)
        forest,components=identify_corner(normal,selected)
        corners.append(dict(black_pair=[x,y],forest=forest,vertex=flookup[tuple(forest)],components=components))
    assert [x['forest']for x in corners]==[
        ['A1','A1','A1','A4','H3','H3'],['A1','A2','A4','H3','H3'],
        ['A4','H3','H3','H3'],['A1','A1','A1','A4','H3','H3']]
    assert corners[2]['vertex']==5572 and forest_source[5572]['forest']==['A4','H3','H3','H3']
    pointed=next(c for c in corners[2]['components']if 14 in c['roots'])
    assert pointed['roots']==[3,14,15] and sorted(pointed['edges'])==[(3,14,3),(14,15,5)]
    # The new minimal core contains two dotted pairs and its native 5-bond.
    local=[5,6,14,15,3];L=sub(G,local);kernel=[2,-1,4,-2,2]
    assert all(sum((L[i][j]*kernel[j]for j in range(5)),K())==0 for i in range(5))
    order=[0,1,4,2,3];piv=diagonal_pivots(sub(L,order))
    assert piv==[2*ta,2*tb,K(2),-PHI/2,K()]
    assert inertia(sub(L,order))==[3,1,1] and inertia(conjugate(sub(L,order)))==[4,0,1]
    assert all(M.integral(2*L[i][j]/L[i][i])for i in range(5)for j in range(5))
    labels=[2,0,2,2,2,0,2,5,3,2]
    normalized_square={(i,j):L[i][j]*L[i][j]*4/(L[i][i]*L[j][j])for i,j in combinations(range(5),2)}
    for (i,j),m in zip(combinations(range(5),2),labels):
        if m:assert normalized_square[i,j]==V.squared(V.edge(m))
        else:assert normalized_square[i,j]==(ta if (i,j)==(0,2)else tb)
    assert [m for (i,j),m in zip(combinations(range(5),2),labels)if i<4 and j<4 and m==5]==[5]
    assert sum(m==5 for c in corners[2]['components']for i,j,m in c['edges'])==3
    # Baseline: right disconnected pairs cost0; framing5 H3 pair costs2.
    assert all(not G[x][y]for x,y in [(5,6),(6,14),(15,5)]) and G[14][15]==-PHI
    a,b,s=alpha[26410],alpha[26880],alpha[580];vy=mu.get(5572,F())
    lower=2*a+b-2*s;gap=lower-3*vy
    assert lower==F(3363761,10**9) and vy==0 and gap>0
    # Exhaust every possible incidence count under the proven cap restrictions.
    table=[]
    for na in (0,1):
        for nb in range(3-na):
            mixed=na*nb;assert mixed>=2*na+nb-2
            table.append([na,nb,mixed,2*na+nb-2])
    paths=[Path(__file__),Path(M.__file__),Path(V.__file__)]
    paths += [ROOT/f'results/{tag}.json'for key,tag in M.TAGS.items()if key!='model']
    paths += [ROOT/f'results/{tag}.json'for tag in (MODEL,FORESTS,FULL)]
    paths += [ROOT/'results'/p for p in ['fresh_simplex_disjoint_ports_20261017.md','fresh_singular_quad_local_realization_20261017.md']]
    return dict(status='PASS_EXACT_CORE23_25_FORCED_QUADRILATERAL_INCIDENCE',dimension=13,
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        common_simplex=dict(face=580,base_type=55301,core=1,core_gram=out(Gc),core_inverse=out(J),
            automorphisms=[list(range(5))],spectator_forest=typ['spectator_forest'],rank=14),
        selected_profiles=records,projection_products=out([[J[3][3],J[3][2]],[J[2][3],J[2][2]]]),
        same_coface_mixed_product=dict(rational_part=(-PHI).out(),radical_coefficient='1/2',radicand=(ta*tb).out(),sign=1),
        combined_gram=out(G),combined_root_roles=dict(core=list(range(5)),qA=5,qB=6,A4=list(range(7,11)),H3=list(range(11,14)),A=14,B=15),
        full_star_projection_products=out(projection),all_common_vertex_residuals=residuals,
        forced_quad=dict(white=W,black_cyclic=black,minimal_root_order=local,labels=labels,
            parameters={'0,2':ta.out(),'1,3':tb.out()},gram=out(L),kernel=kernel,
            pivot_order=order,pivots=[x.out()for x in piv],physical_inertia=[3,1,1],conjugate_inertia=[4,0,1],
            spectator_forest=['A4','H3','H3'],corners=corners,baseline_corner_weights=[0,0,2,0],baseline_total=2,
            native_framing_label=5,induced_angle_cos_squared=((PHI*PHI)/3).out()),
        cap_count_table=table,inequality=dict(positive_faces=[[26410,2],[26880,1]],negative_faces=[[580,2]],
            vertex_capacity=[[5572,3]],face_lower_bound=str(lower),vertex_upper_bound=str(3*vy),gap=str(gap)),
        input_bound=model['bound'],alpha_changed=False,mu_changed=False,optimizer_invocations=0,
        scope='One exact adjacent-face compatibility inequality for selected actual raw23/25 prism profiles at S580. The forced baseline2 quadrilateral is outside the old B1-only domain. Conditional on actual simultaneous occurrences; no ambient realization, complete forced-quad census, repaired witness, optimization or dimension theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    data=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['inequality'],flush=True)
