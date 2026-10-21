#!/usr/bin/env python3
"""A pointed forced quadrilateral couples two finite-vertex budgets.

Reconstructs all84 active six-root supports and all12 actual corner-pair
orientations. No optimizer, pair-outcome changes or frozen-file mutations.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import networkx as nx
import fresh_restricted_290_cyclic_cover_20261017 as Z

A=Z.A;R=Z.R;N=A.N;W=Z.W;V=Z.V;K=Z.K;M=Z.M;ROOT=Z.ROOT
TAG='fresh_forest1325_cross_vertex_incidence_20261017'
INC='fresh_sqrt5_core23_25_spectator_incidence_20261017'
ACTIVE=[0,1,2,6,7,8,13,14,15]
LABELS=['qB','qA','w','h5','hcenter','h3','B','A','C']
WEIGHTS=[3,1,2,4,4,4,0,0]


def normal(x):return json.loads(json.dumps(x,sort_keys=True))
def out(g):return [[x.out()for x in row]for row in g]
def sub(g,ids):return [[g[i][j]for j in ids]for i in ids]


def congruence(g):
    """Exact symmetric elimination, including zero-diagonal 2x2 pivots."""
    a=[row[:]for row in g];labels=list(range(len(a)));pivots=[];sig=[0,0,0]
    while a:
        n=len(a);i=next((i for i in range(n)if a[i][i]),None)
        if i is not None:
            q=a[i][i];sig[0 if q.sign()>0 else 1]+=1
            pivots.append(dict(indices=[labels[i]],value=q.out()))
            keep=[j for j in range(n)if j!=i]
            a=[[a[j][k]-a[j][i]*a[i][k]/q for k in keep]for j in keep]
        else:
            pair=next(((i,j)for i,j in combinations(range(n),2)if a[i][j]),None)
            if pair is None:
                sig[2]+=n;pivots.append(dict(zero_indices=labels));break
            i,j=pair;q=a[i][j];sig[0]+=1;sig[1]+=1
            pivots.append(dict(indices=[labels[i],labels[j]],off_diagonal=q.out()))
            keep=[k for k in range(n)if k not in pair]
            a=[[a[k][l]-(a[k][i]*a[j][l]+a[k][j]*a[i][l])/q for l in keep]for k in keep]
        labels=[labels[i]for i in keep]
    return dict(inertia=sig,pivots=pivots)


def kernel(g):
    a=[row[:]for row in g];n=len(a);piv=[];k=0
    for j in range(n):
        i=next((i for i in range(k,n)if a[i][j]),None)
        if i is None:continue
        a[k],a[i]=a[i],a[k];q=a[k][j];a[k]=[x/q for x in a[k]]
        for i in range(n):
            if i!=k and a[i][j]:
                q=a[i][j];a[i]=[x-q*y for x,y in zip(a[i],a[k])]
        piv.append(j);k+=1
    free=[i for i in range(n)if i not in piv];assert len(free)==1
    v=[K()for _ in range(n)];v[free[0]]=K(1)
    for i,j in enumerate(piv):v[j]=-a[i][free[0]]
    assert all(not sum((x*y for x,y in zip(row,v)),K())for row in g)
    return v


def schur(g,white,black):
    j,det=V.inverse(sub(g,white));assert det.sign()>0
    return [[g[a][b]-sum((g[a][u]*j[i][k]*g[v][b]for i,u in enumerate(white)for k,v in enumerate(white)),K())for b in black]for a in black]


def graph(g):
    q=nx.Graph();q.add_nodes_from(g)
    q.add_edges_from((i,j,dict(label=m))for i,ns in g.items()for j,m in ns.items()if i<j)
    return q


def transported(star,mp):
    contacts=sorted(tuple(sorted((mp[i],tuple(e))for i,e in p))for p in star['contacts'])
    ports=sorted(tuple(sorted((mp[i],mp[j])))for i,j in star['occupied_pairs'])
    return contacts,ports


def run(verify=False):
    start=time.monotonic();p,changed,rows,patches,base=A.prepare()
    old,g,owner,unowned,rdemand,groups,forest,applied=A.current(1325,p,changed,rows,patches)
    nodes=A.T.current(1325,p,changed)[5]
    assert not applied and forest==['A1','A1','A1','H3','H3','H4']
    assert [x['face']for x in groups]==[45769,26918,26918,26362,26362,26448,45573,26683]
    assert all(x['positive_face_coefficients']==[[x['face'],1]]for x in groups)
    assert nodes[12]['origin']['provenance']==541 and nodes[46]['origin']['provenance']==490
    assert nodes[12]['origin']['raw_family']==25 and nodes[46]['origin']['raw_family']==23
    assert nodes[12]['origin']['omitted']==[2,3]and nodes[46]['origin']['omitted']==[2,5]
    cyc=W.read(Z.TAG);cr=next(r for r in cyc['forest_records']if r['vertex']==1325)
    term=next(t for t in cr['cover_terms']if t['literal_stars']==[12,46]);assert term['walls']==[1,5,8]
    # Repeat the target frame and its physical contact/height/cross arithmetic.
    diag=next(r for r in W.read(A.TAG)['forest_records']if r['vertex']==1325)
    assert Z.forest_check(1325,diag,p,changed,rows,patches)==cr
    target=R.P.forest_graph(forest,dict(R.C.models()))
    covectors=[N.native(s)for s in nodes]
    walls=cr['native_wall_contacts'];ids={tuple((i,tuple(e))for i,e in wall):j for j,wall in enumerate(walls)}
    assert [sorted(ids[q]for q in covectors[i])for i in(12,46)]==[[1,8],[5,8]]
    sqs={(i,j):K(*map(F,s))for i,j,s in term['nonzero_squared_edges']}
    tau=K(F(7,2),F(1,2));scales=[K(1),tau]+[K(1)]*7
    G=[[K(2)*scales[i]if i==j else K()for j in range(9)]for i in range(9)]
    for i,j in combinations(range(9),2):
        s=sqs.get(tuple(sorted((ACTIVE[i],ACTIVE[j]))),K())
        q=M.ksqrt(s*scales[i]*scales[j])if s else K();assert q is not None
        G[i][j]=G[j][i]=-q
    full_inertia=congruence(G);assert full_inertia['inertia']==[6,1,2]
    models=dict(R.C.models());lookup={R.C.code(q):name for name,q in models.items()}
    vlookup={tuple(v['forest']):v['index']for v in W.read(R.P.FORESTS)['vertices']}
    # The two true prism faces fix respectively qA,w and qB, in addition to
    # the seven passive roots. Their dimensions are4 and5, not wall counts.
    q25fixed=[1,2];q23fixed=[0]
    assert set(range(13))-set(sum(nodes[12]['occupied_pairs'],[]))=={1,2,3,4,5,9,10,11,12}
    assert set(range(13))-set(sum(nodes[46]['occupied_pairs'],[]))=={0,3,4,5,9,10,11,12}
    supports=[]
    for inds in combinations(range(9),6):
        q=sub(G,inds);cert=congruence(q);pd=cert['inertia']==[6,0,0]
        assert pd==V.pd(q)
        rec=dict(active_support=inds,congruence=cert,positive_definite=pd)
        if pd:
            native={i:{}for i in inds}
            for i,j in combinations(inds,2):
                s=sqs.get(tuple(sorted((ACTIVE[i],ACTIVE[j]))),K())
                if s:
                    m=next((m for m in(3,4,5,6,10)if V.squared(V.edge(m))==s),None);assert m
                    native[i][j]=native[j][i]=m
            fs=sorted(lookup[R.C.code(q)]for q in R.C.pieces(native));whole=sorted(fs+['H3','H4'])
            J,det=V.inverse(q);assert det.sign()>0 and all(x.sign()>=0 for row in J for x in row)
            rec.update(active_forest=fs,full_forest=whole,vertex=vlookup[tuple(whole)],
                inverse=out(J),in_raw25_prism=7 not in inds and all(i in inds for i in q25fixed),
                in_raw23_prism=6 not in inds and all(i in inds for i in q23fixed))
        supports.append(rec)
    pd=[r for r in supports if r['positive_definite']];assert len(supports)==84 and len(pd)==15
    assert sum(r['in_raw25_prism']for r in pd)==8 and sum(r['in_raw23_prism']for r in pd)==10
    outside=[r for r in pd if not(r['in_raw25_prism']or r['in_raw23_prism'])]
    assert len(outside)==1 and outside[0]['vertex']==6537 and outside[0]['active_support']==(2,3,4,5,6,7)
    pds={tuple(r['active_support']):r for r in pd}
    # Genuine source triangle on the fixed eleven walls (0,1,3,4)+passives.
    triangle_white=[0,1,3,4];triangle_sides=[2,5,8]
    triangle_corners=[tuple(sorted(triangle_white+[i,j]))for i,j in combinations(triangle_sides,2)]
    assert all(c in pds for c in triangle_corners)and triangle_corners[0]==tuple(range(6))
    TS=schur(G,triangle_white,triangle_sides);TJ,_=V.inverse(TS)
    assert congruence(TS)['inertia']==[2,1,0]and all(x.sign()<0 for row in TJ for x in row)
    # Forced quadrilateral: w and both H3/H4 spectator components stay fixed.
    quad_white=[2,3,4,5];quad_pairs=[[0,1],[0,7],[6,7],[1,6]]
    quad_corners=[tuple(sorted(quad_white+q))for q in quad_pairs]
    assert all(q in pds for q in quad_corners)
    minimal=[0,1,2,6,7];Q=sub(G,minimal);nu=kernel(Q)
    assert congruence(Q)['inertia']==[3,1,1]
    signs=[nu[minimal.index(i)].sign()for i in(0,1,6,7)]
    assert all(signs[i]*signs[(i+1)%4]<0 for i in range(4))
    assert G[2][7]==-1 and G[6][7]==-K(F(1,2),F(1,2))and not G[2][6]
    # Every compatible selected corner-pair is an automorphic copy, not a
    # free orientation multiplier in the eventual whole-face coefficient.
    pairs=sorted((i,j)for i in groups[1]['literal_stars']for j in groups[5]['literal_stars']if g.has_edge(i,j))
    assert len(pairs)==12
    matcher=nx.algorithms.isomorphism.GraphMatcher(graph(target),graph(target),edge_match=lambda a,b:a['label']==b['label'])
    automorphisms=list(matcher.isomorphisms_iter());assert len(automorphisms)==12
    mapped=[]
    for mp in automorphisms:
        pair=[]
        for s in(12,46):
            key=transported(nodes[s],mp)
            candidates=[i for i in range(len(nodes))if transported(nodes[i],{j:j for j in range(13)})==key]
            assert len(candidates)==1;pair.append(candidates[0])
        mapped.append(dict(pair=pair,base_permutation=[mp[i]for i in range(13)],
            triangle_source_ports=[mp[2],mp[8]],native5_recovery_port=mp[2]))
    assert sorted(tuple(r['pair'])for r in mapped)==pairs
    # Re-enumerate the complete current compatibility graph. The marked
    # configuration count is an actual-corner event, not a mixture amount.
    clique_checks=[];patterns=Counter()
    for cl in nx.find_cliques(g):
        cl=sorted(cl);vec=[Counter(owner[i]for i in cl).get(j,0)for j in range(8)]
        x=sum(i in cl and j in cl for i,j in pairs);charge=sum(a*b for a,b in zip(WEIGHTS,vec))
        assert x<=1 and charge<=4+x
        clique_checks.append(dict(stars=cl,count_vector=vec,mixed_pairs=x,charge=charge,slack=4+x-charge))
        patterns[tuple(vec),x]+=1
    assert len(clique_checks)==522
    # Downward closure: deleting a vertex loses at most one marked pair;
    # every endpoint of such a pair has weight at least1.
    assert WEIGHTS[1]>=1 and WEIGHTS[5]>=1
    assert all(sum(i==a or j==a for i,j in pairs)==1 for a in groups[1]['literal_stars']+groups[5]['literal_stars'])
    coef=Counter()
    for w,q in zip(WEIGHTS,groups):
        for f,c in q['positive_face_coefficients']:coef[f]+=w*c
    terms=[[f,c]for f,c in sorted(coef.items())if c]
    assert terms==[[26362,8],[26448,4],[26918,3],[45769,3]]
    incidence=W.read(INC);oldrow=next(r for r in incidence['constraints']if r['index']==50)
    assert oldrow['positive_faces']==[[26448,2],[26918,1]]and oldrow['negative_faces']==[[311,2]]
    assert oldrow['vertex_capacity']==[[6537,3]]
    alphas=[dict((i,F(x))for i,x in base[k])for k in('face_multipliers','restricted_face_multipliers')]
    repairs=[];loads=[]
    for alpha,old in zip(alphas,base['repairs']):
        mu={i:F(x)for i,x in old['vertex_multipliers']};load=sum((c*alpha.get(i,F())for i,c in terms),F())
        oldrhs=4*mu[1325]+3*mu[6537];gap=load-oldrhs;assert gap>0
        new=dict(mu);new[1325]=max(mu[1325],(load-3*mu[6537])/4)
        assert new[1325]>mu[1325]and 4*new[1325]+3*new[6537]==load
        oldgap=2*alpha[26448]+alpha[26918]-2*alpha[311]-3*mu[6537];assert oldgap<0
        num=F(old['objective_numerator']);den=sum(new.values(),F());assert sum(mu.values(),F())==F(old['normalizer'])
        repairs.append(dict(domain=old['domain'],baseline_bound=old['bound'],objective_numerator=str(num),
            baseline_load=str(load),baseline_capacity=str(oldrhs),gap=str(gap),old_incidence_row50_gap=str(oldgap),
            normalizer=str(den),bound=str(num/den),vertex_updates=[[1325,str(mu[1325]),str(new[1325])]],
            vertex_multipliers=[[v,str(x)]for v,x in sorted(new.items())if x]))
        loads.append(str(load))
    assert F(repairs[0]['bound'])<12<F(repairs[1]['bound'])
    sources=A.bindings()
    tags=[Z.TAG,A.TAG,INC,'fresh_sqrt5_core23_25_forced_quad_20261017','fresh_quad_product_aggregation_20261017','fresh_simplex_disjoint_ports_20261017']
    for tag in tags:
        doc=W.read(tag)
        for path,h in doc.get('sources',{}).items():assert W.digest(path)==h,(tag,path)
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for mod in(R,M,V):
        path=str(Path(mod.__file__).relative_to(ROOT));sources[path]=W.digest(path)
    sources[f'scripts/{TAG}.py']=W.digest(f'scripts/{TAG}.py')
    data=normal(dict(status='PASS_FOREST1325_CROSS_VERTEX_FORCED_INCIDENCE',sources=sources,
        baseline_source=A.S.TAG,vertex=1325,forest=forest,target_vertex=6537,target_forest=['H3','H3','H3','H4'],
        corner_orbits=groups,literal_stars=nodes,selected_stars=[12,46],selected_native_wall_indices=[1,5,8],
        active_full_indices=ACTIVE,active_labels=LABELS,passive_full_indices=[3,4,5,9,10,11,12],
        active_gram=out(G),squared_root_scales=[x.out()for x in scales],active_congruence=full_inertia,
        support_records=supports,forced_new_vertices=outside,
        triangle=dict(white_active=triangle_white,side_active=triangle_sides,corners=triangle_corners,
            source_corner=list(range(6)),source_ports=[2,5],opposite_external_wall=8,
            residual_gram=out(TS),residual_inverse=out(TJ)),
        quadrilateral=dict(white_active=quad_white,side_pairs=quad_pairs,corners=quad_corners,
            corner_forests=[pds[q]['full_forest']for q in quad_corners],minimal_active=minimal,
            minimal_gram=out(Q),kernel=[x.out()for x in nu],native5_roles=dict(center=7,leaf5=6,leaf3=2)),
        compatible_corner_pairs=pairs,pair_automorphisms=sorted(mapped,key=lambda x:x['pair']),
        graph_edges=sorted(tuple(sorted(e))for e in g.edges()),local_orbit_weights=WEIGHTS,
        maximal_clique_checks=sorted(clique_checks,key=lambda x:x['stars']),
        count_patterns=[dict(vector=list(v),mixed_pairs=x,multiplicity=n)for(v,x),n in sorted(patterns.items())],
        pointed_incidence_capacity=3,old_incidence_row=oldrow,
        columns=[dict(index=0,positive_face_coefficients=terms,vertex_capacities=[[1325,4],[6537,3]],loads=loads)],
        repairs=repairs,face_multipliers=base['face_multipliers'],restricted_face_multipliers=base['restricted_face_multipliers'],
        counts=dict(active_supports=84,positive_definite_supports=15,raw25_prism_corners=8,raw23_prism_corners=10,
            vertices_outside_both_prisms=1,compatible_selected_corner_pairs=12,complete_vertex_automorphisms=12,
            maximal_literal_cliques=len(clique_checks),count_patterns=len(patterns),new_columns=1),
        optimizer_invocations=0,new_pair_outcomes=0,
        preservation='All face masses and all historical pair owners/unknowns stay fixed. Onlymu1325 increases in each witness, preserving every439410 inherited nonnegative-variable model inequality including52 signed-face positive-cost columns, all earlier local weighted capacities, and both currently owned1763 unknowns. No old-alpha implicit coverage is claimed.',
        scope='One new cross-vertex necessary incidence column. Complete84-support and12-compatible-orientation domain at forest1325; actual vertices/triangle/quad are justified by strict-center projection and convexity. The pointed triangle/maximal-simplex recovery proof supplies the3-per-target count. No global optimization, complete physical realization, target feasibility, full face/field coverage or dimension theorem.'))
    path=ROOT/f'results/{TAG}.json'
    if verify:assert data==W.read(TAG)
    else:path.write_text(json.dumps(data,separators=(',',':'))+'\n')
    print(data['status'],data['counts'],'bounds',[float(F(x['bound']))for x in repairs],
        'seconds',round(time.monotonic()-start,2),flush=True)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');args=a.parse_args();run(args.verify)
