#!/usr/bin/env python3
"""Exact isometric six-cap insertion at the remote mirrored gap."""
import argparse,json,sys
from pathlib import Path
from itertools import combinations
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_face_bridge_20261013_verify as B
V=B.V;K=B.K;C=B.C;P=B.P
PREFIX='data/enumerations/fresh_h9_prism_completion_20261013'
PERM=[0,2,1,3,4,6,5,7,8,9]
SOURCE_NAMES=C.NEW+B.NEW
NEW=['mirror'+x[0].upper()+x[1:]for x in SOURCE_NAMES]

def table(data,G):
    names,pairs=B.table(G)
    for new,row in zip(B.NEW,data['caps'],strict=True):
        finite=dict(zip(row['source_walls'],row['labels'],strict=True))
        for n in names:pairs['|'.join(sorted((new,n))) ]={'kind':'finite','label':finite[n]}if n in finite else{'kind':'dotted'}
    pairs['|'.join(sorted(B.NEW))]={'kind':'dotted'}
    return names+B.NEW,pairs

def reconstruct(d):
    assert d['source']==B.PREFIX+'_caps.json'and d['source_sha256']==V.sha(d['source'])
    assert d['permutation']==PERM and d['source_normals']==SOURCE_NAMES
    olddata=V.read(d['source']);prior=B.check(olddata)
    G,I,h,g,old,names,caps,center,_,_=B.reconstruct(olddata);old+=caps;names+=B.NEW
    assert len(old)==96
    assert all(G[PERM[i]][PERM[j]]==G[i][j]for i in range(10)for j in range(10))
    expect,pairs=table(olddata,G);assert expect==names
    # Both transpositions are reflections in differences of norm-two seed
    # normals. G(e_i-e_j)=2(e_i-e_j) in coordinates for these pairs.
    for i,j in((1,2),(5,6)):
        assert [K.sub(G[k][i],G[k][j])for k in range(10)]==[K.scale(K.O,2 if k==i else -2 if k==j else 0)for k in range(10)]
    for fixed in('wall17','nextH4'):
        assert all(pairs['|'.join(sorted((fixed,'seed'+str(i))))]=={'kind':'finite','label':2}for i in(1,2,5,6))
    roots=[[old[names.index(n)][i]for i in PERM]for n in SOURCE_NAMES]
    # Finite equalities below are transported by this exact involution;
    # wall17 and nextH4 are fixed because of the four orthogonalities above.
    return G,I,h,g,old,names,roots,center,pairs,prior

def check(d):
    G,I,h,g,old,names,roots,center,pairs,prior=reconstruct(d);hi=list(map(V.golden,h))
    profiles=[];records=[];counts={'right':0,'pi3':0,'dotted':0}
    mapping={'seed1':'seed2','seed2':'seed1','seed5':'seed6','wall17':'wall17','nextH4':'nextH4'}
    for si,r,new in zip(SOURCE_NAMES,roots,NEW,strict=True):
        finite={n:pairs['|'.join(sorted((si,pre)))] for n,pre in mapping.items() if pairs['|'.join(sorted((si,pre)))]['kind']=='finite'}
        cp=V.dot(r,hi);assert cp[0]>0
        assert V.pair(r,g,r)[0]<=2*V.SCALE<=V.pair(r,g,r)[1]
        dotted=[]
        for a,n in zip(old,names,strict=True):
            pp=V.pair(r,g,a)
            if n in finite:
                lab=finite[n]['label'];assert lab in(2,3);want=0 if lab==2 else -V.SCALE
                assert pp[0]<=want<=pp[1];counts['right'if lab==2 else'pi3']+=1
                pairs['|'.join(sorted((new,n)))]=finite[n]
            else:
                assert pp[1]<-2*V.SCALE,(new,n,pp);counts['dotted']+=1
                dotted.append({'old':n,'pairing':V.enc(pp)});pairs['|'.join(sorted((new,n))) ]={'kind':'dotted'}
        profiles.append(finite)
        records.append({'name':new,'source':si,'norm':2,'finite_profile':finite,'center_pairing_dual':V.enc(cp),'dotted_checks':dotted})
    cross=[]
    for i,j in combinations(range(6),2):
        rel=pairs['|'.join(sorted((SOURCE_NAMES[i],SOURCE_NAMES[j])))];pairs['|'.join(sorted((NEW[i],NEW[j])))]=rel
        pp=V.pair(roots[i],g,roots[j]);row={'walls':[NEW[i],NEW[j]],'relation':rel}
        if rel['kind']=='dotted':assert pp[1]<-2*V.SCALE;counts['dotted']+=1;row['pairing']=V.enc(pp)
        else:
            assert rel['label']in(2,3);want=0 if rel['label']==2 else -V.SCALE;assert pp[0]<=want<=pp[1]
            counts['right'if rel['label']==2 else'pi3']+=1
        cross.append(row)
    assert sum(counts.values())==591 and len(pairs)==102*101//2
    allroots=old+roots;allnames=names+NEW;positive=[];bad=0;dottedcount=0;checked=0
    for count in range(1,7):
        for inds in combinations(range(6),count):
            news=[NEW[i]for i in inds]
            if any(pairs['|'.join(sorted((a,b)))]['kind']=='dotted'for a,b in combinations(news,2)):continue
            common=sorted(set.intersection(*(set(profiles[i])for i in inds)))
            for mask in range(1<<len(common)):
                ss=news+[n for i,n in enumerate(common)if mask>>i&1];checked+=1
                if any(pairs['|'.join(sorted((a,b)))]['kind']=='dotted'for a,b in combinations(ss,2)):dottedcount+=1;continue
                typ=P.T.psd(C.gram(ss,pairs));assert typ!='positive_semidefinite',ss
                if typ=='positive_definite':positive.append(ss)
                else:bad+=1
    vertices=[]
    for ss in positive:
        if len(ss)!=9:continue
        M=C.gram(ss,pairs);Mi=K.inverse(M);rs=[allroots[allnames.index(n)]for n in ss]
        lam=V.mv([[V.golden(x)for x in row]for row in Mi],[V.pair(center,g,r)for r in rs]);assert all(x[1]<0 for x in lam)
        point=[V.sub(center[i],V.total(V.mul(t,r[i])for t,r in zip(lam,rs,strict=True)))for i in range(10)]
        q=V.pair(point,g,point);assert q[1]<0 and V.dot(point,hi)[0]>0
        pp=[]
        for n,r in zip(allnames,allroots,strict=True):
            if n in ss:continue
            p=V.pair(point,g,r);assert p[1]<0;pp.append({'wall':n,'pairing':V.enc(p)})
        vertices.append({'incident':ss,'gram':[[K.enc(x)for x in row]for row in M],'point':list(map(V.enc,point)),'norm':V.enc(q),'strict_nonincident_pairings':pp})
    # Every finite or ideal endpoint supplies a ninth wall with a PSD Gram.
    edge_sets=sorted({tuple(sorted(ss))for v in vertices for ss in combinations(v['incident'],8)});frontier=[]
    for ss in edge_sets:
        ends=[]
        for n in allnames:
            if n in ss or any(pairs['|'.join(sorted((n,a)))]['kind']=='dotted'for a in ss):continue
            typ=P.T.psd(C.gram(list(ss)+[n],pairs));assert typ!='positive_semidefinite'
            if typ=='positive_definite':ends.append(n)
        assert len(ends)in(1,2)
        assert all(frozenset(ss+(n,))in{frozenset(v['incident'])for v in vertices}for n in ends)
        frontier.append({'common_facets':list(ss),'finite_end_walls':ends,'status':'finite'if len(ends)==2 else'open_to_light_cone'})
    vertex_sets={frozenset(v['incident'])for v in vertices};reached={next(iter(vertex_sets))}
    while True:
        before=len(reached)
        for edge in frontier:
            if len(edge['finite_end_walls'])!=2:continue
            ends={frozenset(edge['common_facets']+[n])for n in edge['finite_end_walls']}
            if reached&ends:reached.update(ends)
        if len(reached)==before:break
    assert reached==vertex_sets and len(vertices)==3
    assert max(map(len,profiles))==5
    target=tuple(map(K.parse,d['target']));assert d['target']==V.read(d['source'])['escape']
    assert K.sign(K.pair(target,I,target))>0 and K.sign(K.pair(target,I,h))<0
    ti=list(map(V.golden,target));assert all(V.dot(r,ti)[0]>0 for r in old)
    cuts=[V.dot(r,ti)for r in roots];assert cuts[-1][1]<0
    assert prior['affine_census']['new_affine_circuits']==prior['affine_census']['inherited_affine_circuits']==0
    out={'status':'PASS_EXACT_ISOMETRIC_REMOTE_COMPONENT','source_sha256':V.sha(PREFIX+'_caps.json'),'retained_facets':96,'facets':102,
         'permutation':PERM,'fixed_old_anchor_proof':'wall17 and nextH4 are orthogonal to seeds1,2,5,6; the exact Gram difference identities force their swapped coordinates to agree.',
         'new_normals':records,'normal_enclosures':[list(map(V.enc,r))for r in roots],'new_cross_pairs':cross,'new_pair_counts':counts,
         'affine_census':{'new_candidate_subsets':checked,'dotted_excluded':dottedcount,'positive_definite':len(positive),'indefinite':bad,'new_affine_circuits':0,'inherited_affine_circuits':0},
         'actual_finite_vertices':vertices,'new_component_frontier':{'vertices':len(vertices),'connected':True,'finite_edges':sum(e['status']=='finite'for e in frontier),'open_rays':sum(e['status']=='open_to_light_cone'for e in frontier),'edges':frontier,
         'inherited_edge_effect':'No old rank-eight spherical edge acquires an endpoint: each new wall has at most five finite old neighbors. The new three-vertex component is separate from every inherited finite component.'},
         'target_cut_pairings':list(map(V.enc,cuts)),
         'scope':'Six exact isometric images retain every old96 wall and cut their saved remote gap. All new pairs and every PSD subset involving new walls are checked. Strict center and nonpositive pairings prove all102 inequalities are actual facets. This is not a finite-volume claim.'}
    if 'escape'in d:
        y=tuple(map(K.parse,d['escape']));q=K.pair(y,I,y);f=K.pair(y,I,h);assert K.sign(q)>0 and K.sign(f)<0
        pp=[V.dot(r,list(map(V.golden,y)))for r in allroots];assert min(x[0]for x in pp)>0
        out['remaining_escape']={'norm':K.enc(q),'future_center_pairing':K.enc(f),'strict_pairings':list(map(V.enc,pp))};out['scope']+=' An exact strict future spacelike escape proves infinite volume.'
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=V.read(PREFIX+'_caps.json');out=check(d)
    if a.verify:assert out==V.read(PREFIX+'_verify.json')
    else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_pair_counts'],out['new_component_frontier']['vertices'])
