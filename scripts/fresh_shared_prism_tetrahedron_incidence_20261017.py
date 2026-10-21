#!/usr/bin/env python3
"""Actual tetrahedron transports one prism-pair event across four forests.

Two small joint local LPs propose exact covers. Fresh verification uses saved
certificates only; no global model, pair ownership or frozen file is changed.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import fresh_raw23_raw25_cross_vertex_family_20261017 as T

P,A,R,N,W,V,K,M,ROOT=T.P,T.A,T.R,T.N,T.W,T.V,T.K,T.M,T.ROOT
TAG='fresh_shared_prism_tetrahedron_incidence_20261017'
PERM=[0,1,5,8,2,3,6,7,4]
SUPPORTS={1844:[0,1,2,3,4,5],1325:[0,1,2,3,5,8],1640:[0,1,2,4,5,8],982:[0,1,3,4,5,8]}

def normal(x):return json.loads(json.dumps(x,sort_keys=True))

def geometry(reference):
    gr=[[K(*map(F,x))for x in row]for row in reference['active_gram']]
    G=[[K()for _ in range(9)]for _ in range(9)]
    for i in range(9):
        for j in range(9):G[PERM[i]][PERM[j]]=gr[i][j]
    assert P.out(P.sub(G,PERM))==reference['active_gram']
    scales=[K(1),K(F(7,2),F(1,2))]+[K(1)]*7;native={}
    for i,j in combinations(range(9),2):
        sq=G[i][j]*G[i][j]/(scales[i]*scales[j])
        label=next((m for m in(2,3,4,5,6,10)if V.squared(V.edge(m))==sq),None)
        native[i,j]=native[j,i]=V.edge(label)if label else(0,str(sq.a),str(sq.b))
    models={R.C.code(g):name for name,g in R.C.models()}
    lookup={tuple(v['forest']):v['index']for v in W.read(R.P.FORESTS)['vertices']}
    pd=[]
    for r in reference['support_records']:
        if not r['positive_definite']:continue
        ids=sorted(PERM[i]for i in r['active_support'])
        graph={i:{j:native[i,j][0]for j in ids if j!=i and native[i,j][0]!=2}for i in ids}
        fs=sorted([models[R.C.code(g)]for g in R.C.pieces(graph)]+['H3','H4'])
        p25=7 not in ids and {1,5}<=set(ids);p23=6 not in ids and 0 in ids
        pd.append(dict(support=ids,forest=fs,vertex=lookup[tuple(fs)],in_raw25_prism=p25,in_raw23_prism=p23))
    assert len(pd)==15 and sum(r['in_raw25_prism']for r in pd)==8 and sum(r['in_raw23_prism']for r in pd)==10
    common=[r for r in pd if r['in_raw25_prism']and r['in_raw23_prism']]
    assert {r['vertex']:r['support']for r in common}==SUPPORTS
    outside=[r for r in pd if not(r['in_raw25_prism']or r['in_raw23_prism'])]
    assert len(outside)==1 and outside[0]['vertex']==6537
    S=P.schur(G,[0,1,5],[2,3,4,8]);inv,det=V.inverse(S)
    assert P.congruence(S)['inertia']==[3,1,0]and all(x.sign()<0 for row in inv for x in row)
    assert {tuple(sorted([0,1,5]+list(c)))for c in combinations([2,3,4,8],3)}=={tuple(x)for x in SUPPORTS.values()}
    return G,native,dict(reference_to_current_permutation=PERM,scaled_gram=P.out(G),positive_supports=pd,
        common_vertices=common,extra_vertex=outside[0],tetrahedron_fixed_active=[0,1,5],tetrahedron_sides=[2,3,4,8],
        tetrahedron_residual=P.out(S),tetrahedron_inverse=P.out(inv),tetrahedron_inertia=[3,1,0],passive_forest=['H3','H4'])


def physical_map(ids,mp,native,nodes,starids,target):
    J,active,components=R.projection_data(target);external=sorted(set(range(9))-set(ids))
    covectors=[tuple(sorted((mp[j],tuple(native[i,j]))for j in ids if native[i,j][0]!=2))for i in external]
    heights=[]
    for q in covectors:
        field=N.pair_field((q,));projection=R.projection(q,q,J,field)
        assert not any(projection[1:]);height=projection[0]-2;assert height.sign()>0;heights.append(height)
    crosses=[]
    for i,j in combinations(range(3),2):
        edge=tuple(native[external[i],external[j]])
        field=N.pair_field((covectors[i],covectors[j]),(edge,))
        projection=R.projection(covectors[i],covectors[j],J,field);expected=field.negative_contact(edge)
        assert field.signs_difference(tuple(x-y for x,y in zip(projection,expected)),heights[i]*heights[j])==0
        crosses.append(dict(walls=[external[i],external[j]],projection=R.qout(projection),height_product=(heights[i]*heights[j]).out(),native_product=list(edge)))
    for star,absent,face,family,provenance in zip(starids,(6,7),(26448,26918),(23,25),(490,541)):
        s=nodes[star];expected=tuple(sorted(covectors[j]for j,i in enumerate(external)if i!=absent))
        assert tuple(sorted(N.native(s)))==expected
        assert s['origin']['face']==face and s['origin']['raw_family']==family and s['origin']['provenance']==provenance
        hs=[heights[j].out()for j,i in enumerate(external)if i!=absent]
        assert sorted(hs)==sorted(s['squared_timelike_heights'])
        omitted=[i for i in external if i!=absent]
        assert list(native[omitted[0],omitted[1]])==s['native_omitted_product']
        fixed=({0}if family==23 else{1,5})
        passive=set(target)-set(mp.values());assert len(passive)==7
        assert set(target)-set(sum(s['occupied_pairs'],[]))==passive|{mp[i]for i in fixed}
    return dict(base_embedding=sorted(mp.items()),stars=starids,external_roots=external,external_contacts=normal(covectors),
        squared_heights=[h.out()for h in heights],cross_checks=crosses)


def domains(p,changed,rows,patches,native,base):
    result=[]
    for v,ids in SUPPORTS.items():
        old,g,owner,unowned,d,groups,forest,applied=A.current(v,p,changed,rows,patches)
        nodes=A.T.current(v,p,changed)[5];target=R.P.forest_graph(forest,dict(R.C.models()))
        frame={i:{j:native[i,j][0]for j in ids if j!=i and native[i,j][0]!=2}for i in ids}
        lookup={tuple(sorted(N.native(s))):i for i,s in enumerate(nodes)};assert len(lookup)==len(nodes)
        pairs={};checks=[];embeddings=R.active_embeddings(frame,target)
        for mp in embeddings:
            starids=[]
            for absent,face in((6,26448),(7,26918)):
                ext=sorted(set(range(9))-set(ids)-{absent})
                key=tuple(sorted(tuple(sorted((mp[j],tuple(native[i,j]))for j in ids if native[i,j][0]!=2))for i in ext))
                assert key in lookup;star=lookup[key];assert nodes[star]['origin']['face']==face;starids.append(star)
            a,b=starids;assert g.has_edge(a,b)
            q=physical_map(ids,mp,native,nodes,starids,target);checks.append(q)
            pairs.setdefault((a,b),[]).append(len(checks)-1)
        assert len(pairs)==len(embeddings)
        cliques,items=T.enumerate_items(g,owner,len(groups),pairs)
        assert all(q['events']<=1 for q in items)
        demands=[]
        for key in('face_multipliers','restricted_face_multipliers'):
            alpha={i:F(x)for i,x in base[key]}
            demands.append([str(sum((alpha.get(f,F())*m for f,m in q['positive_face_coefficients']),F()))for q in groups])
        result.append(dict(vertex=v,forest=forest,source_support=ids,corner_orbits=groups,literal_stars=nodes,
            graph_edges=sorted(tuple(sorted(e))for e in g.edges()),pair_domain=[list(p)for p in sorted(pairs)],physical_maps=checks,
            maximal_cliques=cliques,items=items,demands=demands))
        print('DOMAIN',v,'physical_pairs',len(pairs),'cliques',len(cliques),'items',len(items),flush=True)
    assert [len(r['pair_domain'])for r in result]==[4,12,8,12]
    return result


def check_joint(records,base,di,c):
    mu={v:F(x)for v,x in base['repairs'][di]['vertex_multipliers']};q=3*mu[6537]
    weights=[list(map(F,w))for w in c['weights']];a=list(map(F,c['source_costs']));z=list(map(F,c['event_coefficients']))
    cap=F(c['quota_cost']);X=F(c['event_mass']);increments=list(map(F,c['increases']))
    assert len(weights)==len(a)==len(z)==len(increments)==len(records)==4
    assert all(0<=x<=1 for x in a)and cap>=max(F(),sum(z,F()))and 0<=X<=q
    assert all(x>=0 for x in increments)
    for r,w,aa,zz,inc,mix in zip(records,weights,a,z,increments,c['mixtures']):
        d=list(map(F,r['demands'][di]));assert len(w)==len(d)and all(x>=0 for x in w)
        assert all(sum(x*y for x,y in zip(w,it['vector']))<=aa+zz*it['events']for it in r['items'])
        terms=[(j,F(x))for j,x in mix];assert len({j for j,x in terms})==len(terms)
        assert all(0<=j<len(r['items'])and x>0 for j,x in terms)
        assert all(sum((x*r['items'][j]['vector'][i]for j,x in terms),F())>=d[i]for i in range(len(d)))
        assert sum((x for j,x in terms),F())<=mu[r['vertex']]+inc
        assert sum((x*r['items'][j]['events']for j,x in terms),F())==X
    load=sum((sum((w*F(d)for w,d in zip(ww,r['demands'][di])),F())for ww,r in zip(weights,records)),F())
    oldcost=sum((aa*mu[r['vertex']]for aa,r in zip(a,records)),F())+cap*q
    gain=sum(increments,F());assert load-oldcost==gain==F(c['gain'])
    num=F(base['repairs'][di]['objective_numerator']);den=F(base['repairs'][di]['normalizer'])+gain
    assert F(c['bound'])==num/den
    return gain


def joint(records,base,di,saved=None):
    if saved is not None:check_joint(records,base,di,saved);return saved
    import numpy as np
    from scipy.optimize import linprog
    mu={v:F(x)for v,x in base['repairs'][di]['vertex_multipliers']};quota=3*mu[6537]
    starts=[];n=0
    for r in records:starts.append(n);n+=len(r['items'])
    delta=list(range(n,n+4));event=n+4;nv=n+5;obj=[0]*nv
    for j in delta:obj[j]=1
    ub=[];rhs=[];labels=[];eq=[]
    for f,r in enumerate(records):
        for i,d in enumerate(r['demands'][di]):
            rr=[0]*nv
            for j,it in enumerate(r['items']):rr[starts[f]+j]=-it['vector'][i]
            ub.append(rr);rhs.append(-F(d));labels.append(('demand',f,i))
        rr=[0]*nv
        for j in range(len(r['items'])):rr[starts[f]+j]=1
        rr[delta[f]]=-1;ub.append(rr);rhs.append(mu[r['vertex']]);labels.append(('cost',f))
        rr=[0]*nv
        for j,it in enumerate(r['items']):rr[starts[f]+j]=it['events']
        rr[event]=-1;eq.append(rr)
    rr=[0]*nv;rr[event]=1;ub.append(rr);rhs.append(quota);labels.append(('quota',))
    scale=max(map(abs,rhs));res=linprog(obj,A_ub=np.array(ub),b_ub=np.array([float(x/scale)for x in rhs]),
        A_eq=np.array(eq),b_eq=np.zeros(4),bounds=(0,None),method='highs',
        options={'time_limit':20,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    assert res.success,res.message
    ids=[j for j,x in enumerate(res.x)if x>1e-10];tight=[i for i,x in enumerate(res.ineqlin.residual)if abs(x)<1e-8]
    vals=W.solve_exact([[r[j]for j in ids]for r in eq]+[[ub[i][j]for j in ids]for i in tight],
        [F()]*4+[rhs[i]for i in tight],len(ids));xx=[F()]*nv
    for j,x in zip(ids,vals):xx[j]=x
    weights=[[F()]*len(r['demands'][di])for r in records];aa=[F()]*4;cap=F()
    for label,x in zip(labels,res.ineqlin.marginals):
        value=F(float(-x)).limit_denominator(1000000)
        if label[0]=='demand':weights[label[1]][label[2]]=value
        elif label[0]=='cost':aa[label[1]]=value
        else:cap=value
    zz=[F(float(-x)).limit_denominator(1000000)for x in res.eqlin.marginals]
    gain=sum((xx[j]for j in delta),F());num=F(base['repairs'][di]['objective_numerator']);den=F(base['repairs'][di]['normalizer'])+gain
    c=dict(weights=[list(map(str,w))for w in weights],source_costs=list(map(str,aa)),event_coefficients=list(map(str,zz)),quota_cost=str(cap),
        event_mass=str(xx[event]),increases=[str(xx[j])for j in delta],gain=str(gain),bound=str(num/den),
        mixtures=[[[j,str(xx[starts[f]+j])]for j in range(len(r['items']))if xx[starts[f]+j]]for f,r in enumerate(records)])
    check_joint(records,base,di,c);return c


def run(verify=False):
    start=time.monotonic();p,changed,rows,patches,old=A.prepare();base=W.read(T.TAG);reference=W.read(P.TAG)
    saved=W.read(TAG)if verify else None
    G,native,geo=geometry(reference);records=domains(p,changed,rows,patches,native,base)
    certificates=[joint(records,base,di,saved['joint_certificates'][di]if saved else None)for di in range(2)]
    columns=[];repairs=[]
    for di,c in enumerate(certificates):
        coefficient=Counter()
        for r,w in zip(records,c['weights']):
            for g,weight in zip(r['corner_orbits'],map(F,w)):
                for f,m in g['positive_face_coefficients']:coefficient[f]+=weight*m
        costs=[[r['vertex'],x]for r,x in zip(records,c['source_costs'])if F(x)]
        if F(c['quota_cost']):costs.append([6537,str(3*F(c['quota_cost']))])
        column=dict(positive_face_coefficients=[[f,str(x)]for f,x in sorted(coefficient.items())if x],vertex_capacities=costs)
        if column not in columns:columns.append(column)
        prior=base['repairs'][di];mu={v:F(x)for v,x in prior['vertex_multipliers']};new=dict(mu)
        for r,inc in zip(records,c['increases']):new[r['vertex']]+=F(inc)
        den=sum(new.values(),F());num=F(prior['objective_numerator']);assert str(num/den)==c['bound']
        repairs.append(dict(domain=prior['domain'],objective_numerator=str(num),normalizer=str(den),bound=str(num/den),
            vertex_updates=[[v,str(mu[v]),str(x)]for v,x in sorted(new.items())if x>mu[v]],vertex_multipliers=[[v,str(x)]for v,x in sorted(new.items())if x]))
    for di,r in enumerate(repairs):
        mu={v:F(x)for v,x in r['vertex_multipliers']};alpha={f:F(x)for f,x in base[('face_multipliers','restricted_face_multipliers')[di]]}
        assert all(sum((F(c)*alpha.get(f,F())for f,c in col['positive_face_coefficients']),F())<=sum((F(c)*mu.get(v,F())for v,c in col['vertex_capacities']),F())for col in columns)
    sources={}
    for tag in(T.TAG,P.TAG):
        doc=W.read(tag)
        for path,h in doc['sources'].items():assert W.digest(path)==h
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for mod in(A,N,R,W):sources[str(__import__('pathlib').Path(mod.__file__).relative_to(ROOT))]=W.digest(str(__import__('pathlib').Path(mod.__file__).relative_to(ROOT)))
    sources[f'scripts/{TAG}.py']=W.digest(f'scripts/{TAG}.py')
    counts=dict(source_forests=4,physical_pair_maps=sum(len(r['physical_maps'])for r in records),
        maximal_cliques=sum(len(r['maximal_cliques'])for r in records),count_event_vectors=sum(len(r['items'])for r in records),
        columns=len(columns),joint_certificates=2)
    assert counts['physical_pair_maps']==36 and counts['maximal_cliques']==2634 and counts['count_event_vectors']==194
    data=normal(dict(status='PASS_SHARED_PRISM_TETRAHEDRON_INCIDENCE',sources=sources,baseline_source=T.TAG,geometry=geo,
        forest_records=records,joint_certificates=certificates,columns=columns,repairs=repairs,counts=counts,
        face_multipliers=base['face_multipliers'],restricted_face_multipliers=base['restricted_face_multipliers'],
        local_joint_proposal_count=2,global_optimizer_invocations=0,new_pair_outcomes=0,
        scope='Compact arithmetic H13 over K=Q(sqrt5). Same actual raw23/raw25 prism pair meets in a compact tetrahedron with one vertex of each of four forest types; their event counts are equal and bounded by3v6537. Exact common-event graph covers and charges certify minimal monotone source-cost repairs at fixed target budget. No global realization, full-model optimization, dimension theorem or change to pair ownership.'))
    if saved:assert data==saved
    else:(ROOT/f'results/{TAG}.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    print(data['status'],counts,'bounds',[float(F(r['bound']))for r in repairs],'seconds',round(time.monotonic()-start,2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args();run(args.verify)
