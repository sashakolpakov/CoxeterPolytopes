#!/usr/bin/env python3
"""Raw88/raw25 forced target, with exact zero-event covers at both old witnesses.

Two small changed local covering LPs propose rational certificates. Verification
uses saved certificates and never invokes an optimizer. No ownership changes.
"""
import argparse,json,time
from fractions import Fraction as F
from itertools import combinations
import fresh_shared_prism_tetrahedron_incidence_20261017 as B
T=B.T
P,A,R,N,W,V,K,M,ROOT=T.P,T.A,T.R,T.N,T.W,T.V,T.K,T.M,T.ROOT
TAG='fresh_raw88_raw25_forced_vertex_frontier_20261018'
SUPPORTS={904:[0,1,2,3,4,5],534:[0,1,2,3,4,8],1022:[0,1,2,3,5,8],690:[0,1,2,4,5,8]}
TARGET=2252
TARGET_SUPPORT=[2,4,5,6,7,8]

def normal(x):return json.loads(json.dumps(x,sort_keys=True))

def geometry():
    scales=[K(5,2)]+[K(1)]*8
    native={(i,j):V.edge(2)for i in range(9)for j in range(9)if i!=j}
    for i,j,e in[(0,6,(0,'5','2')),(1,7,(0,'6','2')),(1,6,V.edge(3)),
        (2,3,V.edge(3)),(3,6,V.edge(3)),(3,7,V.edge(3)),(3,8,V.edge(3)),
        (4,5,V.edge(5)),(5,8,V.edge(3)),(6,7,V.edge(5))]:native[i,j]=native[j,i]=e
    G=[[2*scales[i]if i==j else K()for j in range(9)]for i in range(9)]
    for i,j in combinations(range(9),2):
        sq=V.squared(native[i,j]);z=M.ksqrt(sq*scales[i]*scales[j])if sq else K()
        assert z is not None;G[i][j]=G[j][i]=-z
    assert P.congruence(G)['inertia']==[6,1,2]
    models={R.C.code(g):name for name,g in R.C.models()}
    lookup={tuple(v['forest']):v['index']for v in W.read(R.P.FORESTS)['vertices']}
    records=[];pd=[]
    for ids in combinations(range(9),6):
        inertia=P.congruence(P.sub(G,ids))['inertia'];rec=dict(support=list(ids),inertia=inertia)
        if inertia==[6,0,0]:
            graph={i:{j:native[i,j][0]for j in ids if j!=i and native[i,j][0]!=2}for i in ids}
            assert all(m in(3,4,5,6,10)for row in graph.values()for m in row.values())
            fs=sorted([models[R.C.code(g)]for g in R.C.pieces(graph)]+['A1','H3','H3'])
            rec.update(forest=fs,vertex=lookup[tuple(fs)],in_raw88_prism=7 not in ids and {1,2}<=set(ids),
                in_raw25_prism=6 not in ids and {0,2}<=set(ids))
            pd.append(rec)
        records.append(rec)
    assert len(records)==84 and len(pd)==14
    assert sum(r['in_raw88_prism']for r in pd)==sum(r['in_raw25_prism']for r in pd)==8
    common=[r for r in pd if r['in_raw88_prism']and r['in_raw25_prism']]
    assert {r['vertex']:r['support']for r in common}==SUPPORTS
    outside=[r for r in pd if not(r['in_raw88_prism']or r['in_raw25_prism'])]
    assert [r['vertex']for r in outside]==[1325,TARGET]
    assert next(r['support']for r in outside if r['vertex']==TARGET)==TARGET_SUPPORT
    residual=P.schur(G,[0,1,2],[3,4,5,8]);inv,det=V.inverse(residual)
    assert P.congruence(residual)['inertia']==[3,1,0]and all(x.sign()<0 for row in inv for x in row)
    assert {tuple(sorted([0,1,2]+list(c)))for c in combinations([3,4,5,8],3)}=={tuple(x)for x in SUPPORTS.values()}
    return G,native,dict(squared_scales=[x.out()for x in scales],scaled_gram=P.out(G),
        native_edges=[[i,j,list(native[i,j])]for i,j in combinations(range(9),2)if native[i,j][0]!=2],
        support_records=records,positive_supports=pd,common_vertices=common,outside_both=outside,
        tetrahedron_fixed_active=[0,1,2],tetrahedron_sides=[3,4,5,8],tetrahedron_residual=P.out(residual),
        tetrahedron_inverse=P.out(inv),tetrahedron_inertia=[3,1,0],passive_forest=['A1','H3','H3'])


def frame_recovery(ids,mp,native,target):
    J,active,components=R.projection_data(target);external=sorted(set(range(9))-set(ids))
    covectors=[tuple(sorted((mp[j],tuple(native[i,j]))for j in ids if native[i,j][0]!=2))for i in external]
    assert len(set(covectors))==len(external)==3
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
    return dict(base_embedding=sorted(mp.items()),external_roots=external,external_contacts=normal(covectors),
        squared_heights=[h.out()for h in heights],cross_checks=crosses)


def target_recovery(native,geo):
    rec=next(r for r in geo['outside_both']if r['vertex']==TARGET)
    target=R.P.forest_graph(rec['forest'],dict(R.C.models()));ids=TARGET_SUPPORT
    frame={i:{j:native[i,j][0]for j in ids if j!=i and native[i,j][0]!=2}for i in ids}
    embeddings=R.active_embeddings(frame,target);assert len(embeddings)==12
    records=[];keys=set()
    for mp in embeddings:
        q=frame_recovery(ids,mp,native,target);key=tuple(tuple((i,tuple(e))for i,e in c)for c in q['external_contacts'])
        assert key not in keys;keys.add(key)
        assert q['external_roots']==[0,1,3]
        assert q['squared_heights']==[['5','3'],['11','5'],['11/4','5/4']]
        assert all(c['native_product']==[2,'0','0']for c in q['cross_checks'])
        q['marker']={'chosen_A1':mp[2],'chosen_H3_native5_end':mp[4],
            'chosen_H3_center':mp[5],'chosen_H3_simple_end':mp[8],
            'ordered_I25_ends':[mp[6],mp[7]]}
        assert len(set(target)-set(mp.values()))==7
        records.append(q)
    return dict(vertex=TARGET,forest=rec['forest'],support=ids,capacity=12,
        marker_count_factors=[2,3,2],records=records)


def physical_map(ids,mp,native,nodes,starids,target):
    q=frame_recovery(ids,mp,native,target);external=q['external_roots'];covectors=q['external_contacts'];heights=q['squared_heights']
    for star,absent,face,family in zip(starids,(7,6),(45683,26832),(88,25)):
        s=nodes[star];expected=tuple(sorted(tuple((i,tuple(e))for i,e in covectors[j])for j,i in enumerate(external)if i!=absent))
        assert tuple(sorted(N.native(s)))==expected
        assert s['origin']['face']==face and s['origin']['raw_family']==family and s['origin']['provenance']is None
        hs=[heights[j]for j,i in enumerate(external)if i!=absent]
        assert sorted(hs)==sorted(s['squared_timelike_heights'])
        omitted=[i for i in external if i!=absent]
        assert list(native[omitted[0],omitted[1]])==s['native_omitted_product']
        fixed=({1,2}if family==88 else{0,2})
        passive=set(target)-set(mp.values());assert len(passive)==7
        assert set(target)-set(sum(s['occupied_pairs'],[]))==passive|{mp[i]for i in fixed}
    q['stars']=starids;return q

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
            for absent,face in((7,45683),(6,26832)):
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
    assert [len(r['pair_domain'])for r in result]==[24,48,12,72]
    return result


def check_joint(records,base,di,c):
    mu={v:F(x)for v,x in base['repairs'][di]['vertex_multipliers']};q=12*mu[TARGET]
    weights=[list(map(F,w))for w in c['weights']];a=list(map(F,c['source_costs']));z=list(map(F,c['event_coefficients']))
    cap=F(c['quota_cost']);X=F(c['event_mass']);increments=list(map(F,c['increases']))
    assert len(weights)==len(a)==len(z)==len(increments)==len(c['mixtures'])==len(records)==4
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
    mu={v:F(x)for v,x in base['repairs'][di]['vertex_multipliers']};quota=12*mu[TARGET]
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
    ids=[j for j,x in enumerate(res.x)if x>1e-12];tight=[i for i,x in enumerate(res.ineqlin.residual)if abs(x)<1e-11]
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
    start=time.monotonic();p,changed,rows,patches,old=A.prepare();base=W.read(B.TAG)
    saved=W.read(TAG)if verify else None
    G,native,geo=geometry();recovery=target_recovery(native,geo)
    records=domains(p,changed,rows,patches,native,base)
    certificates=[joint(records,base,di,saved['joint_certificates'][di]if saved else None)for di in range(2)]
    assert all(F(c['gain'])==F(c['event_mass'])==0 and all(F(x)==0 for x in c['increases'])for c in certificates)
    assert all(not any(map(F,c['source_costs']+c['event_coefficients']))and not F(c['quota_cost'])and not any(F(x)for w in c['weights']for x in w)for c in certificates)
    sources={}
    for tag in(B.TAG,P.TAG,N.TAG):
        doc=W.read(tag)
        for path,h in doc['sources'].items():assert W.digest(path)==h
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for mod in(A,N,R,W):
        path=str(__import__('pathlib').Path(mod.__file__).relative_to(ROOT));sources[path]=W.digest(path)
    sources[f'scripts/{TAG}.py']=W.digest(f'scripts/{TAG}.py')
    counts=dict(source_forests=4,active_supports=84,positive_supports=14,extra_vertices=2,target_markers=12,
        physical_pair_maps=sum(len(r['physical_maps'])for r in records),
        maximal_cliques=sum(len(r['maximal_cliques'])for r in records),
        count_event_vectors=sum(len(r['items'])for r in records),columns=0,joint_certificates=2)
    assert counts['physical_pair_maps']==156 and counts['maximal_cliques']==14440 and counts['count_event_vectors']==252
    data=normal(dict(status='PASS_RAW88_RAW25_FORCED_VERTEX_ZERO_EVENT_FRONTIER',sources=sources,baseline_source=B.TAG,
        geometry=geo,target_recovery=recovery,forest_records=records,joint_certificates=certificates,columns=[],repairs=base['repairs'],counts=counts,
        face_multipliers=base['face_multipliers'],restricted_face_multipliers=base['restricted_face_multipliers'],
        local_joint_proposal_count=2,global_optimizer_invocations=0,new_pair_outcomes=0,
        scope='Compact arithmetic H13 over Q(sqrt5). A specific raw88/raw25 prism pair forces two extra finite vertices; complete target-frame recovery gives12v2252 capacity. Four distinct tetrahedral source types count the same pair once. Exact zero-event local covers preserve both current witnesses even under the stronger prohibition of this entire event family. No global realization, dimension theorem, changed ownership or unchanged global optimization.'))
    if saved:assert data==saved
    else:(ROOT/f'results/{TAG}.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    print(data['status'],counts,'bounds',[float(F(r['bound']))for r in data['repairs']],'seconds',round(time.monotonic()-start,2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args();run(args.verify)
