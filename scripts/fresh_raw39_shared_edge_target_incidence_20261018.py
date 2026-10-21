#!/usr/bin/env python3
"""Actual raw39 edge events exclude all represented Q stars at target3324.

Two changed local LPs propose exact joint cover/charge certificates. Verification
rebuilds geometry and graphs and checks saved rational certificates only.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import networkx as nx
import fresh_raw88_raw25_forced_vertex_frontier_20261018 as B
T=B.T
P,A,R,N,W,V,K,M,ROOT=T.P,T.A,T.R,T.N,T.W,T.V,T.K,T.M,T.ROOT
TAG='fresh_raw39_shared_edge_target_incidence_20261018'
SOURCE=5587
TARGET=3324
SOURCE_SUPPORTS=[[0,1,2,3,4],[0,1,2,3,6]]
TARGET_SUPPORTS=[[0,1,4,6,7],[2,3,4,5,6]]
REFLECTION={0:3,1:2,2:1,3:0,4:6,6:4,5:7,7:5}

def normal(x):return json.loads(json.dumps(x,sort_keys=True))

def geometry():
    tau=K(F(41,2),F(17,2));scales=[K(1)]*5+[tau,K(1),tau]
    native={(i,j):V.edge(2)for i in range(8)for j in range(8)if i!=j}
    for i,j,e in [(0,1,V.edge(3)),(1,2,V.edge(3)),(2,3,V.edge(3)),(3,4,V.edge(3)),
        (0,6,V.edge(3)),(4,6,V.edge(5)),(0,5,(0,'41/2','17/2')),
        (3,7,(0,'41/2','17/2')),(5,7,(0,'349','156'))]:native[i,j]=native[j,i]=e
    G=[[2*scales[i]if i==j else K()for j in range(8)]for i in range(8)]
    for i,j in combinations(range(8),2):
        sq=V.squared(native[i,j]);z=M.ksqrt(sq*scales[i]*scales[j])if sq else K()
        assert z is not None;G[i][j]=G[j][i]=-z
    assert P.congruence(G)['inertia']==[5,1,2]
    assert all(native[i,j]==native[REFLECTION[i],REFLECTION[j]]for i,j in combinations(range(8),2))
    models={R.C.code(g):name for name,g in R.C.models()}
    lookup={tuple(v['forest']):v['index']for v in W.read(R.P.FORESTS)['vertices']}
    records=[];pd=[]
    for ids in combinations(range(8),5):
        inertia=P.congruence(P.sub(G,ids))['inertia'];rec=dict(support=list(ids),inertia=inertia)
        if inertia==[5,0,0]:
            graph={i:{j:native[i,j][0]for j in ids if j!=i and native[i,j][0]!=2}for i in ids}
            assert all(m in(3,4,5,6,10)for row in graph.values()for m in row.values())
            fs=sorted([models[R.C.code(g)]for g in R.C.pieces(graph)]+['A5','H3'])
            rec.update(forest=fs,vertex=lookup[tuple(fs)],in_first_quad=7 not in ids and {1,2,3}<=set(ids),
                in_second_quad=5 not in ids and {0,1,2}<=set(ids));pd.append(rec)
        records.append(rec)
    assert len(records)==56 and len(pd)==12
    assert sum(r['in_first_quad']for r in pd)==sum(r['in_second_quad']for r in pd)==4
    shared=[r for r in pd if r['in_first_quad']and r['in_second_quad']]
    assert [r['support']for r in shared]==SOURCE_SUPPORTS and all(r['vertex']==SOURCE for r in shared)
    outside=[r for r in pd if not(r['in_first_quad']or r['in_second_quad'])]
    assert Counter(r['vertex']for r in outside)=={3324:2,1880:2,2622:2}
    assert [r['support']for r in outside if r['vertex']==TARGET]==TARGET_SUPPORTS
    edge=P.schur(G,[0,1,2,3],[4,6]);edge_inv,_=V.inverse(edge)
    assert P.congruence(edge)['inertia']==[1,1,0]and all(x.sign()<0 for row in edge_inv for x in row)
    simplex=[]
    for fixed,sides in[(7,[0,1,2,4,6]),(5,[1,2,3,4,6])]:
        supports=[sorted([fixed]+list(c))for c in combinations(sides,4)]
        assert all(any(r['support']==s and r['inertia']==[5,0,0]for r in pd)for s in supports)
        res=P.schur(G,[fixed],sides);inv,_=V.inverse(res)
        assert P.congruence(res)['inertia']==[4,1,0]and all(x.sign()<0 for row in inv for x in row)
        simplex.append(dict(fixed=[fixed],sides=sides,supports=supports,residual=P.out(res),inverse=P.out(inv)))
    return G,native,dict(squared_scales=[z.out()for z in scales],scaled_gram=P.out(G),
        native_edges=[[i,j,list(native[i,j])]for i,j in combinations(range(8),2)if native[i,j][0]!=2],
        inertia=[5,1,2],support_records=records,positive_supports=pd,shared_edge_vertices=shared,
        outside_both=outside,passive_forest=['A5','H3'],reflection=sorted(REFLECTION.items()),
        edge_fixed=[0,1,2,3],edge_sides=[4,6],edge_residual=P.out(edge),edge_inverse=P.out(edge_inv),
        actual_simplices=simplex)


def recover(ids,mp,native,target):
    J,active,components=R.projection_data(target);external=sorted(set(range(8))-set(ids))
    cov={i:tuple(sorted((mp[j],tuple(native[i,j]))for j in ids if native[i,j][0]!=2))for i in external}
    assert len(set(cov.values()))==3
    heights={};selfchecks=[]
    for i,c in cov.items():
        field=N.pair_field((c,));pr=R.projection(c,c,J,field)
        assert not any(pr[1:]);height=pr[0]-2;assert height.sign()>0;heights[i]=height
        selfchecks.append(dict(wall=i,projection=R.qout(pr),height=height.out()))
    cross=[]
    for i,j in combinations(external,2):
        edge=tuple(native[i,j]);field=N.pair_field((cov[i],cov[j]),(edge,))
        pr=R.projection(cov[i],cov[j],J,field);expected=field.negative_contact(edge)
        assert field.signs_difference(tuple(x-y for x,y in zip(pr,expected)),heights[i]*heights[j])==0
        cross.append(dict(walls=[i,j],projection=R.qout(pr),height_product=(heights[i]*heights[j]).out(),native_product=list(edge)))
    assert len(set(target)-set(mp.values()))==8
    return cov,heights,dict(support=ids,mapping=sorted(mp.items()),external_roots=external,
        contacts=normal([cov[i]for i in external]),heights=[heights[i].out()for i in external],self_checks=selfchecks,cross_checks=cross)


def source_maps(native,target,nodes,owner,g):
    lookup={tuple(sorted(N.native(s))):i for i,s in enumerate(nodes)};assert len(lookup)==len(nodes)
    records=[];domains=[]
    for ids in SOURCE_SUPPORTS:
        frame={i:{j:native[i,j][0]for j in ids if j!=i and native[i,j][0]!=2}for i in ids}
        embeddings=R.active_embeddings(frame,target);assert len(embeddings)==4;pairs=set()
        for mp in embeddings:
            cov,hs,row=recover(ids,mp,native,target);stars=[]
            for absent,fixed in[(7,{1,2,3}),(5,{0,1,2})]:
                ext=sorted(set(range(8))-set(ids)-{absent});key=tuple(sorted(cov[i]for i in ext))
                assert key in lookup;si=lookup[key];s=nodes[si]
                assert s['origin']['face']==28827 and s['origin']['raw_family']==39 and s['origin']['provenance']is None
                assert sorted(hs[i].out()for i in ext)==sorted(s['squared_timelike_heights'])
                assert list(native[ext[0],ext[1]])==s['native_omitted_product']
                assert set(target)-set(sum(s['occupied_pairs'],[]))==(set(target)-set(mp.values()))|{mp[i]for i in fixed}
                stars.append(si)
            pair=tuple(sorted(stars));assert g.has_edge(*pair)and {owner[i]for i in pair}=={1,3}
            assert sorted(nodes[i]['origin']['omitted']for i in pair)==[[0,1],[0,3]]
            assert pair not in pairs;pairs.add(pair);row['stars']=stars;records.append(row)
        domains.append(pairs)
    assert domains[0]==domains[1]=={(4,13),(5,12),(6,15),(7,14)}
    return sorted(domains[0]),records


def target_maps(native,target):
    records=[];parts=[]
    for ids in TARGET_SUPPORTS:
        frame={i:{j:native[i,j][0]for j in ids if j!=i and native[i,j][0]!=2}for i in ids}
        embeddings=R.active_embeddings(frame,target);assert len(embeddings)==1
        mp=embeddings[0];cov,hs,row=recover(ids,mp,native,target);records.append(row);parts.append((mp,cov,hs))
    assert all(parts[1][0][REFLECTION[i]]==j for i,j in parts[0][0].items())
    assert sorted(parts[0][1].values())==sorted(parts[1][1].values())
    assert sorted(h.out()for h in parts[0][2].values())==sorted(h.out()for h in parts[1][2].values())
    mp,cov,hs=parts[0];ports={mp[i]for i in(0,1,4,6)}
    return cov,hs,ports,records


def exclude_target_stars(target,nodes,cov,hs,ports):
    J,active,components=R.projection_data(target);bad=[];allowed=[]
    for si,s in enumerate(nodes):
        port=next((pair for pair in s['occupied_pairs']if set(pair)<=ports),None)
        if port is not None:bad.append(dict(star=si,kind='actual_simplex_triangle_port',ports=port));continue
        failure=None
        for qi,(c,h)in enumerate(zip(N.native(s),s['squared_timelike_heights'])):
            h=K(*map(F,h))
            for i,ec in cov.items():
                if c==ec:continue
                field=N.pair_field((c,ec));pr=R.projection(c,ec,J,field);hp=h*hs[i]
                detail=dict(star=si,star_external_index=qi,event_wall=i,projection=R.qout(pr),height_product=hp.out())
                if field.signs_difference(pr,hp)>0:
                    square=field.qm(pr,pr);assert not any(square[1:])and (square[0]-hp).sign()>0
                    failure=dict(**detail,kind='positive_actual_product',projection_sign=1,projection_square=square[0].out());break
                code,cert=field.arithmetic_wall_gate(tuple((z.a,z.b)for z in pr),hp.a,hp.b)
                if code is not None:
                    assert code==3 and not any(pr)and not M.integral(hp)
                    failure=dict(**detail,kind='nonintegral_two_cycle',two_cycle=hp.out(),trace=str(2*hp.a),norm=str(hp.a*hp.a-5*hp.b*hp.b));break
            if failure:break
        if failure:bad.append(failure)
        else:allowed.append(si)
    assert len(nodes)==15 and not allowed
    assert Counter(r['kind']for r in bad)=={'positive_actual_product':2,'nonintegral_two_cycle':12,'actual_simplex_triangle_port':1}
    return bad,allowed


def enumerate_events(g,owner,n,events,marker):
    cliques=[];items={};start=time.monotonic()
    for cc in nx.find_cliques(g):
        cl=tuple(sorted(cc));cliques.append(cl)
        assert len(cliques)<100000 and time.monotonic()-start<30
        for k in range(1,len(cl)+1):
            for sub in combinations(cl,k):
                ct=Counter(owner[i]for i in sub if i!=marker);vec=tuple(ct.get(j,0)for j in range(n))
                x=sum(set(e)<=set(sub)for e in events);assert x<=1;key=vec,x
                if key not in items or sub<items[key]:items[key]=sub
    return sorted(cliques),[dict(vector=list(v),events=x,clique=list(items[v,x]))for v,x in sorted(items)]


def domains(p,changed,rows,patches,native,base):
    records=[];target_data=None
    for v in(SOURCE,TARGET):
        old,g,owner,unowned,d,groups,forest,applied=A.current(v,p,changed,rows,patches)
        nodes=A.T.current(v,p,changed)[5];target=R.P.forest_graph(forest,dict(R.C.models()))
        inherited_edges=sorted(tuple(sorted(e))for e in g.edges());marker=None
        if v==SOURCE:events,maps=source_maps(native,target,nodes,owner,g)
        else:
            cov,hs,ports,maps=target_maps(native,target);bad,allowed=exclude_target_stars(target,nodes,cov,hs,ports)
            marker=len(nodes);g.add_node(marker);g.add_edges_from((marker,i)for i in allowed);events=[(marker,)]
            target_data=dict(physical_maps=maps,simplex_ports=sorted(ports),exclusions=bad,compatible_stars=allowed,
                event_normal_key=normal(sorted(cov.values())),same_event_marker_count=1)
        assert sorted(tuple(sorted(e))for e in g.edges()if marker not in e)==inherited_edges
        cliques,items=enumerate_events(g,owner,len(groups),events,marker);demands=[]
        for key in('face_multipliers','restricted_face_multipliers'):
            alpha={i:F(x)for i,x in base[key]}
            demands.append([str(sum((alpha.get(f,F())*m for f,m in gr['positive_face_coefficients']),F()))for gr in groups])
        records.append(dict(vertex=v,forest=forest,finite_frame=[[i,sorted(row.items())]for i,row in sorted(target.items())],corner_orbits=groups,literal_stars=nodes,
            inherited_graph_edges=inherited_edges,graph_edges=sorted(tuple(sorted(e))for e in g.edges()),
            inherited_unowned_pairs=normal(sorted(unowned)),marker=marker,event_domain=events,physical_maps=maps,
            maximal_cliques=cliques,items=items,demands=demands))
        print('DOMAIN',v,'cliques',len(cliques),'items',len(items),'eventmaps',len(maps),flush=True)
    return records,target_data

def check_joint(records,base,di,c):
    mu={v:F(x)for v,x in base['repairs'][di]['vertex_multipliers']}
    weights=[list(map(F,w))for w in c['weights']];a=list(map(F,c['source_costs']));z=list(map(F,c['event_coefficients']))
    cap=F(c['quota_cost']);X=F(c['event_mass']);increments=list(map(F,c['increases']))
    assert len(weights)==len(a)==len(z)==len(increments)==len(c['mixtures'])==len(records)==2
    assert all(0<=x<=1 for x in a)and cap==0 and sum(z,F())<=0 and X>=0
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
    oldcost=sum((aa*mu[r['vertex']]for aa,r in zip(a,records)),F())
    gain=sum(increments,F());assert load-oldcost==gain==F(c['gain'])
    num=F(base['repairs'][di]['objective_numerator']);den=F(base['repairs'][di]['normalizer'])+gain
    assert F(c['bound'])==num/den
    return gain


def joint(records,base,di,saved=None):
    if saved is not None:check_joint(records,base,di,saved);return saved
    import numpy as np
    from scipy.optimize import linprog
    mu={v:F(x)for v,x in base['repairs'][di]['vertex_multipliers']}
    starts=[];n=0
    for r in records:starts.append(n);n+=len(r['items'])
    delta=list(range(n,n+2));event=n+2;nv=n+3;obj=[0]*nv
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
    scale=max(map(abs,rhs));res=linprog(obj,A_ub=np.array(ub),b_ub=np.array([float(x/scale)for x in rhs]),
        A_eq=np.array(eq),b_eq=np.zeros(2),bounds=(0,None),method='highs',
        options={'time_limit':20,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    assert res.success,res.message
    ids=[j for j,x in enumerate(res.x)if x>1e-12];tight=[i for i,x in enumerate(res.ineqlin.residual)if abs(x)<1e-11]
    vals=W.solve_exact([[r[j]for j in ids]for r in eq]+[[ub[i][j]for j in ids]for i in tight],
        [F()]*2+[rhs[i]for i in tight],len(ids));xx=[F()]*nv
    for j,x in zip(ids,vals):xx[j]=x
    weights=[[F()]*len(r['demands'][di])for r in records];aa=[F()]*2;cap=F()
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
    start=time.monotonic();p,changed,rows,patches,old=A.prepare();base=W.read(B.TAG);saved=W.read(TAG)if verify else None
    G,native,geo=geometry();records,target_data=domains(p,changed,rows,patches,native,base)
    certificates=[joint(records,base,di,saved['joint_certificates'][di]if saved else None)for di in range(2)]
    assert all(F(c['gain'])==F(934391,2000000000)and c['source_costs']==['1','1']and c['event_coefficients']==['1','-1']for c in certificates)
    columns=[];repairs=[]
    for di,c in enumerate(certificates):
        coefficient=Counter()
        for r,w in zip(records,c['weights']):
            for group,weight in zip(r['corner_orbits'],map(F,w)):
                for f,m in group['positive_face_coefficients']:coefficient[f]+=weight*m
        col=dict(positive_face_coefficients=[[f,str(x)]for f,x in sorted(coefficient.items())if x],
            vertex_capacities=[[r['vertex'],x]for r,x in zip(records,c['source_costs'])if F(x)])
        if col not in columns:columns.append(col)
        prior=base['repairs'][di];mu={v:F(x)for v,x in prior['vertex_multipliers']};new=dict(mu)
        for r,inc in zip(records,c['increases']):new[r['vertex']]+=F(inc)
        den=sum(new.values(),F());num=F(prior['objective_numerator']);assert num/den==F(c['bound'])
        repairs.append(dict(domain=prior['domain'],objective_numerator=str(num),normalizer=str(den),bound=str(num/den),
            vertex_updates=[[v,str(mu[v]),str(x)]for v,x in sorted(new.items())if x>mu[v]],
            vertex_multipliers=[[v,str(x)]for v,x in sorted(new.items())if x]))
    for di,r in enumerate(repairs):
        mu={v:F(x)for v,x in r['vertex_multipliers']};alpha={f:F(x)for f,x in base[('face_multipliers','restricted_face_multipliers')[di]]}
        assert all(sum((F(c)*alpha.get(f,F())for f,c in col['positive_face_coefficients']),F())<=sum((F(c)*mu.get(v,F())for v,c in col['vertex_capacities']),F())for col in columns)
    sources={}
    for tag in(B.TAG,P.TAG,N.TAG):
        doc=W.read(tag)
        for path,h in doc['sources'].items():assert W.digest(path)==h
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for mod in(A,N,R,W):
        path=str(__import__('pathlib').Path(mod.__file__).relative_to(ROOT));sources[path]=W.digest(path)
    sources[f'scripts/{TAG}.py']=W.digest(f'scripts/{TAG}.py')
    counts=dict(active_supports=56,positive_supports=12,outside_vertices=6,source_maps=8,source_events=4,target_maps=2,
        target_event_types=1,new_event_star_exclusions=15,columns=len(columns),joint_certificates=2,
        maximal_cliques=sum(len(r['maximal_cliques'])for r in records),count_event_vectors=sum(len(r['items'])for r in records))
    assert counts['maximal_cliques']==29 and counts['count_event_vectors']==23 and len(columns)==2
    data=normal(dict(status='PASS_RAW39_SHARED_EDGE_TARGET_INCIDENCE',sources=sources,baseline_source=B.TAG,
        geometry=geo,forest_records=records,target_compatibility=target_data,joint_certificates=certificates,columns=columns,repairs=repairs,counts=counts,
        face_multipliers=base['face_multipliers'],restricted_face_multipliers=base['restricted_face_multipliers'],
        source_occurrences_per_pair=2,target_occurrences_per_pair=2,local_joint_proposal_count=2,
        global_optimizer_invocations=0,old_pair_outcome_changes=0,
        scope='Compact arithmetic H13 over Q(sqrt5). The selected raw39 quadrilateral pair has two actual source-edge endpoints and two unique target events. Full-frame recovery and a native reflection identify the same unlabelled pair. All15 represented target Q stars are excluded by necessary actual-wall or triangle-port conditions. Exact coupled local covers/charges produce two new positive cross-vertex columns and monotone fixed-alpha repairs. No global optimization, global physical realization or dimension theorem.'))
    if saved:assert data==saved
    else:(ROOT/f'results/{TAG}.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    print(data['status'],counts,'bounds',[float(F(r['bound']))for r in repairs],'seconds',round(time.monotonic()-start,2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args();run(args.verify)
