#!/usr/bin/env python3
"""Complete raw108 extra-vertex marker family and relative native tests.

Frozen prior geometry/ownership is retained; no optimizer is invoked.
"""
import argparse,json,time
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import networkx as nx
import fresh_raw14_raw112_raw77_raw108_core_incidence_pilot_20261018 as C
import fresh_raw56_raw57_raw76_joint_marker_family_20261018 as J3
T=C.T;U=T.U;D=T.D;A=T.A;R=T.R;M=T.M;H=T.H;K=T.K;ROOT=T.ROOT
TAG='fresh_raw108_complete_marker_family_20261018'
normal=T.normal
KINDS=(1084,1085)
SUPPORTS={1084:(0,1,3,5,6,7),1085:(0,1,3,4,6,7)}
OMITTED={1084:(2,4),1085:(2,5)}
ROLES={1084:{'H4':(1,3,6,5),'I2(5)':(0,7)},1085:{'A1':(4,),'H3':(1,3,6),'I2(5)':(0,7)},
       56:{'A1':(0,5),'I2(10)':(1,3)},57:{'A2':(0,5),'I2(10)':(1,3)},76:{'A1':(0,4),'H3':(1,3,6)}}


def maps(forest,kind,g):
    ids=SUPPORTS[kind];om=OMITTED[kind];target=R.P.forest_graph(forest,dict(R.C.models()))
    active=C.active_graph(g,ids);J,ac,blocks,tens=T.projection_data(target)
    em=R.active_embeddings(active,target);expected=C.expected_maps(C.forest_type(g,ids),forest)
    assert len(em)==expected;markers=[];seen=set()
    ports=[(i,j)for i,j in combinations(sorted(set((1,3,4,5,6))&set(ids)),2)]
    assert len(ports)==6
    for mp in em:
        assert set(mp)==set(ids)and len(set(mp.values()))==6
        cov=[tuple(sorted((mp[j],tuple(g[i][j]))for j in ids if g[i][j]!=A.RIGHT))for i in om]
        assert len(set(cov))==2 and all(i in ac for c in cov for i,e in c)
        hs=[];selfs=[]
        for c in cov:
            f=T.field_for((c,));p=T.projection(c,c,J,tens,f);assert not any(p[1:]);h=p[0]-2
            assert h.sign()>0;hs.append(h);selfs.append(R.qout(p))
        f=T.field_for(cov,(g[om[0]][om[1]],));pr=T.projection(cov[0],cov[1],J,tens,f)
        prescribed=f.negative_contact(g[om[0]][om[1]])
        assert f.signs_difference(tuple(x-y for x,y in zip(pr,prescribed)),hs[0]*hs[1])==0
        order=sorted(range(2),key=lambda i:cov[i]);native=dict(contacts=[[[j,list(e)]for j,e in cov[i]]for i in order],
            squared_timelike_heights=[hs[i].out()for i in order],native_omitted_product=list(g[om[0]][om[1]]),occupied_pairs=[])
        key=H.frozen(native['contacts']);assert key not in seen;seen.add(key)
        markers.append(dict(index=len(markers),kind=kind,native=native,active_embedding=sorted(mp.items()),
            triangle_port_pairs=sorted(sorted((mp[i],mp[j]))for i,j in ports),
            self_projections=[selfs[i]for i in order],mutual_projection=R.qout(pr),height_product=(hs[0]*hs[1]).out()))
    return target,J,tens,markers


def role(kind,marker,name):
    mp=dict(marker['active_embedding']);return tuple(mp[i]for i in ROLES[kind][name])


def relative(ka,a,kb,b):
    assert ka in KINDS
    ten=(0 if role(ka,a,'I2(5)')==role(kb,b,'I2(5)')else 1 if role(ka,a,'I2(5)')==tuple(reversed(role(kb,b,'I2(5)')))else 2)if kb in KINDS else None
    if ka==kb==1084:return(int(role(ka,a,'H4')!=role(kb,b,'H4')),ten)
    if ka==kb==1085:return(int(role(ka,a,'A1')!=role(kb,b,'A1')),int(role(ka,a,'H3')!=role(kb,b,'H3')),ten)
    if ka==1084 and kb==1085:return(ten,)
    if ka==1084:return()
    assert ka==1085
    if kb==57:return()
    ones=role(kb,b,'A1');x=role(ka,a,'A1')[0];same=ones.index(x)if x in ones else 2
    if kb==56:return(same,)
    assert kb==76
    return(same,int(role(ka,a,'H3')!=role(kb,b,'H3')))


def old_maps(forest,kind,old56,old76):
    data=old76 if kind==76 else old56
    r=next(r for r in data['forest_records']if Counter(r['forest'])==Counter(forest))
    assert normal(r['finite_frame'])==normal(R.P.forest_graph(forest,dict(R.C.models())))
    if kind==76:return r['markers']
    return next(x['markers']for x in r['marker_families']if x['raw_family']==kind)


def relative_tests(g,old56,old76):
    specs=[(1084,1084,['A1']+['H4']*2+['I2(5)']*2,5),
           (1085,1085,['A1']*3+['H3']*2+['I2(5)']*2,11),
           (1084,1085,['A1']*2+['H4','H3']+['I2(5)']*2,3),
           (1084,56,['A1']*5+['H4','I2(5)','I2(10)'],1),
           (1084,57,['A1']*3+['A2','H4','I2(5)','I2(10)'],1),
           (1084,76,['A1']*4+['H4','H3','I2(5)'],1),
           (1085,56,['A1']*6+['H3','I2(5)','I2(10)'],3),
           (1085,57,['A1']*4+['A2','H3','I2(5)','I2(10)'],1),
           (1085,76,['A1']*5+['H3']*2+['I2(5)'],6)]
    out=[]
    for ka,kb,forest,n in specs:
        forest=sorted(forest)
        assert sum(len(dict(R.C.models())[name])for name in forest)==13
        target,J,tens,ma=maps(forest,ka,g);anchor=ma[0]
        mb=maps(forest,kb,g)[3]if kb in KINDS else old_maps(forest,kb,old56,old76)
        reps={}
        for b in mb:
            if ka==kb and b['index']==anchor['index']:continue
            reps.setdefault(relative(ka,anchor,kb,b),b)
        assert len(reps)==n,(ka,kb,len(reps),n)
        records=[]
        for pattern,b in sorted(reps.items()):
            result=T.pair_test(J,tens,anchor['native'],b['native'])
            records.append(dict(pattern=list(pattern),other_marker=b,result=result))
            print('RELATIVE',ka,kb,pattern,result['allowed'],[r['code']for r in result['expressions']],flush=True)
        out.append(dict(kinds=[ka,kb],host_forest=forest,finite_frame=target,anchor=anchor,records=records))
    assert sum(len(t['records'])for t in out)==32
    return out


def pilot():
    geo=M.read(C.TAG);gg=next(x for x in geo['geometries']if x['raw_family']==108)
    g=[[tuple(e)for e in row]for row in gg['native_gram']]
    old56=M.read(T.TAG);old76=M.read(U.TAG)
    tests=relative_tests(g,old56,old76)
    return normal(dict(status='PASS_RAW108_COMPLETE_RELATIVE_PILOT',geometry=gg,relative_tests=tests,
        counts=dict(relative_cases=32,retained_cases=sum(r['result']['allowed']for t in tests for r in t['records'])),
        sources={f'{folder}/{tag}.{ext}':M.digest(f'{folder}/{tag}.{ext}')for tag in(C.TAG,T.TAG,U.TAG,J3.TAG)for folder,ext in(('scripts','py'),('results','json'),('results','md'))},
        optimizer_invocations=0,scope='Complete relative marker patterns only; full-family ordinary transport, capacities and charges remain separate.'))


def allowed(ka,a,kb,b):
    if ka not in KINDS:return allowed(kb,b,ka,a)
    if kb not in KINDS:
        if kb==57:return False
        if ka==1084:return True
        return kb==56 and role(ka,a,'A1')[0]not in role(kb,b,'A1')
    if ka==kb==1085:return role(ka,a,'A1')==role(kb,b,'A1')and role(ka,a,'H3')==role(kb,b,'H3')
    return True


def automorphism(target,base,other):
    graph,components=U.graph_components(target);a=dict(base['active_embedding']);b=dict(other['active_embedding'])
    perm={a[i]:b[i]for i in a};used=set(perm);dest=set(perm.values())
    left=[c for c in components if not set(c)&used];right=[c for c in components if not set(c)&dest]
    assert sum(map(len,left))==sum(map(len,right))==7
    for c in left:
        for k,e in enumerate(right):
            iso=nx.algorithms.isomorphism.GraphMatcher(graph.subgraph(c),graph.subgraph(e),edge_match=lambda a,b:a['label']==b['label'])
            if iso.is_isomorphic():perm.update(next(iso.isomorphisms_iter()));right.pop(k);break
        else:raise AssertionError('unmatched passive component')
    assert not right and set(perm)==set(perm.values())==set(range(13))
    assert all(target[i].get(j,2)==target[perm[i]].get(perm[j],2)for i,j in combinations(range(13),2))
    return perm


def process(rows,g,p,domains,local,patches,oldby):
    row=rows[0];v=row['target_vertex'];forest=row['target_forest']
    target=R.P.forest_graph(forest,dict(R.C.models()));J,ac,blocks,tens=T.projection_data(target)
    nodes,groups,graph,owner,got,history=U.complete_current(v,row,p,domains,local,patches)
    assert got==forest
    for old in oldby.get(v,[]):
        assert old['forest']==forest and normal(old['finite_frame'])==normal(target)
        assert (old['old_nodes'],old['old_groups'],old['old_graph_edges'],old['prior_history'])==(normal(nodes),normal(groups),normal(sorted(tuple(sorted(e))for e in graph.edges())),normal(history))
    canonical=[U.transport_native(s,dict(zip(range(13),range(13))))for s in nodes]
    oldkeys={H.frozen(s['contacts']):i for i,s in enumerate(canonical)};assert len(oldkeys)==len(nodes)
    newf=[]
    for row in rows:
        kind=row['kind'];_,_,_,markers=maps(forest,kind,g);assert len(markers)==row['active_role_orientations']
        keys={H.frozen(m['native']['contacts']):i for i,m in enumerate(markers)}
        assert len(keys)==len(markers)and not set(oldkeys)&set(keys)
        anchor=[]
        for i,s in enumerate(nodes):
            ports=[pair for pair in markers[0]['triangle_port_pairs']if pair in s['occupied_pairs']]
            result=dict(allowed=False,triangle_port_exclusion=ports,expressions=[])if ports else T.pair_test(J,tens,s,markers[0]['native'])
            anchor.append(dict(ordinary=i,result=result))
        good={r['ordinary']for r in anchor if r['result']['allowed']};transports=[];edges=[]
        for marker in markers:
            perm=automorphism(target,markers[0],marker);on=[];mn=[]
            for i,s in enumerate(nodes):
                t=U.transport_native(s,perm);j=oldkeys[H.frozen(t['contacts'])]
                assert all(normal(t[k])==normal(canonical[j][k])for k in H.NATIVE_FIELDS)and owner[i]==owner[j];on.append(j)
            for m in markers:
                t=U.transport_native(m['native'],perm);j=keys[H.frozen(t['contacts'])]
                assert all(normal(t[k])==normal(markers[j]['native'][k])for k in H.NATIVE_FIELDS);mn.append(j)
            assert len(set(on))==len(nodes)and len(set(mn))==len(markers)and mn[0]==marker['index']
            edges.extend([on[i],marker['index']]for i in sorted(good))
            transports.append(dict(frame_permutation=sorted(perm.items()),ordinary_permutation=on,marker_permutation=mn))
        co=Counter(forest);cap=2*co['I2(5)']*(co['H4']if kind==1084 else 1)
        if kind==1084:chosen=list(range(len(markers)))
        else:
            chosen=[m['index']for m in markers if role(kind,m,'A1')==role(kind,markers[0],'A1')and role(kind,m,'H3')==role(kind,markers[0],'H3')]
        assert len(chosen)==cap and all(allowed(kind,markers[i],kind,markers[j])for i,j in combinations(chosen,2))
        newf.append(dict(face=row['face'],kind=kind,omitted=row['omitted'],actual_extra_occurrences=1,
            markers=markers,anchor_decisions=anchor,transports=transports,marker_ordinary_edges=sorted(edges),
            capacity=cap,attaining_marker_clique=chosen))
    keys=[H.frozen(m['native']['contacts'])for f in newf for m in f['markers']]
    prior=[]
    for old in oldby.get(v,[]):
        prior.extend(old['markers']if 'markers'in old else[m for f in old['marker_families']for m in f['markers']])
    assert len(keys)==len(set(keys))and not set(keys)&{H.frozen(m['native']['contacts'])for m in prior}
    counts=dict(ordinary_stars=len(nodes),markers=len(keys),anchor_tests=len(nodes)*len(rows),
        ordinary_pair_positions=len(nodes)*len(keys),old_marker_pair_positions=len(prior)*len(keys),
        new_marker_pair_positions=len(keys)*(len(keys)-1)//2,
        anchor_unknowns=sum(z['result']['allowed']for f in newf for z in f['anchor_decisions']),
        retained_ordinary_pairs=sum(len(f['marker_ordinary_edges'])for f in newf))
    return normal(dict(vertex=v,forest=forest,finite_frame=target,old_nodes=nodes,old_groups=groups,
        old_graph_edges=sorted(tuple(sorted(e))for e in graph.edges()),prior_history=history,
        marker_families=newf,counts=counts))


def solve(a,b):
    a=[list(map(F,row))+[F(x)]for row,x in zip(a,b)];n=len(a)
    for i in range(n):
        piv=next((k for k in range(i,n)if a[k][i]),None)
        if piv is None:return None
        a[i],a[piv]=a[piv],a[i];z=a[i][i];a[i]=[x/z for x in a[i]]
        for k in range(n):
            if k!=i:
                z=a[k][i];a[k]=[x-z*y for x,y in zip(a[k],a[i])]
    return tuple(r[-1]for r in a)


def envelope(row,old76,old56):
    forest=row['forest'];co=Counter(forest);a,h,t=co['A1'],co['H3'],co['I2(10)']
    fam={f['kind']:f['markers']for f in row['marker_families']}
    if old76:fam[76]=old76['markers']
    if old56:fam.update({f['raw_family']:f['markers']for f in old56['marker_families']})
    cA=2*co['H4']*co['I2(5)']if 1084 in fam else 0;cB=2*co['I2(5)']if 1085 in fam else 0
    c56=(a//2)*t if 56 in fam else 0;c76=min(a-1,h)if 76 in fam else 0;c57=t if 57 in fam else 0
    boundary=[(c56,0,0,cA,0)]if c56 or cA else[]
    for k in range(1,c76+1):boundary.append((((a-k-1)//2)*t if 56 in fam else 0,0,k,cA,0))
    if cB:boundary.append((((a-1)//2)*t if 56 in fam else 0,0,0,cA,cB))
    if c57:boundary.append((0,c57,0,0,0))
    boundary=sorted(set(p for p in boundary if any(p)))
    assert boundary
    order=(56,57,76,1084,1085);active=[i for i,k in enumerate(order)if k in fam and k!=57];n=len(active)
    constraints=[(tuple(F(row[i])for i in active),F(1))for row in boundary if any(row[i]for i in active)]
    constraints +=[(tuple(F(int(i==j))for i in range(n)),F())for j in range(n)]
    vertices=set()
    for cc in combinations(constraints,n):
        x=solve([p for p,b in cc],[b for p,b in cc])
        if x is not None and min(x)>=0 and all(sum(F(row[i])*x[j]for j,i in enumerate(active))<=1 for row in boundary):vertices.add(x)
    maximal=sorted(x for x in vertices if not any(x!=y and all(a<=b for a,b in zip(x,y))for y in vertices))
    weights=[]
    for x in maximal:
        w=[F()]*5
        for i,z in zip(active,x):w[i]=z
        if c57:w[1]=F(1,c57)
        assert all(sum(n*z for n,z in zip(p,w))<=1 for p in boundary);weights.append(w)
    # Literal graph witnesses show that the count upper rows are jointly
    # attainable in the conservative graph; no physical realization follows.
    target={int(i):r for i,r in row['finite_frame'].items()};ones=sorted(i for i,r in target.items()if not r)
    witnesses=[]
    for q,r,k,ca,cb in boundary:
        chosen=[];used=set()
        if k:
            blocks=sorted({role(76,m,'H3')for m in fam[76]})
            for j in range(k):
                m=next(m for m in fam[76]if role(76,m,'A1')==(ones[0],ones[j+1])and role(76,m,'H3')==blocks[j])
                chosen.append((76,m['index']))
            used.update(ones[:k+1])
        if cb:
            base=fam[1085][0];used.add(role(1085,base,'A1')[0])
            chosen.extend((1085,m['index'])for m in fam[1085]if allowed(1085,base,1085,m))
        if ca:chosen.extend((1084,m['index'])for m in fam[1084])
        if q:
            avail=[i for i in ones if i not in used];pairs=set(zip(avail[::2],avail[1::2]))
            chosen.extend((56,m['index'])for m in fam[56]if role(56,m,'A1')in pairs and role(56,m,'I2(10)')[0]<role(56,m,'I2(10)')[1])
        if r:
            base=role(57,fam[57][0],'A2')
            chosen.extend((57,m['index'])for m in fam[57]if role(57,m,'A2')==base and role(57,m,'I2(10)')[0]<role(57,m,'I2(10)')[1])
        assert Counter(kind for kind,i in chosen)==Counter({kind:z for kind,z in zip(order,(q,r,k,ca,cb))if z})
        for(ka,i),(kb,j)in combinations(chosen,2):
            assert allowed(ka,fam[ka][i],kb,fam[kb][j])if ka in KINDS or kb in KINDS else J3.compatible(ka,fam[ka][i],kb,fam[kb][j])
        witnesses.append(dict(counts=[q,r,k,ca,cb],marker_clique=chosen))
    return dict(kind_order=list(order),boundary_counts=boundary,coefficient_vertices=weights,attaining_marker_cliques=witnesses)


def export(rows,old56,old76,joint,model,domains):
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']};baseline={v:F(x)for v,x in joint['exported']['repair']['vertex_multipliers']}
    assert dual['face_multipliers']==joint['exported']['repair']['face_multipliers']
    newmu=dict(baseline);oldjoint={r['vertex']:r for r in joint['forest_records']}
    older={r['vertex']:r for r in old76['forest_records']};oldest={r['vertex']:r for r in old56['forest_records']}
    candidates=defaultdict(list);vectors={v:d['domain']['maximal_count_vectors']for v,d in domains.items()}
    for tag in(D.L.TAG,D.B.TAG,D.THIRD):
        for r in M.read(tag)['forest_records']:
            candidates[r['vertex']].append((tag,r['certificate']['weights']))
            if tag!=D.L.TAG:vectors[r['vertex']]=r['maximal_count_vectors']
    columns=[];quotas=[]
    for r in rows:
        v=r['vertex'];groups=r['old_groups'];owner={i:z['index']for z in groups for i in z['literal_stars']}
        omitted={owner[i]for f in r['marker_families']for i,k in f['marker_ordinary_edges']}
        oldterms={}
        if v in oldjoint:
            oj=oldjoint[v];omitted.update(oj['ordinary_omitted_orbits'])
            if v in older:oldterms[76]=older[v]['face']
            if v in oldest:oldterms.update({f['raw_family']:f['face']for f in oldest[v]['marker_families']})
        omitted=sorted(omitted);newfamilies={f['kind']:f for f in r['marker_families']}
        env=envelope(r,older.get(v),oldest.get(v))
        originalc=candidates[v]if groups else[('marker_only',[])];vv=vectors[v]if groups else[[]];assert originalc
        unique={}
        for tag,ww in originalc:
            original=list(map(F,ww))+[F()]*(len(groups)-len(ww));assert len(original)==len(groups)and all(x>=0 for x in original)
            assert all(len(x)==len(original)and sum(w*n for w,n in zip(original,x))<=1 for x in vv)
            w=[F()if i in omitted else x for i,x in enumerate(original)]
            for mw in env['coefficient_vertices']:
                terms=[[k,oldterms[k]if k in oldterms else newfamilies[k]['face'],str(weight)]for k,weight in zip(env['kind_order'],mw)if weight]
                cf=Counter()
                for weight,z in zip(w,groups):
                    for face,n in z['positive_face_coefficients']:cf[face]+=weight*n
                for k,face,weight in terms:cf[face]+=F(weight)
                key=tuple((f,x)for f,x in sorted(cf.items())if x);load=sum((alpha.get(f,F())*x for f,x in key),F())
                col=dict(vertex=v,source_charge=tag,other_charge_sources=[],old_orbit_weights=list(map(str,original)),
                    omitted_orbits=omitted,retained_orbit_weights=list(map(str,w)),marker_terms=terms,
                    positive_face_coefficients=[[f,str(x)]for f,x in key],vertex_cost='1',load=str(load),
                    baseline_cost=str(baseline.get(v,F())),violation=str(load-baseline.get(v,F())),
                    old_vector_count=len(vv),old_vectors_sha256=D.L.jhash(vv))
                if key in unique:unique[key]['other_charge_sources'].append(tag)
                else:unique[key]=col
        cs=list(unique.values());columns.extend(cs);newmu[v]=max([baseline.get(v,F())]+[F(c['load'])for c in cs])
        quotas.append(dict(vertex=v,baseline_cost=str(baseline.get(v,F())),new_cost=str(newmu[v]),increase=str(newmu[v]-baseline.get(v,F())),
            omitted_orbits=omitted,envelope=env))
    preserved=[]
    for tag,data in((U.TAG,old76),(T.TAG,old56),(J3.TAG,joint)):
        for c in data['exported']['columns']:
            load=sum((alpha.get(f,F())*F(x)for f,x in c['positive_face_coefficients']),F());v=c['vertex']
            assert c['vertex_cost']=='1'and load<=baseline.get(v,F())<=newmu.get(v,F())
            preserved.append(dict(source=tag,vertex=v,load=str(load),cost=str(newmu.get(v,F()))))
    multi=[]
    for c in model['prepared']['inherited_multi_vertex_columns']:
        load=sum((alpha.get(f,F())*n for f,n in c['integer_face_coefficients']),F())
        oldcost=sum((baseline.get(v,F())*n for v,n in c['integer_vertex_capacities']),F())
        cost=sum((newmu.get(v,F())*n for v,n in c['integer_vertex_capacities']),F());assert load<=oldcost<=cost
        multi.append(dict(variable=c['variable'],all_vertex_costs=c['integer_vertex_capacities'],load=str(load),baseline_cost=str(oldcost),new_cost=str(cost)))
    assert len(preserved)==2884 and len(multi)==14 and all(newmu.get(v,F())>=x for v,x in baseline.items())
    num=F(dual['objective_numerator']);den=sum(newmu.values(),F());gap=num/12-den
    return dict(columns=columns,target_allowances=quotas,preserved_columns=preserved,preserved_multi_owner_costs=multi,
        repair=dict(face_multipliers=dual['face_multipliers'],vertex_multipliers=[[v,str(x)]for v,x in sorted(newmu.items())if x],
            objective_numerator=str(num),normalizer=str(den),bound=str(num/den),baseline_bound=joint['exported']['repair']['bound'],
            remaining_normalizer_gap_to12=str(gap),scope='Least monotone fixed-alpha repair for exported and prior columns, not an optimum of enlarged cones or a global optimization.'),
        counts=dict(exported_columns=len(columns),strictly_violated_columns=sum(F(c['violation'])>0 for c in columns),
            raised_allowances=sum(F(q['increase'])>0 for q in quotas)))


def run():
    start=time.monotonic();inp=D.Q.Inputs();geo=inp.read(C.TAG);old56=inp.read(T.TAG);old76=inp.read(U.TAG);joint=inp.read(J3.TAG);model=inp.read(T.MODEL)
    for data in(geo,old56,old76,joint):
        for path,h in data['sources'].items():assert M.digest(path)==h,path
    gg=next(x for x in geo['geometries']if x['raw_family']==108);g=[[tuple(e)for e in row]for row in gg['native_gram']]
    assert gg['required_quads']==3 and gg['native_automorphisms']==[list(range(8))]
    assert [x['omitted']for x in gg['extra_vertices']]==[[2,4],[2,5]]
    tests=relative_tests(g,old56,old76)
    for t in tests:
        ka,kb=t['kinds']
        for z in t['records']:assert z['result']['allowed']==allowed(ka,t['anchor'],kb,z['other_marker'])
    family=[r for r in geo['complete_family_inventory']if r['raw_family']==108]
    rows=[dict(r,kind=1084 if r['omitted']==[2,4]else 1085)for r in geo['extra_marker_inventory']if r['raw_family']==108]
    assert len(family)==196 and len(rows)==392 and sum(r['active_role_orientations']for r in rows)==1676
    assert all(r['actual_extra_occurrences']==1 for r in rows)
    for kind in KINDS:assert len({r['target_vertex']for r in rows if r['kind']==kind})==196
    domains,represented=D.L.current_domains();p=D.Z.prepare();local={r['vertex']:r for r in inp.read(D.L.TAG)['forest_records']}
    patches={r['vertex']:r for r in inp.read(D.Z.TAG)['forest_records']};grouped=defaultdict(list);oldby=defaultdict(list)
    for r in rows:grouped[r['target_vertex']].append(r)
    for data in(old56,old76):
        for r in data['forest_records']:oldby[r['vertex']].append(r)
    print('PREFLIGHT',len(grouped),len(rows),sum(r['active_role_orientations']for r in rows),flush=True)
    records=[]
    for i,(v,rr)in enumerate(sorted(grouped.items()),1):
        record=process(rr,g,p,domains,local,patches,oldby);records.append(record)
        if i%25==0 or record['counts']['anchor_unknowns']:
            print('FOREST',i,'/',len(grouped),v,record['counts'],'seconds',round(time.monotonic()-start,2),flush=True)
    exported=export(records,old56,old76,joint,model,domains)
    paths=set(inp.paths)|set(p['sources'])
    for mod in(C,J3,T,U,D,D.L,R,A,H,M):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    for tag in(C.TAG,J3.TAG,T.TAG,U.TAG,D.L.TAG,D.B.TAG,D.THIRD,T.MODEL):
        paths.update(f'{folder}/{tag}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    paths.add(str(Path(__file__).relative_to(ROOT)))
    counts=dict(classes=196,pointed_marker_domains=392,target_forests=len(records),relative_cases=32,
        retained_relative_cases=sum(z['result']['allowed']for t in tests for z in t['records']),
        **{key:sum(r['counts'][key]for r in records)for key in records[0]['counts']},**exported['counts'])
    return normal(M.enc(dict(status='PASS_RAW108_COMPLETE_EXTRA_VERTEX_MARKER_FAMILY',sources={p:M.digest(p)for p in sorted(paths)},
        geometry=gg,family=family,pointed_family=rows,relative_tests=tests,forest_records=records,exported=exported,counts=counts,
        optimizer_invocations=0,old_outcomes_changed=False,scope='Complete raw108 family with two once-per-aggregate extra vertices, one allowance shared with ordinary and raw56/57/76 markers. Retained compatibility is unknown, not realization. Only exported columns and fixed-mass repair are asserted.')))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    if args.pilot:
        out=pilot();path=ROOT/f'results/{TAG}_pilot.json';path.write_text(json.dumps(out,indent=2)+'\n')
    else:
        out=run();path=ROOT/f'results/{TAG}.json'
        if args.verify:assert out==json.loads(path.read_text())
        else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],out.get('exported',{}).get('repair',{}).get('bound'),flush=True)
