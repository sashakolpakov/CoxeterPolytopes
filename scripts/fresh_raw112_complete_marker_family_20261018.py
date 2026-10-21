#!/usr/bin/env python3
"""Complete raw112 relative-role pilot against the frozen five-family rule.

No old compatibility outcome is changed and no optimizer is invoked.
"""
import argparse,json,time
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_raw108_complete_marker_family_20261018 as N
C=N.C;T=N.T;U=N.U;D=N.D;A=N.A;R=N.R;M=N.M;H=N.H;ROOT=N.ROOT
TAG='fresh_raw112_complete_marker_family_20261018'
KINDS=(1124,1125)
ROLES={**N.ROLES,1124:{'A1':(0,7),'H4':(1,3,6,5)},1125:{'A1':(0,4,7),'H3':(1,3,6)}}
normal=N.normal
ORDER=(56,57,76,1084,1085,1124,1125)
ALLOWED={(1124,1124):{(0,1,1),(2,3,0),(2,3,1)},(1125,1125):{(3,1,4,0)},
    (1124,1125):{(0,2,1),(2,3,4)},(1124,56):{(2,3)},(1124,57):set(),(1124,76):{(2,3)},
    (1125,56):{(3,4)},(1125,57):set(),(1125,76):set(),
    (1124,1084):{(0,),(1,)},(1124,1085):{(2,)},(1125,1084):{()},(1125,1085):{(1,0)}}


def field_for(covectors,tens):
    # Include the non-K I2(10) inverse only when a contacted block needs it.
    # An unrelated passive I2(10) squareclass must not enlarge a wall field.
    contacts={i for c in covectors for i,e in c}
    extra=((0,str(T.TAU.a),str(T.TAU.b)),)if any(contacts&set(t)for t in tens)else()
    return A.pair_field(tuple(covectors),extra)


def projection(c,d,J,tens,field):
    contacts={i for p in(c,d)for i,e in p}
    active=[t for t in tens if contacts&set(t)]
    return T.projection(c,d,J,active,field)if active else R.projection(c,d,J,field)


def pair_test(J,tens,a,b):
    expressions=[]
    for i,c in enumerate(A.native(a)):
        for j,d in enumerate(A.native(b)):
            if c==d:continue
            field=field_for((c,d),tens);hs=[]
            for p in(c,d):
                z=projection(p,p,J,tens,field);assert not any(z[1:]);h=z[0]-2;assert h.sign()>0;hs.append(h)
            z=projection(c,d,J,tens,field);hp=hs[0]*hs[1]
            code,cert=field.arithmetic_wall_gate(tuple((x.a,x.b)for x in z),hp.a,hp.b)
            sign=field.signs_difference(z,hp)
            if sign>0:code='positive'
            elif sign<0 and field.signs_difference(tuple(x+(T.K(1)if k==0 else T.K())for k,x in enumerate(z)),hp)>0:code='gap'
            expressions.append(dict(walls=[i,j],contacts=normal((c,d)),projection=R.qout(z),
                heights=[h.out()for h in hs],height_product=hp.out(),code=code,gate_certificate=cert))
            if code is not None:return dict(allowed=False,expressions=expressions)
    return dict(allowed=True,expressions=expressions)


def maps(forest,kind,g):
    # The two omission supports and simplex-side labels coincide with raw108.
    # N.maps derives all products, active types and expected map counts from g.
    assert kind in KINDS
    target,J,tens,markers=N.maps(forest,1084 if kind==1124 else 1085,g)
    for m in markers:m['kind']=kind
    return target,J,tens,markers


def role(kind,m,name):
    p=dict(m['active_embedding']);return tuple(p[i]for i in ROLES[kind].get(name,()))


def relative(ka,a,kb,b):
    assert ka in KINDS
    left=role(ka,a,'A1');right=role(kb,b,'A1');names={i:j for j,i in enumerate(left)}
    for i in right:
        if i not in names:names[i]=len(names)
    signature=tuple(names[i]for i in right)
    for kind in('H3','H4'):
        if kind in ROLES[ka]and kind in ROLES[kb]:signature+=(int(role(ka,a,kind)!=role(kb,b,kind)),)
    return signature


def compatible(ka,a,kb,b):
    if ka not in KINDS and kb not in KINDS:
        return N.allowed(ka,a,kb,b)if ka in N.KINDS or kb in N.KINDS else N.J3.compatible(ka,a,kb,b)
    if ka not in KINDS or(ka==1125 and kb==1124):return compatible(kb,b,ka,a)
    return relative(ka,a,kb,b)in ALLOWED[ka,kb]


def old_maps(forest,kind,old56,old76,old108):
    if kind in(56,57,76):return N.old_maps(forest,kind,old56,old76)
    assert kind in N.KINDS
    row=next(r for r in old108['forest_records']if r['forest']==forest)
    target=R.P.forest_graph(forest,dict(R.C.models()))
    assert normal(target)==normal(row['finite_frame'])
    return next(f['markers']for f in row['marker_families']if f['kind']==kind)


def specs():
    return[(1124,1124,['A1']*5+['H4']*2,13),
        (1125,1125,['A1']*7+['H3']*2,67),
        (1124,1125,['A1']*6+['H4','H3'],13),
        (1124,56,['A1']*7+['H4','I2(10)'],7),
        (1124,57,['A1']*5+['H4','A2','I2(10)'],1),
        (1124,76,['A1']*6+['H4','H3'],7),
        (1125,56,['A1']*8+['H3','I2(10)'],13),
        (1125,57,['A1']*6+['H3','A2','I2(10)'],1),
        (1125,76,['A1']*7+['H3']*2,26),
        (1124,1084,['A1']*3+['H4']*2+['I2(5)'],2),
        (1124,1085,['A1']*4+['H4','H3','I2(5)'],3),
        (1125,1084,['A1']*4+['H4','H3','I2(5)'],1),
        (1125,1085,['A1']*5+['H3']*2+['I2(5)'],8)]


def relative_tests(g,old56,old76,old108):
    out=[];cache={};models=dict(R.C.models())
    for ka,kb,forest,expected in specs():
        forest=sorted(forest);assert sum(len(models[x])for x in forest)==13
        def mm(kind):
            key=(tuple(forest),kind)
            if key not in cache:cache[key]=maps(forest,kind,g)
            return cache[key]
        target,J,tens,aa=mm(ka);anchor=aa[0]
        bb=mm(kb)[3]if kb in KINDS else old_maps(forest,kb,old56,old76,old108)
        reps={}
        for b in bb:
            if ka==kb and b['index']==anchor['index']:continue
            reps.setdefault(relative(ka,anchor,kb,b),b)
        assert len(reps)==expected,(ka,kb,len(reps),expected)
        records=[]
        for pattern,b in sorted(reps.items()):
            result=pair_test(J,tens,anchor['native'],b['native'])
            records.append(dict(pattern=list(pattern),other_marker=b,result=result))
            print('RELATIVE',ka,kb,pattern,result['allowed'],[x['code']for x in result['expressions']],flush=True)
        out.append(dict(kinds=[ka,kb],host_forest=forest,finite_frame=target,anchor=anchor,records=records))
    assert sum(len(t['records'])for t in out)==162
    return out


def process(rows,g,p,domains,local,patches,oldby):
    row=rows[0];v=row['target_vertex'];forest=row['target_forest']
    target=R.P.forest_graph(forest,dict(R.C.models()));J,ac,blocks,tens=T.projection_data(target)
    nodes,groups,graph,owner,got,history=U.complete_current(v,row,p,domains,local,patches)
    assert got==forest
    for old in oldby.get(v,[]):
        assert old['forest']==forest and normal(old['finite_frame'])==normal(target)
        assert(old['old_nodes'],old['old_groups'],old['old_graph_edges'],old['prior_history'])==(normal(nodes),normal(groups),normal(sorted(tuple(sorted(e))for e in graph.edges())),normal(history))
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
            result=dict(allowed=False,triangle_port_exclusion=ports,expressions=[])if ports else pair_test(J,tens,s,markers[0]['native'])
            anchor.append(dict(ordinary=i,result=result))
        good={r['ordinary']for r in anchor if r['result']['allowed']};transports=[];edges=[]
        for marker in markers:
            perm=N.automorphism(target,markers[0],marker);on=[]
            for i,s in enumerate(nodes):
                moved=U.transport_native(s,perm);j=oldkeys[H.frozen(moved['contacts'])]
                assert all(normal(moved[k])==normal(canonical[j][k])for k in H.NATIVE_FIELDS)and owner[i]==owner[j];on.append(j)
            moved=U.transport_native(markers[0]['native'],perm)
            assert all(normal(moved[k])==normal(marker['native'][k])for k in H.NATIVE_FIELDS)
            assert len(set(on))==len(nodes)
            edges.extend([on[i],marker['index']]for i in sorted(good))
            transports.append(dict(frame_permutation=sorted(perm.items()),ordinary_permutation=on))
        newf.append(dict(face=row['face'],kind=kind,omitted=row['omitted'],actual_extra_occurrences=1,
            markers=markers,anchor_decisions=anchor,transports=transports,marker_ordinary_edges=sorted(edges)))
    keys=[H.frozen(m['native']['contacts'])for f in newf for m in f['markers']];prior=[]
    for old in oldby.get(v,[]):prior.extend(old['markers']if'markers'in old else[m for f in old['marker_families']for m in f['markers']])
    assert len(keys)==len(set(keys))and not set(keys)&{H.frozen(m['native']['contacts'])for m in prior}
    counts=dict(ordinary_stars=len(nodes),markers=len(keys),anchor_tests=len(nodes)*len(rows),
        ordinary_pair_positions=len(nodes)*len(keys),old_marker_pair_positions=len(prior)*len(keys),
        new_marker_pair_positions=len(keys)*(len(keys)-1)//2,
        anchor_unknowns=sum(z['result']['allowed']for f in newf for z in f['anchor_decisions']),
        retained_ordinary_pairs=sum(len(f['marker_ordinary_edges'])for f in newf))
    return normal(dict(vertex=v,forest=forest,finite_frame=target,old_nodes=nodes,old_groups=groups,
        old_graph_edges=sorted(tuple(sorted(e))for e in graph.edges()),prior_history=history,marker_families=newf,counts=counts))


def envelope(row,old76,old56,old108):
    forest=row['forest'];co=Counter(forest);a,h,h4,t,s=co['A1'],co['H3'],co['H4'],co['I2(10)'],co['I2(5)']
    fam={f['kind']:f['markers']for f in row['marker_families']}
    if old76:fam[76]=old76['markers']
    if old56:fam.update({f['raw_family']:f['markers']for f in old56['marker_families']})
    if old108:fam.update({f['kind']:f['markers']for f in old108['marker_families']})
    t=t if 56 in fam else 0;c57=co['I2(10)']if 57 in fam else 0
    c76=min(a-1,h)if 76 in fam else 0;c108A=2*h4*s if 1084 in fam else 0;c108B=2*s if 1085 in fam else 0
    x=h4 if 1124 in fam else 0;y=int(1125 in fam);boundary=[];recipes=[]
    for k in range(c76+1):
        used=k+1 if k else 0;g=(a-used)//2
        for m in range(g+1)if x else(0,):
            q=(g-m)*t;counts=(q,0,k,c108A,0,m*x,0)
            if any(counts):boundary.append(counts);recipes.append(dict(mode='old76',k=k,m=m,counts=counts))
    if c108B or y:
        g=(a-1)//2
        for m in range(g+1)if x or y else(0,):
            counts=((g-m)*t,0,0,c108A,c108B,m*x,m*y)
            if any(counts):boundary.append(counts);recipes.append(dict(mode='centered',k=0,m=m,counts=counts))
    if c57:
        counts=(0,c57,0,0,0,0,0);boundary.append(counts);recipes.append(dict(mode='raw57',k=0,m=0,counts=counts))
    boundary=sorted(set(boundary));assert boundary
    active=[i for i,k in enumerate(ORDER)if k in fam and k!=57];dim=len(active)
    constraints=[(tuple(F(p[i])for i in active),F(1))for p in boundary if any(p[i]for i in active)]
    constraints +=[(tuple(F(int(i==j))for i in range(dim)),F())for j in range(dim)]
    vertices=set()
    for cc in combinations(constraints,dim):
        z=N.solve([p for p,b in cc],[b for p,b in cc])
        if z is not None and min(z)>=0 and all(sum(F(p[i])*z[j]for j,i in enumerate(active))<=1 for p in boundary):vertices.add(z)
    maximal=sorted(z for z in vertices if not any(z!=w and all(x<=y for x,y in zip(z,w))for w in vertices))
    weights=[]
    for z in maximal:
        w=[F()]*7
        for i,coefficient in zip(active,z):w[i]=coefficient
        if c57:w[1]=F(1,c57)
        assert all(sum(n*x for n,x in zip(p,w))<=1 for p in boundary);weights.append(w)
    # Every boundary row has a literal clique in the conservative role graph.
    target={int(i):r for i,r in row['finite_frame'].items()};ones=sorted(i for i,r in target.items()if not r)
    witnesses=[];seen=set()
    for recipe in recipes:
        counts=recipe['counts']
        if counts in seen:continue
        seen.add(counts);mode=recipe['mode'];k=recipe['k'];m=recipe['m'];chosen=[]
        if mode=='raw57':
            anchor=role(57,fam[57][0],'A2')
            chosen=[(57,z['index'])for z in fam[57]if role(57,z,'A2')==anchor and role(57,z,'I2(10)')[0]<role(57,z,'I2(10)')[1]]
        else:
            used=set();central=None;block=None
            if k:
                hs=sorted({role(76,z,'H3')for z in fam[76]})
                for j in range(k):
                    z=next(z for z in fam[76]if role(76,z,'A1')==(ones[0],ones[j+1])and role(76,z,'H3')==hs[j])
                    chosen.append((76,z['index']))
                used.update(ones[:k+1])
            if mode=='centered':
                central=ones[0];used.add(central)
                kk=1125 if y else 1085;block=role(kk,fam[kk][0],'H3')
                if c108B:chosen.extend((1085,z['index'])for z in fam[1085]if role(1085,z,'A1')==(central,)and role(1085,z,'H3')==block)
            avail=[i for i in ones if i not in used];pairs=list(zip(avail[::2],avail[1::2]));own=set(pairs[:m]);other=set(pairs[m:])
            if x:chosen.extend((1124,z['index'])for z in fam[1124]if role(1124,z,'A1')in own)
            if mode=='centered'and y:chosen.extend((1125,z['index'])for z in fam[1125]if role(1125,z,'A1')in{(p,central,q)for p,q in own}and role(1125,z,'H3')==block)
            if t:chosen.extend((56,z['index'])for z in fam[56]if role(56,z,'A1')in other and role(56,z,'I2(10)')[0]<role(56,z,'I2(10)')[1])
            if c108A:chosen.extend((1084,z['index'])for z in fam[1084])
        assert Counter(kind for kind,i in chosen)==Counter({kind:z for kind,z in zip(ORDER,counts)if z}),(forest,recipe,Counter(kind for kind,i in chosen))
        assert all(compatible(ka,fam[ka][i],kb,fam[kb][j])for(ka,i),(kb,j)in combinations(chosen,2))
        witnesses.append(dict(counts=counts,marker_clique=chosen))
    assert set(tuple(w['counts'])for w in witnesses)==set(boundary)
    return dict(kind_order=list(ORDER),boundary_counts=boundary,coefficient_vertices=weights,attaining_marker_cliques=witnesses)


def export(rows,old56,old76,joint,old108,model,domains):
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']};baseline={v:F(x)for v,x in old108['exported']['repair']['vertex_multipliers']}
    assert dual['face_multipliers']==old108['exported']['repair']['face_multipliers']
    newmu=dict(baseline);oldjoint={r['vertex']:r for r in joint['forest_records']}
    oldquota={r['vertex']:r for r in old108['exported']['target_allowances']}
    older={r['vertex']:r for r in old76['forest_records']};oldest={r['vertex']:r for r in old56['forest_records']}
    previous={r['vertex']:r for r in old108['forest_records']}
    candidates=defaultdict(list);vectors={v:d['domain']['maximal_count_vectors']for v,d in domains.items()}
    for tag in(D.L.TAG,D.B.TAG,D.THIRD):
        for r in M.read(tag)['forest_records']:
            candidates[r['vertex']].append((tag,r['certificate']['weights']))
            if tag!=D.L.TAG:vectors[r['vertex']]=r['maximal_count_vectors']
    columns=[];quotas=[]
    for row in rows:
        v=row['vertex'];groups=row['old_groups'];owner={i:z['index']for z in groups for i in z['literal_stars']}
        omitted={owner[i]for f in row['marker_families']for i,k in f['marker_ordinary_edges']}
        if v in oldquota:omitted.update(oldquota[v]['omitted_orbits'])
        elif v in oldjoint:omitted.update(oldjoint[v]['ordinary_omitted_orbits'])
        omitted=sorted(omitted);faces={f['kind']:f['face']for f in row['marker_families']}
        if v in older:faces[76]=older[v]['face']
        if v in oldest:faces.update({f['raw_family']:f['face']for f in oldest[v]['marker_families']})
        if v in previous:faces.update({f['kind']:f['face']for f in previous[v]['marker_families']})
        env=envelope(row,older.get(v),oldest.get(v),previous.get(v))
        cc=candidates[v]if groups else[('marker_only',[])];vv=vectors[v]if groups else[[]];assert cc
        unique={}
        for tag,ww in cc:
            original=list(map(F,ww))+[F()]*(len(groups)-len(ww));assert len(original)==len(groups)and all(w>=0 for w in original)
            assert all(len(p)==len(original)and sum(w*n for w,n in zip(original,p))<=1 for p in vv)
            w=[F()if i in omitted else x for i,x in enumerate(original)]
            for mw in env['coefficient_vertices']:
                terms=[[kind,faces[kind],str(weight)]for kind,weight in zip(ORDER,mw)if weight];cf=Counter()
                for weight,z in zip(w,groups):
                    for face,n in z['positive_face_coefficients']:cf[face]+=weight*n
                for kind,face,weight in terms:cf[face]+=F(weight)
                key=tuple((f,x)for f,x in sorted(cf.items())if x);load=sum((alpha.get(f,F())*x for f,x in key),F())
                column=dict(vertex=v,source_charge=tag,other_charge_sources=[],old_orbit_weights=list(map(str,original)),
                    omitted_orbits=omitted,retained_orbit_weights=list(map(str,w)),marker_terms=terms,
                    positive_face_coefficients=[[f,str(x)]for f,x in key],vertex_cost='1',load=str(load),
                    baseline_cost=str(baseline.get(v,F())),violation=str(load-baseline.get(v,F())),
                    old_vector_count=len(vv),old_vectors_sha256=D.L.jhash(vv))
                if key in unique:unique[key]['other_charge_sources'].append(tag)
                else:unique[key]=column
        cs=list(unique.values());columns.extend(cs);newmu[v]=max([baseline.get(v,F())]+[F(c['load'])for c in cs])
        quotas.append(dict(vertex=v,baseline_cost=str(baseline.get(v,F())),new_cost=str(newmu[v]),increase=str(newmu[v]-baseline.get(v,F())),
            omitted_orbits=omitted,envelope=env))
    preserved=[]
    for tag,data in((U.TAG,old76),(T.TAG,old56),(N.J3.TAG,joint),(N.TAG,old108)):
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
    assert len(preserved)==3299 and len(multi)==14 and all(newmu.get(v,F())>=x for v,x in baseline.items())
    num=F(dual['objective_numerator']);den=sum(newmu.values(),F());gap=num/12-den
    return dict(columns=columns,target_allowances=quotas,preserved_columns=preserved,preserved_multi_owner_costs=multi,
        repair=dict(face_multipliers=dual['face_multipliers'],vertex_multipliers=[[v,str(x)]for v,x in sorted(newmu.items())if x],
            objective_numerator=str(num),normalizer=str(den),bound=str(num/den),baseline_bound=old108['exported']['repair']['bound'],
            remaining_normalizer_gap_to12=str(gap),scope='Least monotone fixed-alpha repair for exported and prior columns, not enlarged-cone optimality or a global optimization.'),
        counts=dict(exported_columns=len(columns),strictly_violated_columns=sum(F(c['violation'])>0 for c in columns),raised_allowances=sum(F(q['increase'])>0 for q in quotas)))


def run():
    start=time.monotonic();inp=D.Q.Inputs();geo=inp.read(C.TAG);old56=inp.read(T.TAG);old76=inp.read(U.TAG)
    joint=inp.read(N.J3.TAG);old108=inp.read(N.TAG);model=inp.read(T.MODEL)
    for d in(geo,old56,old76,joint,old108):
        for p,h in d['sources'].items():assert M.digest(p)==h,p
    gg=next(r for r in geo['geometries']if r['raw_family']==112);g=[[tuple(e)for e in row]for row in gg['native_gram']]
    assert gg['required_quads']==3 and gg['native_automorphisms']==[list(range(8))]
    assert[r['omitted']for r in gg['extra_vertices']]==[[2,4],[2,5]]
    tests=relative_tests(g,old56,old76,old108)
    for t in tests:
        for z in t['records']:assert z['result']['allowed']==compatible(t['kinds'][0],t['anchor'],t['kinds'][1],z['other_marker'])
    family=[r for r in geo['complete_family_inventory']if r['raw_family']==112]
    rows=[dict(r,kind=1124 if r['omitted']==[2,4]else 1125)for r in geo['extra_marker_inventory']if r['raw_family']==112]
    assert len(family)==196 and len(rows)==392 and sum(r['active_role_orientations']for r in rows)==11012
    assert all(r['actual_extra_occurrences']==1 for r in rows)
    for kind in KINDS:assert len({r['target_vertex']for r in rows if r['kind']==kind})==196
    domains,represented=D.L.current_domains();p=D.Z.prepare();local={r['vertex']:r for r in inp.read(D.L.TAG)['forest_records']}
    patches={r['vertex']:r for r in inp.read(D.Z.TAG)['forest_records']};grouped=defaultdict(list);oldby=defaultdict(list)
    for row in rows:grouped[row['target_vertex']].append(row)
    for data in(old56,old76,old108):
        for row in data['forest_records']:oldby[row['vertex']].append(row)
    print('PREFLIGHT',len(grouped),len(rows),sum(r['active_role_orientations']for r in rows),flush=True)
    records=[]
    for i,(v,rr)in enumerate(sorted(grouped.items()),1):
        record=process(rr,g,p,domains,local,patches,oldby);records.append(record)
        if i%25==0 or record['counts']['anchor_unknowns']:print('FOREST',i,'/',len(grouped),v,record['counts'],'seconds',round(time.monotonic()-start,2),flush=True)
    exported=export(records,old56,old76,joint,old108,model,domains)
    paths=set(inp.paths)|set(p['sources'])
    for mod in(C,N,N.J3,T,U,D,D.L,R,A,H,M):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    for tag in(C.TAG,N.TAG,N.J3.TAG,T.TAG,U.TAG,D.L.TAG,D.B.TAG,D.THIRD,T.MODEL):paths.update(f'{folder}/{tag}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    paths.add(str(Path(__file__).relative_to(ROOT)))
    counts=dict(classes=196,pointed_marker_domains=392,target_forests=len(records),relative_cases=162,
        retained_relative_cases=sum(z['result']['allowed']for t in tests for z in t['records']),
        **{key:sum(row['counts'][key]for row in records)for key in records[0]['counts']},**exported['counts'])
    return normal(M.enc(dict(status='PASS_RAW112_COMPLETE_EXTRA_VERTEX_MARKER_FAMILY',sources={p:M.digest(p)for p in sorted(paths)},
        geometry=gg,family=family,pointed_family=rows,relative_tests=tests,forest_records=records,exported=exported,counts=counts,
        optimizer_invocations=0,old_outcomes_changed=False,scope='Complete raw112 family, once per aggregate at each extra vertex. Ordinary and five prior marker types share one allowance; unknowns retained. Exported fixed-mass repair only.')))


def pilot():
    geo=M.read(C.TAG);old56=M.read(T.TAG);old76=M.read(U.TAG);old108=M.read(N.TAG)
    for d in(geo,old56,old76,old108):
        for p,h in d['sources'].items():assert M.digest(p)==h,p
    gg=next(r for r in geo['geometries']if r['raw_family']==112)
    assert gg['required_quads']==3 and gg['native_automorphisms']==[list(range(8))]
    assert[r['omitted']for r in gg['extra_vertices']]==[[2,4],[2,5]]
    g=[[tuple(e)for e in row]for row in gg['native_gram']]
    tests=relative_tests(g,old56,old76,old108)
    family=[r for r in geo['complete_family_inventory']if r['raw_family']==112]
    pointed=[dict(r,kind=1124 if r['omitted']==[2,4]else 1125)for r in geo['extra_marker_inventory']if r['raw_family']==112]
    assert len(family)==196 and len(pointed)==392 and sum(r['active_role_orientations']for r in pointed)==11012
    assert all(r['actual_extra_occurrences']==1 for r in pointed)
    paths={str(Path(__file__).relative_to(ROOT))}
    for tag in(C.TAG,T.TAG,U.TAG,N.TAG):paths.update(f'{folder}/{tag}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return normal(dict(status='PASS_RAW112_COMPLETE_RELATIVE_PILOT',sources={p:M.digest(p)for p in sorted(paths)},
        geometry=gg,family=family,pointed_family=pointed,relative_tests=tests,
        counts=dict(classes=196,pointed_domains=392,role_orientations=11012,relative_cases=162,
            retained_cases=sum(r['result']['allowed']for t in tests for r in t['records'])),
        optimizer_invocations=0,old_outcomes_changed=False,
        scope='Complete relative native pair domain only. Shared capacity, ordinary transport, columns and fixed-mass repair are not yet asserted.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--pilot',action='store_true');args=ap.parse_args()
    out=pilot()if args.pilot else run();path=ROOT/f'results/{TAG}{"_pilot"if args.pilot else""}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
