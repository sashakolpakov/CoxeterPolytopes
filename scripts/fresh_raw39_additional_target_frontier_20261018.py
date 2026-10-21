#!/usr/bin/env python3
"""Explore the four previously unused outside vertices of the raw39 pair.

Exact necessary gates only; conservative compatibility is not realization.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import time
import networkx as nx
import fresh_raw39_shared_edge_target_incidence_20261018 as B

P,A,R,N,W,V,K,M,ROOT=B.P,B.A,B.R,B.N,B.W,B.V,B.K,B.M,B.ROOT
TAG='fresh_raw39_additional_target_frontier_20261018'
TARGETS=(1880,2622)


def normal(x):
    return json.loads(json.dumps(x,sort_keys=True))


def covkey(cov):
    return tuple(sorted(cov.values()))


def necessary_gate(covs,heights,other_covs,other_heights,J):
    for i,c in enumerate(covs):
        for j,d in enumerate(other_covs):
            if c==d:
                assert heights[i]==other_heights[j]
                continue
            field=N.pair_field((c,d))
            pr=R.projection(c,d,J,field);hp=heights[i]*other_heights[j]
            rec=dict(first=i,second=j,projection=R.qout(pr),height_product=hp.out())
            if field.signs_difference(pr,hp)>0:
                return dict(**rec,kind='positive_actual_product')
            if not any(pr) and not M.integral(hp):
                return dict(**rec,kind='nonintegral_two_cycle',trace=str(2*hp.a),
                            norm=str(hp.a*hp.a-5*hp.b*hp.b))
    return None


def event_domain(vertex,target,native,geo):
    records=[];by_support=[]
    supports=[r['support'] for r in geo['outside_both'] if r['vertex']==vertex]
    assert len(supports)==2
    for ids in supports:
        frame={i:{j:native[i,j][0] for j in ids if j!=i and native[i,j][0]!=2} for i in ids}
        embeddings=R.active_embeddings(frame,target);keys=set()
        for mp in embeddings:
            chosen=set(mp.values());passive=set(target)-chosen
            assert len(chosen)==5 and len(passive)==8
            assert all(not(set(target[i])&passive) for i in chosen)
            models={R.C.code(g):name for name,g in R.C.models()}
            rest={i:{j:m for j,m in target[i].items() if j in passive} for i in passive}
            assert sorted(models[R.C.code(g)] for g in R.C.pieces(rest))==['A5','H3']
            cov,hs,rec=B.recover(ids,mp,native,target)
            key=covkey(cov);keys.add(key)
            # An outside vertex lies on exactly one of the two forced compact
            # four-simplices. Those four ports therefore have triangular faces.
            simplices=[s for s in geo['actual_simplices'] if set(s['fixed'])<=set(ids)
                       and len(set(ids)&set(s['sides']))==4]
            assert len(simplices)==1
            ports=sorted(mp[i] for i in set(ids)&set(simplices[0]['sides']))
            records.append(dict(**rec,simplex_ports=ports,event_key=normal(key)))
        by_support.append(keys)
    assert by_support[0]==by_support[1]
    events=[]
    for key in sorted(by_support[0]):
        same=[r for r in records if normal(key)==r['event_key']]
        assert len(same)==2 and same[0]['simplex_ports']==same[1]['simplex_ports']
        # Their two presentations exchange the original unlabelled Q pair.
        one,two=same
        assert dict(two['mapping'])=={B.REFLECTION[i]:j for i,j in one['mapping']}
        heights={tuple((j,tuple(e)) for j,e in c):K(*map(F,h))
                 for c,h in zip(one['contacts'],one['heights'])}
        assert set(heights)==set(key)
        events.append(dict(key=key,heights=[heights[c] for c in key],
                           simplex_ports=one['simplex_ports'],presentations=same))
    return events,records


def target_record(vertex,p,changed,rows,patches,native,geo,base):
    old,g,owner,unowned,domain,groups,forest,applied=A.current(vertex,p,changed,rows,patches)
    nodes=A.T.current(vertex,p,changed)[5]
    target=R.P.forest_graph(forest,dict(R.C.models()))
    J,active,components=R.projection_data(target)
    inherited_edges=sorted(tuple(sorted(e)) for e in g.edges())
    events,maps=event_domain(vertex,target,native,geo)
    gates=[];markers=[]
    for ei,event in enumerate(events):
        marker=len(nodes)+ei;markers.append(marker);g.add_node(marker)
        for si,star in enumerate(nodes):
            port=next((pair for pair in star['occupied_pairs']
                       if set(pair)<=set(event['simplex_ports'])),None)
            if port is not None:gate=dict(kind='actual_simplex_triangle_port',ports=port)
            else:
                hs=[K(*map(F,h)) for h in star['squared_timelike_heights']]
                gate=necessary_gate(event['key'],event['heights'],N.native(star),hs,J)
            if gate is None:g.add_edge(marker,si)
            gates.append(dict(event=ei,star=si,gate=gate))
    event_gates=[]
    for i,j in combinations(range(len(events)),2):
        e,f=events[i],events[j]
        gate=necessary_gate(e['key'],e['heights'],f['key'],f['heights'],J)
        if gate is None:g.add_edge(markers[i],markers[j])
        event_gates.append(dict(events=[i,j],gate=gate))
    assert sorted(tuple(sorted(e)) for e in g.edges() if all(x<len(nodes) for x in e))==inherited_edges
    cliques=[];items={};start=time.monotonic();maxevent=0
    for cc in nx.find_cliques(g):
        cl=tuple(sorted(cc));cliques.append(cl)
        assert len(cliques)<100000 and time.monotonic()-start<45
        for k in range(1,len(cl)+1):
            for sub in combinations(cl,k):
                ct=Counter(owner[i] for i in sub if i<len(nodes))
                vec=tuple(ct.get(j,0) for j in range(len(groups)))
                ev=sum(i>=len(nodes) for i in sub);maxevent=max(maxevent,ev)
                key=(vec,ev)
                if key not in items or sub<items[key]:items[key]=sub
    demands=[]
    for field in ('face_multipliers','restricted_face_multipliers'):
        alpha={i:F(x) for i,x in base[field]}
        demands.append([str(sum((alpha.get(f,F())*m for f,m in gr['positive_face_coefficients']),F()))
                        for gr in groups])
    record=dict(vertex=vertex,forest=forest,finite_frame=[[i,sorted(row.items())] for i,row in sorted(target.items())],
        corner_orbits=groups,literal_stars=nodes,physical_maps=maps,
        event_records=[dict(key=normal(e['key']),heights=[h.out() for h in e['heights']],
                            simplex_ports=e['simplex_ports']) for e in events],
        event_star_gates=gates,event_event_gates=event_gates,
        inherited_graph_edges=inherited_edges,graph_edges=sorted(tuple(sorted(e)) for e in g.edges()),
        inherited_unowned_pairs=normal(sorted(unowned)),markers=markers,
        maximal_cliques=sorted(cliques),
        items=[dict(vector=list(v),events=e,clique=list(items[v,e])) for v,e in sorted(items)],
        event_capacity=maxevent,demands=demands)
    print('TARGET',vertex,'stars',len(nodes),'events',len(events),'capacity',maxevent,
          'gates',Counter('compatible' if r['gate'] is None else r['gate']['kind'] for r in gates),
          'cliques',len(cliques),'vectors',len(items),flush=True)
    return record


def check_joint(records,base,di,c):
    mu={v:F(x) for v,x in base['repairs'][di]['vertex_multipliers']}
    weights=[list(map(F,w)) for w in c['weights']]
    aa=list(map(F,c['source_costs']));zz=list(map(F,c['event_coefficients']))
    inc=list(map(F,c['increases']));X=F(c['event_mass'])
    assert len(weights)==len(aa)==len(zz)==len(inc)==len(c['mixtures'])==len(records)
    assert all(0<=x<=1 for x in aa) and sum(zz,F())<=0 and X>=0
    assert all(x>=0 for x in inc)
    for r,w,a,z,d,mix in zip(records,weights,aa,zz,inc,c['mixtures']):
        demand=list(map(F,r['demands'][di]))
        assert len(w)==len(demand) and all(x>=0 for x in w)
        assert all(sum(x*y for x,y in zip(w,it['vector']))<=a+z*it['events'] for it in r['items'])
        terms=[(i,F(v)) for i,v in mix]
        assert len({i for i,v in terms})==len(terms)
        assert all(0<=i<len(r['items']) and v>0 for i,v in terms)
        assert all(sum((v*r['items'][i]['vector'][j] for i,v in terms),F())>=demand[j]
                   for j in range(len(demand)))
        assert sum((v for i,v in terms),F())<=mu.get(r['vertex'],F())+d
        assert sum((v*r['items'][i]['events'] for i,v in terms),F())==X
    load=sum((sum((w*F(d) for w,d in zip(ww,r['demands'][di])),F())
              for ww,r in zip(weights,records)),F())
    cost=sum((a*mu.get(r['vertex'],F()) for a,r in zip(aa,records)),F())
    gain=sum(inc,F());assert load-cost==gain==F(c['gain'])
    num=F(base['repairs'][di]['objective_numerator']);den=F(base['repairs'][di]['normalizer'])+gain
    assert F(c['bound'])==num/den


def joint(records,base,di,saved=None):
    if saved is not None:
        check_joint(records,base,di,saved);return saved
    import numpy as np
    from scipy.optimize import linprog
    mu={v:F(x) for v,x in base['repairs'][di]['vertex_multipliers']}
    starts=[];n=0
    for r in records:starts.append(n);n+=len(r['items'])
    count=len(records);deltas=list(range(n,n+count));event=n+count;nv=event+1
    objective=[0]*nv
    for i in deltas:objective[i]=1
    ub=[];rhs=[];labels=[];eq=[]
    for f,r in enumerate(records):
        for i,d in enumerate(r['demands'][di]):
            row=[0]*nv
            for j,it in enumerate(r['items']):row[starts[f]+j]=-it['vector'][i]
            ub.append(row);rhs.append(-F(d));labels.append(('demand',f,i))
        row=[0]*nv
        for j in range(len(r['items'])):row[starts[f]+j]=1
        row[deltas[f]]=-1;ub.append(row);rhs.append(mu.get(r['vertex'],F()));labels.append(('cost',f))
        row=[0]*nv
        for j,it in enumerate(r['items']):row[starts[f]+j]=it['events']
        row[event]=-1;eq.append(row)
    scale=max(map(abs,rhs))
    res=linprog(objective,A_ub=np.array(ub),b_ub=np.array([float(x/scale) for x in rhs]),
                A_eq=np.array(eq),b_eq=np.zeros(count),bounds=(0,None),method='highs',
                options={'time_limit':20,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    assert res.success,res.message
    ids=[j for j,x in enumerate(res.x) if x>1e-12]
    tight=[i for i,x in enumerate(res.ineqlin.residual) if abs(x)<1e-11]
    vals=W.solve_exact([[r[j] for j in ids] for r in eq]+
                       [[ub[i][j] for j in ids] for i in tight],
                       [F()]*count+[rhs[i] for i in tight],len(ids))
    xx=[F()]*nv
    for j,x in zip(ids,vals):xx[j]=x
    weights=[[F()]*len(r['demands'][di]) for r in records];aa=[F()]*count
    for label,x in zip(labels,res.ineqlin.marginals):
        value=F(float(-x)).limit_denominator(1000000)
        if label[0]=='demand':weights[label[1]][label[2]]=value
        else:aa[label[1]]=value
    zz=[F(float(-x)).limit_denominator(1000000) for x in res.eqlin.marginals]
    gain=sum((xx[j] for j in deltas),F())
    num=F(base['repairs'][di]['objective_numerator']);den=F(base['repairs'][di]['normalizer'])+gain
    out=dict(weights=[list(map(str,w)) for w in weights],source_costs=list(map(str,aa)),
        event_coefficients=list(map(str,zz)),event_mass=str(xx[event]),
        increases=[str(xx[j]) for j in deltas],gain=str(gain),bound=str(num/den),
        mixtures=[[[j,str(xx[starts[f]+j])] for j in range(len(r['items'])) if xx[starts[f]+j]]
                  for f,r in enumerate(records)])
    check_joint(records,base,di,out)
    print('JOINT',di,'gain',out['gain'],'X',out['event_mass'],'increases',out['increases'],
          'bound',float(F(out['bound'])),flush=True)
    return out


def run(verify=False):
    started=time.monotonic();base=W.read(B.TAG)
    saved=W.read(TAG) if verify else None
    for path,h in base['sources'].items():assert W.digest(path)==h
    p,changed,rows,patches,_=A.prepare()
    G,native,geo=B.geometry()
    assert normal(geo)==base['geometry']
    graph=nx.Graph();graph.add_nodes_from(range(8))
    for i,j in combinations(range(8),2):
        if native[i,j][0]!=2:graph.add_edge(i,j,label=tuple(native[i,j]))
    autos=list(nx.isomorphism.GraphMatcher(graph,graph,
               edge_match=lambda a,b:a['label']==b['label']).isomorphisms_iter())
    assert len(autos)==2 and {tuple(m[i] for i in range(8)) for m in autos}=={
        tuple(range(8)),tuple(B.REFLECTION[i] for i in range(8))}
    records=[target_record(v,p,changed,rows,patches,native,geo,base) for v in TARGETS]
    assert [len(r['physical_maps']) for r in records]==[8,8]
    assert [len(r['event_records']) for r in records]==[4,4]
    assert all(r['event_capacity']==1 for r in records)
    combined=base['forest_records']+records
    certs=[joint(combined,base,i,saved['joint_certificates'][i] if saved else None) for i in range(2)]
    assert [c['gain'] for c in certs]==['934391/2000000000','456567/2500000000']
    columns=[];repairs=[]
    for di,c in enumerate(certs):
        coefficients=Counter()
        for r,weights in zip(combined,c['weights']):
            for group,w in zip(r['corner_orbits'],map(F,weights)):
                for f,m in group['positive_face_coefficients']:coefficients[f]+=w*m
        col=dict(positive_face_coefficients=[[f,str(x)] for f,x in sorted(coefficients.items()) if x],
                 vertex_capacities=[[r['vertex'],a] for r,a in zip(combined,c['source_costs']) if F(a)])
        if col not in columns:columns.append(col)
        prior=base['repairs'][di];mu={v:F(x) for v,x in prior['vertex_multipliers']};new=dict(mu)
        for r,x in zip(combined,c['increases']):new[r['vertex']]=new.get(r['vertex'],F())+F(x)
        den=sum(new.values(),F());num=F(prior['objective_numerator'])
        assert num/den==F(c['bound'])
        repairs.append(dict(domain=prior['domain'],objective_numerator=str(num),normalizer=str(den),
            bound=str(num/den),vertex_updates=[[v,str(mu.get(v,F())),str(x)] for v,x in sorted(new.items()) if x>mu.get(v,F())],
            vertex_multipliers=[[v,str(x)] for v,x in sorted(new.items()) if x]))
    for di,r in enumerate(repairs):
        mu={v:F(x) for v,x in r['vertex_multipliers']}
        alpha={f:F(x) for f,x in base[('face_multipliers','restricted_face_multipliers')[di]]}
        for col in columns:
            assert sum((F(x)*alpha.get(f,F()) for f,x in col['positive_face_coefficients']),F())<=sum(
                (F(x)*mu.get(v,F()) for v,x in col['vertex_capacities']),F())
    paths=[Path(__file__),Path(B.__file__),ROOT/f'results/{B.TAG}.json',ROOT/f'results/{B.TAG}.md']
    paths += [Path(m.__file__) for m in (A,R,N,W)]
    paths += [ROOT/'scripts/check_raw39_shared_edge_target_saved_20261018.py']
    output=normal(dict(status='PASS_RAW39_ADDITIONAL_TARGET_INCIDENCE',
        sources={str(q.relative_to(ROOT)):sha256(q.read_bytes()).hexdigest() for q in paths},
        target_records=records,joint_record_vertices=[r['vertex'] for r in combined],
        joint_certificates=certs,source_occurrences_per_pair=2,
        native_configuration_automorphisms=[sorted(a.items()) for a in autos],
        columns=columns,repairs=repairs,face_multipliers=base['face_multipliers'],
        restricted_face_multipliers=base['restricted_face_multipliers'],
        target_occurrences_per_pair={str(v):2 for v in TARGETS},
        local_joint_proposals=2,global_optimizer_invocations=0,
        scope='Complete frame maps and necessary compatibility gates at two additional raw39 target types. The written recovery proof identifies each unlabelled source pair at its two targets of each type. Exact four-vertex coupled local covers/charges give two new positive columns and monotone fixed-face-mass repairs. Conservative compatible records do not assert physical realization. No global optimization, new dimension bound or example.'))
    path=ROOT/f'results/{TAG}.json'
    if verify:assert output==json.loads(path.read_text())
    else:path.write_text(json.dumps(output,separators=(',',':'))+'\n')
    print(output['status'],'seconds',round(time.monotonic()-started,2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    run(args.verify)
