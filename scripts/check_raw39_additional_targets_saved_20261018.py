#!/usr/bin/env python3
"""Independent saved-map, gate, graph and four-vertex certificate checks.

Uses the older independent checker's Q(sqrt5) arithmetic, never a producer or
primary project arithmetic module. Positive radical comparisons are enclosed
by exact rational intervals, not floating point.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json
from math import isqrt
from pathlib import Path
import check_raw39_shared_edge_target_saved_20261018 as S

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw39_additional_target_frontier_20261018'
Q=S.Q


def freeze(x):
    return tuple(map(freeze,x)) if isinstance(x,list) else x


def sqrt_interval(x,bits):
    assert x>=0
    scale=1<<bits;k=isqrt(x.numerator*scale*scale//x.denominator)
    lo,hi=F(k,scale),F(k+1,scale)
    assert lo*lo<=x<hi*hi
    return lo,hi


def q_interval(x,bits):
    lo,hi=sqrt_interval(F(5),bits)
    return (x.a+x.b*lo,x.a+x.b*hi) if x.b>=0 else (x.a+x.b*hi,x.a+x.b*lo)


def positive_gate(c,d,J,gate):
    pr,ps,hp=S.projections(c,d,J)
    assert hp==Q(*gate['height_product'])
    basis=[]
    for i,e in c+d:
        t=S.squared(e)
        if t!=0 and S.root(t) is None and all(S.root(t/u) is None for u in basis):basis.append(t)
    assert len(pr)==1<<len(basis)
    assert all(x.sign()>=0 for x in pr) and any(x.sign()>0 for x in pr)
    powers=[]
    for mask in range(len(pr)):
        p=Q(1)
        for j,b in enumerate(basis):
            if mask>>j&1:p=p*b
        powers.append(p)
    for bits in (16,32,64,128,256):
        lower=F();usable=True
        for x,p in zip(pr,powers):
            xl,xh=q_interval(x,bits);pl,ph=q_interval(p,bits)
            if xl<0 or pl<0:usable=False;break
            lower+=xl*sqrt_interval(pl,bits)[0]
        if usable and lower*lower>q_interval(hp,bits)[1]:return bits
    raise AssertionError('positive-product gate not separated by exact intervals')


def complete_embeddings(ids,native,target):
    # Whole active components are enforced by matching the full incident
    # native-label multiset at each image, then every pair product.
    order=sorted(ids,key=lambda i:-sum(native[i,j][0]!=2 for j in ids if i!=j))
    sig={i:sorted(native[i,j] for j in ids if i!=j and native[i,j][0]!=2) for i in ids}
    tsig={i:sorted(e for j,e in row.items() if e[0]!=2) for i,row in target.items()}
    candidates={i:[j for j in target if sig[i]==tsig[j]] for i in ids}
    out=[]
    def visit(mp):
        if len(mp)==len(order):out.append(tuple(sorted(mp.items())));return
        i=order[len(mp)]
        for j in candidates[i]:
            if j in mp.values():continue
            if all(native[i,k]==target[j].get(v,(2,'0','0')) for k,v in mp.items()):
                visit(mp|{i:j})
    visit({})
    return sorted(out)


def verify_joint(records,base,data,di,counts):
    field=('face_multipliers','restricted_face_multipliers')[di]
    assert data[field]==base[field]
    alpha={i:F(x) for i,x in data[field]};prior=base['repairs'][di]
    mu={v:F(x) for v,x in prior['vertex_multipliers']}
    c=data['joint_certificates'][di];weights=[list(map(F,w)) for w in c['weights']]
    aa=list(map(F,c['source_costs']));zz=list(map(F,c['event_coefficients']))
    delta=list(map(F,c['increases']));X=F(c['event_mass'])
    assert len(weights)==len(aa)==len(zz)==len(delta)==len(c['mixtures'])==len(records)==4
    assert all(0<=x<=1 for x in aa) and all(x>=0 for x in delta) and X>=0 and sum(zz,F())<=0
    new=dict(mu);load=F();coefficient=Counter()
    for r,w,a,z,inc,mix in zip(records,weights,aa,zz,delta,c['mixtures']):
        demand=[sum((alpha.get(f,F())*m for f,m in g['positive_face_coefficients']),F()) for g in r['corner_orbits']]
        assert list(map(str,demand))==r['demands'][di] and len(w)==len(demand) and all(x>=0 for x in w)
        assert all(sum(x*y for x,y in zip(w,it['vector']))<=a+z*it['events'] for it in r['items'])
        terms=[(j,F(x)) for j,x in mix]
        assert len({j for j,x in terms})==len(terms) and all(0<=j<len(r['items']) and x>0 for j,x in terms)
        assert all(sum((x*r['items'][j]['vector'][i] for j,x in terms),F())>=demand[i] for i in range(len(demand)))
        assert sum((x*r['items'][j]['events'] for j,x in terms),F())==X
        assert sum((x for j,x in terms),F())<=mu.get(r['vertex'],F())+inc
        load+=sum((x*y for x,y in zip(w,demand)),F());new[r['vertex']]=new.get(r['vertex'],F())+inc
        for weight,g in zip(w,r['corner_orbits']):
            for f,m in g['positive_face_coefficients']:coefficient[f]+=weight*m
        counts.update(charge_inequalities=len(r['items']),mixture_terms=len(terms))
    gain=load-sum((a*mu.get(r['vertex'],F()) for a,r in zip(aa,records)),F())
    assert gain==sum(delta,F())==F(c['gain'])
    col=dict(positive_face_coefficients=[[f,str(x)] for f,x in sorted(coefficient.items()) if x],
             vertex_capacities=[[r['vertex'],str(a)] for r,a in zip(records,aa) if a])
    repair=data['repairs'][di];den=sum(new.values(),F());num=F(prior['objective_numerator'])
    assert repair['vertex_multipliers']==[[v,str(x)] for v,x in sorted(new.items()) if x]
    assert repair['vertex_updates']==[[v,str(mu.get(v,F())),str(x)] for v,x in sorted(new.items()) if x>mu.get(v,F())]
    assert F(repair['normalizer'])==den and F(repair['objective_numerator'])==num
    assert F(repair['bound'])==F(c['bound'])==num/den
    for column in data['columns']:
        assert sum((F(x)*alpha.get(f,F()) for f,x in column['positive_face_coefficients']),F())<=sum(
            (F(x)*new.get(v,F()) for v,x in column['vertex_capacities']),F())
    return col


def run():
    data=json.loads((ROOT/f'results/{TAG}.json').read_text())
    assert data['status']=='PASS_RAW39_ADDITIONAL_TARGET_INCIDENCE'
    for name,h in data['sources'].items():assert S.digest(ROOT/name)==h,name
    base=json.loads((ROOT/f'results/{S.TAG}.json').read_text())
    geo=base['geometry'];G=S.mat(geo['scaled_gram']);scales=[Q(*x) for x in geo['squared_scales']]
    native={(i,j):(2,'0','0') for i in range(8) for j in range(8) if i!=j}
    for i,j,e in geo['native_edges']:native[i,j]=native[j,i]=tuple(e)
    graph={i:{j:native[i,j] for j in range(8) if j!=i and native[i,j][0]!=2} for i in range(8)}
    autos=complete_embeddings(list(range(8)),native,graph)
    assert list(map(lambda x:list(map(list,x)),autos))==sorted(data['native_configuration_automorphisms']) and len(autos)==2
    counts=Counter();precision=Counter()
    for rec in data['target_records']:
        frame,J=S.frame_matrix(rec);nodes=rec['literal_stars'];events=rec['event_records'];n=len(nodes)
        expected_keys={freeze(sorted(row['contacts'])) for row in rec['physical_maps']}
        assert len(events)==len(expected_keys)==4
        assert {freeze(e['key']) for e in events}==expected_keys
        assert rec['markers']==list(range(n,n+len(events)))
        assert sorted((row['event'],row['star']) for row in rec['event_star_gates'])==[
            (e,s) for e in range(len(events)) for s in range(n)]
        assert sorted(tuple(row['events']) for row in rec['event_event_gates'])==list(combinations(range(len(events)),2))
        target={i:{j:(m,'0','0') for j,m in row} for i,row in rec['finite_frame']}
        supports=sorted({tuple(r['support']) for r in rec['physical_maps']})
        assert supports==sorted(tuple(r['support']) for r in geo['outside_both'] if r['vertex']==rec['vertex'])
        for ids in supports:
            actual=complete_embeddings(list(ids),native,target)
            saved=sorted(tuple(map(tuple,r['mapping'])) for r in rec['physical_maps'] if tuple(r['support'])==ids)
            assert actual==saved and len(actual)==4
            counts['complete_support_embedding_domains']+=1
        for row in rec['physical_maps']:
            S.check_map(G,scales,native,row);mp=dict(row['mapping']);chosen=set(mp.values())
            assert all(frame[i][j]==0 for i in chosen for j in set(range(13))-chosen)
            assert all(frame[mp[i]][mp[j]]==-S.root(S.squared(native[i,j])) for i,j in combinations(row['support'],2))
            counts['physical_maps']+=1
        for e in events:
            maps=[r for r in rec['physical_maps'] if sorted(r['contacts'])==e['key']]
            assert len(maps)==2 and maps[0]['simplex_ports']==maps[1]['simplex_ports']==e['simplex_ports']
            m0,m1=(dict(r['mapping']) for r in maps);reflection=dict(geo['reflection'])
            assert all(m1[reflection[i]]==j for i,j in m0.items())
            for row in maps:
                simplex=next(s for s in geo['actual_simplices'] if set(s['fixed'])<=set(row['support']))
                mp=dict(row['mapping']);assert sorted(mp[i] for i in set(row['support'])&set(simplex['sides']))==e['simplex_ports']
            for c,h in zip(e['key'],e['heights']):
                pr,ps,hp=S.projections(c,c,J);assert all(x==0 for x in pr[1:]) and pr[0]-2==Q(*h)
        edges=set(map(tuple,rec['inherited_graph_edges']))
        for row in rec['event_star_gates']+rec['event_event_gates']:
            gate=row['gate']
            if 'event' in row:
                ei=row['event'];si=row['star'];a=events[ei]['key'];b=nodes[si]['contacts'];pair=tuple(sorted((n+ei,si)))
            else:
                i,j=row['events'];a=events[i]['key'];b=events[j]['key'];pair=(n+i,n+j)
            if gate is None:edges.add(pair);counts['conservative_compatible_pairs']+=1;continue
            if gate['kind']=='actual_simplex_triangle_port':
                assert 'event' in row and gate['ports'] in nodes[si]['occupied_pairs']
                assert set(gate['ports'])<=set(events[ei]['simplex_ports'])
            else:
                c,d=a[gate['first']],b[gate['second']];assert c!=d
                if gate['kind']=='positive_actual_product':precision[positive_gate(c,d,J,gate)]+=1
                else:
                    assert gate['kind']=='nonintegral_two_cycle'
                    pr,ps,hp=S.projections(c,d,J);assert all(x==0 for x in pr) and hp==Q(*gate['height_product'])
                    trace,norm=2*hp.a,hp.a*hp.a-5*hp.b*hp.b
                    assert str(trace)==gate['trace'] and str(norm)==gate['norm']
                    assert trace.denominator!=1 or norm.denominator!=1
            counts[gate['kind']]+=1
        assert sorted(edges)==list(map(tuple,rec['graph_edges']))
        assert sorted(e for e in edges if e[1]<n)==list(map(tuple,rec['inherited_graph_edges']))
        adj={i:set() for i in range(n+len(events))}
        for i,j in edges:assert i<j;adj[i].add(j);adj[j].add(i)
        actual=S.cliques(adj);assert list(map(list,actual))==rec['maximal_cliques']
        owner={i:j for j,g in enumerate(rec['corner_orbits']) for i in g['literal_stars']}
        assert sorted(owner)==list(range(n));items={}
        for cc in actual:
            for k in range(1,len(cc)+1):
                for sub in combinations(cc,k):
                    ct=Counter(owner[i] for i in sub if i<n)
                    vec=tuple(ct.get(j,0) for j in range(len(rec['corner_orbits'])))
                    event=sum(i>=n for i in sub);assert event<=1;key=vec,event
                    if key not in items or sub<items[key]:items[key]=sub
        assert rec['items']==[dict(vector=list(v),events=e,clique=list(items[v,e])) for v,e in sorted(items)]
        assert rec['event_capacity']==1
        counts.update(maximal_cliques=len(actual),count_event_vectors=len(items))
    records=base['forest_records']+data['target_records'];columns=[]
    assert [r['vertex'] for r in records]==data['joint_record_vertices']==[5587,3324,1880,2622]
    for di in range(2):
        col=verify_joint(records,base,data,di,counts)
        if col not in columns:columns.append(col)
    assert columns==data['columns']
    print(json.dumps(dict(status='PASS_INDEPENDENT_RAW39_ADDITIONAL_TARGET_SAVED_CHECK',
        counts=dict(counts),positive_gate_interval_bits=dict(precision),bindings=len(data['sources']),
        scope='Independent saved physical maps and complete active embeddings, rationally enclosed radical gates, graph/subset enumeration and exact four-vertex cover/charge/column repairs. Uses the bound earlier independent arithmetic helper; imports no producer. Historical star ownership/provenance and actual-face existence theorem are inherited. No global optimization or physical realization.'),sort_keys=True))


if __name__=='__main__':run()
