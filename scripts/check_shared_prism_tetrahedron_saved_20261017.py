#!/usr/bin/env python3
"""No-project-import saved graph/event/joint certificate arithmetic."""
import ast,json,re
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_shared_prism_tetrahedron_incidence_20261017'

def digest(p):
    h=sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def cliques(adj):
    def visit(s,p,x):
        if not p and not x:yield tuple(sorted(s));return
        pivot=max(p|x,key=lambda i:len(p&adj[i]))
        for i in sorted(p-adj[pivot]):
            yield from visit(s|{i},p&adj[i],x&adj[i]);p.remove(i);x.add(i)
    return sorted(visit(set(),set(adj),set()))

def run():
    path=ROOT/f'results/{TAG}.json';data=json.loads(path.read_text())
    assert data['status']=='PASS_SHARED_PRISM_TETRAHEDRON_INCIDENCE'
    for name,h in data['sources'].items():assert digest(ROOT/name)==h,name
    base=json.loads((ROOT/f"results/{data['baseline_source']}.json").read_text())
    ref=json.loads((ROOT/'results/fresh_forest1325_cross_vertex_incidence_20261017.json').read_text())
    geo=data['geometry'];perm=geo['reference_to_current_permutation']
    assert len(set(perm))==len(perm)==9
    assert [[geo['scaled_gram'][i][j]for j in perm]for i in perm]==ref['active_gram']
    supports={tuple(r['support'])for r in geo['positive_supports']}
    assert supports=={tuple(sorted(perm[i]for i in r['active_support']))for r in ref['support_records']if r['positive_definite']}
    assert len(supports)==15
    common=[r for r in geo['positive_supports']if r['in_raw23_prism']and r['in_raw25_prism']]
    assert common==geo['common_vertices']and len(common)==4
    assert {r['vertex']for r in common}=={1844,1325,1640,982}
    assert {tuple(r['support'])for r in common}=={tuple(sorted([0,1,5]+list(c)))for c in combinations([2,3,4,8],3)}
    records=data['forest_records'];assert [r['vertex']for r in records]==[1844,1325,1640,982]
    counts=Counter();columns=[]
    for r in records:
        groups=r['corner_orbits'];nodes=r['literal_stars'];owner={i:j for j,g in enumerate(groups)for i in g['literal_stars']}
        assert len(owner)==sum(len(g['literal_stars'])for g in groups)==len(nodes)and sorted(owner)==list(range(len(nodes)))
        edges=set(map(tuple,r['graph_edges']));assert len(edges)==len(r['graph_edges'])
        adj={i:set()for i in owner}
        for i,j in edges:assert 0<=i<j<len(nodes);adj[i].add(j);adj[j].add(i)
        pairs=set(map(tuple,r['pair_domain']));assert len(pairs)==len(r['pair_domain'])==len(r['physical_maps'])
        assert pairs=={tuple(x['stars'])for x in r['physical_maps']}
        for a,b in pairs:
            assert b in adj[a]and nodes[a]['origin']['face']==26448 and nodes[b]['origin']['face']==26918
        actual=cliques(adj);assert list(map(list,actual))==r['maximal_cliques'];vectors={}
        for clique in actual:
            for size in range(1,len(clique)+1):
                for sub in combinations(clique,size):
                    vv=tuple(Counter(owner[i]for i in sub).get(j,0)for j in range(len(groups)))
                    x=sum(a in sub and b in sub for a,b in pairs);assert x<=1;key=vv,x
                    if key not in vectors or sub<vectors[key]:vectors[key]=sub
        assert r['items']==[dict(vector=list(v),events=x,clique=list(vectors[v,x]))for v,x in sorted(vectors)]
        counts.update(forests=1,pairs=len(pairs),maximal_cliques=len(actual),count_event_vectors=len(vectors))
    summaries=[]
    for di,key in enumerate(('face_multipliers','restricted_face_multipliers')):
        assert data[key]==base[key];alpha={f:F(x)for f,x in data[key]}
        prior=base['repairs'][di];mu={v:F(x)for v,x in prior['vertex_multipliers']};quota=3*mu[6537]
        c=data['joint_certificates'][di];a=list(map(F,c['source_costs']));z=list(map(F,c['event_coefficients']))
        cap=F(c['quota_cost']);X=F(c['event_mass']);delta=list(map(F,c['increases']));weights=[list(map(F,w))for w in c['weights']]
        assert all(0<=x<=1 for x in a)and all(x>=0 for x in delta)and 0<=X<=quota and cap>=max(F(),sum(z,F()))
        coefficient=Counter();load=F();new=dict(mu)
        for r,w,aa,zz,inc,mix in zip(records,weights,a,z,delta,c['mixtures']):
            demand=[sum((alpha.get(f,F())*m for f,m in g['positive_face_coefficients']),F())for g in r['corner_orbits']]
            assert list(map(str,demand))==r['demands'][di]and len(w)==len(demand)and all(x>=0 for x in w)
            assert all(sum(x*y for x,y in zip(w,q['vector']))<=aa+zz*q['events']for q in r['items'])
            terms=[(j,F(x))for j,x in mix];assert len(terms)==len({j for j,x in terms})
            assert all(0<=j<len(r['items'])and x>0 for j,x in terms)
            assert all(sum((x*r['items'][j]['vector'][i]for j,x in terms),F())>=demand[i]for i in range(len(demand)))
            assert sum((x*r['items'][j]['events']for j,x in terms),F())==X
            assert sum((x for j,x in terms),F())<=mu[r['vertex']]+inc
            new[r['vertex']]+=inc;load+=sum((x*y for x,y in zip(w,demand)),F())
            for weight,g in zip(w,r['corner_orbits']):
                for f,m in g['positive_face_coefficients']:coefficient[f]+=weight*m
            counts.update(charge_inequalities=len(r['items']),mixture_terms=len(terms),event_equalities=1)
        lower=load-sum((aa*mu[r['vertex']]for aa,r in zip(a,records)),F())-cap*quota
        assert lower==sum(delta,F())==F(c['gain'])
        costs=[[r['vertex'],str(aa)]for r,aa in zip(records,a)if aa]
        if cap:costs.append([6537,str(3*cap)])
        col=dict(positive_face_coefficients=[[f,str(x)]for f,x in sorted(coefficient.items())if x],vertex_capacities=costs)
        if col not in columns:columns.append(col)
        repair=data['repairs'][di];den=sum(new.values(),F());num=F(prior['objective_numerator'])
        assert repair['vertex_multipliers']==[[v,str(x)]for v,x in sorted(new.items())if x]
        assert repair['vertex_updates']==[[v,str(mu[v]),str(x)]for v,x in sorted(new.items())if x>mu[v]]
        assert F(repair['normalizer'])==den and F(repair['objective_numerator'])==num and F(repair['bound'])==num/den==F(c['bound'])
        for col in data['columns']:
            assert sum((F(x)*alpha.get(f,F())for f,x in col['positive_face_coefficients']),F())<=sum((F(x)*new[v]for v,x in col['vertex_capacities']),F())
        summaries.append(dict(domain=repair['domain'],gain=str(lower),bound=str(num/den),common_event=str(X)))
    assert columns==data['columns']and counts['pairs']==36 and counts['maximal_cliques']==2634 and counts['count_event_vectors']==194
    assert F(data['repairs'][0]['bound'])<12<F(data['repairs'][1]['bound'])
    ast.parse((ROOT/f'scripts/{TAG}.py').read_text());ast.parse(Path(__file__).read_text())
    md=ROOT/f'results/{TAG}.md';links=[]
    for ref in re.findall(r'\]\(([^)]+)\)',md.read_text()):
        if not ref.startswith(('http:','https:','#')):assert(md.parent/ref.split('#')[0]).exists(),ref;links.append(ref)
    print(json.dumps(dict(status='PASS_NO_PRODUCER_SHARED_TETRAHEDRON_SAVED_CHECK',counts=dict(counts),repairs=summaries,
        source_bindings=len(data['sources']),links=len(links),json_sha256=digest(path),checker_sha256=digest(Path(__file__)),
        scope='Independent saved graph enumeration/subset events, support-permutation and saved pair-domain consistency, actual coefficient aggregation, all rational joint charges/covers/equal-event masses and two repairs. No independent physical projection, geometric tetrahedron/injection proof or inherited pair-owner reconstruction.'),sort_keys=True))

if __name__=='__main__':run()
