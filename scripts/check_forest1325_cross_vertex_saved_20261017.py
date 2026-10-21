#!/usr/bin/env python3
"""Independent standard-library checks of saved1325 graph/count arithmetic.

No project imports, physical wall reconstruction or geometric injection proof.
Re-enumerates maximal cliques from the saved graph, then verifies the marked
charge and both multi-vertex repairs against their frozen baseline.
"""
import ast,json,re
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_forest1325_cross_vertex_incidence_20261017'


def digest(path):
    h=sha256()
    with path.open('rb')as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def maximal_cliques(adj):
    """Bron--Kerbosch recursion with a complete pivot branch partition."""
    def visit(selected,possible,earlier):
        if not possible and not earlier:
            yield tuple(sorted(selected));return
        pivot=max(possible|earlier,key=lambda i:len(possible&adj[i]))
        for i in sorted(possible-adj[pivot]):
            yield from visit(selected|{i},possible&adj[i],earlier&adj[i])
            possible.remove(i);earlier.add(i)
    return sorted(visit(set(),set(adj),set()))


def run():
    path=ROOT/f'results/{TAG}.json';data=json.loads(path.read_text())
    assert data['status']=='PASS_FOREST1325_CROSS_VERTEX_FORCED_INCIDENCE'
    for name,h in data['sources'].items():assert digest(ROOT/name)==h,name
    basepath=ROOT/f"results/{data['baseline_source']}.json"
    assert data['sources'][str(basepath.relative_to(ROOT))]==digest(basepath)
    base=json.loads(basepath.read_text())
    groups=data['corner_orbits'];nodes=data['literal_stars'];n=len(nodes)
    owner={i:j for j,q in enumerate(groups)for i in q['literal_stars']}
    assert len(owner)==sum(len(q['literal_stars'])for q in groups)==n
    assert sorted(owner)==list(range(n)) and len(groups)==8
    edges=set(map(tuple,data['graph_edges']))
    assert len(edges)==len(data['graph_edges']) and all(0<=i<j<n for i,j in edges)
    adj={i:set()for i in range(n)}
    for i,j in edges:adj[i].add(j);adj[j].add(i)
    pairs=set(map(tuple,data['compatible_corner_pairs']))
    assert len(pairs)==12 and pairs=={(i,j)for i in groups[1]['literal_stars']for j in groups[5]['literal_stars']if j in adj[i]}
    endpoint_degrees=Counter(i for p in pairs for i in p)
    assert set(endpoint_degrees)==set(groups[1]['literal_stars']+groups[5]['literal_stars'])
    assert all(x==1 for x in endpoint_degrees.values())
    weights=data['local_orbit_weights'];assert weights==[3,1,2,4,4,4,0,0]
    assert all(x>=0 for x in weights)and all(weights[owner[i]]>=1 for i in endpoint_degrees)
    actual=maximal_cliques(adj);saved={tuple(r['stars']):r for r in data['maximal_clique_checks']}
    assert len(actual)==len(saved)==522 and actual==sorted(saved)
    patterns=Counter()
    for cl in actual:
        assert all(j in adj[i]for i,j in combinations(cl,2))
        vv=[Counter(owner[i]for i in cl).get(j,0)for j in range(8)]
        mixed=sum(i in cl and j in cl for i,j in pairs)
        charge=sum(x*y for x,y in zip(weights,vv));slack=4+mixed-charge
        assert mixed<=1 and slack>=0
        assert saved[cl]==dict(stars=list(cl),count_vector=vv,mixed_pairs=mixed,charge=charge,slack=slack)
        patterns[tuple(vv),mixed]+=1
    expected=[dict(vector=list(v),mixed_pairs=x,multiplicity=m)for(v,x),m in sorted(patterns.items())]
    assert len(patterns)==8 and expected==data['count_patterns']
    coefficient=Counter()
    for weight,q in zip(weights,groups):
        assert q['positive_face_coefficients']==[[q['face'],1]]
        for i,c in q['positive_face_coefficients']:coefficient[i]+=weight*F(c)
    terms=[[i,int(c)]for i,c in sorted(coefficient.items())if c]
    assert terms==[[26362,8],[26448,4],[26918,3],[45769,3]]
    column=data['columns'][0];assert len(data['columns'])==1
    assert column['positive_face_coefficients']==terms and column['vertex_capacities']==[[1325,4],[6537,3]]
    assert data['pointed_incidence_capacity']==3
    assert len(data['repairs'])==len(base['repairs'])==2
    summaries=[]
    for k,key in enumerate(('face_multipliers','restricted_face_multipliers')):
        assert data[key]==base[key]
        alpha={i:F(x)for i,x in data[key]};old=base['repairs'][k];new=data['repairs'][k]
        mu={i:F(x)for i,x in old['vertex_multipliers']};nm={i:F(x)for i,x in new['vertex_multipliers']}
        load=sum((F(c)*alpha.get(i,F())for i,c in terms),F())
        budget=sum((F(c)*mu.get(i,F())for i,c in column['vertex_capacities']),F())
        finalbudget=sum((F(c)*nm.get(i,F())for i,c in column['vertex_capacities']),F())
        gap=load-budget;assert gap>0 and finalbudget==load
        assert load==F(column['loads'][k])==F(new['baseline_load'])
        assert budget==F(new['baseline_capacity'])and gap==F(new['gap'])
        expectedmu=dict(mu);expectedmu[1325]+=gap/4;assert nm==expectedmu
        assert new['vertex_updates']==[[1325,str(mu[1325]),str(nm[1325])]]
        assert sum(mu.values(),F())==F(old['normalizer'])
        numerator=F(old['objective_numerator']);den=sum(nm.values(),F());ratio=numerator/den
        assert new['objective_numerator']==old['objective_numerator']and new['baseline_bound']==old['bound']
        assert ratio==F(new['bound'])and den==F(new['normalizer'])
        assert new['domain']==old['domain']
        oldinc=2*alpha[26448]+alpha[26918]-2*alpha[311]-3*mu[6537]
        assert oldinc==F(new['old_incidence_row50_gap'])<0
        summaries.append(dict(domain=new['domain'],gap=str(gap),increase=str(gap/4),bound=str(ratio)))
    assert F(summaries[0]['bound'])<12<F(summaries[1]['bound'])
    ast.parse((ROOT/f'scripts/{TAG}.py').read_text());ast.parse(Path(__file__).read_text())
    report=ROOT/f'results/{TAG}.md';links=[]
    for ref in re.findall(r'\]\(([^)]+)\)',report.read_text()):
        if not ref.startswith(('http:','https:','#')):
            assert(report.parent/ref.split('#')[0]).exists(),ref;links.append(ref)
    print(json.dumps(dict(status='PASS_NO_PRODUCER_FOREST1325_SAVED_GRAPH_AND_REPAIRS',
        producer_json_sha256=digest(path),checker_sha256=digest(Path(__file__)),
        sources=len(data['sources']),maximal_cliques=len(actual),count_patterns=len(patterns),
        matching_pairs=len(pairs),source_face_coefficients=len(terms),vertex_capacities=2,repairs=summaries,
        local_links=len(links),
        scope='Independent saved-graph maximal-clique enumeration and charge/matching/count arithmetic, actual coefficient aggregation and both exact monotone multi-vertex repairs. No producer imports, physical projection/support/orientation reconstruction, pair-outcome verification or proof of the geometric pointed injection.'),sort_keys=True))


if __name__=='__main__':run()
