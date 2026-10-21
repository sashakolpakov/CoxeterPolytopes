#!/usr/bin/env python3
"""Independent saved-graph and Fraction checks; no project arithmetic imports."""
import ast,json,re
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw23_raw25_cross_vertex_family_20261017'

def digest(path):
    h=sha256()
    with path.open('rb')as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()

def maximal_cliques(adj):
    def visit(s,p,x):
        if not p and not x:yield tuple(sorted(s));return
        pivot=max(p|x,key=lambda i:len(p&adj[i]))
        for i in sorted(p-adj[pivot]):
            yield from visit(s|{i},p&adj[i],x&adj[i])
            p.remove(i);x.add(i)
    return sorted(visit(set(),set(adj),set()))

def check(items,d,quota,c,event=True):
    w=list(map(F,c['weights']));z=F(c['event_penalty'])if event else F()
    assert len(w)==len(d)and min(w)>=0 and z>=0
    assert all(sum(a*b for a,b in zip(w,it['vector']))<=1+z*it['events']for it in items)
    mix=[(j,F(x))for j,x in c['mixture']]
    assert len(mix)==len({j for j,x in mix})and all(x>0 and 0<=j<len(items)for j,x in mix)
    assert all(sum((x*items[j]['vector'][i]for j,x in mix),F())>=d[i]for i in range(len(d)))
    events=sum((x*items[j]['events']for j,x in mix),F())
    assert events<=quota
    mass=sum((x for j,x in mix),F());assert mass==sum((a*b for a,b in zip(w,d)),F())-z*quota==F(c['optimal_mass'])
    if event:assert events==F(c['event_mass'])
    return len(items),len(mix)

def column(groups,c,r,event=True):
    coeff=Counter()
    for w,g in zip(map(F,c['weights']),groups):
        for f,m in g['positive_face_coefficients']:coeff[f]+=w*m
    return dict(source_vertex=r['vertex'],target_vertex=r['target_vertex'],
        target_vertex_cost=str(F(c['event_penalty'])*r['pointed_capacity'])if event else '0',
        positive_face_coefficients=[[f,str(x)]for f,x in sorted(coeff.items())if x])

def run():
    file=ROOT/f'results/{TAG}.json';q=json.loads(file.read_text())
    assert q['status']=='PASS_RAW23_RAW25_CROSS_VERTEX_FAMILY'
    for path,h in q['sources'].items():assert digest(ROOT/path)==h,path
    b=json.loads((ROOT/f"results/{q['baseline_source']}.json").read_text())
    alphas=[];mus=[]
    for i,k in enumerate(('face_multipliers','restricted_face_multipliers')):
        assert q[k]==b[k];alphas.append({f:F(x)for f,x in q[k]})
        mus.append({v:F(x)for v,x in b['repairs'][i]['vertex_multipliers']})
    sources={r['vertex']for r in q['forest_records']};targets={r['target_vertex']for r in q['forest_records']}
    assert len(sources)==len(targets)==len(q['forest_records'])==7 and not sources&targets
    assert sources=={453,806,986,1022,1062,1325,1329}
    counts=Counter();cols=[];hcols=[]
    for r in q['forest_records']:
        groups=r['corner_orbits'];nodes=r['literal_stars'];n=len(nodes)
        owner={i:j for j,g in enumerate(groups)for i in g['literal_stars']}
        assert len(owner)==sum(len(g['literal_stars'])for g in groups)==n and sorted(owner)==list(range(n))
        adj={i:set()for i in owner};edges=set(map(tuple,r['graph_edges']));assert len(edges)==len(r['graph_edges'])
        for i,j in edges:assert 0<=i<j<n;adj[i].add(j);adj[j].add(i)
        af,bf=r['faces']
        aa=[i for i,s in enumerate(nodes)if s['origin']['face']==af and s['origin']['omitted']==[2,5]]
        bb=[i for i,s in enumerate(nodes)if s['origin']['face']==bf and s['origin']['omitted']==[2,3]]
        pairs=sorted((a,b)for a in aa for b in bb if b in adj[a]);assert pairs==list(map(tuple,r['compatible_selected_pairs']))
        assert [x['stars']for x in r['native_maps']]==r['compatible_selected_pairs']
        assert all(x['target_forest'].count('H3')==r['pointed_capacity']for x in r['native_maps'])
        cl=maximal_cliques(adj);assert list(map(list,cl))==r['maximal_cliques'];items={}
        for clique in cl:
            for size in range(1,len(clique)+1):
                for sub in combinations(clique,size):
                    v=tuple(Counter(owner[i]for i in sub).get(j,0)for j in range(len(groups)))
                    x=sum(a in sub and b in sub for a,b in pairs);key=v,x
                    if key not in items or sub<items[key]:items[key]=sub
        expected=[dict(vector=list(v),events=x,clique=list(items[v,x]))for v,x in sorted(items)]
        assert expected==r['items'];zi=[x for x in expected if x['events']==0];assert zi==r['zero_event_items']
        for di,alpha in enumerate(alphas):
            demand=[sum((alpha.get(f,F())*m for f,m in g['positive_face_coefficients']),F())for g in groups]
            quota=r['pointed_capacity']*mus[di].get(r['target_vertex'],F())
            assert list(map(str,demand))==r['demands'][di]and str(quota)==r['quotas'][di]
            for event,domain,c in((True,expected,r['certificates'][di]),(False,zi,r['zero_event_certificates'][di])):
                ni,nm=check(domain,demand,quota if event else F(),c,event)
                counts.update(inequalities=ni,mixture_terms=nm,certificates=1)
                co=column(groups,c,r,event);target=cols if event else hcols
                if co not in target:target.append(co)
        counts.update(forests=1,complete_pairs=len(pairs),maximal_cliques=len(cl),event_vectors=len(expected))
    assert cols==q['columns']and hcols==q['hypothetical_zero_event_columns']
    for r in q['zero_endpoint_contexts']:
        for di,alpha in enumerate(alphas):
            assert any(not alpha.get(f,F())for f in r['faces'])
            demand=[sum((alpha.get(f,F())*m for f,m in g['positive_face_coefficients']),F())for g in r['corner_orbits']]
            cover=r['projected_covers'][di];assert list(map(str,demand))==cover['demands']
            assert all(sum((F(t['amount'])*t['vector'][i]for t in cover['terms']),F())>=demand[i]for i in range(len(demand)))
            assert sum((F(t['amount'])for t in cover['terms']),F())==F(cover['cover_mass'])<=mus[di].get(r['vertex'],F())
    summaries=[]
    for name,columns in(('repairs',cols),('hypothetical_zero_event_repairs',hcols)):
        for di,alpha in enumerate(alphas):
            mu=dict(mus[di]);loads=[]
            for co in columns:
                load=sum((F(x)*alpha.get(f,F())for f,x in co['positive_face_coefficients']),F());loads.append(str(load))
                v=co['source_vertex'];mu[v]=max(mu.get(v,F()),load-F(co['target_vertex_cost'])*mus[di].get(co['target_vertex'],F()))
            r=q[name][di];assert r['loads']==loads and r['vertex_multipliers']==[[v,str(x)]for v,x in sorted(mu.items())if x]
            assert r['vertex_updates']==[[v,str(mus[di].get(v,F())),str(x)]for v,x in sorted(mu.items())if x>mus[di].get(v,F())]
            numerator=F(b['repairs'][di]['objective_numerator']);den=sum(mu.values(),F())
            assert F(r['normalizer'])==den and F(r['bound'])==numerator/den and F(r['objective_numerator'])==numerator
            for co in cols:
                assert sum((F(x)*alpha.get(f,F())for f,x in co['positive_face_coefficients']),F())<=mu[co['source_vertex']]+F(co['target_vertex_cost'])*mu.get(co['target_vertex'],F())
            for fr in q['forest_records']:
                key='certificates'if name=='repairs'else'zero_event_certificates'
                assert mu[fr['vertex']]==max(mus[di].get(fr['vertex'],F()),F(fr[key][di]['optimal_mass']))
            summaries.append(dict(mode=name,domain=di,bound=r['bound'],updates=len(r['vertex_updates'])))
    assert len(q['zero_endpoint_contexts'])==6 and counts['forests']==7 and counts['complete_pairs']==78
    assert F(q['hypothetical_zero_event_repairs'][1]['bound'])>12
    ast.parse((ROOT/f'scripts/{TAG}.py').read_text());ast.parse(Path(__file__).read_text())
    links=[];md=ROOT/f'results/{TAG}.md'
    if md.exists():
        for ref in re.findall(r'\]\(([^)]+)\)',md.read_text()):
            if not ref.startswith(('http:','https:','#')):assert(md.parent/ref.split('#')[0]).exists(),ref;links.append(ref)
    print(json.dumps(dict(status='PASS_NO_PRODUCER_CROSS_VERTEX_FAMILY_SAVED_CHECK',counts=dict(counts),sources=len(q['sources']),
        columns=len(cols),hypothetical_columns=len(hcols),repairs=summaries,links=len(links),json_sha256=digest(file),checker_sha256=digest(Path(__file__)),
        scope='Independent saved-graph maximal cliques and every subset event/count vector; all exact charges, covers, actual occurrence coefficients and four monotone repairs. Native contact maps and actual geometric injection remain source/proof-reviewed inputs; no independent physical projection or frozen ownership reconstruction.'),sort_keys=True))

if __name__=='__main__':run()
