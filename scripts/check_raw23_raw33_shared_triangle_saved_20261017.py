#!/usr/bin/env python3
"""No-project-import saved graph/event/joint certificate arithmetic."""
import ast,json,re
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw23_raw33_shared_triangle_frontier_20261017'

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

# Independent ordered arithmetic in Q(sqrt5), used only on saved matrices.
class Q:
    def __init__(self,a=0,b=0):self.a,self.b=F(a),F(b)
    def __add__(self,o):
        if not isinstance(o,Q):o=Q(o)
        return Q(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return Q(-self.a,-self.b)
    def __sub__(self,o):return self+-o
    def __mul__(self,o):
        if not isinstance(o,Q):o=Q(o)
        return Q(self.a*o.a+5*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        if not isinstance(o,Q):o=Q(o)
        n=o.a*o.a-5*o.b*o.b;assert n
        return self*Q(o.a/n,-o.b/n)
    def sign(self):
        if not self.b:return (self.a>0)-(self.a<0)
        if not self.a:return (self.b>0)-(self.b<0)
        if self.a*self.b>0:return (self.a>0)-(self.a<0)
        n=self.a*self.a-5*self.b*self.b
        return ((self.a>0)-(self.a<0))*((n>0)-(n<0))
    def __eq__(self,o):
        if not isinstance(o,Q):o=Q(o)
        return self.a==o.a and self.b==o.b

def mat(x):return [[Q(*y)for y in row]for row in x]
def principal(g,ids):return [[g[i][j]for j in ids]for i in ids]
def positive(g):
    a=[r[:]for r in g]
    for k in range(len(a)):
        if a[k][k].sign()<=0:return False
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):a[i][j]=a[i][j]-a[i][k]*a[k][j]/a[k][k]
    return True

def inverse(g):
    n=len(g);a=[row[:]+[Q(int(i==j))for j in range(n)]for i,row in enumerate(g)]
    for i in range(n):
        k=next(k for k in range(i,n)if a[k][i]!=0);a[i],a[k]=a[k],a[i]
        z=a[i][i];a[i]=[x/z for x in a[i]]
        for k in range(n):
            if k!=i:
                z=a[k][i];a[k]=[x-z*y for x,y in zip(a[k],a[i])]
    return [row[n:]for row in a]

def run():
    path=ROOT/f'results/{TAG}.json';data=json.loads(path.read_text())
    assert data['status']=='PASS_RAW23_RAW33_SHARED_TRIANGLE_FRONTIER'
    for name,h in data['sources'].items():assert digest(ROOT/name)==h,name
    base=json.loads((ROOT/f"results/{data['baseline_source']}.json").read_text())
    geo=data['geometry'];G=mat(geo['scaled_gram']);records0=geo['support_records']
    assert [r['support']for r in records0]==list(map(list,combinations(range(9),6)))
    for r in records0:assert positive(principal(G,r['support']))==(r['inertia']==[6,0,0])
    pd=[r for r in records0 if r['inertia']==[6,0,0]];assert pd==geo['positive_supports']and len(pd)==13
    assert sum(r['in_raw23_prism']for r in pd)==10 and sum(r['in_raw33_prism']for r in pd)==6
    for r in pd:
        ids=set(r['support'])
        assert r['in_raw23_prism']==(7 not in ids and 1 in ids)
        assert r['in_raw33_prism']==(6 not in ids and {0,4,5}<=ids)
        assert r['in_raw23_prism']or r['in_raw33_prism']
    assert not geo['outside_both']
    common=[r for r in pd if r['in_raw23_prism']and r['in_raw33_prism']]
    assert common==geo['common_vertices']and {r['vertex']for r in common}=={1806,1505,944}
    fixed,sides=[0,1,4,5],[2,3,8]
    assert {tuple(r['support'])for r in common}=={tuple(sorted(fixed+list(c)))for c in combinations(sides,2)}
    J=inverse(principal(G,fixed))
    residual=[[G[i][j]-sum((G[i][a]*J[k][l]*G[b][j]for k,a in enumerate(fixed)for l,b in enumerate(fixed)),Q())for j in sides]for i in sides]
    assert residual==mat(geo['triangle_residual'])
    inv=inverse(residual);assert inv==mat(geo['triangle_inverse'])and all(x.sign()<0 for row in inv for x in row)
    assert G[6][7]==Q(-2,-1)and (G[6][7]+2).sign()<0
    records=data['forest_records'];assert [r['vertex']for r in records]==[1806,1505,944]
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
            assert b in adj[a]and nodes[a]['origin']['face']==26410 and nodes[b]['origin']['face']==27350
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
        prior=base['repairs'][di];mu={v:F(x)for v,x in prior['vertex_multipliers']}
        c=data['joint_certificates'][di];a=list(map(F,c['source_costs']));z=list(map(F,c['event_coefficients']))
        cap=F(c['quota_cost']);X=F(c['event_mass']);delta=list(map(F,c['increases']));weights=[list(map(F,w))for w in c['weights']]
        assert all(0<=x<=1 for x in a)and all(x>=0 for x in delta)and X>=0 and cap==0 and sum(z,F())<=0
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
        lower=load-sum((aa*mu[r['vertex']]for aa,r in zip(a,records)),F())
        assert lower==sum(delta,F())==F(c['gain'])
        costs=[[r['vertex'],str(aa)]for r,aa in zip(records,a)if aa]
        col=dict(positive_face_coefficients=[[f,str(x)]for f,x in sorted(coefficient.items())if x],vertex_capacities=costs)
        assert col==dict(positive_face_coefficients=[],vertex_capacities=[])
        assert lower==0 and all(x==0 for x in delta)
        repair=data['repairs'][di];den=sum(new.values(),F());num=F(prior['objective_numerator'])
        assert repair['vertex_multipliers']==[[v,str(x)]for v,x in sorted(new.items())if x]
        assert new==mu and repair==prior
        assert F(repair['normalizer'])==den and F(repair['objective_numerator'])==num and F(repair['bound'])==num/den==F(c['bound'])
        for col in data['columns']:
            assert sum((F(x)*alpha.get(f,F())for f,x in col['positive_face_coefficients']),F())<=sum((F(x)*new[v]for v,x in col['vertex_capacities']),F())
        summaries.append(dict(domain=repair['domain'],gain=str(lower),bound=str(num/den),common_event=str(X)))
    assert columns==data['columns']and counts['pairs']==28 and counts['maximal_cliques']==7766 and counts['count_event_vectors']==317
    assert data['repairs']==base['repairs']and F(data['repairs'][0]['bound'])<12<F(data['repairs'][1]['bound'])
    ast.parse((ROOT/f'scripts/{TAG}.py').read_text());ast.parse(Path(__file__).read_text())
    md=ROOT/f'results/{TAG}.md';links=[]
    for ref in re.findall(r'\]\(([^)]+)\)',md.read_text()):
        if not ref.startswith(('http:','https:','#')):assert(md.parent/ref.split('#')[0]).exists(),ref;links.append(ref)
    print(json.dumps(dict(status='PASS_NO_PRODUCER_SHARED_TRIANGLE_SAVED_CHECK',counts=dict(counts),repairs=summaries,
        source_bindings=len(data['sources']),links=len(links),json_sha256=digest(path),checker_sha256=digest(Path(__file__)),
        scope='Independent saved84-support positivity and triangle inverse in Q(sqrt5), graph enumeration/subset events, saved pair-domain consistency, actual coefficient aggregation, all rational joint charges/covers/equal-event masses and unchanged witnesses. No independent native projection, actual-face proof or inherited pair-owner reconstruction.'),sort_keys=True))

if __name__=='__main__':run()
