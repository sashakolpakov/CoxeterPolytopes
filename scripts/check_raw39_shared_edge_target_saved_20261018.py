#!/usr/bin/env python3
"""No-project-import saved graph/event/joint certificate arithmetic."""
import ast,json,re
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw39_shared_edge_target_incidence_20261018'

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

def squared(e):
    m,a,b=e
    if m==0:return Q(a,b)
    return {2:Q(),3:Q(1),4:Q(2),5:Q('3/2','1/2'),6:Q(3),10:Q('5/2','1/2')}[m]

def rational_root(x):
    from math import isqrt
    if x<0:return None
    a,b=isqrt(x.numerator),isqrt(x.denominator)
    return F(a,b)if a*a==x.numerator and b*b==x.denominator else None

def root(q):
    if q.sign()<0:return None
    if q==0:return Q()
    if q.b==0:
        d=rational_root(q.a/5)
        if d is not None:return Q(0,d)
    z=rational_root(q.a*q.a-5*q.b*q.b)
    if z is None:return None
    for a in((q.a+z)/2,(q.a-z)/2):
        c=rational_root(a)
        if c is None or not c:continue
        r=Q(c,q.b/(2*c))
        if r*r==q:return r if r.sign()>0 else -r
    return None

def frame_matrix(record):
    rows={i:dict(row)for i,row in record['finite_frame']};assert set(rows)==set(range(13))
    g=[[Q(2)if i==j else Q()for j in range(13)]for i in range(13)]
    for i,row in rows.items():
        for j,m in row.items():
            assert rows[j][i]==m;z=root(squared([m,'0','0']));assert z is not None;g[i][j]=-z
    assert positive(g)
    return g,inverse(g)

def projections(c,d,J):
    # Build only the radical classes present in these two native covectors.
    basis=[]
    for i,e in c+d:
        t=squared(e)
        if t!=0 and root(t)is None and all(root(t/u)is None for u in basis):basis.append(t)
    assert len(basis)<=2
    size=1<<len(basis);powers=[Q(1)]
    for mask in range(1,size):
        x=Q(1)
        for k,t in enumerate(basis):
            if mask>>k&1:x=x*t
        powers.append(x)
    def mul(a,b):
        out=[Q()for _ in range(size)]
        for i,x in enumerate(a):
            for j,y in enumerate(b):out[i^j]=out[i^j]+x*y*powers[i&j]
        return out
    def contact(e):
        t=squared(e);a=[Q()for _ in range(size)];z=root(t)
        if z is not None:a[0]=-z;return a
        for k,u in enumerate(basis):
            z=root(t/u)
            if z is not None:a[1<<k]=-z;return a
        raise AssertionError('missing native squareclass')
    def product(x,y):
        out=[Q()for _ in range(size)]
        for i,e in x:
            for j,f in y:
                p=mul(contact(e),contact(f));out=[a+J[i][j]*b for a,b in zip(out,p)]
        return out
    cc,dd,p=product(c,c),product(d,d),product(c,d)
    assert all(x==0 for x in cc[1:]+dd[1:])
    hp=(cc[0]-2)*(dd[0]-2);assert (cc[0]-2).sign()>0 and (dd[0]-2).sign()>0
    return p,mul(p,p),hp

def check_map(G,scales,native,row):
    ids=row['support'];mp=dict(row['mapping']);assert len(ids)==len(mp)==len(set(mp.values()))==5 and set(mp)==set(ids)
    ext=sorted(set(range(8))-set(ids));assert row['external_roots']==ext
    expected=[[[mp[j],list(native[i,j])]for j in sorted(ids,key=lambda j:mp[j])if native[i,j][0]!=2]for i in ext]
    assert row['contacts']==expected
    J=inverse(principal(G,ids));proj={}
    for i in ext:
        for j in ext:proj[i,j]=sum((G[i][a]*J[k][l]*G[b][j]for k,a in enumerate(ids)for l,b in enumerate(ids)),Q())
    h={i:proj[i,i]-G[i][i]for i in ext}
    for i,x in zip(ext,row['heights']):assert h[i].sign()>0 and h[i]/scales[i]==Q(*x)
    assert len(row['cross_checks'])==3
    for (i,j),r in zip(combinations(ext,2),row['cross_checks']):
        assert r['walls']==[i,j]and tuple(r['native_product'])==native[i,j]
        d=proj[i,j]-G[i][j];assert d.sign()>0 and d*d==h[i]*h[j]
        assert Q(*r['height_product'])==h[i]*h[j]/(scales[i]*scales[j])

def residual(G,fixed,sides):
    J=inverse(principal(G,fixed))
    return [[G[i][j]-sum((G[i][a]*J[k][l]*G[b][j]for k,a in enumerate(fixed)for l,b in enumerate(fixed)),Q())for j in sides]for i in sides]

def run():
    path=ROOT/f'results/{TAG}.json';data=json.loads(path.read_text());assert data['status']=='PASS_RAW39_SHARED_EDGE_TARGET_INCIDENCE'
    for name,h in data['sources'].items():assert digest(ROOT/name)==h,name
    base=json.loads((ROOT/f"results/{data['baseline_source']}.json").read_text());geo=data['geometry'];G=mat(geo['scaled_gram']);scales=[Q(*x)for x in geo['squared_scales']]
    native={(i,j):(2,'0','0')for i in range(8)for j in range(8)if i!=j}
    for i,j,e in geo['native_edges']:native[i,j]=native[j,i]=tuple(e)
    for i in range(8):
        assert G[i][i]==2*scales[i]
        for j in range(i+1,8):assert G[i][j].sign()<=0 and G[i][j]*G[i][j]==squared(native[i,j])*scales[i]*scales[j]
    reflection=dict(geo['reflection']);assert set(reflection)==set(reflection.values())==set(range(8))
    assert all(native[i,j]==native[reflection[i],reflection[j]]for i,j in combinations(range(8),2))
    supports=geo['support_records'];assert [r['support']for r in supports]==list(map(list,combinations(range(8),5)))
    for r in supports:assert positive(principal(G,r['support']))==(r['inertia']==[5,0,0])
    pd=[r for r in supports if r['inertia']==[5,0,0]];assert pd==geo['positive_supports']and len(pd)==12
    for r in pd:
        ids=set(r['support']);assert r['in_first_quad']==(7 not in ids and {1,2,3}<=ids)
        assert r['in_second_quad']==(5 not in ids and {0,1,2}<=ids)
    assert sum(r['in_first_quad']for r in pd)==sum(r['in_second_quad']for r in pd)==4
    assert geo['shared_edge_vertices']==[r for r in pd if r['in_first_quad']and r['in_second_quad']]
    assert geo['outside_both']==[r for r in pd if not(r['in_first_quad']or r['in_second_quad'])]
    assert Counter(r['vertex']for r in geo['outside_both'])=={3324:2,1880:2,2622:2}
    edge=residual(G,geo['edge_fixed'],geo['edge_sides']);assert edge==mat(geo['edge_residual'])
    assert inverse(edge)==mat(geo['edge_inverse'])and all(x.sign()<0 for row in inverse(edge)for x in row)
    for s in geo['actual_simplices']:
        assert s['supports']==[sorted(s['fixed']+list(c))for c in combinations(s['sides'],4)]
        assert all(any(r['support']==q for r in pd)for q in s['supports'])
        res=residual(G,s['fixed'],s['sides']);inv=inverse(res)
        assert res==mat(s['residual'])and inv==mat(s['inverse'])and all(x.sign()<0 for row in inv for x in row)
    records=data['forest_records'];assert [r['vertex']for r in records]==[5587,3324]
    counts=Counter(active_supports=56,positive_supports=12);target=data['target_compatibility']
    assert target['physical_maps']==records[1]['physical_maps']and len(target['physical_maps'])==2
    maps=target['physical_maps'];assert sorted(maps[0]['contacts'])==sorted(maps[1]['contacts'])
    m0,m1=(dict(r['mapping'])for r in maps);assert all(m1[reflection[i]]==j for i,j in m0.items())
    assert target['simplex_ports']==sorted(m0[i]for i in(0,1,4,6))
    frame,J=frame_matrix(records[1]);cov=dict(zip(maps[0]['external_roots'],maps[0]['contacts']));nodes=records[1]['literal_stars']
    assert sorted(r['star']for r in target['exclusions'])==list(range(15))and not target['compatible_stars']
    for r in target['exclusions']:
        s=nodes[r['star']]
        if r['kind']=='actual_simplex_triangle_port':
            assert r['ports']in s['occupied_pairs']and set(r['ports'])<=set(target['simplex_ports']);counts['triangle_port_exclusions']+=1;continue
        c=s['contacts'][r['star_external_index']];d=cov[r['event_wall']];assert c!=d
        pr,ps,hp=projections(c,d,J);assert hp==Q(*r['height_product'])
        if r['kind']=='positive_actual_product':
            assert all(x.sign()>=0 for x in pr)and any(x.sign()>0 for x in pr)
            assert all(x==0 for x in ps[1:])and ps[0]==Q(*r['projection_square'])and (ps[0]-hp).sign()>0
            counts['positive_product_exclusions']+=1
        else:
            assert r['kind']=='nonintegral_two_cycle'and all(x==0 for x in pr)
            assert Q(*r['two_cycle'])==hp and F(r['trace'])==2*hp.a and F(r['norm'])==hp.a*hp.a-5*hp.b*hp.b
            assert F(r['trace']).denominator!=1 or F(r['norm']).denominator!=1;counts['nonintegral_exclusions']+=1
    for rec in records:
        frame0,J0=frame_matrix(rec);groups=rec['corner_orbits'];nodes=rec['literal_stars'];marker=rec['marker']
        owner={i:j for j,g in enumerate(groups)for i in g['literal_stars']};assert sorted(owner)==list(range(len(nodes)))
        for row in rec['physical_maps']:
            check_map(G,scales,native,row);mp=dict(row['mapping'])
            assert all(frame0[i][j]==0 for i in mp.values()for j in set(range(13))-set(mp.values()))
            for i,j in combinations(row['support'],2):
                assert frame0[mp[i]][mp[j]]==-root(squared(native[i,j]))
            counts['physical_saved_maps']+=1
        if marker is None:
            domains=[]
            for support in([0,1,2,3,4],[0,1,2,3,6]):
                rows=[r for r in rec['physical_maps']if r['support']==support];assert len(rows)==4
                keys={tuple(sorted(r['stars']))for r in rows};assert len(keys)==4;domains.append(keys)
                for r in rows:
                    ext=r['external_roots'];covs=dict(zip(ext,r['contacts']));mp=dict(r['mapping'])
                    for star,absent,fixed in zip(r['stars'],(7,5),({1,2,3},{0,1,2})):
                        s=nodes[star];assert s['origin']['face']==28827 and s['origin']['raw_family']==39
                        assert sorted(s['contacts'])==sorted(covs[i]for i in ext if i!=absent)
                        assert set(range(13))-set(sum(s['occupied_pairs'],[]))==(set(range(13))-set(mp.values()))|{mp[i]for i in fixed}
            assert domains[0]==domains[1]==set(map(tuple,rec['event_domain']))
        else:assert marker==len(nodes)==15 and rec['event_domain']==[[15]]
        edges=set(map(tuple,rec['graph_edges']));assert len(edges)==len(rec['graph_edges'])
        assert sorted(e for e in edges if marker not in e)==sorted(map(tuple,rec['inherited_graph_edges']))
        adj={i:set()for i in range(len(nodes)+(marker is not None))}
        for i,j in edges:assert i<j;adj[i].add(j);adj[j].add(i)
        actual=cliques(adj);assert list(map(list,actual))==rec['maximal_cliques'];items={}
        for cc in actual:
            for k in range(1,len(cc)+1):
                for sub in combinations(cc,k):
                    ct=Counter(owner[i]for i in sub if i!=marker);vec=tuple(ct.get(j,0)for j in range(len(groups)))
                    e=sum(set(q)<=set(sub)for q in rec['event_domain']);assert e<=1;key=vec,e
                    if key not in items or sub<items[key]:items[key]=sub
        assert rec['items']==[dict(vector=list(v),events=e,clique=list(items[v,e]))for v,e in sorted(items)]
        counts.update(maximal_cliques=len(actual),event_vectors=len(items))
    columns=[];summaries=[]
    for di,key in enumerate(('face_multipliers','restricted_face_multipliers')):
        assert data[key]==base[key];alpha={i:F(x)for i,x in data[key]};prior=base['repairs'][di];mu={v:F(x)for v,x in prior['vertex_multipliers']}
        c=data['joint_certificates'][di];weights=[list(map(F,w))for w in c['weights']];a=list(map(F,c['source_costs']));z=list(map(F,c['event_coefficients']));delta=list(map(F,c['increases']));X=F(c['event_mass'])
        assert len(weights)==len(a)==len(z)==len(delta)==len(c['mixtures'])==len(records)==2
        assert all(0<=x<=1 for x in a)and all(x>=0 for x in delta)and X>=0 and F(c['quota_cost'])==0 and sum(z,F())<=0
        load=F();coefficient=Counter();new=dict(mu)
        for rec,w,aa,zz,inc,mix in zip(records,weights,a,z,delta,c['mixtures']):
            d=[sum((alpha.get(f,F())*m for f,m in g['positive_face_coefficients']),F())for g in rec['corner_orbits']]
            assert list(map(str,d))==rec['demands'][di]and len(w)==len(d)and all(x>=0 for x in w)
            assert all(sum(x*y for x,y in zip(w,it['vector']))<=aa+zz*it['events']for it in rec['items'])
            terms=[(j,F(x))for j,x in mix];assert len(terms)==len({j for j,x in terms})and all(0<=j<len(rec['items'])and x>0 for j,x in terms)
            assert all(sum((x*rec['items'][j]['vector'][i]for j,x in terms),F())>=d[i]for i in range(len(d)))
            assert sum((x*rec['items'][j]['events']for j,x in terms),F())==X
            assert sum((x for j,x in terms),F())<=mu[rec['vertex']]+inc
            load+=sum((x*y for x,y in zip(w,d)),F());new[rec['vertex']]+=inc
            for weight,g in zip(w,rec['corner_orbits']):
                for f,m in g['positive_face_coefficients']:coefficient[f]+=weight*m
            counts.update(charge_inequalities=len(rec['items']),mixture_terms=len(terms))
        gain=load-sum((aa*mu[r['vertex']]for aa,r in zip(a,records)),F())
        assert gain==sum(delta,F())==F(c['gain'])==F(934391,2000000000)
        col=dict(positive_face_coefficients=[[f,str(x)]for f,x in sorted(coefficient.items())if x],vertex_capacities=[[r['vertex'],str(aa)]for r,aa in zip(records,a)if aa])
        if col not in columns:columns.append(col)
        repair=data['repairs'][di];den=sum(new.values(),F());num=F(prior['objective_numerator'])
        assert repair['vertex_multipliers']==[[v,str(x)]for v,x in sorted(new.items())if x]
        assert repair['vertex_updates']==[[v,str(mu[v]),str(x)]for v,x in sorted(new.items())if x>mu[v]]
        assert all(new[v]>=x for v,x in mu.items())and F(repair['normalizer'])==den and F(repair['objective_numerator'])==num and F(repair['bound'])==num/den==F(c['bound'])
        for col in data['columns']:assert sum((F(x)*alpha.get(f,F())for f,x in col['positive_face_coefficients']),F())<=sum((F(x)*new[v]for v,x in col['vertex_capacities']),F())
        summaries.append(dict(gain=str(gain),event_mass=str(X),bound=str(num/den)))
    assert columns==data['columns']and len(columns)==2 and counts['physical_saved_maps']==10
    assert counts['maximal_cliques']==29 and counts['event_vectors']==23 and counts['charge_inequalities']==46
    for p in [ROOT/f'scripts/{TAG}.py',Path(__file__)]:ast.parse(p.read_text())
    md=ROOT/f'results/{TAG}.md';links=re.findall(r'\]\(([^)]+)\)',md.read_text())
    for ref in links:assert(md.parent/ref.split('#')[0]).exists(),ref
    print(json.dumps(dict(status='PASS_NO_PRODUCER_RAW39_SHARED_EDGE_SAVED_CHECK',counts=dict(counts),repairs=summaries,
        source_bindings=len(data['sources']),links=len(links),json_sha256=digest(path),checker_sha256=digest(Path(__file__)),
        scope='Independent saved Gram/support/inverse/frame reconstruction,15 event-star gate tests from saved native contacts, graph/subset enumeration, source-map consistency, exact charge/cover/column aggregation and both repairs. No historical provenance/owner or complete embedding recensus, no proof of actual face existence or global realization.'),sort_keys=True))

if __name__=='__main__':run()
