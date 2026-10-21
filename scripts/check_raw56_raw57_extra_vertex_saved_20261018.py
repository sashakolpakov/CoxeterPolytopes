#!/usr/bin/env python3
"""Stdlib saved checker for raw56/57 geometry, maps, gates and charges.

Uses a positively scaled Q(sqrt5) inverse on contacted components, with
untouched orthogonal blocks replaced by harmless diagonal blocks. This
is independent of the producer's unscaled I2(10) formula. Historical
owner/graph census and catalogue completeness remain inherited.
"""
import ast,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw56_raw57_extra_vertex_family_20261018'
MODEL='fresh_raw76_complete_family_full_model_20261018'
LOCAL='fresh_current324_new5073_local_covers_20261018'
PATCH='fresh_new5073_raw76_high_mass_local_covers_20261018'
THIRD='fresh_new5073_raw76_third_class_local_covers_20261018'
OLDMARKER='fresh_raw76_extra_vertex_complete_family_20261018'
IDS=(0,1,3,5)
def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
def digest(p):
 h=sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def key(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def jhash(x):return sha256(key(x).encode()).hexdigest()
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
    rows={int(i):{int(k):m for k,m in row.items()}for i,row in record['finite_frame'].items()};assert set(rows)==set(range(13))
    g=[[Q(2)if i==j else Q()for j in range(13)]for i in range(13)]
    for i,row in rows.items():
        for j,m in row.items():
            assert rows[j][i]==m;z=root(squared([m,'0','0']));assert z is not None;g[i][j]=-z
    assert positive(g)
    return g,inverse(g)

def residual(g,fixed,sides):
    j=inverse(principal(g,fixed))
    return [[g[a][b]-sum((g[a][x]*j[i][k]*g[y][b]
        for i,x in enumerate(fixed)for k,y in enumerate(fixed)),Q())for b in sides]for a in sides]

class Field:
    """Independent multiquadratic construction, including the height root."""
    def __init__(self,radicands):
        self.basis=[];self.powers=[Q(1)]
        for r in radicands:
            assert r.sign()>0
            if any(root(r/s)is not None for s in self.powers):continue
            self.basis.append(r);self.powers +=[s*r for s in self.powers]
        assert len(self.basis)<=3
        self.n=len(self.powers)
    def scalar(self,q):return [q if isinstance(q,Q)else Q(q)]+[Q()for _ in range(self.n-1)]
    def sqrt(self,q):
        if q==0:return self.scalar(0)
        for i,s in enumerate(self.powers):
            z=root(q/s)
            if z is not None:
                a=self.scalar(0);a[i]=z;return a
        raise AssertionError('missing squareclass')
    def mul(self,a,b):
        c=self.scalar(0)
        for i,x in enumerate(a):
            for j,y in enumerate(b):c[i^j]=c[i^j]+x*y*self.powers[i&j]
        return c
    def sign(self,a):
        def rec(z,k):
            if k==0:return z[0].sign()
            m=1<<(k-1);x=z[:m];y=z[m:];sx,sy=rec(x,k-1),rec(y,k-1)
            if not sx:return sy
            if not sy or sx==sy:return sx
            xx=self.mul(x+[Q()]* (self.n-m),x+[Q()]* (self.n-m))[:m]
            yy=self.mul(y+[Q()]* (self.n-m),y+[Q()]* (self.n-m))[:m]
            norm=[a-self.basis[k-1]*b for a,b in zip(xx,yy)]
            return sx*rec(norm,k-1)
        return rec(a,len(self.basis))


TAU=Q('5/2','1/2')
def frame(data):
    adj={int(i):{int(j):m for j,m in r.items()}for i,r in data.items()}
    assert set(adj)==set(range(13))and all(adj[j][i]==m for i,r in adj.items()for j,m in r.items())
    unseen=set(adj);components=[]
    while unseen:
        todo=[min(unseen)];seen=set()
        while todo:
            i=todo.pop()
            if i in seen:continue
            seen.add(i);todo.extend(adj[i])
        unseen-=seen;components.append(seen)
    scales=[Q(1)]*13;g=[[Q(2)if i==j else Q()for j in range(13)]for i in range(13)];active=set()
    for c in components:
        labels={m for i in c for m in adj[i].values()}
        if labels<={3,5}:
            active.update(c)
            for i in c:
                for j,m in adj[i].items():g[i][j]=-root(squared([m,'0','0']))
        elif len(c)==2 and labels=={10}:
            i,j=sorted(c);active.update(c);scales[i]=TAU;g[i][i]=2*TAU;g[i][j]=g[j][i]=-TAU
    assert all(g[i][j]==g[j][i]for i,j in product(range(13),repeat=2))and positive(g)
    return adj,scales,inverse(g),active,components


def project(c,d,j,scales,f):
    z=f.scalar(0)
    for i,e in c:
        for k,g in d:
            p=f.mul(f.sqrt(squared(e)*scales[i]),f.sqrt(squared(g)*scales[k]))
            z=[a+j[i][k]*b for a,b in zip(z,p)]
    return z


def physical(c,d,j,scales):
    rad=[squared(e)*scales[i]for i,e in c+d if squared(e)!=0];f=Field(rad)
    cc=project(c,c,j,scales,f);dd=project(d,d,j,scales,f)
    assert all(x==0 for x in cc[1:]+dd[1:]);hs=[cc[0]-2,dd[0]-2]
    assert all(h.sign()>0 for h in hs);hp=hs[0]*hs[1];f=Field(rad+[hp])
    p=project(c,d,j,scales,f);z=[a-b for a,b in zip(p,f.sqrt(hp))];zz=f.mul(z,z)
    sign=f.sign(z);gap=sign<0 and f.sign([z[0]+1]+z[1:])>0
    ink=all(x==0 for x in zz[1:])
    integral=ink and(2*zz[0].a).denominator==1 and(zz[0].a**2-5*zz[0].b**2).denominator==1
    code='positive'if sign>0 else'gap'if gap else 4 if not ink else 3 if not integral else None
    return code,hs,hp,f,p,z


def expected_maps(raw,adj,components):
    tens=[sorted(c)for c in components if len(c)==2 and all(list(adj[i].values())==[10]for i in c)]
    if raw==56:
        ones=sorted(i for i in adj if not adj[i]);roles=list(permutations(ones,2))
    else:
        twos=[sorted(c)for c in components if len(c)==2 and all(list(adj[i].values())==[3]for i in c)]
        roles=[p for c in twos for p in(c,c[::-1])]
    return {tuple(sorted({0:x,5:y,1:u,3:v}.items()))for x,y in roles for c in tens for u,v in(c,c[::-1])}


def normalize(s,perm):
    co=[sorted([[perm[i],e]for i,e in c])for c in s['contacts']];order=sorted(range(2),key=lambda i:co[i])
    return dict(contacts=[co[i]for i in order],occupied_pairs=sorted(sorted(perm[i]for i in p)for p in s['occupied_pairs']),
        squared_timelike_heights=[s['squared_timelike_heights'][i]for i in order],native_omitted_product=s['native_omitted_product'])


def marker(m,raw,native,adj,scales,j,active):
    mp=dict(m['active_embedding']);assert set(mp)==set(IDS)and len(set(mp.values()))==4
    assert all(adj[mp[i]].get(mp[k],2)==native[i][k][0]for i,k in combinations(IDS,2))
    assert all(not(set(adj[i])-set(mp.values()))for i in mp.values())
    cov=[sorted([[mp[k],native[i][k]]for k in IDS if native[i][k][0]!=2])for i in(2,4)]
    cov.sort();assert cov==m['native']['contacts']and cov[0]!=cov[1]
    assert all(i in active for c in cov for i,e in c)
    code,hs,hp,f,p,z=physical(cov[0],cov[1],j,scales)
    assert all(x==0 for x in z)and native[2][4][0]==2 and code is None
    assert hs==[Q(*h)for h in m['native']['squared_timelike_heights']]
    assert hp==Q(*m['height_product'])and m['triangle_ports']==sorted([mp[1],mp[3]])
    assert m['native']['native_omitted_product']==native[2][4]and m['native']['occupied_pairs']==[]


def pair_result(a,b,result,j,scales):
    if result.get('triangle_port_exclusion')is not None:
        assert not result['allowed']and result['triangle_port_exclusion']in a['occupied_pairs'];return 0
    expected=[];bad=None
    for i,c in enumerate(a['contacts']):
        for k,d in enumerate(b['contacts']):
            if c==d:continue
            code,hs,hp,f,p,z=physical(c,d,j,scales);expected.append(([i,k],c,d,code,hs,hp))
            if code is not None:bad=code;break
        if bad is not None:break
    assert len(result['expressions'])==len(expected)
    for rec,(walls,c,d,code,hs,hp)in zip(result['expressions'],expected):
        assert rec['walls']==walls and rec['contacts']==[c,d]and rec['code']==code
        assert [Q(*h)for h in rec['heights']]==hs and Q(*rec['height_product'])==hp
    assert result['allowed']==(bad is None);return len(expected)


def relative(ra,a,rb,b):
    a=dict(a['active_embedding']);b=dict(b['active_embedding'])
    ten=0 if(a[1],a[3])==(b[1],b[3])else 1 if(a[1],a[3])==(b[3],b[1])else 2
    if ra==rb==56:
        names={a[0]:0,a[5]:1}
        for i in(b[0],b[5]):
            if i not in names:names[i]=len(names)
        return(names[b[0]],names[b[5]],ten)
    if ra==rb==57:
        two=0 if(a[0],a[5])==(b[0],b[5])else 1 if(a[0],a[5])==(b[5],b[0])else 2
        return(two,ten)
    assert(ra,rb)==(56,57);return(ten,)


def permitted(ra,rb,p):
    return False if ra!=rb else tuple(p)==(0,2)if ra==57 else tuple(p)in((0,1,2),(2,3,0),(2,3,2))


def check():
    data=read(TAG);assert data['status']=='PASS_RAW56_RAW57_COMPLETE_EXTRA_VERTEX_FAMILY'
    for p,h in data['sources'].items():assert digest(ROOT/p)==h,p
    ast.parse((ROOT/f'scripts/{TAG}.py').read_text());native={}
    for ri,s in data['geometries'].items():
        ri=int(ri);g=mat(s['scaled_gram']);ng=s['native_gram'];native[ri]=ng;scales=[Q(*x)for x in s['squared_scales']]
        assert len(g)==len(ng)==6 and all(len(row)==6 for row in g+ng)
        assert all(g[i][k]==g[k][i]and ng[i][k]==ng[k][i]for i,k in product(range(6),repeat=2))
        for i in range(6):
            assert g[i][i]==2*scales[i]and scales[i].sign()>0
            for k in range(i):assert g[i][k].sign()<=0 and g[i][k]*g[i][k]==squared(ng[i][k])*scales[i]*scales[k]
        support={tuple(r['omitted']):r for r in s['support_records']};assert set(support)==set(combinations(range(6),2))
        pd=[]
        for omitted,r in support.items():
            ids=sorted(set(range(6))-set(omitted));assert ids==r['support']
            p=positive(principal(g,ids));assert p==(r['inertia']==[4,0,0])
            if p:pd.append(omitted)
        assert set(pd)=={(0,1),(0,3),(1,2),(2,3),(2,4)}
        assert s['extra_support']==list(IDS)and s['extra_omission']==[2,4]
        rr=residual(g,[0,5],[1,3,4]);assert rr==mat(s['triangle_residual'])
        inv=inverse(rr);assert inv==mat(s['triangle_inverse'])and all(x.sign()<0 for row in inv for x in row)
        assert s['triangle_supports']==[sorted([0,5]+list(q))for q in combinations([1,3,4],2)]
        assert [list(p)for p in permutations(range(6))if all(ng[i][k]==ng[p[i]][p[k]]for i,k in combinations(range(6),2))]==s['native_automorphisms']==[list(range(6))]
    assert not any(all(native[56][i][j]==native[57][p[i]][p[j]]for i,j in combinations(range(6),2))for p in permutations(range(6)))
    cases=0;expressions=0
    for th in data['marker_theorem']['tests']:
        ra,rb=th['raw_families'];adj,scales,j,active,components=frame(th['finite_frame'])
        marker(th['anchor'],ra,native[ra],adj,scales,j,active)
        seen=set()
        for r in th['records']:
            marker(r['other_marker'],rb,native[rb],adj,scales,j,active)
            p=relative(ra,th['anchor'],rb,r['other_marker']);assert list(p)==r['pattern']and p not in seen;seen.add(p)
            expressions+=pair_result(th['anchor']['native'],r['other_marker']['native'],r['result'],j,scales)
            assert r['result']['allowed']==permitted(ra,rb,p);cases+=1
        # Enumerate relative patterns directly from all component-role maps.
        full=[]
        for mp in expected_maps(rb,adj,components):
            m={'active_embedding':mp}
            if ra==rb and dict(mp)==dict(th['anchor']['active_embedding']):continue
            full.append(relative(ra,th['anchor'],rb,m))
        assert seen==set(full)
    assert cases==31
    rows={r['vertex']:r for r in data['forest_records']};assert len(data['forest_records'])==len(rows)==1212
    family={(r['target_vertex'],r['raw_family']):r for r in data['family']};assert len(data['family'])==len(family)==1408
    total=Counter();all_sources=[];compatible_orbits={};vectors={};coefficients={}
    for tag in(LOCAL,PATCH,THIRD):
        for r in read(tag)['forest_records']:
            if 'maximal_count_vectors'in r:vectors[r['vertex']]=r['maximal_count_vectors']
    # The base600 count-vector domains are bound historical inputs. Import
    # their serialized rows only, never producer geometry or optimization.
    local=read(LOCAL)
    cache={}
    for r in local['forest_records']:
        if r['vertex']in vectors:continue
        tag=r['geometry_source']
        if tag not in cache:cache[tag]={z['vertex']:z for z in read(tag)['forest_records']}
        rr=cache[tag][r['vertex']];vectors[r['vertex']]=rr['maximal_count_vectors']
    for v,r in rows.items():
        adj,scales,j,active,components=frame(r['finite_frame']);nodes=r['old_nodes'];groups=r['old_groups'];n=len(nodes)
        owner={i:g['index']for g in groups for i in g['literal_stars']};assert set(owner)==set(range(n))
        assert all(i in active for s in nodes for c in s['contacts']for i,e in c)
        for g in groups:assert Counter(o['face']for o in g['origins'])==Counter(dict(g['positive_face_coefficients']))
        oldkeys={key(normalize(s,dict(zip(range(13),range(13))))['contacts']):i for i,s in enumerate(nodes)};assert len(oldkeys)==n
        allkeys=set(oldkeys);omitted=set();mcount=0;crosscount=0;unk=0;tests=0
        for f in r['marker_families']:
            ri=f['raw_family'];fr=family[v,ri];assert fr['face']==f['face']and fr['actual_extra_occurrences']==1
            ms=f['markers'];mcount+=len(ms);assert len(ms)==fr['marker_orientations']
            maps={tuple(map(tuple,m['active_embedding']))for m in ms};assert maps==expected_maps(ri,adj,components)
            keys={key(m['native']['contacts']):k for k,m in enumerate(ms)};assert len(keys)==len(ms)and not allkeys&set(keys);allkeys.update(keys)
            for k,m in enumerate(ms):assert m['index']==k;marker(m,ri,native[ri],adj,scales,j,active)
            assert len(f['anchor_decisions'])==n and len(f['transports'])==len(ms)
            allowed=[]
            for i,z in enumerate(f['anchor_decisions']):
                assert z['ordinary']==i
                if 'triangle_port_exclusion'in z['result']:assert z['result']['triangle_port_exclusion']==ms[0]['triangle_ports']
                expressions+=pair_result(nodes[i],ms[0]['native'],z['result'],j,scales);tests+=1
                if z['result']['allowed']:allowed.append(i);unk+=1
            edges=[]
            for k,t in enumerate(f['transports']):
                perm=dict(t['frame_permutation']);assert set(perm)==set(perm.values())==set(range(13))
                assert all(adj[i].get(z,2)==adj[perm[i]].get(perm[z],2)for i,z in combinations(range(13),2))
                om=t['ordinary_permutation'];mm=t['marker_permutation']
                assert len(om)==n and set(om)==set(range(n))and len(mm)==len(ms)and set(mm)==set(range(len(ms)))and mm[0]==k
                for i,s in enumerate(nodes):assert normalize(s,perm)==normalize(nodes[om[i]],dict(zip(range(13),range(13))))and owner[i]==owner[om[i]]
                for i,m in enumerate(ms):assert normalize(m['native'],perm)==ms[mm[i]]['native']
                edges.extend([om[i],k]for i in allowed)
            assert sorted(edges)==f['marker_ordinary_edges'];crosscount+=len(edges);omitted.update(owner[i]for i,k in edges)
            c=Counter(r['forest']);cap=(c['A1']//2)*c['I2(10)']if ri==56 else c['I2(10)']
            cl=f['attaining_marker_clique'];assert f['marker_capacity']==cap==len(cl)==len(set(cl))
            assert all(permitted(ri,ri,relative(ri,ms[i],ri,ms[k]))for i,k in combinations(cl,2))
            all_sources.append((v,ri,f['face']))
        counts=dict(ordinary_stars=n,markers=mcount,anchor_tests=tests,cross_pairs=n*mcount,
            marker_pairs=mcount*(mcount-1)//2,retained_cross_pairs=crosscount,anchor_unknowns=unk)
        assert r['counts']==counts;total.update(counts);compatible_orbits[v]=omitted;coefficients[v]=[g['positive_face_coefficients']for g in groups]
    assert {(v,ri)for v,ri,f in all_sources}==set(family)
    for k,x in total.items():assert data['counts'][k]==x
    model=read(MODEL);dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={i:F(x)for i,x in dual['face_multipliers']};oldmu={i:F(x)for i,x in dual['vertex_multipliers']};mu=dict(oldmu)
    charges=data['exported']['columns'];assert len(charges)==len({c['vertex']for c in charges})==1212
    assert {c['vertex']for c in charges}==set(rows)
    allowances=data['exported']['target_allowances']
    assert len(allowances)==len({q['vertex']for q in allowances})==1212
    assert {q['vertex']for q in allowances}==set(rows)
    inequalities=0
    for c in charges:
        v=c['vertex'];r=rows[v];w=list(map(F,c['old_orbit_weights']));kept=list(map(F,c['retained_orbit_weights']))
        assert len(w)==len(kept)==len(r['old_groups'])and all(x>=0 for x in w)
        vv=vectors[v]if w else[[]];assert len(vv)==c['old_vector_count']and jhash(vv)==c['old_vectors_sha256']
        for x in vv:assert len(x)==len(w)and sum((a*b for a,b in zip(w,x)),F())<=1;inequalities+=1
        assert set(c['omitted_orbits'])==compatible_orbits[v]
        assert kept==[F()if i in compatible_orbits[v]else x for i,x in enumerate(w)]
        cf=Counter()
        for a,terms in zip(kept,coefficients[v]):
            for f,m in terms:cf[f]+=a*m
        assert c['marker_terms']==[[z['raw_family'],z['face'],str(F(1,z['marker_capacity']))]for z in r['marker_families']]
        for ri,f,x in c['marker_terms']:cf[f]+=F(x)
        cf={f:x for f,x in cf.items()if x};assert cf=={f:F(x)for f,x in c['positive_face_coefficients']}and c['vertex_cost']=='1'
        load=sum((alpha.get(f,F())*x for f,x in cf.items()),F());assert load==F(c['load'])and F(c['violation'])==load-oldmu.get(v,F())
        mu[v]=max(mu.get(v,F()),load)
    export=data['exported'];repair=export['repair']
    assert repair['face_multipliers']==dual['face_multipliers']and {v:F(x)for v,x in repair['vertex_multipliers']}=={v:x for v,x in mu.items()if x}
    for q in export['target_allowances']:
        v=q['vertex'];assert F(q['old_vertex_mass'])==oldmu.get(v,F())and F(q['new_vertex_mass'])==mu.get(v,F())
        assert F(q['increase'])==mu.get(v,F())-oldmu.get(v,F())and set(q['ordinary_omitted_orbits'])==compatible_orbits[v]
    assert F(repair['objective_numerator'])==F(dual['objective_numerator'])and F(repair['normalizer'])==sum(mu.values(),F())
    assert F(repair['bound'])==F(repair['objective_numerator'])/F(repair['normalizer'])and F(repair['bound'])>12
    oldcols=read(OLDMARKER)['exported']['columns'];assert len(oldcols)==len(export['preserved_raw76_columns'])==242
    for c,r in zip(oldcols,export['preserved_raw76_columns']):
        v=c['vertex'];load=sum((alpha.get(f,F())*F(x)for f,x in c['positive_face_coefficients']),F())
        assert r['vertex']==v and r['source_face']==c['source_face']and F(r['load'])==load<=oldmu.get(v,F())<=mu.get(v,F())
    mult=model['prepared']['inherited_multi_vertex_columns'];assert len(mult)==len(export['preserved_multi_owner_costs'])==14
    for c,r in zip(mult,export['preserved_multi_owner_costs']):
        load=sum((alpha.get(f,F())*F(x)for f,x in c['integer_face_coefficients']),F())
        old=sum((oldmu.get(v,F())*F(x)for v,x in c['integer_vertex_capacities']),F());new=sum((mu.get(v,F())*F(x)for v,x in c['integer_vertex_capacities']),F())
        assert r['variable']==c['variable']and r['all_vertex_costs']==c['integer_vertex_capacities']and load<=old<=new
        assert list(map(F,(r['load'],r['old_cost'],r['new_cost'])))==[load,old,new]
    print('PASS_RAW56_RAW57_SAVED_CHECK',dict(relative_patterns=cases,marker_maps=total['markers'],anchor_tests=total['anchor_tests'],
        checked_wall_expressions=expressions,ordinary_charge_inequalities=inequalities,columns=len(charges),preserved_raw76=242,multi_owner=14),flush=True)
    print('EXACT_EXPORTED_REPAIR',repair['bound'],float(F(repair['bound'])),flush=True)


if __name__=='__main__':check()

