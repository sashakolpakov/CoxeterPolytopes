#!/usr/bin/env python3
"""Stdlib-only complete raw76 family transport/field/charge checker.

Independently implements ordered Q(sqrt5) and square-root fields. Historical
owner reconstruction and catalogue completeness remain bound producer inputs.
"""
import ast,json,re
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw76_extra_vertex_complete_family_20261018'
PILOT='fresh_raw76_forced_extra_vertex_family_20261018'
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

def native_projection(c,d,j,field):
    out=field.scalar(0)
    for i,e in c:
        for k,f in d:
            # Two negative native contacts have positive product.
            p=field.mul(field.sqrt(squared(e)),field.sqrt(squared(f)))
            out=[a+j[i][k]*b for a,b in zip(out,p)]
    return out

def wall_pair(c,d,j):
    rad=[squared(e)for _,e in c+d if squared(e)!=0]
    f=Field(rad);cc=native_projection(c,c,j,f);dd=native_projection(d,d,j,f)
    assert all(x==0 for x in cc[1:]+dd[1:])
    hs=[cc[0]-2,dd[0]-2];assert all(h.sign()>0 for h in hs)
    hp=hs[0]*hs[1];f=Field(rad+[hp]);p=native_projection(c,d,j,f)
    z=[a-b for a,b in zip(p,f.sqrt(hp))];zz=f.mul(z,z)
    s=f.sign(z);gap=s<0 and f.sign([z[0]+1]+z[1:])>0
    ink=all(x==0 for x in zz[1:])
    integral=ink and (2*zz[0].a).denominator==1 and (zz[0].a**2-5*zz[0].b**2).denominator==1
    code='positive'if s>0 else 'gap'if gap else 4 if not ink else 3 if not integral else None
    return code,hs,hp,f,p,z

def jhash(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())

def normalized(star,perm):
    co=[sorted([[perm[i],e]for i,e in c])for c in star['contacts']];order=sorted(range(2),key=lambda i:co[i])
    return dict(contacts=[co[i]for i in order],occupied_pairs=sorted(sorted(perm[i]for i in p)for p in star['occupied_pairs']),
        squared_timelike_heights=[star['squared_timelike_heights'][i]for i in order],native_omitted_product=star['native_omitted_product'])

def key(x):return json.dumps(x,sort_keys=True,separators=(',',':'))

def role(m):
    p=dict(m['active_embedding']);return p[0],p[4],tuple(p[i]for i in(1,3,6))

def relative(a,b):
    x,y,h=role(a);u,v,k=role(b);names={x:0,y:1}
    for i in(u,v):
        if i not in names:names[i]=len(names)
    return names[u],names[v],int(h!=k)

def allowed(a,b):
    x,y,h=role(a);u,v,k=role(b);return x==u and y!=v and h!=k

def finite_data(data):
    adj={int(i):{int(j):m for j,m in r.items()}for i,r in data.items()};assert set(adj)==set(range(13))
    assert all(adj[j][i]==m for i,r in adj.items()for j,m in r.items())
    unseen=set(adj);components=[]
    while unseen:
        todo=[min(unseen)];c=set()
        while todo:
            i=todo.pop()
            if i in c:continue
            c.add(i);todo.extend(adj[i])
        unseen-=c;components.append(c)
    active=set();g=[[Q(2)if i==j else Q()for j in range(13)]for i in range(13)]
    for c in components:
        if any(m not in(3,5)for i in c for m in adj[i].values()):continue
        active.update(c)
        for i in c:
            for j,m in adj[i].items():
                z=root(squared([m,'0','0']));assert z is not None;g[i][j]=-z
    # Non-K finite components are wholly passive for every checked covector.
    # The dummy identity blocks below are never used in any projection.
    assert positive(g);return adj,g,inverse(g),active,components

def check_marker(m,native,g,j,active):
    mp=dict(m['active_embedding']);assert set(mp)=={0,1,3,4,6}and len(set(mp.values()))==5
    assert set(mp.values())<=active
    assert all(g[mp[i]][mp[k]]==-root(squared(native[i][k]))for i,k in combinations(mp,2))
    assert all(g[i][k]==0 for i in mp.values()for k in set(range(13))-set(mp.values()))
    cov=sorted([[[mp[k],native[i][k]]for k in sorted(mp,key=lambda k:mp[k])if native[i][k][0]!=2]for i in(2,5)])
    assert cov==m['native']['contacts']and m['triangle_ports']==sorted([mp[1],mp[3]])
    code,hs,hp,field,p,z=wall_pair(cov[0],cov[1],j)
    assert code is None and all(x==0 for x in z)
    assert hs==[Q(*x)for x in m['native']['squared_timelike_heights']]and hp==Q(*m['height_product'])

def check_decisions(dec,expressions,nodes,anchor,j,counts):
    checked=[]
    for e in expressions:
        c,d=e['contacts'];code,hs,hp,_,_,_=wall_pair(c,d,j)
        assert code==e['code']and hs==[Q(*x)for x in e['heights']]and hp==Q(*e['height_product'])
        checked.append((c,d,code));counts['new_wall_expressions']+=1
    assert [x['pair']for x in dec]==[[i,len(nodes)]for i in range(len(nodes))]
    for x in dec:
        i,k=x['pair'];bad=x['exclusion']
        if bad is None:
            # Recheck every distinct cross wall for a genuinely retained pair.
            for c in nodes[i]['contacts']:
                for d in anchor['native']['contacts']:
                    if c!=d:assert wall_pair(c,d,j)[0]is None
            continue
        if bad['kind']=='actual_triangle_port':
            assert bad['ports']==anchor['triangle_ports']and bad['ports']in nodes[i]['occupied_pairs']
        else:
            assert bad['kind']=='physical_product';a,b=bad['wall_indices'];c,d,code=checked[bad['expression']]
            assert sorted([nodes[i]['contacts'][a],anchor['native']['contacts'][b]])==[c,d]
            assert code==bad['code']and code is not None
        counts['direct_pair_exclusions']+=1

def run():
    d=read(TAG);assert d['status']=='PASS_RAW76_COMPLETE_FAMILY_COUPLED_CHARGES'
    for p,h in d['sources'].items():assert digest(ROOT/p)==h,p
    pilot=read(PILOT);assert d['core_geometry']==pilot['core_geometry']and d['family']==pilot['family']
    native=d['core_geometry']['native_gram'];model=read('fresh_current327_full_model_20261018')
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']};mu={v:F(x)for v,x in dual['vertex_multipliers']}
    th=d['marker_theorem'];adj,g,j,active,comps=finite_data(th['finite_frame']);marks=th['markers']
    assert Counter(th['host']['target_forest'])=={'A1':7,'H3':2}and len(marks)==84 and th['anchor']==0
    for m in marks:check_marker(m,native,g,j,active)
    patterns={relative(marks[0],m)for m in marks[1:]};assert len(patterns)==13
    assert {tuple(r['relative'])for r in th['relative_patterns']}==patterns
    counts=Counter(relative_patterns=13,host_marker_maps=84)
    for r in th['relative_patterns']:
        m=marks[r['other_marker']];assert relative(marks[0],m)==tuple(r['relative'])and r['allowed']==allowed(marks[0],m)
        check_decisions([r['decision']],r['wall_expressions'],[m['native']],marks[0],j,counts)
        assert r['allowed']==(r['decision']['exclusion']is None)
    records=d['forest_records'];assert len(records)==235
    assert [(r['face'],r['vertex'])for r in records]==[(r['face'],r['target_vertex'])for r in d['family']]
    frozen={r['vertex']:r for r in pilot['pilots']};family={r['face']:r for r in d['family']};byvertex={}
    for r in records:
        v=r['vertex'];byvertex[v]=r;nodes=r['old_nodes'];groups=r['old_groups'];marks=r['markers'];n=len(nodes);m=len(marks)
        adj,g,j,active,components=finite_data(r['finite_frame'])
        assert all(i in active for s in nodes for c in s['contacts']for i,e in c)
        owner={i:z['index']for z in groups for i in z['literal_stars']};assert set(owner)==set(range(n))
        assert all(Counter(o['face']for o in z['origins'])==Counter(dict(z['positive_face_coefficients']))for z in groups)
        c=Counter(r['forest']);assert m==c['A1']*(c['A1']-1)*c['H3']==family[r['face']]['marker_orientations']
        assert len(r['transports'])==m and [x['index']for x in marks]==list(range(m))
        canon=[normalized(s,dict(zip(range(13),range(13))))for s in nodes]
        assert len({key(s['contacts'])for s in canon})==n
        assert len({key(s['native']['contacts'])for s in marks})==m
        for marker in marks:check_marker(marker,native,g,j,active)
        for k,t in enumerate(r['transports']):
            perm=dict(t['frame_permutation']);on=t['ordinary_permutation'];mn=t['marker_permutation']
            assert set(perm)==set(perm.values())==set(range(13))
            assert all(adj[i].get(z,2)==adj[perm[i]].get(perm[z],2)for i,z in combinations(range(13),2))
            assert sorted(on)==list(range(n))and sorted(mn)==list(range(m))and mn[0]==k
            for i,z in enumerate(on):assert normalized(nodes[i],perm)==canon[z]and owner[i]==owner[z]
            for i,z in enumerate(mn):assert normalized(marks[i]['native'],perm)==marks[z]['native']
            counts['full_frame_transports']+=1
        if v in frozen:
            old=frozen[v];assert r['reused_frozen_pilot']and nodes==old['old_nodes']and groups==old['old_groups']and marks==old['markers']
            assert r['old_graph_edges']==old['old_graph_edges']
            assert r['anchor_decisions']==[x for x in old['new_pair_decisions']if x['pair'][1]==n]
            assert r['wall_expressions']==old['wall_expressions'];counts['reused_anchor_tests']+=n
        else:
            assert not r['reused_frozen_pilot'];check_decisions(r['anchor_decisions'],r['wall_expressions'],nodes,marks[0],j,counts)
            counts['new_anchor_tests']+=n
        anchor={x['pair'][0]for x in r['anchor_decisions']if x['exclusion']is None}
        cross=[[i,k]for k,t in enumerate(r['transports'])for i in sorted(t['ordinary_permutation'][i]for i in anchor)]
        assert cross==r['marker_ordinary_edges']
        edges=[list(z)for z in combinations(range(m),2)if allowed(*(marks[i]for i in z))]
        assert edges==r['marker_edges']
        assert all(relative(marks[i],marks[z])in patterns for i,z in combinations(range(m),2))
        cap=min(c['A1']-1,c['H3']);assert r['marker_capacity']==cap
        cl=r['attaining_marker_clique'];assert len(cl)==len(set(cl))==cap and all(list(sorted(z))in edges for z in combinations(cl,2))
        counts['cross_positions']+=n*m;counts['marker_pair_positions']+=m*(m-1)//2
        counts['retained_cross_pairs']+=len(cross);counts['retained_marker_pairs']+=len(edges)
    # Obtain unchanged complete old count vectors by their frozen geometry owner.
    base=read('fresh_current324_new5073_local_covers_20261018');rows={r['vertex']:r for r in base['forest_records']}
    docs={};vectors={}
    for v,r in rows.items():
        tag=r['geometry_source']
        if tag not in docs:docs[tag]={x['vertex']:x for x in read(tag)['forest_records']}
        rr=docs[tag][v];vectors[v]=rr['maximal_count_vectors']
        assert jhash(vectors[v])==r['count_vectors_sha256']
    for tag in('fresh_new5073_raw76_high_mass_local_covers_20261018','fresh_new5073_raw76_third_class_local_covers_20261018'):
        for r in read(tag)['forest_records']:vectors[r['vertex']]=r['maximal_count_vectors']
    exported=d['exported'];newmu=dict(mu);loads={};columns=exported['columns'];assert len(columns)==242
    for col in columns:
        r=byvertex[col['vertex']];groups=r['old_groups'];owner={i:z['index']for z in groups for i in z['literal_stars']}
        original=list(map(F,col['old_orbit_weights']));w=list(map(F,col['retained_orbit_weights']));assert len(original)==len(w)==len(groups)
        assert min(original,default=F())>=0
        vv=vectors[r['vertex']]if groups else [[]]
        assert len(vv)==col['old_vector_count']and jhash(vv)==col['old_vectors_sha256']
        assert all(len(x)==len(w)and sum((a*b for a,b in zip(original,x)),F())<=1 for x in vv)
        counts['old_charge_vector_inequalities']+=len(vv)
        omitted=sorted({owner[i]for i,k in r['marker_ordinary_edges']})
        assert omitted==col['omitted_orbits']and w==[F()if i in omitted else z for i,z in enumerate(original)]
        assert not omitted or(r['vertex']==1304 and omitted==[0]and groups[0]['positive_face_coefficients']==[[60808,1]])
        cf=Counter()
        for a,z in zip(w,groups):
            for f,n in z['positive_face_coefficients']:cf[f]+=a*n
        oldload=sum((alpha.get(f,F())*a for f,a in cf.items()),F())
        weight=F(1,r['marker_capacity']);assert weight==F(col['marker_weight'])
        cf[r['face']]+=weight
        assert [[f,str(a)]for f,a in sorted(cf.items())if a]==col['positive_face_coefficients']
        load=sum((alpha.get(f,F())*a for f,a in cf.items()),F())
        assert oldload==F(col['old_charge_load'])and load==F(col['load'])and col['vertex_cost']=='1'
        assert F(col['baseline_cost'])==mu.get(r['vertex'],F())and F(col['violation'])==load-mu.get(r['vertex'],F())
        newmu[r['vertex']]=max(newmu.get(r['vertex'],F()),load)
        loads.setdefault(r['vertex'],[]).append(col)
    assert set(loads)==set(byvertex)
    for q in exported['target_allowances']:
        v=q['vertex'];assert F(q['old_vertex_mass'])==mu.get(v,F())and F(q['new_vertex_mass'])==newmu[v]
        assert F(q['increase'])==newmu[v]-mu.get(v,F())
    rep=exported['repair'];assert rep['face_multipliers']==dual['face_multipliers']
    assert rep['vertex_multipliers']==[[v,str(x)]for v,x in sorted(newmu.items())if x]
    num=F(dual['objective_numerator']);den=sum(newmu.values(),F());assert num==F(rep['objective_numerator'])and den==F(rep['normalizer'])and num/den==F(rep['bound'])
    assert len(exported['preserved_multi_owner_costs'])==14
    for out,col in zip(exported['preserved_multi_owner_costs'],model['prepared']['inherited_multi_vertex_columns']):
        assert out['variable']==col['variable']and out['all_vertex_costs']==col['integer_vertex_capacities']
        load=sum((alpha.get(f,F())*n for f,n in col['integer_face_coefficients']),F())
        old=sum((mu.get(v,F())*n for v,n in col['integer_vertex_capacities']),F());new=sum((newmu.get(v,F())*n for v,n in col['integer_vertex_capacities']),F())
        assert load==F(out['load'])<=old==F(out['old_cost'])<=new==F(out['new_cost'])
    assert counts['new_anchor_tests']==10181 and counts['reused_anchor_tests']==464
    assert counts['cross_positions']==267150 and counts['marker_pair_positions']==39966
    assert counts['retained_cross_pairs']==6 and counts['retained_marker_pairs']==750
    assert sum(F(c['violation'])>0 for c in columns)==10
    assert sum(newmu.get(v,F())>x for v,x in {v:mu.get(v,F())for v in byvertex}.items())==8
    assert d['optimizer_invocations']==0
    counts.update(classes=235,columns=242,multi_owner_columns=14)
    print('PASS_RAW76_COMPLETE_FAMILY_SAVED_CHECK',dict(counts))

if __name__=='__main__':run()

