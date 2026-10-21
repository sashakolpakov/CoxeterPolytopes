#!/usr/bin/env python3
"""Independent stdlib checker for the saved raw108 extension.

Reconstructs the new cross-product gates with ordered multiquadratic
arithmetic and a scaled contacted-block inverse; unused orthogonal blocks
are harmless diagonals after zero-contact guards. Rechecks every saved
count envelope, graph witness, charge and repair. The frozen geometry,
whole-orbit/owner and historical ordinary graph/vector completeness are
inputs, not recounted here. No producer imports or optimization.
"""
import ast,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw108_complete_marker_family_20261018'
NEW='fresh_raw56_raw57_extra_vertex_family_20261018'
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

JOINT='fresh_raw56_raw57_raw76_joint_marker_family_20261018'
PILOT='fresh_raw14_raw112_raw77_raw108_core_incidence_pilot_20261018'
NEWKINDS=(1084,1085)
SUPPORT={1084:(0,1,3,5,6,7),1085:(0,1,3,4,6,7),56:(0,1,3,5),57:(0,1,3,5),76:(0,1,3,4,6)}
OM={1084:(2,4),1085:(2,5),56:(2,4),57:(2,4),76:(2,5)}
ROLES={1084:((1,3,6,5),(0,7)),1085:((4,),(1,3,6),(0,7)),56:((0,5),(1,3)),57:((0,5),(1,3)),76:((0,4),(1,3,6))}

def normalize(s,p):
    co=[sorted([[p[i],e]for i,e in c])for c in s['contacts']];order=sorted(range(2),key=lambda i:co[i])
    return dict(contacts=[co[i]for i in order],occupied_pairs=sorted(sorted(p[i]for i in x)for x in s['occupied_pairs']),
        squared_timelike_heights=[s['squared_timelike_heights'][i]for i in order],native_omitted_product=s['native_omitted_product'])

def parts(adj):
    unseen=set(adj);out=[]
    while unseen:
        seen=set();todo=[min(unseen)]
        while todo:
            i=todo.pop()
            if i in seen:continue
            seen.add(i);todo.extend(adj[i])
        unseen-=seen;out.append(tuple(sorted(seen)))
    return out

def embeddings(kind,ng,adj):
    ids=SUPPORT[kind];small={i:{j:ng[i][j][0]for j in ids if i!=j and ng[i][j][0]!=2}for i in ids}
    blocks=parts(small);targets=parts(adj);choices=[]
    for block in blocks:
        cc=[]
        for target in targets:
            if len(target)!=len(block):continue
            for order in permutations(target):
                p=dict(zip(block,order))
                if all(small[i].get(j,2)==adj[p[i]].get(p[j],2)for i,j in combinations(block,2)):
                    cc.append(p)
        choices.append(cc)
    out=set()
    for choiceset in product(*choices):
        p={i:j for m in choiceset for i,j in m.items()}
        if len(set(p.values()))==len(ids):out.add(tuple(sorted(p.items())))
    return out

def signature(ka,a,kb,b):
    pa,pb=dict(a['active_embedding']),dict(b['active_embedding'])
    aa=[tuple(pa[i]for i in r)for r in ROLES[ka]];bb=[tuple(pb[i]for i in r)for r in ROLES[kb]]
    if kb in NEWKINDS:
        t=0 if aa[-1]==bb[-1]else 1 if aa[-1]==bb[-1][::-1]else 2
        if ka==kb==1084:return(int(aa[0]!=bb[0]),t)
        if ka==kb==1085:return(int(aa[0]!=bb[0]),int(aa[1]!=bb[1]),t)
        assert(ka,kb)==(1084,1085);return(t,)
    if ka==1084 or kb==57:return()
    i=bb[0].index(aa[0][0])if aa[0][0]in bb[0]else 2
    return(i,)if kb==56 else(i,int(aa[1]!=bb[1]))

def compatible(ka,a,kb,b):
    if kb not in NEWKINDS:
        if kb==57:return False
        if ka==1084:return True
        p,q=dict(a['active_embedding']),dict(b['active_embedding'])
        return kb==56 and p[4]not in(q[0],q[5])
    if ka!=1085 or kb!=1085:return True
    p,q=dict(a['active_embedding']),dict(b['active_embedding'])
    return p[4]==q[4]and all(p[i]==q[i]for i in(1,3,6))

def check_marker(m,kind,ng,adj,scales,j,active):
    p=dict(m['active_embedding']);ids=SUPPORT[kind];om=OM[kind]
    assert set(p)==set(ids)and len(set(p.values()))==len(ids)
    assert all(adj[p[i]].get(p[k],2)==ng[i][k][0]for i,k in combinations(ids,2))
    assert all(not(set(adj[i])-set(p.values()))for i in p.values())
    cov=[sorted([[p[k],ng[i][k]]for k in ids if ng[i][k][0]!=2])for i in om];cov.sort()
    assert cov==m['native']['contacts']and cov[0]!=cov[1]and all(i in active for c in cov for i,e in c)
    code,hs,hp,f,projection,z=physical(cov[0],cov[1],j,scales)
    prescribed=f.sqrt(squared(ng[om[0]][om[1]]))
    assert all(a==-b for a,b in zip(z,prescribed))
    assert hs==[Q(*h)for h in m['native']['squared_timelike_heights']]and hp==Q(*m['height_product'])
    assert m['native']['native_omitted_product']==ng[om[0]][om[1]]and m['native']['occupied_pairs']==[]
    if kind in NEWKINDS:
        expected=sorted(sorted((p[i],p[k]))for i,k in combinations(sorted(set(ids)&{1,3,4,5,6}),2))
        assert m['triangle_port_pairs']==expected and len(expected)==6

def edge(ka,a,kb,b):
    if ka in NEWKINDS:return compatible(ka,a,kb,b)
    if kb in NEWKINDS:return compatible(kb,b,ka,a)
    if ka>kb:return edge(kb,b,ka,a)
    p,q=dict(a['active_embedding']),dict(b['active_embedding'])
    if ka==kb==76:return p[0]==q[0]and p[4]!=q[4]and tuple(p[i]for i in(1,3,6))!=tuple(q[i]for i in(1,3,6))
    if kb==76:return ka==56 and not({p[0],p[5]}&{q[0],q[4]})
    if ka!=kb:return False
    x,y=(p[0],p[5]),(q[0],q[5]);u,v=(p[1],p[3]),(q[1],q[3])
    if ka==57:return x==y and set(u)!=set(v)
    return set(u)!=set(v)if x==y else not(set(x)&set(y))and(set(u)!=set(v)or u==v)


def coefficient_vertices(boundary,have):
    order=(56,57,76,1084,1085);ids=[i for i,k in enumerate(order)if k in have and k!=57];n=len(ids)
    planes={(tuple(F(p[i])for i in ids),F(1))for p in boundary if any(p[i]for i in ids)}
    planes.update((tuple(F(int(i==j))for i in range(n)),F())for j in range(n));verts=set()
    for cc in combinations(sorted(planes),n):
        mat=[list(a)+[b]for a,b in cc];ok=True
        for i in range(n):
            pivot=next((k for k in range(i,n)if mat[k][i]),None)
            if pivot is None:ok=False;break
            mat[i],mat[pivot]=mat[pivot],mat[i];z=mat[i][i];mat[i]=[x/z for x in mat[i]]
            for k in range(n):
                if k==i:continue
                z=mat[k][i];mat[k]=[x-z*y for x,y in zip(mat[k],mat[i])]
        if not ok:continue
        x=tuple(row[-1]for row in mat)
        if min(x)>=0 and all(sum(F(p[i])*x[j]for j,i in enumerate(ids))<=1 for p in boundary):verts.add(x)
    maximal=[x for x in verts if not any(y!=x and all(a<=b for a,b in zip(x,y))for y in verts)];out=[]
    for x in maximal:
        w=[F()]*5
        for i,z in zip(ids,x):w[i]=z
        if 57 in have:w[1]=F(1,max(p[1]for p in boundary))
        out.append(tuple(w))
    return sorted(out)


def check_envelope(r,prior,q):
    fam={f['kind']:f['markers']for f in r['marker_families']}
    for old in prior:
        if 'markers'in old:fam[76]=old['markers']
        else:fam.update({f['raw_family']:f['markers']for f in old['marker_families']})
    c=Counter(r['forest']);a,h,t=c['A1'],c['H3'],c['I2(10)']
    capA=2*c['H4']*c['I2(5)']if 1084 in fam else 0;capB=2*c['I2(5)']if 1085 in fam else 0
    cap56=(a//2)*t if 56 in fam else 0;cap76=min(a-1,h)if 76 in fam else 0
    boundary=[]
    if capA or cap56:boundary.append((cap56,0,0,capA,0))
    for k in range(1,cap76+1):boundary.append((((a-k-1)//2)*t if 56 in fam else 0,0,k,capA,0))
    if capB:boundary.append((((a-1)//2)*t if 56 in fam else 0,0,0,capA,capB))
    if 57 in fam:boundary.append((0,t,0,0,0))
    boundary=sorted(set(boundary));env=q['envelope'];assert env['kind_order']==[56,57,76,1084,1085]
    assert env['boundary_counts']==list(map(list,boundary))
    assert sorted(tuple(map(F,w))for w in env['coefficient_vertices'])==coefficient_vertices(boundary,fam)
    witnesses=env['attaining_marker_cliques'];assert len(witnesses)==len(boundary)
    for row,p in zip(witnesses,boundary):
        assert row['counts']==list(p);chosen=row['marker_clique']
        assert len(chosen)==len(set(map(tuple,chosen)))and Counter(k for k,i in chosen)==Counter({k:z for k,z in zip(env['kind_order'],p)if z})
        assert all(k in fam and 0<=i<len(fam[k])for k,i in chosen)
        assert all(edge(ka,fam[ka][i],kb,fam[kb][j])for(ka,i),(kb,j)in combinations(chosen,2))
    return env


def check():
    data=read(TAG);assert data['status']=='PASS_RAW108_COMPLETE_EXTRA_VERTEX_MARKER_FAMILY'
    for p,h in data['sources'].items():assert digest(ROOT/p)==h,p
    ast.parse((ROOT/f'scripts/{TAG}.py').read_text())
    pilot=read(PILOT);old76=read(OLDMARKER);old56=read(NEW);joint=read(JOINT);model=read(MODEL)
    geo=next(g for g in pilot['geometries']if g['raw_family']==108);assert data['geometry']==geo
    ng=geo['native_gram'];g=mat(geo['scaled_gram']);sc=[Q(*x)for x in geo['squared_scales']]
    assert len(g)==len(ng)==8 and all(len(r)==8 for r in g+ng)
    assert all(g[i][k]==g[k][i]and ng[i][k]==ng[k][i]for i,k in product(range(8),repeat=2))
    for i in range(8):
        assert sc[i].sign()>0 and g[i][i]==2*sc[i]
        for k in range(i):assert g[i][k].sign()<=0 and g[i][k]*g[i][k]==squared(ng[i][k])*sc[i]*sc[k]
    supports={tuple(r['omitted']):r for r in geo['support_records']};assert len(supports)==28 and set(supports)==set(combinations(range(8),2))
    pd=[]
    for om,r in supports.items():
        ids=sorted(set(range(8))-set(om));assert r['support']==ids
        yes=positive(principal(g,ids));assert yes==(r['inertia']==[6,0,0])
        if yes:pd.append(om)
    assert set(pd)-set(map(tuple,geo['counted_vertex_omissions']))=={(2,4),(2,5)}and geo['required_quads']==3
    for sx in geo['compact_simplex_faces']:
        r=residual(g,sx['fixed'],sx['sides']);assert r==mat(sx['residual_gram'])
        inv=inverse(r);assert inv==mat(sx['inverse'])and all(x.sign()<0 for row in inv for x in row)
    native={1084:ng,1085:ng,56:old56['geometries']['56']['native_gram'],57:old56['geometries']['57']['native_gram'],76:old76['core_geometry']['native_gram']}
    cases=expressions=0;seenpairs=set()
    for test in data['relative_tests']:
        ka,kb=test['kinds'];assert(ka,kb)not in seenpairs;seenpairs.add((ka,kb))
        adj,scales,j,active,components=frame(test['finite_frame']);anchor=test['anchor']
        check_marker(anchor,ka,native[ka],adj,scales,j,active)
        ea=embeddings(ka,native[ka],adj);eb=embeddings(kb,native[kb],adj)
        assert tuple(map(tuple,anchor['active_embedding']))in ea
        expected={signature(ka,anchor,kb,{'active_embedding':mp})for mp in eb if ka!=kb or dict(mp)!=dict(anchor['active_embedding'])}
        seen=set()
        for rec in test['records']:
            m=rec['other_marker'];check_marker(m,kb,native[kb],adj,scales,j,active)
            assert tuple(map(tuple,m['active_embedding']))in eb
            pattern=signature(ka,anchor,kb,m);assert pattern not in seen and list(pattern)==rec['pattern'];seen.add(pattern)
            expressions+=pair_result(anchor['native'],m['native'],rec['result'],j,scales)
            assert rec['result']['allowed']==compatible(ka,anchor,kb,m);cases+=1
        assert seen==expected
    assert cases==32 and len(seenpairs)==9
    family=[r for r in pilot['complete_family_inventory']if r['raw_family']==108]
    pointed=[dict(r,kind=1084 if r['omitted']==[2,4]else 1085)for r in pilot['extra_marker_inventory']if r['raw_family']==108]
    assert data['family']==family and data['pointed_family']==pointed and len(family)==196 and len(pointed)==392
    fr={(r['target_vertex'],r['kind']):r for r in pointed};assert len(fr)==392
    rows={r['vertex']:r for r in data['forest_records']};assert len(rows)==len(data['forest_records'])==383 and set(rows)=={r['target_vertex']for r in pointed}
    prior=defaultdict(list)
    for old in(old56,old76):
        for r in old['forest_records']:prior[r['vertex']].append(r)
    total=Counter();omitted={};rowfamily={}
    for v,r in rows.items():
        adj,scales,j,active,components=frame(r['finite_frame']);nodes=r['old_nodes'];groups=r['old_groups'];owner={i:g['index']for g in groups for i in g['literal_stars']}
        assert set(owner)==set(range(len(nodes)))and all(i in active for s in nodes for c in s['contacts']for i,e in c)
        for old in prior[v]:
            for field in('forest','finite_frame','old_nodes','old_groups','old_graph_edges','prior_history'):assert r[field]==old[field]
        for z in groups:assert Counter(o['face']for o in z['origins'])==Counter(dict(z['positive_face_coefficients']))
        oldkeys={key(normalize(s,dict(zip(range(13),range(13))))['contacts']):i for i,s in enumerate(nodes)};assert len(oldkeys)==len(nodes)
        allkeys=set(oldkeys);omit=set();rowfamily[v]={}
        for ff in r['marker_families']:
            kind=ff['kind'];assert kind not in rowfamily[v];rowfamily[v][kind]=ff
            fp=fr[v,kind];assert ff['face']==fp['face']and ff['omitted']==fp['omitted']and ff['actual_extra_occurrences']==1
            markers=ff['markers'];expected=embeddings(kind,ng,adj)
            assert len(markers)==fp['active_role_orientations']==len(expected)
            assert {tuple(map(tuple,m['active_embedding']))for m in markers}==expected
            mk={key(m['native']['contacts']):i for i,m in enumerate(markers)};assert len(mk)==len(markers)and not allkeys&set(mk);allkeys.update(mk)
            for i,m in enumerate(markers):assert m['index']==i;check_marker(m,kind,ng,adj,scales,j,active)
            assert len(ff['anchor_decisions'])==len(nodes)and [x['ordinary']for x in ff['anchor_decisions']]==list(range(len(nodes)))
            good=[]
            for rec in ff['anchor_decisions']:
                s=nodes[rec['ordinary']];ports=[p for p in markers[0]['triangle_port_pairs']if p in s['occupied_pairs']]
                result=rec['result']
                if ports:assert result==dict(allowed=False,triangle_port_exclusion=ports,expressions=[])
                else:expressions+=pair_result(s,markers[0]['native'],result,j,scales)
                if result['allowed']:good.append(rec['ordinary']);omit.add(owner[rec['ordinary']])
            assert len(ff['transports'])==len(markers);edges=[]
            for mi,tr in enumerate(ff['transports']):
                perm=dict(tr['frame_permutation']);assert set(perm)==set(perm.values())==set(range(13))
                assert all(adj[i].get(k,2)==adj[perm[i]].get(perm[k],2)for i,k in combinations(range(13),2))
                on=[oldkeys[key(normalize(s,perm)['contacts'])]for s in nodes]
                mn=[mk[key(normalize(m['native'],perm)['contacts'])]for m in markers]
                assert tr['ordinary_permutation']==on and tr['marker_permutation']==mn and mn[0]==mi
                assert len(set(on))==len(nodes)and len(set(mn))==len(markers)
                for i,k in enumerate(on):assert owner[i]==owner[k]and normalize(nodes[i],perm)==normalize(nodes[k],dict(zip(range(13),range(13))))
                for i,k in enumerate(mn):assert normalize(markers[i]['native'],perm)==markers[k]['native']
                edges.extend([on[i],mi]for i in good)
            assert ff['marker_ordinary_edges']==sorted(edges)
            co=Counter(r['forest']);cap=2*co['I2(5)']*(co['H4']if kind==1084 else 1)
            chosen=ff['attaining_marker_clique'];assert ff['capacity']==len(chosen)==cap
            assert len(set(chosen))==len(chosen)and all(0<=i<len(markers)for i in chosen)
            assert all(compatible(kind,markers[i],kind,markers[k])for i,k in combinations(chosen,2))
            total['markers']+=len(markers);total['anchor_tests']+=len(nodes);total['anchor_unknowns']+=len(good);total['retained_ordinary_pairs']+=len(edges)
        assert set(rowfamily[v])=={kind for vv,kind in fr if vv==v};omitted[v]=omit
        oldmarkers=[m for old in prior[v]for m in(old['markers']if'markers'in old else[m for ff in old['marker_families']for m in ff['markers']])]
        assert not{key(m['native']['contacts'])for m in oldmarkers}&(allkeys-set(oldkeys))
        n=sum(len(ff['markers'])for ff in r['marker_families'])
        total['ordinary_stars']+=len(nodes);total['ordinary_pair_positions']+=len(nodes)*n
        total['old_marker_pair_positions']+=len(oldmarkers)*n;total['new_marker_pair_positions']+=n*(n-1)//2
    assert total['markers']==1676
    for k,v in total.items():assert data['counts'][k]==v,(k,v)
    # Frozen ordinary count-vector domains and geometrical provenance remain
    # inherited. Their exact saved charge arithmetic is checked independently.
    vectors={};weights=defaultdict(dict);cache={}
    for tag in(LOCAL,PATCH,THIRD):
        for r in read(tag)['forest_records']:
            weights[r['vertex']][tag]=r['certificate']['weights']
            if 'maximal_count_vectors'in r:vectors[r['vertex']]=r['maximal_count_vectors']
    for r in read(LOCAL)['forest_records']:
        if r['vertex']in vectors:continue
        tag=r['geometry_source']
        if tag not in cache:cache[tag]={z['vertex']:z for z in read(tag)['forest_records']}
        vectors[r['vertex']]=cache[tag][r['vertex']]['maximal_count_vectors']
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']};base={v:F(x)for v,x in joint['exported']['repair']['vertex_multipliers']};mu=dict(base)
    assert data['exported']['repair']['face_multipliers']==dual['face_multipliers']==joint['exported']['repair']['face_multipliers']
    oldjoint={r['vertex']:r for r in joint['forest_records']};oldterms=defaultdict(dict)
    for r in old76['forest_records']:oldterms[r['vertex']][76]=r['face']
    for r in old56['forest_records']:
        for ff in r['marker_families']:oldterms[r['vertex']][ff['raw_family']]=ff['face']
    quotas={q['vertex']:q for q in data['exported']['target_allowances']};assert len(quotas)==len(data['exported']['target_allowances'])==383 and set(quotas)==set(rows)
    columns=data['exported']['columns'];assert len(columns)==415 and {c['vertex']for c in columns}==set(rows)
    envs={v:check_envelope(rows[v],prior[v],q)for v,q in quotas.items()}
    for c in columns:
        v=c['vertex'];r=rows[v];groups=r['old_groups'];ww=list(map(F,c['old_orbit_weights']));w=list(map(F,c['retained_orbit_weights']));tag=c['source_charge']
        raw=weights[v][tag]if groups else[];assert ww==list(map(F,raw))+[F()]*(len(groups)-len(raw))
        vv=vectors[v]if groups else[[]];assert len(vv)==c['old_vector_count']and jhash(vv)==c['old_vectors_sha256']
        assert len(ww)==len(w)==len(groups)and all(x>=0 for x in ww)
        assert all(len(x)==len(ww)and sum(n*z for n,z in zip(x,ww))<=1 for x in vv)
        omit=omitted[v]|set(oldjoint[v]['ordinary_omitted_orbits']if v in oldjoint else[])
        assert c['omitted_orbits']==sorted(omit)and w==[F()if i in omit else x for i,x in enumerate(ww)]
        terms={kind:(face,F(x))for kind,face,x in c['marker_terms']};assert len(terms)==len(c['marker_terms'])
        mw=tuple(terms.get(kind,(None,F()))[1]for kind in(56,57,76,1084,1085))
        assert mw in[tuple(map(F,w))for w in envs[v]['coefficient_vertices']]
        for k,(f,x)in terms.items():
            if k in NEWKINDS:assert f==rowfamily[v][k]['face']
            else:assert f==oldterms[v][k]
        cf=Counter()
        for z,group in zip(w,groups):
            for f,n in group['positive_face_coefficients']:cf[f]+=n*z
        for f,x in terms.values():cf[f]+=x
        cf=[[f,str(x)]for f,x in sorted(cf.items())if x];assert cf==c['positive_face_coefficients']and c['vertex_cost']=='1'
        load=sum((alpha.get(f,F())*F(x)for f,x in cf),F());assert load==F(c['load'])and F(c['baseline_cost'])==base.get(v,F())and F(c['violation'])==load-base.get(v,F())
        mu[v]=max(mu.get(v,F()),load)
    for v,q in quotas.items():
        assert F(q['baseline_cost'])==base.get(v,F())and F(q['new_cost'])==mu[v]and F(q['increase'])==mu[v]-base.get(v,F())
    prior=[]
    for tag,old in((OLDMARKER,old76),(NEW,old56),(JOINT,joint)):
        for c in old['exported']['columns']:
            v=c['vertex'];load=sum((alpha.get(f,F())*F(x)for f,x in c['positive_face_coefficients']),F())
            assert load<=base.get(v,F())<=mu.get(v,F());prior.append(dict(source=tag,vertex=v,load=str(load),cost=str(mu.get(v,F()))))
    assert prior==data['exported']['preserved_columns']and len(prior)==2884
    actual=data['exported']['preserved_multi_owner_costs'];original=model['prepared']['inherited_multi_vertex_columns'];assert len(actual)==len(original)==14
    for saved,c in zip(actual,original):
        load=sum((alpha.get(f,F())*n for f,n in c['integer_face_coefficients']),F());old=sum((base.get(v,F())*n for v,n in c['integer_vertex_capacities']),F());cost=sum((mu.get(v,F())*n for v,n in c['integer_vertex_capacities']),F())
        assert saved==dict(variable=c['variable'],all_vertex_costs=c['integer_vertex_capacities'],load=str(load),baseline_cost=str(old),new_cost=str(cost))and load<=old<=cost
    rep=data['exported']['repair'];num=F(dual['objective_numerator']);den=sum(mu.values(),F())
    assert num==F(rep['objective_numerator'])and den==F(rep['normalizer'])and num/den==F(rep['bound'])and num/12-den==F(rep['remaining_normalizer_gap_to12'])
    assert rep['vertex_multipliers']==[[v,str(x)]for v,x in sorted(mu.items())if x]
    assert data['counts']['strictly_violated_columns']==sum(F(c['violation'])>0 for c in columns)==8
    assert data['counts']['raised_allowances']==sum(F(q['increase'])>0 for q in quotas.values())==8
    assert data['optimizer_invocations']==0 and data['old_outcomes_changed']is False
    print('PASS_RAW108_SAVED_CHECK',dict(relative_cases=cases,expressions=expressions,**total,columns=len(columns),bound=str(num/den)),flush=True)

if __name__=='__main__':check()
