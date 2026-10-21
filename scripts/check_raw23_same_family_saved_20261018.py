#!/usr/bin/env python3
"""Independent saved raw23 native relative and twelve-wall support checks.

No producer imports. Whole-component maps and relative signatures, exact
native gates, scaled union Gram matrices and every nine-support decision
are reconstructed. Historical raw23 catalogue/source-family and actual
finite-vertex incidence are inherited; no capacity or matrix claim.
"""
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,combinations_with_replacement,permutations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
REL='fresh_raw23_same_family_relative_pilot_20261018'
UNION='fresh_raw23_same_family_union_support_pilot_20261018'
def read(s):return json.loads((ROOT/f'results/{s}.json').read_text())
def key(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def bind(d):
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
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


def components(adj):
    unseen=set(adj);out=[]
    while unseen:
        seen=set();todo=[min(unseen)]
        while todo:
            i=todo.pop()
            if i in seen:continue
            seen.add(i);todo.extend(adj[i])
        unseen-=seen;out.append(sorted(seen))
    return out


def maps(active,target):
    ac=components(active);tc=components(target);ans=[]
    def rec(k,used,mp):
        if k==len(ac):ans.append(tuple(sorted(mp.items())));return
        a=ac[k]
        for b in tc:
            if used.intersection(b)or len(a)!=len(b):continue
            for q in permutations(b):
                m=dict(zip(a,q))
                if all(active[i].get(j,2)==target[m[i]].get(m[j],2)for i,j in combinations(a,2)):
                    rec(k+1,used|set(b),mp|m)
    rec(0,set(),{});return set(ans)


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
        assert not result['allowed']and result['triangle_port_exclusion']and all(p in a['occupied_pairs']for p in result['triangle_port_exclusion']);return 0
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


def active_parts(adj,m):
    mp=dict(m['active_embedding']);back={v:k for k,v in mp.items()}
    assert len(back)==len(mp)
    out=[]
    for block in components(adj):
        if set(block)&set(back):
            assert set(block)<=set(back)
            out.append(tuple(sorted(block,key=back.get)))
    return sorted(out,key=lambda c:min(back[i]for i in c))


def relative(adj,ka,a,kb,b):
    left=active_parts(adj,a);right=active_parts(adj,b);mp=dict(a['active_embedding'])
    labels={mp[i]:j for j,i in enumerate(sorted(mp))};pattern=[]
    for block in right:
        match=[j for j,part in enumerate(left)if set(part)==set(block)]
        if match:
            assert len(match)==1;j=match[0];pattern.append([j,[left[j].index(i)for i in block]])
        else:
            assert not set(block)&set(labels);pattern.append([-1,[]])
            labels.update({i:len(labels)+j for j,i in enumerate(block)})
    ids=sorted(labels,key=labels.get);cov=[];hs=[]
    for marker in(a,b):
        cs=[sorted([[labels[i],e]for i,e in c])for c in marker['native']['contacts']]
        order=sorted(range(2),key=lambda i:cs[i]);cov.append([cs[i]for i in order]);hs.append([marker['native']['squared_timelike_heights'][i]for i in order])
    payload=dict(contacted_native_gram=[[0 if i==j else adj[i].get(j,2)for j in ids]for i in ids],
        full_contacted_covectors=cov,heights=hs,omitted_products=[a['native']['native_omitted_product'],b['native']['native_omitted_product']])
    return key([ka,kb,pattern]),pattern,payload


def inertia(g):
    a=[r[:]for r in g];out=[0,0,0]
    while a:
        n=len(a);p=next((i for i in range(n)if a[i][i]!=0),None)
        if p is not None:
            ids=[p]+[i for i in range(n)if i!=p];a=principal(a,ids);z=a[0][0]
            out[0 if z.sign()>0 else 1]+=1
            a=[[a[i][j]-a[i][0]*a[0][j]/z for j in range(1,n)]for i in range(1,n)]
            continue
        pair=next(((i,j)for i,j in combinations(range(n),2)if a[i][j]!=0),None)
        if pair is None:out[2]+=n;break
        p,q=pair;ids=[p,q]+[i for i in range(n)if i not in pair];a=principal(a,ids);z=a[0][1]
        out[0]+=1;out[1]+=1
        a=[[a[i][j]-(a[i][0]*a[1][j]+a[i][1]*a[0][j])/z for j in range(2,n)]for i in range(2,n)]
    return out


def native(g,s,mp,j,scales):
    ids=s['support'];om=s['omitted'];cv=[]
    for i in om:cv.append(sorted([[mp[k],g[i][k]]for k in ids if g[i][k][0]!=2]))
    cv.sort();assert cv[0]!=cv[1]
    code,hs,hp,f,pr,z=physical(cv[0],cv[1],j,scales)
    assert all(a==-b for a,b in zip(z,f.sqrt(squared(g[om[0]][om[1]]))))
    return dict(contacts=cv,squared_timelike_heights=[[str(x.a),str(x.b)]for x in hs],
        native_omitted_product=g[om[0]][om[1]],occupied_pairs=[])


def main():
    d,u=read(REL),read(UNION);bind(d);bind(u)
    assert d['status']=='PASS_RAW23_SAME_FAMILY_RELATIVE_NATIVE_PILOT'
    assert u['status']=='PASS_RAW23_THREE_RETAINED_UNION_SUPPORT_PILOT'
    for x in(d,u):assert x['optimizer_invocations']==0 and x['exported_columns']==[]
    assert not d['old_graph_outcomes_changed']
    ng=d['geometry']['native_gram'];sg=mat(d['geometry']['scaled_gram']);sc=[Q(*x)for x in d['geometry']['squared_scales']]
    assert len(ng)==len(sg)==len(sc)==7 and all(len(r)==7 for r in ng+sg)
    assert all(sc[i].sign()>0 and sg[i][i]==2*sc[i]for i in range(7))
    for i,k in combinations(range(7),2):
        assert ng[i][k]==ng[k][i]and sg[i][k]==sg[k][i]
        assert sg[i][k]*sg[i][k]==squared(ng[i][k])*sc[i]*sc[k]and sg[i][k].sign()<=0
    assert inertia(sg)==[5,1,1]
    ss=d['geometry']['finite_vertex_types'];assert len(ss)==10 and [r['kind']for r in ss]==list(range(10))
    assert [list(om)for om in combinations(range(7),2)if positive(principal(sg,sorted(set(range(7))-set(om))))]==[s['omitted']for s in ss]
    for s in ss:assert s['support']==sorted(set(range(7))-set(s['omitted']))
    quad=read('fresh_sqrt5_rank13_packing_quad_separation_20261017')['canonical_cores'][19]
    assert quad['source_core']==23 and quad['required']==10 and quad['quad_partitions']==d['geometry']['quad_partitions']
    occ=Counter(tuple(sorted(x['vertex_omitted']))for x in quad['aggregate_corner_occurrences'])
    assert occ==Counter({tuple(s['omitted']):4 for s in ss})
    assert d['geometry']['quadrilateral_corner_occurrences']==40 and d['geometry']['actual_vertex_occurrences']==10
    # Family catalogue/provenance is inherited but exact saved coverage and masses are compared.
    inc=read('fresh_sqrt5_prism_simplex_base_incidence_20261017')
    ap=read('fresh_sqrt5_actual_provenance_dual_20261017')['provenance_masses']
    p=read('fresh_all_marker_families_dual_continuation_20261018')['prepared']['provenance_to_face']
    explicit={ap[i]['source_context']:f for i,f in p}
    expected=[(ctx,explicit.get(ctx,f),bs)for ctx,f,c,k,bs in inc['contexts']if c==19]
    fam=d['complete_source_family'];assert len(fam)==len({r['face']for r in fam})==235
    assert [(r['source_context'],r['face'],r['bases'])for r in fam]==expected
    a={f:F(x)for f,x in read('fresh_raw23_surviving_coface_witness_20261018')['candidate']['face_multipliers']}
    assert all(r['required_quads']==10 and F(r['current_mass'])==a.get(r['face'],F())for r in fam)
    assert sum(F(r['current_mass'])>0 for r in fam)==50
    cases={key([*r['kinds'],r['pattern']]):r for r in d['relative_cases']}
    assert len(cases)==len(d['relative_cases'])==175
    hosts={tuple(r['kinds']):r for r in d['hosts']};assert len(hosts)==len(d['hosts'])==55
    assert set(hosts)==set(combinations_with_replacement(range(10),2))
    used=Counter();positions=identical=0
    for kinds,host in hosts.items():
        saved=next(c for c in cases.values()if c['kinds']==list(kinds))
        assert saved['forest']==host['forest']==sorted(ss[kinds[0]]['forest']+ss[kinds[1]]['forest']+['A1']*3)
        adj,scales,j,active,parts=frame(saved['finite_frame']);allmaps=[]
        for kind in kinds:
            s=ss[kind];ids=s['support'];ag={i:{k:ng[i][k][0]for k in ids if i!=k and ng[i][k][0]!=2}for i in ids}
            allmaps.append(maps(ag,adj))
        assert [len(x)for x in allmaps]==[host['left_orientations'],host['right_orientations']]
        anchor=saved['anchor'];amp=tuple(map(tuple,anchor['active_embedding']));assert amp in allmaps[0]
        assert anchor['native']==native(ng,ss[kinds[0]],dict(amp),j,scales)
        local=0
        for mp in allmaps[1]:
            b=dict(active_embedding=[list(x)for x in mp],native=native(ng,ss[kinds[1]],dict(mp),j,scales))
            if anchor['native']['contacts']==b['native']['contacts']:identical+=1;continue
            k,pat,payload=relative(adj,kinds[0],anchor,kinds[1],b)
            assert k in cases and payload==cases[k]['payload'];used[k]+=1;positions+=1;local+=1
        assert local==host['anchor_pairs']
    assert positions==1482 and identical==10 and set(used)==set(cases)
    wallchecks=0
    for k,c in cases.items():
        adj,scales,j,active,parts=frame(c['finite_frame'])
        assert relative(adj,c['kinds'][0],c['anchor'],c['kinds'][1],c['other'])==(k,c['pattern'],c['payload'])
        assert all(i in active for m in(c['anchor'],c['other'])for cv in m['native']['contacts']for i,e in cv)
        wallchecks+=pair_result(c['anchor']['native'],c['other']['native'],c['result'],j,scales)
    retained=[i for i,c in enumerate(d['relative_cases'])if c['result']['allowed']]
    assert retained==[12,22,103]and len(u['cases'])==3
    assert [c['relative_index']for c in u['cases']]==retained
    supports=0
    for z in u['cases']:
        c=d['relative_cases'][z['relative_index']];assert z['kinds']==c['kinds']and z['full_frame']==c['finite_frame']
        g=mat(z['scaled_gram']);sc=[Q(*x)for x in z['squared_scales']]
        assert len(g)==len(sc)==12 and all(len(r)==12 for r in g)
        assert all(sc[i].sign()>0 and g[i][i]==2*sc[i]for i in range(12))
        assert all(g[i][j]==g[j][i]and g[i][j].sign()<=0 for i,j in combinations(range(12),2))
        assert inertia(g)==z['inertia']==[9,1,2]
        assert inertia([[Q(x.a,-x.b)for x in r]for r in g])==z['conjugate_inertia']==[10,0,2]
        vv=[list(v)for v in combinations(range(12),9)if positive(principal(g,v))];supports+=220
        assert vv==z['all_positive_definite_nine_supports']and len(vv)==18
        known=[];assert len(z['prism_maps'])==2
        for side,pm in enumerate(z['prism_maps']):
            mp=dict(pm['raw_map']);assert set(mp)==set(range(7))and len(set(mp.values()))==7
            for i,k in combinations(range(7),2):assert g[mp[i]][mp[k]]*g[mp[i]][mp[k]]==squared(ng[i][k])*sc[mp[i]]*sc[mp[k]]
            passive=sorted(set(range(9))-set(mp.values()));assert passive==pm['passive']and len(passive)==4
            assert all(g[i][mp[k]]==0 for i in passive for k in range(7))
            v=[sorted(passive+[mp[i]for i in s['support']])for s in ss]
            assert v==pm['vertices']and all(x in vv for x in v);known+=v
        shared=sorted(set(map(tuple,z['prism_maps'][0]['vertices']))&set(map(tuple,z['prism_maps'][1]['vertices'])))
        assert len(shared)==2 and len(set(shared[0])&set(shared[1]))==8
        endpoint_kinds=[[pm['vertices'].index(list(v))for pm in z['prism_maps']]for v in shared]
        assert sorted(endpoint_kinds)==sorted([z['kinds'],[7-k for k in z['kinds']]])
        known=set(map(tuple,known));assert len(known)==z['inherited_union_vertex_count']==18
        assert [v for v in vv if tuple(v)not in known]==z['extra_positive_definite_supports']==[]
        # Recompute all exterior products from the actual thirteen-root native frame.
        adj,scale,j,active,parts=frame(z['full_frame'])
        covs=c['anchor']['native']['contacts']+c['other']['native']['contacts']
        common=[cv for cv in covs[:2]if cv in covs[2:]];assert len(common)==1
        cvs=common+[next(cv for cv in covs[2*i:2*i+2]if cv!=common[0])for i in range(2)]
        assert len(z['cross_product_checks'])==3
        for (i,k),rec in zip(combinations(range(3),2),z['cross_product_checks']):
            code,hs,hp,f,pr,phys=physical(cvs[i],cvs[k],j,scale)
            sq=f.mul(phys,phys);assert all(x==0 for x in sq[1:])and f.sign(phys)<=0
            expected=g[9+i][9+k]*g[9+i][9+k]/(sc[9+i]*sc[9+k])
            assert sq[0]==expected==Q(*rec['physical_squared_product'])
            assert rec['nodes']==[9+i,9+k]and Q(*rec['scaled_product'])==g[9+i][9+k]
        assert z['passive_rank']==4 and z['actual_vertex_argument_inherited']
    assert supports==u['complete_support_decisions']==660
    print('PASS_INDEPENDENT_RAW23_RELATIVE_AND_UNION_SAVED',dict(source_classes=235,host_pairs=55,anchor_positions=positions,
        relative_cases=175,wall_expressions=wallchecks,retained_cases=3,support_decisions=supports,extra_supports=0),flush=True)


if __name__=='__main__':main()

