#!/usr/bin/env python3
"""Independent saved raw357 incidence pilot checker (stdlib only).

Ordered field and contacted-block routines reuse the standalone raw56/57
checker implementation, not producer code. Historical family/catalogue and
ordinary-domain completeness remain inherited; no pair gates or optimization.
"""
import ast,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw357_extra_vertex_incidence_pilot_20261018'
def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
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


def native_autos(g):
    n=len(g);sig=[sorted(tuple(g[i][j])for j in range(n)if j!=i)for i in range(n)]
    choices=[[j for j in range(n)if sig[j]==sig[i]]for i in range(n)];out=[]
    order=sorted(range(n),key=lambda i:(len(choices[i]),i))
    def rec(k,mp):
        if k==n:out.append([mp[i]for i in range(n)]);return
        i=order[k]
        for j in choices[i]:
            if j not in mp.values()and all(g[i][a]==g[j][b]for a,b in mp.items()):rec(k+1,mp|{i:j})
    rec(0,{});return sorted(out)


def main():
    data=read(TAG);assert data['status']=='PASS_BOUNDED_RAW357_EXTRA_VERTEX_INCIDENCE_PILOT'
    for p,h in data['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    raw=read('fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017')
    quad=read('fresh_sqrt5_rank13_packing_quad_separation_20261017')
    forests=read('fresh_sqrt5_rank13_combined_disjoint_ports_20261017')
    spectators=read('fresh_sqrt5_rank13_dotted_b1_spectators_20261017')
    positive_faces={r['face']:r for r in read('fresh_all_marker_families_model_support_20261018')['positive_faces']}
    r=data['geometry'];native=r['native_gram'];g=mat(r['scaled_gram']);sc=[Q(*x)for x in r['squared_scales']]
    assert len(g)==len(native)==len(sc)==10 and all(len(x)==10 for x in g+native)
    assert all(g[i][j]==g[j][i]and native[i][j]==native[j][i]for i,j in product(range(10),repeat=2))
    assert all(sc[i].sign()>0 and g[i][i]==2*sc[i]for i in range(10))
    source=raw['records'][357]
    for(i,j),m in zip(combinations(range(10),2),source['row']['labels']):
        assert native[i][j][0]==m
        if m==0:assert native[i][j][1:]==source['parameters'][f'{i},{j}']
        assert g[i][j].sign()<=0 and g[i][j]*g[i][j]==squared(native[i][j])*sc[i]*sc[j]
    pending=set(combinations(range(10),2));pd=set();assert len(r['support_records'])==45
    for z in r['support_records']:
        om=tuple(z['omitted']);assert om in pending;pending.remove(om)
        ids=sorted(set(range(10))-set(om));assert z['support']==ids
        good=positive(principal(g,ids));assert good==(z['inertia']==[8,0,0])
        if good:pd.add(om)
    assert not pending and len(pd)==7
    c=quad['canonical_cores'][318];assert c['source_core']==357 and c['required']==1
    assert sorted(map(list,pd))==c['pd_deletions']
    used={tuple(sorted(x['vertex_omitted']))for x in c['aggregate_corner_occurrences']}
    assert len(used)==4 and sorted(map(list,used))==r['counted_vertex_omissions']
    assert pd-used=={(2,4),(2,5),(2,6)}and [x['omitted']for x in r['extra_vertices']]==[[2,4],[2,5],[2,6]]
    assert native_autos(native)==r['native_automorphisms']==[list(range(10))]
    # A nonsingular9-root minor with one negative pivot and a full null
    # vector independently certify inertia(8,1,1).
    basis=source['quotient_basis'];a=principal(g,basis);assert len(a)==9;neg=0
    for k in range(9):
        assert a[k][k]!=0;neg+=a[k][k].sign()<0
        for i in range(k+1,9):
            for j in range(k+1,9):a[i][j]=a[i][j]-a[i][k]*a[k][j]/a[k][k]
    assert neg==1
    outside=next(i for i in range(10)if i not in basis);j=inverse(principal(g,basis))
    coeff=[-sum((j[a][b]*g[basis[b]][outside]for b in range(9)),Q())for a in range(9)]
    assert all(g[i][outside]+sum((g[i][b]*z for b,z in zip(basis,coeff)),Q())==0 for i in range(10))
    assert len(r['compact_simplex_faces'])==1;sx=r['compact_simplex_faces'][0]
    assert sx['fixed']==[0,7,8,9]and sx['sides']==[1,3,4,5,6]and sx['excluded_wall']==2
    sg=residual(g,sx['fixed'],sx['sides']);iv=inverse(sg)
    assert sg==mat(sx['residual_gram'])and iv==mat(sx['inverse'])and all(x.sign()<0 for row in iv for x in row)
    supports=[sorted(sx['fixed']+list(set(sx['sides'])-{i}))for i in sx['sides']]
    assert sx['vertex_supports']==supports and all(positive(principal(g,ids))for ids in supports)
    for ex in r['extra_vertices']:
        om=ex['omitted'];ids=ex['support'];j=inverse(principal(g,ids));cov=[[g[i][k]for k in ids]for i in om]
        assert cov==mat(ex['scaled_covectors'])
        pp=[[sum((cov[a][i]*j[i][k]*cov[b][k]for i,k in product(range(8),repeat=2)),Q())for b in range(2)]for a in range(2)]
        hh=[pp[i][i]-g[om[i]][om[i]]for i in range(2)];z=pp[0][1]-g[om[0]][om[1]]
        assert pp==mat(ex['scaled_projection'])and hh==[Q(*x)for x in ex['scaled_heights']]
        assert z==Q(*ex['positive_height_product_root'])and z.sign()>0 and all(x.sign()>0 for x in hh)and z*z==hh[0]*hh[1]
        ports=sorted(set(sx['sides'])&set(ids));assert len(ports)==4
        assert ex['forced_triangle_port_pairs']==list(map(list,combinations(ports,2)))
    family={z['face']:z for z in data['complete_family_inventory']};assert len(family)==len(data['complete_family_inventory'])==47
    expected={f['index']for f in forests['faces']if any(ci==318 for _,ci in f.get('source_contexts',[]))};assert set(family)==expected
    markers={(z['face'],tuple(z['omitted'])):z for z in data['extra_marker_inventory']}
    assert len(markers)==len(data['extra_marker_inventory'])==141
    assert set(markers)=={(f,om)for f in family for om in((2,4),(2,5),(2,6))}
    lookup={tuple(z['forest']):z['index']for z in forests['vertices']};total=Counter();mass=F();npos=0
    for f,z in family.items():
        ctx=z['source_context'];assert forests['faces'][f]['source_contexts']==[[ctx,318]]
        assert z['spectator_forest']==spectators['contexts'][ctx]['spectator_forest']
        mm=F(positive_faces.get(f,{}).get('mass','0'));assert mm==F(z['current_mass']);mass+=mm;npos+=mm>0
        assert z['required_quads']==1 and z['extra_vertex_occurrences']==3
        for ex in r['extra_vertices']:
            t=markers[f,tuple(ex['omitted'])];assert t['active_forest']==ex['active_forest']and t['spectator_forest']==z['spectator_forest']
            assert t['target_forest']==sorted(t['active_forest']+t['spectator_forest'])and lookup[tuple(t['target_forest'])]==t['target_vertex']
            cc=Counter(t['target_forest']);b=cc['A2'];om=ex['omitted']
            n=4*b*(b-1)*cc['H4']if om==[2,4]else 4*cc['A1']*b*(b-1)*cc['H3']if om==[2,5]else 16*b*(b-1)*(b-2)*cc['I2(5)']
            assert n==t['active_role_orientations']and t['actual_extra_occurrences']==1 and F(t['current_mass'])==mm
            total[tuple(om)]+=n
    assert npos==9 and mass==F(data['current_positive_mass'])==F(83284683,10**10)
    assert dict(total)=={(2,4):552,(2,5):1144,(2,6):8928}
    assert len(data['native_recovery_pilots'])==3;pilotmaps=0
    for z in data['native_recovery_pilots']:
        assert z['face']==70408;ex=next(x for x in r['extra_vertices']if x['omitted']==z['omitted']);ids=ex['support']
        adj,sc,j,active,_=frame(z['finite_frame']);ag={i:{k:native[i][k][0]for k in ids if i!=k and native[i][k][0]!=2}for i in ids}
        expected=maps(ag,adj);assert len(expected)==len(z['complete_active_maps'])==z['active_role_orientations'];keys=set()
        for saved in z['complete_active_maps']:
            mp=dict(saved['active_embedding']);key0=tuple(sorted(mp.items()));assert key0 in expected;expected.remove(key0)
            assert all(not adj[i].get(k)for i in mp.values()for k in adj if k not in mp.values())
            cov=[sorted((mp[k],native[i][k])for k in ids if native[i][k][0]!=2)for i in z['omitted']]
            assert [[list(x)for x in cv]for cv in cov]==saved['full_frame_covectors']
            key=tuple(sorted(tuple((i,tuple(e))for i,e in cv)for cv in cov));assert key not in keys;keys.add(key)
            assert all(i in active for cv in cov for i,e in cv)
            sq=[squared(e)*sc[i]for cv in cov for i,e in cv];field=Field(sq);hh=[]
            for cv in cov:
                pp=project(cv,cv,j,sc,field);assert all(x==0 for x in pp[1:]);hh.append(pp[0]-2)
            assert hh==[Q(*x)for x in saved['positive_heights']]and all(x.sign()>0 for x in hh)
            contact=squared(native[z['omitted'][0]][z['omitted'][1]])
            field=Field(sq+[hh[0]*hh[1]]+([contact]if contact!=0 else[]))
            pp=project(cov[0],cov[1],j,sc,field)
            assert [a-b for a,b in zip(pp,field.sqrt(hh[0]*hh[1]))]==[-x for x in field.sqrt(contact)];pilotmaps+=1
        assert not expected
    old=Counter();new=Counter();ordinary={};cache={}
    for tag in('fresh_raw76_extra_vertex_complete_family_20261018','fresh_raw56_raw57_extra_vertex_family_20261018','fresh_raw108_complete_marker_family_20261018','fresh_raw112_complete_marker_family_20261018'):
        for z in read(tag)['forest_records']:old[z['vertex']]+=z['counts']['markers']
    for tag in('fresh_current324_new5073_local_covers_20261018','fresh_new5073_raw76_high_mass_local_covers_20261018','fresh_new5073_raw76_third_class_local_covers_20261018'):
        for z in read(tag)['forest_records']:
            if 'domain'in z:q=z['domain']
            else:
                name=z['group_source']
                if name not in cache:cache[name]={q['vertex']:q for q in read(name)['forest_records']}
                q=cache[name][z['vertex']];q=q.get('domain',q.get('union',q))
            keys=[i for group in q['corner_orbits']for i in group['literal_stars']]
            assert len(keys)==len(set(keys))and set(keys)==set(range(len(keys)));ordinary[z['vertex']]=len(keys)
    assert len(ordinary)==600
    for z in markers.values():new[z['target_vertex']]+=z['active_role_orientations']
    work=data['future_pair_workload'];assert len(work['records'])==len(new)and{z['vertex']for z in work['records']}==set(new)
    for z in work['records']:
        v=z['vertex'];n=new[v];assert(z['new_markers'],z['old_ordinary_stars'],z['old_family_markers'])==(n,ordinary.get(v,0),old[v])
        assert z['new_ordinary_pairs']==n*ordinary.get(v,0)and z['new_old_marker_pairs']==n*old[v]and z['new_marker_pairs']==n*(n-1)//2
    for k in('new_ordinary_pairs','new_old_marker_pairs','new_marker_pairs'):assert work[k]==sum(z[k]for z in work['records'])
    assert data['counts']['role_orientations']==sum(total.values())==10624 and pilotmaps==data['counts']['native_pilot_maps']==152
    assert data['optimizer_invocations']==data['new_pair_tests']==0 and data['exported_columns']==[]and not data['old_outcomes_changed']
    print('PASS_INDEPENDENT_RAW357_SAVED_CHECK',45,'supports',1,'simplex',47,'classes',141,'domains',10624,'orientations',pilotmaps,'native maps',flush=True)


if __name__=='__main__':main()


