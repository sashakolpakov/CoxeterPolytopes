#!/usr/bin/env python3
"""Independent saved four-core incidence pilot checker (stdlib only).

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
TAG='fresh_raw14_raw112_raw77_raw108_core_incidence_pilot_20261018'
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


def main():
    data=read(TAG);assert data['status']=='PASS_BOUNDED_FOUR_RAW_CORE_INCIDENCE_PILOT'
    for p,h in data['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    ast.parse((ROOT/f'scripts/{TAG}.py').read_text())
    raw=read('fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017');quad=read('fresh_sqrt5_rank13_packing_quad_separation_20261017')
    forests=read('fresh_sqrt5_rank13_combined_disjoint_ports_20261017');spect=read('fresh_sqrt5_rank13_dotted_b1_spectators_20261017')
    support=read('fresh_raw76_complete_family_support_20261018');positive_faces={r['face']:r for r in support['positive_faces']}
    geom={r['raw_family']:r for r in data['geometries']};assert len(geom)==len(data['geometries'])==4 and set(geom)=={14,112,77,108}
    tests=simplices=pilotmaps=0
    for ri,r in geom.items():
        g=mat(r['scaled_gram']);native=r['native_gram'];n=r['nodes'];sc=list(map(lambda x:Q(*x),r['squared_scales']))
        assert len(g)==len(native)==len(sc)==n and all(len(x)==n for x in g+native)
        assert all(g[i][j]==g[j][i]and native[i][j]==native[j][i]for i,j in product(range(n),repeat=2))
        assert all(sc[i].sign()>0 and g[i][i]==2*sc[i]for i in range(n))
        source=raw['records'][ri];assert source['nodes']==n
        for (i,j),m in zip(combinations(range(n),2),source['row']['labels']):
            assert native[i][j][0]==m
            if m==0:assert native[i][j][1:]==source['parameters'][f'{i},{j}']
            assert g[i][j].sign()<=0 and g[i][j]*g[i][j]==squared(native[i][j])*sc[i]*sc[j]
        omset=set(combinations(range(n),2));assert len(r['support_records'])==len(omset)
        pd=set()
        for row in r['support_records']:
            om=tuple(row['omitted']);assert om in omset;omset.remove(om)
            ids=sorted(set(range(n))-set(om));assert row['support']==ids
            good=positive(principal(g,ids));assert good==(row['inertia']==[n-2,0,0]);tests+=1
            if good:pd.add(om)
        assert not omset
        c=quad['canonical_cores'][r['canonical_core']];assert c['source_core']==ri
        assert sorted(map(list,pd))==c['pd_deletions']
        used={tuple(sorted(x['vertex_omitted']))for x in c['aggregate_corner_occurrences']}
        assert sorted(map(list,used))==r['counted_vertex_omissions']and len(c['aggregate_corner_occurrences'])==4*r['required_quads']
        assert sorted(pd-used)==list(map(lambda x:tuple(x['omitted']),r['extra_vertices']))
        autom=[]
        for p in permutations(range(n)):
            if all(native[i][j]==native[p[i]][p[j]]for i,j in combinations(range(n),2)):autom.append(list(p))
        assert autom==r['native_automorphisms']==[list(range(n))]
        for sx in r['compact_simplex_faces']:
            fixed=sx['fixed'];sides=sx['sides'];e=sx['excluded_wall']
            assert sides==sorted(next(i for i in om if i!=e)for om in pd if e in om)
            assert fixed==sorted(set(range(n))-set(sides)-{e})and positive(principal(g,fixed))
            sg=residual(g,fixed,sides);inv=inverse(sg)
            assert sg==mat(sx['residual_gram'])and inv==mat(sx['inverse'])
            assert all(x.sign()<0 for row in inv for x in row)
            assert sx['dimension']==len(sides)-1
            assert all(positive(principal(g,ids))for ids in sx['vertex_supports']);simplices+=1
        for ex in r['extra_vertices']:
            om=ex['omitted'];ids=ex['support'];j=inverse(principal(g,ids));cov=[[g[i][k]for k in ids]for i in om]
            assert cov==mat(ex['scaled_covectors'])
            p=[[sum((cov[a][i]*j[i][k]*cov[b][k]for i,k in product(range(n-2),repeat=2)),Q())for b in range(2)]for a in range(2)]
            hs=[p[i][i]-g[om[i]][om[i]]for i in range(2)];z=p[0][1]-g[om[0]][om[1]]
            assert p==mat(ex['scaled_projection'])and hs==[Q(*x)for x in ex['scaled_heights']]
            assert z==Q(*ex['positive_height_product_root'])and z.sign()>0 and all(x.sign()>0 for x in hs)and z*z==hs[0]*hs[1]
            expected=set()
            for sx in r['compact_simplex_faces']:
                if set(ids)==set(sx['fixed'])|(set(sx['sides'])-{om[1]}):expected.update(combinations(sorted(set(sx['sides'])&set(ids)),2))
            assert sorted(expected)==list(map(tuple,ex['forced_triangle_port_pairs']))
    # Complete family selection reconstructed from the inherited catalogue.
    selected={r['face']:r for r in data['complete_family_inventory']};assert len(selected)==len(data['complete_family_inventory'])==1331
    expected={fi for fi,f in enumerate(forests['faces'])if f['kind']=='canonical_dotted_quad_aggregate'and any(ci in(12,92,66,88)for _,ci in f['source_contexts'])}
    assert set(selected)==expected
    markerkeys=set();summary={ri:dict(classes=0,positive=0,mass=F(),markers=0,orient=0,extra_mass=F())for ri in geom}
    lookup={tuple(r['forest']):r['index']for r in forests['vertices']}
    for fi,row in selected.items():
        ri=row['raw_family'];ctx=row['source_context'];assert forests['faces'][fi]['source_contexts']==[[ctx,row['canonical_core']]]
        assert row['spectator_forest']==spect['contexts'][ctx]['spectator_forest']
        mass=F(positive_faces.get(fi,{}).get('mass','0'));assert mass==F(row['current_mass'])
        assert row['required_quads']==geom[ri]['required_quads']and row['extra_vertex_occurrences']==len(geom[ri]['extra_vertices'])
        z=summary[ri];z['classes']+=1;z['positive']+=mass>0;z['mass']+=mass
        markerkeys.update((fi,tuple(ex['omitted']))for ex in geom[ri]['extra_vertices'])
    assert len(data['extra_marker_inventory'])==len(markerkeys)==784
    for row in data['extra_marker_inventory']:
        k=(row['face'],tuple(row['omitted']));assert k in markerkeys;markerkeys.remove(k)
        base=selected[row['face']];ri=row['raw_family'];ex=next(x for x in geom[ri]['extra_vertices']if x['omitted']==row['omitted'])
        assert row['active_forest']==ex['active_forest']and row['spectator_forest']==base['spectator_forest']
        assert row['target_forest']==sorted(row['spectator_forest']+ex['active_forest'])and lookup[tuple(row['target_forest'])]==row['target_vertex']
        c=Counter(row['target_forest']);a=c['A1'];h=c['H3'];v=c['H4'];p=c['I2(5)']
        count=a*(a-1)*v if ri==112 and row['omitted']==[2,4]else a*(a-1)*(a-2)*h if ri==112 else 2*v*p if row['omitted']==[2,4]else 2*a*h*p
        assert count==row['active_role_orientations']and row['actual_extra_occurrences']==1 and row['current_mass']==base['current_mass']
        z=summary[ri];z['markers']+=1;z['orient']+=count;z['extra_mass']+=F(row['current_mass'])
    assert not markerkeys
    for row in data['native_recovery_pilots']:
        ri=row['raw_family'];r=geom[ri];g=r['native_gram'];ex=next(x for x in r['extra_vertices']if x['omitted']==row['omitted']);ids=ex['support']
        adj,sc,j,active,cc=frame(row['finite_frame']);ag={i:{k:g[i][k][0]for k in ids if i!=k and g[i][k][0]!=2}for i in ids}
        expected=maps(ag,adj);keys=set();assert len(expected)==len(row['complete_active_maps'])==row['active_role_orientations']
        for saved in row['complete_active_maps']:
            mp=dict(saved['active_embedding']);mk=tuple(sorted(mp.items()));assert mk in expected;expected.remove(mk)
            cov=[sorted((mp[k],g[i][k])for k in ids if g[i][k][0]!=2)for i in row['omitted']]
            assert [[list(x)for x in c]for c in cov]==saved['full_frame_covectors']
            key=tuple(sorted(tuple((i,tuple(e))for i,e in c)for c in cov));assert key not in keys;keys.add(key)
            assert all(i in active for c in cov for i,e in c)
            sq=[squared(e)*sc[i]for c in cov for i,e in c];f=Field(sq);hh=[]
            for c in cov:
                pp=project(c,c,j,sc,f);assert all(x==0 for x in pp[1:]);hh.append(pp[0]-2)
            assert hh==[Q(*x)for x in saved['positive_heights']]and all(x.sign()>0 for x in hh)
            native=squared(g[row['omitted'][0]][row['omitted'][1]]);f=Field(sq+[hh[0]*hh[1]]+([native]if native!=0 else[]))
            pp=project(cov[0],cov[1],j,sc,f);value=[a-b for a,b in zip(pp,f.sqrt(hh[0]*hh[1]))]
            assert value==[-x for x in f.sqrt(native)];pilotmaps+=1
        assert not expected
    for r in data['summaries']:
        z=summary[r['raw_family']]
        assert (r['complete_classes'],r['positive_classes'],F(r['positive_mass']),r['marker_domains'],r['role_orientation_count'],F(r['positive_extra_occurrence_mass']))==(z['classes'],z['positive'],z['mass'],z['markers'],z['orient'],z['extra_mass'])
    work=data['future_pair_workload'];new=Counter();old=Counter()
    for r in data['extra_marker_inventory']:new[r['target_vertex']]+=r['active_role_orientations']
    for tag in('fresh_raw76_extra_vertex_complete_family_20261018','fresh_raw56_raw57_extra_vertex_family_20261018'):
        for r in read(tag)['forest_records']:old[r['vertex']]+=r['counts']['markers']
    assert len(work['records'])==len(new)and{r['vertex']for r in work['records']}==set(new)
    for r in work['records']:
        n=new[r['vertex']];assert r['new_markers']==n and r['old_family_markers']==old[r['vertex']]
        assert r['new_ordinary_pairs']==n*r['old_ordinary_stars']and r['new_old_marker_pairs']==n*old[r['vertex']]and r['new_marker_pairs']==n*(n-1)//2
    for k in('new_ordinary_pairs','new_old_marker_pairs','new_marker_pairs'):assert work[k]==sum(r[k]for r in work['records'])
    assert tests==92 and simplices==6 and pilotmaps==40
    assert data['new_pair_tests']==data['optimizer_invocations']==0 and data['exported_columns']==[]and not data['old_outcomes_changed']
    print('PASS_INDEPENDENT_SAVED_FOUR_CORE_PILOT',tests,'supports',simplices,'simplexes',pilotmaps,'native_maps',len(selected),'classes',flush=True)

if __name__=='__main__':main()
