#!/usr/bin/env python3
"""Stdlib-only saved raw76 extra-vertex geometry/graph/charge checker.

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
TAG='fresh_raw76_forced_extra_vertex_family_20261018'
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

def run():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    assert d['status']=='PASS_RAW76_FORCED_EXTRA_VERTEX_BOUNDED_FAMILY_PILOT'
    for p,h in d['sources'].items():assert digest(ROOT/p)==h,p
    model=json.loads((ROOT/f"results/{d['current_global_model']}.json").read_text())
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={i:F(x)for i,x in dual['face_multipliers']};mu={i:F(x)for i,x in dual['vertex_multipliers']}
    geo=d['core_geometry'];g=mat(geo['scaled_gram']);native=geo['native_gram'];scale=[Q(*x)for x in geo['squared_scales']]
    assert len(g)==len(native)==len(scale)==7
    assert all(len(row)==7 for row in g+native)
    assert all(g[i][k]==g[k][i]and native[i][k]==native[k][i]for i in range(7)for k in range(7))
    for i in range(7):
        assert g[i][i]==2*scale[i]
        for k in range(i):assert g[i][k].sign()<=0 and g[i][k]*g[i][k]==squared(native[i][k])*scale[i]*scale[k]
    supports=geo['support_records'];assert [r['omitted']for r in supports]==list(map(list,combinations(range(7),2)))
    pd=[]
    for r in supports:
        ids=sorted(set(range(7))-set(r['omitted']));assert r['support']==ids
        p=positive(principal(g,ids));assert p==(r['inertia']==[5,0,0])
        if p:pd.append(r['omitted'])
    assert pd==[[0,1],[0,3],[1,2],[2,3],[2,5]]
    assert geo['quad_fixed']==[4,5,6]and geo['extra_support']==[0,1,3,4,6]
    t=residual(g,geo['triangle_fixed'],geo['triangle_sides']);assert t==mat(geo['triangle_residual'])
    assert inverse(t)==mat(geo['triangle_inverse'])and all(x.sign()<0 for r in inverse(t)for x in r)
    assert geo['triangle_supports']==[sorted(geo['triangle_fixed']+list(s))for s in combinations(geo['triangle_sides'],2)]
    assert all(sorted(set(range(7))-set(s))in pd for s in geo['triangle_supports'])
    autos=[p for p in permutations(range(7))if all(native[i][j]==native[p[i]][p[j]]for i,j in combinations(range(7),2))]
    assert autos==[tuple(range(7))]
    family=d['family'];assert len(family)==len({r['target_vertex']for r in family})==235
    forests=json.loads((ROOT/'results/fresh_sqrt5_rank13_combined_disjoint_ports_20261017.json').read_text())
    catalog={i:f for i,f in enumerate(forests['faces'])if f['kind']=='canonical_dotted_quad_aggregate'and any(ci==65 for _,ci in f['source_contexts'])}
    assert set(catalog)=={r['face']for r in family}
    targets={r['index']:r['forest']for r in forests['vertices']}
    provenance={f:p for p,f in model['prepared']['provenance_to_face']}
    for r in family:
        assert catalog[r['face']]['source_contexts']==[[r['source_context'],65]]
        assert r['provenance']==provenance.get(r['face'])
        assert sorted(r['spectator_forest']+['A1','A1','H3'])==r['target_forest']==targets[r['target_vertex']]
        assert F(r['current_mass'])==alpha.get(r['face'],0)and r['actual_extra_occurrences']==1
        c=Counter(r['target_forest']);assert r['marker_orientations']==c['A1']*(c['A1']-1)*c['H3']
    assert sum(r['marker_orientations']for r in family)==3204
    assert sum(r['marker_orientations']*(r['marker_orientations']-1)//2 for r in family)==39966
    assert sum(F(r['current_mass'])>0 for r in family)==16
    assert sum((F(r['current_mass'])for r in family),F())==F('57192559/5000000000')
    assert [r['face']for r in d['pilots']]==[r['face']for r in sorted(family,key=lambda r:(-F(r['current_mass']),r['face']))[:3]]
    counts=Counter(supports=21,positive_supports=5,family_classes=235)
    for r in d['pilots']:
        frame,j=frame_matrix(r);nodes=r['old_nodes'];groups=r['old_groups'];markers=r['markers'];oldn=len(nodes)
        owner={i:z['index']for z in groups for i in z['literal_stars']};assert set(owner)==set(range(oldn))
        assert [z['index']for z in groups]==list(range(len(groups)))
        expected=next(x['marker_orientations']for x in family if x['face']==r['face']);assert len(markers)==expected
        seen=set()
        for marker in markers:
            mp=dict(marker['active_embedding']);assert set(mp)=={0,1,3,4,6}and len(set(mp.values()))==5
            assert all(frame[mp[i]][mp[k]]==-root(squared(native[i][k]))for i,k in combinations(mp,2))
            assert all(frame[i][k]==0 for i in mp.values()for k in set(range(13))-set(mp.values()))
            cov=[[[mp[k],native[i][k]]for k in sorted(mp,key=lambda k:mp[k])if native[i][k][0]!=2]for i in(2,5)]
            cov=sorted(cov);assert cov==marker['native']['contacts']
            key=json.dumps(cov);assert key not in seen;seen.add(key)
            assert marker['triangle_ports']==sorted([mp[1],mp[3]])
            code,hs,hp,f,p,z=wall_pair(cov[0],cov[1],j)
            assert code is None and all(x==0 for x in z)
            assert hs==[Q(*x)for x in marker['native']['squared_timelike_heights']]
            assert hp==Q(*marker['height_product']);counts['marker_maps']+=1
        expressions=r['wall_expressions'];checked=[]
        for e in expressions:
            c,d0=e['contacts'];code,hs,hp,_,_,_=wall_pair(c,d0,j)
            assert code==e['code']and hs==[Q(*x)for x in e['heights']]and hp==Q(*e['height_product'])
            checked.append((c,d0,code));counts['wall_expressions']+=1
        new=r['new_pair_decisions'];nn=oldn+len(markers)
        assert [x['pair']for x in new]==[list(z)for z in combinations(range(nn),2)if z[1]>=oldn]
        oldedges={tuple(x)for x in r['old_graph_edges']};edges=set(oldedges)
        allnodes=nodes+[z['native']for z in markers]
        for x in new:
            i,k=x['pair'];bad=x['exclusion']
            if bad is None:edges.add((i,k));counts['retained_unknowns']+=1;continue
            if bad['kind']=='actual_triangle_port':
                assert i<oldn and bad['ports']==markers[k-oldn]['triangle_ports']and bad['ports']in nodes[i]['occupied_pairs']
            else:
                assert bad['kind']=='physical_product';a,b=bad['wall_indices']
                c,d0,code=checked[bad['expression']]
                assert sorted([allnodes[i]['contacts'][a],allnodes[k]['contacts'][b]])==[c,d0]
                assert code==bad['code']and code is not None
            counts['pair_exclusions']+=1
        assert edges==set(map(tuple,r['graph_edges']))and {e for e in edges if max(e)<oldn}==oldedges
        adj={i:set()for i in range(nn)}
        for i,k in edges:assert i<k;adj[i].add(k);adj[k].add(i)
        actual=cliques(adj);assert actual==sorted(map(tuple,r['maximal_cliques']))
        counts['maximal_cliques']+=len(actual)
        assert max(sum(i>=oldn for i in cl)for cl in actual)==r['marker_capacity']
        bysource={z['source']:z for z in r['charge_envelopes']}
        for env in r['charge_envelopes']:
            w=list(map(F,env['weights']));assert len(w)==len(groups)and min(w)>=0
            maxima={}
            for cl in actual:
                k=sum(i>=oldn for i in cl);x=sum((w[owner[i]]for i in cl if i<oldn),F());assert x<=1
                maxima[k]=max(maxima.get(k,F()),x)
            assert maxima=={z['markers']:F(z['old_charge'])for z in env['count_maxima']}
        for ch in r['charges']:
            env=bysource[ch['old_charge_source']];w=list(map(F,env['weights']));lam,t=F(ch['lambda_']),F(ch['marker_weight'])
            assert lam>=0 and t>0
            for z in env['count_maxima']:assert lam*F(z['old_charge'])+t*z['markers']<=1;counts['charge_inequalities']+=1
            cf=Counter({r['face']:t})
            for a,group in zip(w,groups):
                for f,n in group['positive_face_coefficients']:cf[f]+=lam*a*n
            assert [[f,str(a)]for f,a in sorted(cf.items())if a]==ch['positive_face_coefficients']
            load=sum((alpha.get(f,F())*a for f,a in cf.items()),F())
            assert load==F(ch['load'])and mu.get(r['vertex'],F())==F(ch['baseline_cost'])
            assert load-mu.get(r['vertex'],F())==F(ch['violation'])and ch['vertex_cost']=='1'
            counts['charges']+=1
        assert r['best_charge']==max(r['charges'],key=lambda c:F(c['violation']))
    assert counts['retained_unknowns']==6 and counts['marker_maps']==22
    assert counts['pair_exclusions']==3809 and counts['maximal_cliques']==12804
    assert [F(r['best_charge']['violation'])for r in d['pilots']]==[F('-4039181/10000000000'),F('241899/80000000'),F('23181393/20000000000')]
    assert d['optimizer_invocations']==0 and d['new_dimension_bound']is None and d['current_graphs_changed']is False
    print('PASS_RAW76_EXTRA_VERTEX_SAVED_CHECK',dict(counts))

if __name__=='__main__':run()

