#!/usr/bin/env python3
"""Independent stdlib checker for the saved three-family envelope.

Reconstructs eight new cross-product gates with ordered multiquadratic
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
TAG='fresh_raw56_raw57_raw76_joint_marker_family_20261018'
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


def cross_pattern(old,new,raw):
    p=dict(old['active_embedding']);q=dict(new['active_embedding'])
    if raw==57:return ()
    names={p[0]:0,p[4]:1}
    for i in (q[0],q[5]):
        if i not in names:names[i]=len(names)
    return names[q[0]],names[q[5]]


def edge(ra,x,rb,y):
    if ra>rb:return edge(rb,y,ra,x)
    a,b=dict(x['active_embedding']),dict(y['active_embedding'])
    if ra==rb==76:
        return a[0]==b[0] and a[4]!=b[4] and tuple(a[i]for i in(1,3,6))!=tuple(b[i]for i in(1,3,6))
    if rb==76:return ra==56 and not({a[0],a[5]}&{b[0],b[4]})
    if ra!=rb:return False
    t1,t2=(a[1],a[3]),(b[1],b[3]);p1,p2=(a[0],a[5]),(b[0],b[5])
    if ra==57:return p1==p2 and set(t1)!=set(t2)
    if p1==p2:return set(t1)!=set(t2)
    return not(set(p1)&set(p2)) and (set(t1)!=set(t2) or t1==t2)


def coefficient_vertices(counts,have):
    # Independent intersection enumeration for the nonnegative 56/76 plane.
    cap={r:max(n[i]for n in counts)for i,r in enumerate((56,57,76))}
    if 56 not in have or 76 not in have:
        return [tuple(F(1,cap[r])if r in have else F()for r in(56,57,76))]
    planes=set((F(n[0]),F(n[2]),F(1))for n in counts if n[0]or n[2])
    planes|={(F(1),F(),F()),(F(),F(1),F())}
    candidates=set()
    for p,q in combinations(sorted(planes),2):
        det=p[0]*q[1]-p[1]*q[0]
        if not det:continue
        x=(p[2]*q[1]-p[1]*q[2])/det;y=(p[0]*q[2]-p[2]*q[0])/det
        if min(x,y)>=0 and all(x*n[0]+y*n[2]<=1 for n in counts):candidates.add((x,y))
    return sorted((x,F(1,cap[57])if 57 in have else F(),y)for x,y in candidates
        if not any(u>=x and v>=y and(u,v)!=(x,y)for u,v in candidates))


def check():
    data=read(TAG);assert data['status']=='PASS_RAW56_RAW57_RAW76_JOINT_COMPLETE_FAMILY'
    for p,h in data['sources'].items():assert digest(ROOT/p)==h,p
    ast.parse((ROOT/f'scripts/{TAG}.py').read_text())
    old,new,model=read(OLDMARKER),read(NEW),read(MODEL)
    oo={r['vertex']:r for r in old['forest_records']};nn={r['vertex']:r for r in new['forest_records']}
    assert len(oo)==len(old['forest_records'])==235 and len(nn)==len(new['forest_records'])==1212
    assert len(oo.keys()&nn.keys())==53
    cases=expressions=0
    for test in data['relative_tests']:
        raw=test['new_raw'];v=test['host_vertex'];o=oo[v];n=nn[v]
        assert test['finite_frame']==o['finite_frame']==n['finite_frame'] and test['host_forest']==o['forest']
        anchor=test['raw76_anchor'];assert anchor==o['markers'][0]
        ms=next(f['markers']for f in n['marker_families']if f['raw_family']==raw)
        adj,scales,j,active,components=frame(test['finite_frame'])
        expected={cross_pattern(anchor,m,raw)for m in ms};seen=set()
        for rec in test['records']:
            m=rec['new_marker'];assert m==ms[m['index']]
            p=cross_pattern(anchor,m,raw);assert list(p)==rec['pattern'] and p not in seen;seen.add(p)
            assert all(i in active for z in(anchor,m)for c in z['native']['contacts']for i,e in c)
            expressions+=pair_result(anchor['native'],m['native'],rec['result'],j,scales)
            assert rec['result']['allowed']==(raw==56 and p==(2,3));cases+=1
        assert seen==expected and len(seen)==(7 if raw==56 else 1)
    assert cases==data['relative_cases']==8 and data['retained_patterns']==1
    # Import saved complete ordinary vectors only. No historical graph census.
    vectors={};weights=defaultdict(dict);cache={}
    for tag in(LOCAL,PATCH,THIRD):
        for r in read(tag)['forest_records']:
            weights[r['vertex']][tag]=list(map(F,r['certificate']['weights']))
            if 'maximal_count_vectors'in r:vectors[r['vertex']]=r['maximal_count_vectors']
    for r in read(LOCAL)['forest_records']:
        if r['vertex']in vectors:continue
        tag=r['geometry_source']
        if tag not in cache:cache[tag]={z['vertex']:z for z in read(tag)['forest_records']}
        vectors[r['vertex']]=cache[tag][r['vertex']]['maximal_count_vectors']
    rows={r['vertex']:r for r in data['forest_records']}
    assert len(rows)==len(data['forest_records'])==1394 and set(rows)==oo.keys()|nn.keys()
    charges=defaultdict(list)
    for c in data['exported']['columns']:charges[c['vertex']].append(c)
    assert set(charges)==set(rows)
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']};mu={v:F(x)for v,x in dual['vertex_multipliers']}
    base={v:F(x)for v,x in new['exported']['repair']['vertex_multipliers']};repaired=dict(base)
    assert dual['face_multipliers']==new['exported']['repair']['face_multipliers']==data['exported']['repair']['face_multipliers']
    totals=Counter();totalcols=0;globalkeys=set()
    for v,rec in rows.items():
        o,n=oo.get(v),nn.get(v);r=n or o
        if o and n:
            for name in('forest','finite_frame','old_nodes','old_groups','old_graph_edges','prior_history'):assert o[name]==n[name]
        assert rec['forest']==r['forest']
        for name,field in(('old_nodes','old_nodes_sha256'),('old_groups','old_groups_sha256'),('old_graph_edges','old_graph_sha256')):assert rec[field]==jhash(r[name])
        fam={f['raw_family']:f['markers']for f in n['marker_families']}if n else{}
        faces={f['raw_family']:f['face']for f in n['marker_families']}if n else{}
        cross=[(i,k)for f in n['marker_families']for i,k in f['marker_ordinary_edges']]if n else[]
        if o:fam[76]=o['markers'];faces[76]=o['face'];cross+=o['marker_ordinary_edges']
        assert rec['raw_families']==sorted(fam) and rec['source_faces']==[list(z)for z in sorted(faces.items())]
        assert rec['source_records']==[[OLDMARKER if raw==76 else NEW,v,raw]for raw in sorted(fam)]
        keys=[key(m['native']['contacts'])for ms in fam.values()for m in ms]
        assert len(set(keys))==len(keys) and not set(keys)&{key(s['contacts'])for s in r['old_nodes']}
        cc=Counter(r['forest']);a,h,t=cc['A1'],cc['H3'],cc['I2(10)']
        states=[]
        if 56 in fam:states.append((a//2*t,0,0))
        if 76 in fam:
            for k in range(1,min(a-1,h)+1):states.append((((a-k-1)//2*t)if 56 in fam else 0,0,k))
        if 57 in fam:states.append((0,t,0))
        env=rec['envelope'];assert env['raw_order']==[56,57,76] and env['role_counts']==dict(A1=a,H3=h,I210=t)
        assert env['boundary_counts']==[list(s)for s in states]
        ww=coefficient_vertices(states,fam);assert [list(map(str,w))for w in ww]==env['coefficient_vertices']
        assert len(env['attaining_graph_cliques'])==len(states)
        for s,witness in zip(states,env['attaining_graph_cliques']):
            assert witness['counts']==list(s);q=witness['marker_clique'];assert len(q)==len(set(map(tuple,q)))
            assert Counter(raw for raw,i in q)==Counter({raw:z for raw,z in zip((56,57,76),s)if z})
            assert all(0<=i<len(fam[raw])for raw,i in q)
            assert all(edge(ra,fam[ra][i],rb,fam[rb][j])for(ra,i),(rb,j)in combinations(q,2))
        pairs=retained=0
        if o and n:
            for raw in fam:
                if raw==76:continue
                for x in fam[76]:
                    for y in fam[raw]:pairs+=1;retained+=edge(76,x,raw,y)
        assert pairs==rec['cross_marker_positions'] and retained==rec['retained_cross_marker_pairs']
        totals.update(cross_marker_positions=pairs,retained_cross_marker_pairs=retained)
        groups=r['old_groups'];owner={i:g['index']for g in groups for i in g['literal_stars']}
        assert len(owner)==len(r['old_nodes'])
        for g in groups:assert Counter(x['face']for x in g['origins'])==Counter(dict(g['positive_face_coefficients']))
        omitted=sorted({owner[i]for i,k in cross});assert omitted==rec['ordinary_omitted_orbits'] and len(cross)==rec['marker_ordinary_unknowns']
        vv=vectors[v]if groups else[[]]
        candidates=weights[v]if groups else{'marker_only':[]};expected={}
        for source,original in candidates.items():
            original=original+[F()]*(len(groups)-len(original));assert len(original)==len(groups)and min(original,default=F())>=0
            assert all(len(z)==len(original)and sum(w*k for w,k in zip(original,z))<=1 for z in vv)
            totals['ordinary_vector_checks']+=len(vv)
            kept=[F()if i in omitted else x for i,x in enumerate(original)]
            for mw in ww:
                cf=Counter()
                for w,g in zip(kept,groups):
                    for f,nc in g['positive_face_coefficients']:cf[f]+=w*nc
                for raw,z in zip((56,57,76),mw):
                    if z:cf[faces[raw]]+=z
                k=tuple((f,z)for f,z in sorted(cf.items())if z)
                expected.setdefault(k,[]).append((source,original,kept,mw))
        assert len(charges[v])==len(expected)
        seen=set();loads=[]
        for c in charges[v]:
            kk=tuple((f,F(x))for f,x in c['positive_face_coefficients']);assert kk in expected and kk not in seen;seen.add(kk)
            source,original,kept,mw=expected[kk][0]
            assert c['source_charge']==source and c['other_charge_sources']==[s for s,*_ in expected[kk][1:]]
            assert list(map(F,c['old_orbit_weights']))==original and list(map(F,c['retained_orbit_weights']))==kept
            assert c['omitted_orbits']==omitted and c['old_vector_count']==len(vv)and c['old_vectors_sha256']==jhash(vv)
            assert c['marker_terms']==[[raw,faces[raw],str(z)]for raw,z in zip((56,57,76),mw)if z]
            assert all(sum(z*n for z,n in zip(mw,s))<=1 for s in states)
            load=sum((alpha.get(f,F())*z for f,z in kk),F());loads.append(load)
            assert c['vertex_cost']=='1' and load==F(c['load']) and F(c['violation'])==load-base.get(v,F())
            assert F(c['old_model_cost'])==mu.get(v,F()) and F(c['baseline_cost'])==base.get(v,F())
            globalkeys.add((v,kk));totalcols+=1
        repaired[v]=max([base.get(v,F())]+loads)
        assert F(rec['baseline_cost'])==base.get(v,F()) and F(rec['new_cost'])==repaired[v]
        assert F(rec['increase'])==repaired[v]-base.get(v,F())
    assert len(globalkeys)==totalcols==1430
    preserved=data['exported']['preserved_columns'];assert len(preserved)==1454
    idx=0
    for tag,obj in((OLDMARKER,old),(NEW,new)):
        for c in obj['exported']['columns']:
            q=preserved[idx];idx+=1;v=c['vertex'];load=sum((alpha.get(f,F())*F(z)for f,z in c['positive_face_coefficients']),F())
            assert q==dict(source=tag,vertex=v,load=str(load),cost=str(repaired.get(v,F())))
            assert load<=base.get(v,F())<=repaired.get(v,F())
    costs=data['exported']['preserved_multi_owner_costs'];inputs=model['prepared']['inherited_multi_vertex_columns'];assert len(costs)==len(inputs)==14
    for c,q in zip(inputs,costs):
        load=sum((alpha.get(f,F())*n for f,n in c['integer_face_coefficients']),F())
        before=sum((mu.get(v,F())*n for v,n in c['integer_vertex_capacities']),F());after=sum((repaired.get(v,F())*n for v,n in c['integer_vertex_capacities']),F())
        assert q==dict(variable=c['variable'],all_vertex_costs=c['integer_vertex_capacities'],load=str(load),old_cost=str(before),new_cost=str(after))and load<=before<=after
    repair=data['exported']['repair'];assert repair['vertex_multipliers']==[[v,str(x)]for v,x in sorted(repaired.items())if x]
    num=F(dual['objective_numerator']);den=sum(repaired.values(),F())
    assert F(repair['objective_numerator'])==num and F(repair['normalizer'])==den and F(repair['bound'])==num/den>12
    assert repair['baseline_bound']==new['exported']['repair']['bound']
    support=read('fresh_raw76_complete_family_support_20261018');remaining=[r for r in support['ranked_unambiguous_missing_families']if r['raw_family']not in(56,57)and r['missing_classes']>0]
    sd=data['remaining_support_diagnostic'];assert len(remaining)==sd['unrepresented_unambiguous_families']==211
    assert sum(r['missing_classes']for r in remaining)==sd['unrepresented_unambiguous_classes']==1213
    assert sd['highest_saved_demand_families']==remaining[:10]and F(sd['added_normalizer_needed_to_reach12'])==num/12-den>0
    count=data['counts'];assert count['cross_marker_positions']==totals['cross_marker_positions']==30336
    assert count['retained_cross_marker_pairs']==totals['retained_cross_marker_pairs']==9888
    assert count['targets']==1394 and count['overlap_targets']==53 and count['family_classes']==1643 and count['marker_orientations']==24788
    assert count['exported_columns']==totalcols and count['strictly_violated_columns']==sum(F(c['violation'])>0 for c in data['exported']['columns'])==3
    assert count['raised_allowances']==sum(F(r['increase'])>0 for r in rows.values())==2
    print('PASS_RAW56_RAW57_RAW76_SAVED_CHECK',dict(relative_cases=cases,wall_expressions=expressions,
        targets=len(rows),columns=totalcols,**totals),repair['bound'],flush=True)

if __name__=='__main__':check()
