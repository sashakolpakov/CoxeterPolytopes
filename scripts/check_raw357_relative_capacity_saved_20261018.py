#!/usr/bin/env python3
"""Independent saved raw357 relative-gate/capacity/ordinary transport checks.

No producer imports. Historical incidence, source-family and ordinary
owner/graph completeness are inherited. Contacted blocks are inverted
exactly; unused passive blocks are harmless diagonal placeholders only
with zero-contact checks. No charge or global-matrix assertion.
"""
import json
from collections import Counter,deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
REL='fresh_raw357_relative_marker_pilot_20261018'
ORD='fresh_raw357_capacity_ordinary_pilot_20261018'
INC='fresh_raw357_extra_vertex_incidence_pilot_20261018'
OLD=('fresh_raw76_extra_vertex_complete_family_20261018','fresh_raw56_raw57_extra_vertex_family_20261018',
 'fresh_raw108_complete_marker_family_20261018','fresh_raw112_complete_marker_family_20261018')
NEW=(3574,3575,3576)
def read(s):return json.loads((ROOT/f'results/{s}.json').read_text())
def key(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
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


def rotate_pair(u,v):
    choices=[]
    for _ in range(4):choices.append((u,v));u,v=v,tuple(reversed(u))
    return min(choices)


def rule(adj,a,b):
    ka,kb=a['kind'],b['kind']
    if ka>kb:return rule(adj,b,a)
    x,y=active_parts(adj,a),active_parts(adj,b)
    same=rotate_pair(x[0],x[-1])==rotate_pair(y[0],y[-1])
    disjoint=set(x[0]+x[-1]).isdisjoint(y[0]+y[-1])
    if ka==kb==3574:return same or(disjoint and x[1]==y[1])
    if ka==kb in(3575,3576):return x[1:3]==y[1:3]and(same or disjoint)
    return same if(ka,kb)==(3574,3575)else False


def normalize(native,p):
    co=[sorted([[p[i],e]for i,e in c])for c in native['contacts']]
    order=sorted(range(2),key=lambda i:co[i])
    return dict(contacts=[co[i]for i in order],squared_timelike_heights=[native['squared_timelike_heights'][i]for i in order],
        native_omitted_product=native['native_omitted_product'],occupied_pairs=sorted(sorted(p[i]for i in c)for c in native['occupied_pairs']))


def main():
    data,ordata,inc=read(REL),read(ORD),read(INC)
    assert data['status']=='PASS_RAW357_COMPLETE_RELATIVE_MARKER_PILOT'and ordata['status']=='PASS_RAW357_COMPLETE_CAPACITY_ORDINARY_PILOT'
    for d in(data,ordata):
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
        assert d['exported_columns']==[]and d['optimizer_invocations']==0 and not d['old_outcomes_changed']
        assert d['source_quad_demand']==1 and d['extra_occurrences_per_source']==3
    assert data['family']==inc['complete_family_inventory']and data['pointed_family']==inc['extra_marker_inventory']
    assert len(data['family'])==47 and len(data['pointed_family'])==141
    rows={r['vertex']:r for r in data['forest_records']};ordinary={r['vertex']:r for r in ordata['forest_records']}
    assert len(rows)==len(data['forest_records'])==len(ordinary)==len(ordata['forest_records'])==138 and set(rows)==set(ordinary)
    expected={(r['target_vertex'],3570+r['omitted'][1]):r for r in inc['extra_marker_inventory']}
    assert len(expected)==len(inc['extra_marker_inventory'])==141
    observed=[(v,f['kind'])for v,r in rows.items()for f in r['marker_families']]
    assert len(observed)==len(set(observed))==141 and set(observed)==set(expected)
    for v,r in rows.items():
        for f in r['marker_families']:
            e=expected[v,f['kind']]
            assert r['forest']==e['target_forest']and f['source_face']==e['face']
            assert len(f['markers'])==e['active_role_orientations']and f['actual_extra_occurrences']==e['actual_extra_occurrences']==1
    previous={v:{}for v in rows}
    for tag in OLD:
        prior=read(tag)
        for r in prior['forest_records']:
            if r['vertex']not in rows:continue
            fs=[{'kind':76,'markers':r['markers']}]if tag==OLD[0]else r['marker_families']
            for f in fs:
                k=f.get('kind',f.get('raw_family'));assert k not in previous[r['vertex']]
                previous[r['vertex']][k]=f['markers']
            assert r['finite_frame']==rows[r['vertex']]['finite_frame']
        del prior
    cases={key([c['kinds'][0],c['kinds'][1],c['pattern']]):c for c in data['relative_cases']}
    assert len(cases)==len(data['relative_cases'])==578
    used=Counter();comparisons=native_count=anchor_count=unknowns=retained=gate_count=0
    ng=inc['geometry']['native_gram'];extra={3570+x['omitted'][1]:x for x in inc['geometry']['extra_vertices']}
    for v,row in rows.items():
        adj,scales,j,active,parts=frame(row['finite_frame']);fs={f['kind']:f for f in row['marker_families']}
        assert len(fs)==len(row['marker_families'])
        allkeys=set()
        for kind,fam in fs.items():
            ex=extra[kind];ids=ex['support'];om=ex['omitted']
            act={i:{k:ng[i][k][0]for k in ids if i!=k and ng[i][k][0]!=2}for i in ids}
            ms=fam['markers'];got={tuple(map(tuple,m['active_embedding']))for m in ms}
            assert len(got)==len(ms)and got==maps(act,adj)
            assert [m['index']for m in ms]==list(range(len(ms)))
            for m in ms:
                assert m['kind']==kind;mp=dict(m['active_embedding'])
                cov=[sorted([[mp[k],ng[i][k]]for k in ids if ng[i][k][0]!=2])for i in om];cov.sort()
                assert cov==m['native']['contacts']and cov[0]!=cov[1]
                assert all(i in active for c in cov for i,e in c)
                code,hs,hp,field,pr,z=physical(cov[0],cov[1],j,scales)
                assert all(a==-b for a,b in zip(z,field.sqrt(squared(ng[om[0]][om[1]]))))
                assert hs==[Q(*h)for h in m['native']['squared_timelike_heights']]
                assert m['native']['native_omitted_product']==ng[om[0]][om[1]]and m['native']['occupied_pairs']==[]
                ports=sorted(sorted((mp[i],mp[k]))for i,k in ex['forced_triangle_port_pairs'])
                assert ports==m['triangle_port_pairs']and len(ports)==6
                k=key(cov);assert k not in allkeys;allkeys.add(k);native_count+=1
            assert fam['actual_extra_occurrences']==1
        assert not allkeys&{key(m['native']['contacts'])for ms in previous[v].values()for m in ms}
        for ka,fam in fs.items():
            anchor=fam['markers'][0];others={k:f['markers']for k,f in fs.items()if k>=ka}|previous[v]
            for kb,ms in others.items():
                for m in ms:
                    if ka==kb and m['index']==0:continue
                    k,pattern,payload=relative(adj,ka,anchor,kb,m);assert k in cases and payload==cases[k]['payload']
                    used[k]+=1;comparisons+=1
        o=ordinary[v];assert o['forest']==row['forest']and o['finite_frame']==row['finite_frame']
        nodes=o['ordinary_nodes'];owner={int(i):k for i,k in o['old_owner'].items()}
        groupowner={i:g['index']for g in o['ordinary_groups']for i in g['literal_stars']}
        assert owner==groupowner and set(owner)==set(range(len(nodes)))
        assert all(0<=i<k<len(nodes)for i,k in o['old_graph_edges'])
        canonical=[normalize(n,dict(zip(range(13),range(13))))for n in nodes]
        mm={f['kind']:f for f in o['marker_families']};assert set(mm)==set(fs)and len(mm)==len(o['marker_families'])
        for z in o['generators']:
            p=dict(z['frame_permutation']);assert set(p)==set(p.values())==set(range(13))
            assert all(adj[i].get(k,2)==adj[p[i]].get(p[k],2)for i,k in combinations(range(13),2))
            action=z['ordinary_permutation'];assert sorted(action)==list(range(len(nodes)))
            assert all(normalize(n,p)==canonical[action[i]]and owner[i]==owner[action[i]]for i,n in enumerate(nodes))
        for kind,fam in fs.items():
            ms=fam['markers'];q=mm[kind];assert q['marker_count']==len(ms)and q['actual_extra_occurrences']==1 and q['source_face']==fam['source_face']
            actions=q['generator_marker_permutations'];assert len(actions)==len(o['generators'])
            for act,gen in zip(actions,o['generators']):
                assert sorted(act)==list(range(len(ms)));p=dict(gen['frame_permutation'])
                for i,m in enumerate(ms):
                    assert normalize(m['native'],p)==ms[act[i]]['native']
                    assert sorted(sorted(p[k]for k in pair)for pair in m['triangle_port_pairs'])==ms[act[i]]['triangle_port_pairs']
            reached={0};todo=[0]
            while todo:
                i=todo.pop()
                for action in actions:
                    k=action[i]
                    if k not in reached:reached.add(k);todo.append(k)
            assert reached==set(range(len(ms)))
            decisions=q['anchor_decisions'];assert [x['ordinary']for x in decisions]==list(range(len(nodes)))
            for rec in decisions:
                i=rec['ordinary'];result=rec['result']
                if result.get('triangle_port_exclusion')is not None:
                    assert all(p in ms[0]['triangle_port_pairs']for p in result['triangle_port_exclusion'])
                assert all(i in active for n in(nodes[i],ms[0]['native'])for c in n['contacts']for i,e in c)
                gate_count+=pair_result(nodes[i],ms[0]['native'],result,j,scales)
            n=sum(x['result']['allowed']for x in decisions)
            assert q['anchor_unknowns']==n and q['retained_ordinary_pairs']==n*len(ms)
            anchor_count+=len(decisions);unknowns+=n;retained+=n*len(ms)
        co=Counter(row['forest']);caps={r['kind']:r for r in o['capacity']['individual']};assert set(caps)==set(fs)
        for kind,c in caps.items():
            b,h=co['A2'],co['H4'];expected=4*max(h,b//2)if kind==3574 else 4*(b//2)if kind==3575 else 4*((b-1)//2)
            chosen=c['witness'];assert c['capacity']==expected==len(chosen)==len(set(chosen))
            assert all(0<=i<len(fs[kind]['markers'])for i in chosen)
            for i,k in combinations(chosen,2):assert rule(adj,fs[kind]['markers'][i],fs[kind]['markers'][k])
        joint=max(c['capacity']for c in caps.values())
        if 3574 in caps and 3575 in caps:
            assert co['A2']==2 and co['H3']==co['H4']==1;joint=max(joint,caps[3574]['capacity']+caps[3575]['capacity'])
        assert joint==o['capacity']['joint_new_marker_upper']<=8
    assert set(used)==set(cases)and comparisons==11557 and native_count==10624
    for k,c in cases.items():
        assert used[k]==c['anchor_comparisons'];adj,sc,j,active,parts=frame(c['finite_frame'])
        assert relative(adj,*[c['kinds'][0],c['anchor'],c['kinds'][1],c['other']])==(k,c['pattern'],c['payload'])
        assert all(i in active for m in(c['anchor'],c['other'])for cv in m['native']['contacts']for i,e in cv)
        gate_count+=pair_result(c['anchor']['native'],c['other']['native'],c['result'],j,sc)
        if c['kinds'][1]in NEW:assert rule(adj,c['anchor'],c['other'])==c['result']['allowed']
    assert sum(c['result']['allowed']for c in cases.values())==27
    assert (anchor_count,unknowns,retained)==(7282,192,20848)
    print('PASS_INDEPENDENT_RAW357_RELATIVE_CAPACITY_ORDINARY',dict(native_maps=native_count,relative_cases=len(cases),anchor_comparisons=comparisons,ordinary_anchors=anchor_count,retained_unknowns=unknowns,retained_pairs=retained,wall_expressions=gate_count),flush=True)

if __name__=='__main__':main()
