#!/usr/bin/env python3
"""Independent saved triangle-domain check; Fraction only, no project imports.

Finite diagrams use rooted-tree canonical codes. Exact PSD tests use reverse
Schur order in Q(sqrt2,sqrt5); the historical ambient geometry is not rerun.
"""
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_compact29_three_triangle_star_constraint_20261018'
def K(x=0):return tuple(map(F,x))if isinstance(x,(tuple,list))else(F(x),F(),F(),F())
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def neg(x):return tuple(-v for v in x)
def sub(x,y):return add(x,neg(y))
def mul(x,y):
    z=[F()]*4
    for i,a in enumerate(x):
        for j,b in enumerate(y):z[i^j]+=a*b*(2 if i&j&1 else 1)*(5 if i&j&2 else 1)
    return tuple(z)
def inv(x):
    y=(x[0],x[1],-x[2],-x[3]);u=mul(x,y);assert not u[2]and not u[3]
    den=u[0]**2-2*u[1]**2;assert den
    return mul(y,(u[0]/den,-u[1]/den,F(),F()))
def sign2(a,b):
    sa,sb=(a>0)-(a<0),(b>0)-(b<0)
    if not sb:return sa
    if not sa or sa==sb:return sb
    z=a*a-2*b*b;return sa*((z>0)-(z<0))
def sign(x):
    a,b,c,d=x;u,v=sign2(a,b),sign2(c,d)
    if not v:return u
    if not u or u==v:return v
    z=sub(mul((a,b,F(),F()),(a,b,F(),F())),mul(K(5),mul((c,d,F(),F()),(c,d,F(),F()))))
    return u*sign2(z[0],z[1])
COS={2:K(),3:K(1),4:K((0,1,0,0)),5:K((F(1,2),0,F(1,2),0))}
def psd(g):
    a=[r[::-1]for r in g[::-1]];zero=False
    for i in range(len(a)):
        s=sign(a[i][i])
        if s<0:return-1
        if not s:
            if any(any(a[i][j])for j in range(i+1,len(a))):return-1
            zero=True;continue
        p=inv(a[i][i])
        for j in range(i+1,len(a)):
            for k in range(i+1,len(a)):a[j][k]=sub(a[j][k],mul(mul(a[j][i],p),a[i][k]))
    return 0 if zero else 1
def rational_inverse(g):
    n=len(g);a=[r[:]+[F(i==j)for j in range(n)]for i,r in enumerate(g)]
    for i in range(n):
        p=a[i][i];assert p>0;a[i]=[v/p for v in a[i]]
        for j in range(n):
            if i!=j:
                z=a[j][i];a[j]=[x-z*y for x,y in zip(a[j],a[i])]
    return[r[n:]for r in a]
def tree_code(nodes,edges):
    adj={i:[]for i in nodes}
    for(i,j),m in edges.items():
        if i in adj and j in adj:adj[i].append((j,m));adj[j].append((i,m))
    if sum(map(len,adj.values()))!=2*(len(nodes)-1):return None
    def rooted(i,p):return tuple(sorted((m,rooted(j,i))for j,m in adj[i]if j!=p))
    return min(rooted(i,None)for i in nodes)
def catalogue():
    out=set()
    def put(n,e):out.add((n,tree_code(list(range(n)),e)))
    for n in range(1,8):put(n,{(i,i+1):3 for i in range(n-1)})
    for n in range(2,8):put(n,{(i,i+1):4 if i==0 else 3 for i in range(n-1)})
    for n in range(4,8):put(n,{**{(i,i+1):3 for i in range(n-3)},(n-3,n-2):3,(n-3,n-1):3})
    for n in(6,7):put(n,{**{(i,i+1):3 for i in range(n-2)},(2,n-1):3})
    put(4,{(0,1):3,(1,2):4,(2,3):3})
    for n in(2,3,4):put(n,{(i,i+1):5 if i==0 else 3 for i in range(n-1)})
    return out
CAT=catalogue()
def finite(nodes,edges):
    todo=set(nodes)
    while todo:
        part={min(todo)};stack=list(part)
        for i in stack:
            for(a,b),m in edges.items():
                j=b if a==i else a if b==i else None
                if j in todo and j not in part:part.add(j);stack.append(j)
        todo-=part
        if(len(part),tree_code(sorted(part),edges))not in CAT:return False
    return True
def main():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    assert d['status']=='PASS_COMPACT29_UNIVERSAL_THREE_TRIANGLE_STAR_CONSTRAINT'
    assert len(d['sources'])==8
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    e={(i,i+1):3 for i in range(5)};e[2,6]=3
    g=[[F(2*(i==j))for j in range(7)]for i in range(7)]
    for i,j in e:g[i][j]=g[j][i]=F(-1)
    J=rational_inverse(g)
    assert [[K(x)for x in row]for row in g]==[[K(x)for x in row]for row in d['E7_Gram']]
    assert [[K(x)for x in row]for row in J]==[[K(x)for x in row]for row in d['E7_inverse']]
    assert d['marked_pairs']==[[2,4],[2,5],[3,5]]and len(d['complete_triangle_domains'])==3
    total=Counter()
    for marks,row in zip(d['marked_pairs'],d['complete_triangle_domains'],strict=True):
        a,b=marks;assert row['marked']==marks and(a,b)not in e
        records={tuple(r['labels']):r for r in row['finite_eligible']};assert len(records)==len(row['finite_eligible'])
        found=set();rejected=0;ledger=sha256();counts=Counter()
        for ms in product(range(2,6),repeat=7):
            ee=e|{(i,7):m for i,m in enumerate(ms)if m>2}
            corners=[finite([i for i in range(8)if i!=omit],ee)for omit in marks]
            ledger.update(json.dumps([ms,corners],separators=(',',':')).encode()+b'\n')
            if not all(corners):rejected+=1;continue
            found.add(ms);r=records[ms];p=[neg(COS[m])for m in ms];h=K(-2)
            for i,j in product(range(7),repeat=2):h=add(h,mul(K(J[i][j]),mul(p[i],p[j])))
            assert h==K(r['height_squared'])and r['corner_finite']==corners
            if sign(h)<=0:assert r['status']=='nonpositive_height'
            else:
                gg=[[K(2*(i==j))for j in range(8)]for i in range(8)]
                for(i,j),m in ee.items():gg[i][j]=gg[j][i]=neg(COS[m])
                if r['status']=='parabolic_obstruction':
                    ns=r['parabolic_subset'];assert len(ns)==len(set(ns))and set(ns)<=set(range(8))
                    assert psd([[gg[i][j]for j in ns]for i in ns])==0
                else:
                    assert r['status']=='unique_necessary_triangle_root'and ms==(2,2,2,2,2,5,2)
                    assert h==K((F(1,4),0,F(3,4),0))
                    for size in range(1,9):
                        for ns in combinations(range(8),size):assert psd([[gg[i][j]for j in ns]for i in ns])!=0
            counts[r['status']]+=1
        assert found==records.keys()and rejected==row['finite_rejections']
        assert len(found)+rejected==row['all_low_vectors']==16384
        assert counts==row['status_counts']and ledger.hexdigest()==row['finite_filter_ledger_sha256']
        total.update(counts)
    assert total=={'nonpositive_height':9,'parabolic_obstruction':14,'unique_necessary_triangle_root':3}
    h=K(d['forced_height_squared']);z=mul(h,h);assert z==K(d['cross_factor_squared_product'])==K((F(23,8),0,F(3,8),0))
    assert sign(sub(z,K(F(37,10))))>0 and sign(sub(K(F(93,25)),z))>0
    gap=d['native_gap'];assert gap['labels']==[11,12]
    def atan(x):
        value=sum((-1)**j*x**(2*j+1)/F(2*j+1)for j in range(30));return value,value+x**61/F(61)
    a,b=atan(F(1,5));c,e0=atan(F(1,239));pi=(16*a-4*e0,16*b-4*c)
    def cos(x,n):
        term=value=F(1)
        for j in range(1,n+1):term=-term*x*x/F((2*j-1)*2*j);value+=term
        return value
    for m in(11,12):
        lo,hi=map(F,gap[f'label{m}_squared']);x,y=2*pi[0]/m,2*pi[1]/m
        assert 0<x<y<1 and lo<=2+2*cos(y,17)<=2+2*cos(x,16)<=hi
    assert F(gap['label11_squared'][1])<F(37,10)and F(93,25)<F(gap['label12_squared'][0])
    assert [d[k]for k in('all_low_vectors','total_finite_eligible','total_nonpositive_heights','total_parabolic_obstructions','total_unique_roots')]==[49152,26,9,14,3]
    assert [d[k]for k in('triangle_factor_capacity','designated_angles','designated_triangle_upper','designated_nonsimplex_faces_lower')]==[1,9,3,6]
    assert not d['geometric_scope']['new_face_charge_saving_proved']and d['geometric_scope']['new_compact_dimension_bound']is None
    print('PASS49,152 exact finite-filter vectors;26 eligible;3 unique triangles;coupled capacity1')
if __name__=='__main__':main()
