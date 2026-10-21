#!/usr/bin/env python3
"""Saved first-root gate checker; Fraction arithmetic, no project imports."""
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_compact29_coupled_first_roots_20261018'
d=json.loads((ROOT/f'results/{TAG}.json').read_text())
assert d['status']=='PASS_FIXED_COMPACT29_COUPLED_FIRST_ROOT_GATES'
for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
assert len(d['sources'])==7

def k(x=0):return tuple(map(F,x))if isinstance(x,(tuple,list))else(F(x),F(),F(),F())
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def scale(x,s):return tuple(a*s for a in x)
def mul(x,y):
    z=[F()]*4
    for i,a in enumerate(x):
        for j,b in enumerate(y):z[i^j]+=a*b*(2 if i&j&1 else 1)*(5 if i&j&2 else 1)
    return tuple(z)
def enclosure(x):
    # Independent rational enclosures with explicit square checks.
    a,b=F(1414213562,10**9),F(1414213563,10**9)
    c,e=F(2236067977,10**9),F(2236067978,10**9)
    assert a*a<2<b*b and c*c<5<e*e
    bs=[(F(1),F(1)),(a,b),(c,e),(a*c,b*e)]
    return(sum(v*(lo if v>=0 else hi)for v,(lo,hi)in zip(x,bs)),
           sum(v*(hi if v>=0 else lo)for v,(lo,hi)in zip(x,bs)))
cos={2:k(),3:k(1),4:k((0,1,0,0)),5:k((F(1,2),0,F(1,2),0))}

def inverse(g):
    n=len(g);a=[list(row)+[F(i==j)for j in range(n)]for i,row in enumerate(g)]
    for j in range(n):
        assert a[j][j]>0 # Symmetric finite E-type elimination pivots.
        z=a[j][j];a[j]=[v/z for v in a[j]]
        for i in range(n):
            if i!=j:
                z=a[i][j];a[i]=[v-z*w for v,w in zip(a[i],a[j])]
    out=[row[n:]for row in a]
    assert all(sum(g[i][j]*out[j][t]for j in range(n))==int(i==t)for i in range(n)for t in range(n))
    return out

old=json.loads((ROOT/'results/fresh_compact29_small_faces_20261015.json').read_text())
expected=[(i,row['recipe'])for i,row in enumerate(old['cases'])if row['recipe']['type']in('E7','E8')]
assert len(expected)==len(d['roots'])==49
heights={};kinds={}
for rid,((idx,recipe),row)in enumerate(zip(expected,d['roots'],strict=True)):
    assert row['id']==rid and row['source_index']==idx and row['recipe']==recipe
    n=int(recipe['type'][1:]);g=[[F(2*int(i==j))for j in range(n)]for i in range(n)]
    for i in range(n-2):g[i][i+1]=g[i+1][i]=F(-1)
    g[2][n-1]=g[n-1][2]=F(-1);j=inverse(g)
    contacts={i:m for i,m in recipe.get('contacts',recipe.get('c_contacts',[]))}
    if recipe['kind']=='quadrilateral':contacts[recipe['marked'][0]]=recipe['ac_label']
    p=[scale(cos[contacts[i]],-1)if i in contacts else k()for i in range(n)]
    assert [k(x)for x in row['contacts']]==p
    h=k(-2)
    for a,b in product(range(n),repeat=2):h=add(h,scale(mul(p[a],p[b]),j[a][b]))
    assert h==k(row['height_squared'])and enclosure(h)[0]>0
    heights[rid]=h;kinds[rid]=recipe['type']

# Build literal placements on four distinguished component slots. Their
# quotient, rather than a hand-picked pair list, must equal the saved domain.
slots=[[i for i in heights if kinds[i]==t]for t in('E7','E7','E7','E8')]
multiplicities=Counter()
for a,b in combinations(range(4),2):
    for i,j in product(slots[a],slots[b]):
        if kinds[i]==kinds[j]:i,j=sorted((i,j))
        multiplicities[kinds[i]+'/'+kinds[j],i,j]+=1
assert len(multiplicities)==819 and sum(multiplicities.values())==3087
assert len(d['pairs'])==819
badheight=k((F(1,4),0,F(3,4),0));badproduct=k((F(23,8),0,F(3,8),0))
assert mul(badheight,badheight)==badproduct
bad=literal_bad=0;seen=set()
for row in d['pairs']:
    key=(row['kind'],row['left'],row['right']);assert key not in seen;seen.add(key)
    assert row['literal_multiplicity']==multiplicities[key]
    z=mul(heights[row['left']],heights[row['right']]);assert z==k(row['squared_mutual_product'])
    if row['status']=='strictly_dotted_first_root_pair':assert enclosure(z)[0]>4
    else:
        assert row['status']=='excluded_native_gap_11_12'
        assert heights[row['left']]==heights[row['right']]==badheight and z==badproduct
        bad+=1;literal_bad+=row['literal_multiplicity']
assert seen==multiplicities.keys()and(bad,literal_bad)==(15,75)

# Separate rational Taylor bounds verify the saved native-angle enclosures
# contain the true values; use more atan terms than the producer.
def atan(x,n=30):
    s=sum((-1)**j*x**(2*j+1)/F(2*j+1)for j in range(n));e=x**(2*n+1)/F(2*n+1)
    return s,s+e
a,b=atan(F(1,5));c,e=atan(F(1,239));pi=(16*a-4*e,16*b-4*c)
def cosine(x,n):
    term=F(1);s=term
    for j in range(1,n+1):term=-term*x*x/F((2*j-1)*2*j);s+=term
    return s
gap=d['native_gap'];assert(gap['lower_label'],gap['upper_label'])==(11,12)
for m in(11,12):
    lo,hi=map(F,gap[f'label{m}_squared']);x,y=2*pi[0]/m,2*pi[1]/m
    assert 0<x<y<1
    true_lo,true_hi=2+2*cosine(y,17),2+2*cosine(x,16)
    assert lo<=true_lo<=true_hi<=hi
assert F(gap['label11_squared'][1])<F(37,10)<enclosure(badproduct)[0]
assert enclosure(badproduct)[1]<F(93,25)<F(gap['label12_squared'][0])
assert d['excluded_E7_marked_pairs']==[[2,4],[2,5],[3,5],[4,5],[5,6]]
assert k(d['excluded_height_squared'])==badheight and k(d['excluded_squared_product'])==badproduct
assert d['forest']==['E7','E7','E7','E8']and d['rank']==29
assert [d[k]for k in('root_count','quotient_pairs','literal_pairs','excluded_quotient_pairs','excluded_literal_pairs','dotted_literal_pairs')]==[49,819,3087,15,75,3012]
assert d['common_vertex_orientation_required']and d['first_root_gates_only']
assert not d['all_star_recipes_classified']and not d['second_root_gates_checked']
assert d['new_compact_dimension_bound']is None and d['new_finite_volume_example']is False
print('PASS49 exact projections,819 quotient/3087 literal gates,75 exclusions; no full-star or global claim')
