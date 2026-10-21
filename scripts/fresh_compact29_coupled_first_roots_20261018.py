#!/usr/bin/env python3
"""Exact necessary cross gates for49 fixed E7/E8 small-face recipes."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement,product
import fresh_compact29_small_faces_20261015 as P
ROOT=P.ROOT
TAG='fresh_compact29_coupled_first_roots_20261018'

def atan_bounds(x,n=24):
    s=sum(((-1)**j*x**(2*j+1)/F(2*j+1)for j in range(n)),F())
    e=x**(2*n+1)/F(2*n+1)
    return(s,s+e)if n%2==0 else(s-e,s)

def cosine_square_bounds(m):
    a,b=atan_bounds(F(1,5));c,d=atan_bounds(F(1,239))
    lo,hi=2*(16*a-4*d)/m,2*(16*b-4*c)/m
    assert 0<lo<hi<1
    def partial(x,n):
        term=F(1);s=term
        for j in range(1,n+1):
            term=-term*x*x/F((2*j-1)*2*j);s+=term
        return s
    # Cosine decreases on this interval; odd/even last degree give bounds.
    return 2+2*partial(hi,13),2+2*partial(lo,12)

def run():
    path=ROOT/'results/fresh_compact29_small_faces_20261015.json'
    old=json.loads(path.read_text())
    assert old['status']=='PASS_EXACT_LOCAL_SMALL_FACE_BARRIER_79_MARKED_PAIRS'
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    roots=[]
    for index,row in enumerate(old['cases']):
        r=row['recipe'];kind=r['type']
        if kind not in('E7','E8'):continue
        n,e=P.diagram(kind);G=P.base_gram(n,e,n);J=P.inverse(G)
        contacts={i:m for i,m in(r['contacts']if r['kind']=='triangle'else r['c_contacts'])}
        if r['kind']=='quadrilateral':contacts[r['marked'][0]]=r['ac_label']
        c=[-P.COS[contacts[i]]if i in contacts else P.K()for i in range(n)]
        height=P.dot(c,c,J)-2;assert height.sign()>0
        roots.append(dict(id=len(roots),source_index=index,recipe=r,contacts=c,height_squared=height))
    e7=[r for r in roots if r['recipe']['type']=='E7'];e8=[r for r in roots if r['recipe']['type']=='E8']
    assert(len(e7),len(e8))==(21,28)
    badheight=P.K((F(1,4),0,F(3,4),0));badproduct=badheight*badheight
    assert badproduct==P.K((F(23,8),0,F(3,8),0))
    low11,up11=cosine_square_bounds(11);low12,up12=cosine_square_bounds(12)
    assert up11<F(37,10) and(P.K(F(37,10))-badproduct).sign()<0
    assert(badproduct-P.K(F(93,25))).sign()<0 and F(93,25)<low12
    pairs=[]
    domain=[('E7/E7',a,b,3*(1 if a['id']==b['id']else 2))for a,b in combinations_with_replacement(e7,2)]
    domain +=[('E7/E8',a,b,3)for a,b in product(e7,e8)]
    for kind,a,b,mult in domain:
        z=a['height_squared']*b['height_squared']
        if(z-4).sign()>0:status='strictly_dotted_first_root_pair'
        else:
            assert a['height_squared']==b['height_squared']==badheight
            assert z==badproduct;status='excluded_native_gap_11_12'
        pairs.append(dict(kind=kind,left=a['id'],right=b['id'],literal_multiplicity=mult,
                          squared_mutual_product=z,status=status))
    assert len(pairs)==819 and sum(r['literal_multiplicity']for r in pairs)==3087
    bad=[r for r in pairs if r['status'].startswith('excluded')]
    assert len(bad)==15 and sum(r['literal_multiplicity']for r in bad)==75
    marked=[r['recipe']['marked']for r in e7 if r['height_squared']==badheight]
    assert marked==[[2,4],[2,5],[3,5],[4,5],[5,6]]
    paths=[f'scripts/{TAG}.py','scripts/fresh_compact29_small_faces_20261015.py',
      'results/fresh_compact29_small_faces_20261015.json','results/fresh_compact29_small_faces_20261015.md',
      'results/fresh_compact29_small_faces_20261015_recipes.json',
      'results/fresh_compact29_savings_frontier_20261018.json','results/fresh_compact29_savings_frontier_20261018.md']
    return P.enc(dict(status='PASS_FIXED_COMPACT29_COUPLED_FIRST_ROOT_GATES',
      sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
      forest=['E7','E7','E7','E8'],rank=29,roots=roots,pairs=pairs,
      root_count=49,quotient_pairs=819,literal_pairs=3087,excluded_quotient_pairs=15,
      excluded_literal_pairs=75,dotted_literal_pairs=3012,
      excluded_E7_marked_pairs=marked,excluded_height_squared=badheight,excluded_squared_product=badproduct,
      native_gap=dict(lower_label=11,upper_label=12,label11_squared=[low11,up11],
        label12_squared=[low12,up12],product_lower='37/10',product_upper='93/25'),
      common_vertex_orientation_required=True,first_root_gates_only=True,
      all_star_recipes_classified=False,second_root_gates_checked=False,
      new_compact_dimension_bound=None,new_finite_volume_example=False,
      conclusion='The five specified frozen E7 recipes cannot be simultaneously used on two distinct E7 components at one finite vertex. Other first-root pairs are dotted; this does not certify complete stars or rule out alternative recipes.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert a.verify or not p.exists();out=run()
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['quotient_pairs'],out['excluded_literal_pairs'],flush=True)
