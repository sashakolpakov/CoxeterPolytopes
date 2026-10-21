#!/usr/bin/env python3
"""All finite four-special corner patterns, with free opposite Gram products."""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_mixed_d4_d5_migration_20261017 as E
import fresh_sqrt5_d4_short_face_incidence_20261017 as P
import fresh_sqrt5_sealed_subcore_closure_20261015_review as J
A=P.A;K,Z,O=P.K,P.Z,P.O
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_mixed_d4_d5_quad_patterns_20261017'
INPUT=f'results/{E.TAG}.json'

def encode(p):return[[list(k),str(v)]for k,v in sorted(p.items())]
def decode(p):return{tuple(k):F(v)for k,v in p}
def swapped(p):return{(j,i):v for(i,j),v in p.items()}

def literal(g,i,j,ai,aj,bond):
    n=len(g);h={x:set(v)for x,v in g.items()};h[n]=set(ai);h[n+1]=set(aj)
    for x in ai:h[x].add(n)
    for x in aj:h[x].add(n+1)
    if bond:h[n].add(n+1);h[n+1].add(n)
    ids=[x for x in h if x not in(i,j)]
    opposite=E.classify({x:h[x]&set(ids)for x in ids})
    if opposite is None:return None
    order=[i,j,n,n+1]+[x for x in g if x not in(i,j)];mp={x:k for k,x in enumerate(order)}
    known={tuple(sorted((mp[x],mp[y])))for x in h for y in h[x]}
    row=dict(nodes=n+2,white=list(range(4,n+2)),labels=[0 if p in((0,2),(1,3))else(3 if p in known else 2)for p in combinations(range(n+2),2)])
    return row,opposite,order

def symbolic(row):
    n=row['nodes'];labels=dict(zip(combinations(range(n),2),row['labels']))
    g=[[K(2*(i==j))for j in range(n)]for i in range(n)]
    for(i,j),m in labels.items():
        if m==3:g[i][j]=g[j][i]=-O
    white=list(range(4,n));inv=A.inverse(A.sub(g,white))
    contacts=[[g[i][w]for w in white]for i in range(4)]
    residual=[]
    for i in range(4):
        line=[]
        for j in range(4):
            q=g[i][j]-A.dot(contacts[i],A.mv(inv,contacts[j]));assert not q.b
            v=P.pc(q.a)
            if tuple(sorted((i,j)))==(0,2):v=P.plus(v,{(1,0):F(1)})
            if tuple(sorted((i,j)))==(1,3):v=P.plus(v,{(0,1):F(1)})
            line.append(v)
        residual.append(line)
    for i,j in[(0,1),(1,2),(2,3),(3,0)]:
        a=residual[i][i].get((0,0),F(0));b=residual[j][j].get((0,0),F(0));c=residual[i][j].get((0,0),F(0))
        assert a>0 and b>0 and a*b-c*c>0
    determinant=P.determinant(residual)
    # The square classes are read from the graph of fixed label-three edges.
    unseen=set(range(n));components=[];which={}
    while unseen:
        group={min(unseen)};stack=list(group)
        while stack:
            i=stack.pop()
            for j in range(n):
                if i!=j and labels.get(tuple(sorted((i,j))))==3 and j not in group:group.add(j);stack.append(j)
        for i in group:which[i]=len(components)
        unseen-=group;components.append(sorted(group))
    u_in=which[0]==which[2];v_in=which[1]==which[3]
    coupled=not u_in and not v_in and sorted([which[0],which[2]])==sorted([which[1],which[3]])
    if u_in and v_in:kind='BOTH_PRODUCTS_IN_K'
    elif u_in or v_in:kind='ONE_PRODUCT_IN_K_ONE_SQUARE_IN_K'
    elif coupled:kind='SQUARES_AND_MUTUAL_PRODUCT_IN_K'
    else:kind='TWO_INDEPENDENT_SQUARES_IN_K'
    if not u_in and not coupled:assert all(a%2==0 for a,b in determinant)
    if not v_in and not coupled:assert all(b%2==0 for a,b in determinant)
    return dict(residual=[[encode(p)for p in r]for r in residual],determinant=encode(determinant),
        known_components=components,field_kind=kind,u_in_K=u_in,v_in_K=v_in,uv_in_K=coupled or(u_in and v_in),
        rank_condition='The displayed determinant must vanish; physical inertia must be(2,1,1) on this4-root residual, and conjugate PSD of rank3 must hold in a valid ground-field gauge. These extra conditions are not solved by this pattern census.')

def run():
    source=json.loads((ROOT/INPUT).read_text())
    for p,h in source['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    targets=[]
    for target in source['targets']:
        names=target['forest'];g,_=E.forest(names);n=len(g);wanted=Counter(names);rows=[];maps=[];buckets=defaultdict(list);pairs=set();counts=Counter()
        for i,j in combinations(range(n),2):
            for ai,aj in product(target['ports'][i]['attachments'],target['ports'][j]['attachments']):
                for bond in(0,1):
                    hit=literal(g,i,j,ai['contacts'],aj['contacts'],bond)
                    if hit is None or Counter(c['type']for c in hit[1])!=wanted:counts['REJECT_FOURTH_VERTEX_TYPE']+=1;continue
                    row,opposite,order=hit;s=symbolic(row);key=J.invariant(J.matrix(row));found=None
                    for idx in buckets[key]:
                        mp=J.isomorphism(J.matrix(row),J.matrix(rows[idx]['row']))
                        if mp is not None:
                            J.check_map(J.matrix(row),J.matrix(rows[idx]['row']),mp)
                            old=decode(rows[idx]['symbolic']['determinant']);now=decode(s['determinant'])
                            swap=tuple(sorted((mp[0],mp[2])))==(1,3)
                            assert (swapped(now)if swap else now)==old
                            found=(idx,mp,swap);break
                    if found is None:
                        idx=len(rows);rows.append(dict(index=idx,row=row,symbolic=s));buckets[key].append(idx)
                        found=(idx,list(range(n+2)),False)
                    pairs.add((i,j));counts[s['field_kind']]+=1
                    maps.append(dict(ports=[i,j],first_contacts=ai['contacts'],second_contacts=aj['contacts'],
                        first_copy=ai['copy'],second_copy=aj['copy'],new_new_label=3 if bond else 2,
                        original_order=order,target=found[0],map=found[1],swap_free_products=found[2]))
            if i==0:print('MIXED_QUAD_PAIR',target['name'],i,j,'classes',len(rows),flush=True)
        targets.append(dict(name=target['name'],dimension=n,forest=names,literal_patterns=len(maps),colored_patterns=len(rows),
            field_counts=dict(counts),possible_corners=[list(p)for p in sorted(pairs)],corner_count=len(pairs),records=rows,maps=maps))
        print('MIXED_QUAD_TARGET',target['name'],len(maps),len(rows),len(pairs),dict(counts),flush=True)
    paths=[f'scripts/{TAG}.py',INPUT,f'scripts/{E.TAG}.py',f'scripts/{P.TAG}.py',f'results/{P.TAG}.json',f'scripts/{J.TAG}.py']
    return dict(status='PASS_COMPLETE_SELECTED_MIXED_D_QUADRILATERAL_FINITE_PATTERNS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},targets=targets,
        variable_convention='u=Gram(black0,black2), v=Gram(black1,black3), norm-two roots. Label0 in the exported rows denotes a free opposite product, not an assertion that it is dotted. The new-new adjacent product is fully enumerated as0 or-1.',
        completeness='Every all-special quadrilateral has one of the already complete same-type attachments at each of two adjacent edges, and its fourth finite corner is the same target forest. Every pair of those literal attachments and both possible adjacent new-new products are tested. Full black-opposite/white marked isomorphisms and parameter swaps are recorded.',
        field_scope='A fixed-edge path makes its framing product a signed cyclic product, hence integral in K. Otherwise the two-cycle only places its square in K; when both free edges close the same contracted cycle their mutual product is also in K. These are necessary field/integrality conditions, not an admissible parameter enumeration.',
        scope='Exact finite corner patterns and symbolic rank equations only, for four specified forests. No framing products have been specialized. Lorentz signature, conjugate positivity, chamber coherence and Coxeter/arithmetic parameter solutions remain to be classified. No mixed-type incidence inequality or dimension bound is claimed.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
