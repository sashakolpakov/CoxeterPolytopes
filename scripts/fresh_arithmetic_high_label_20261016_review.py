#!/usr/bin/env python3
"""Independent exact field, graph and counting checks for the high-bond lemma.

Only Python standard library; no producer or inherited graph imports.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import comb, gcd, lcm
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_arithmetic_high_label_20261016_review'
PRIMARY='results/fresh_arithmetic_high_label_20261016.json'


def sha(path): return sha256((ROOT/path).read_bytes()).hexdigest()
def edge(a,b): return tuple(sorted((a,b)))


def field_checks():
    def c2(a): return 1 if a%8 in (1,7) else -1
    def c3(a): return 1 if a%12 in (1,11) else -1
    def c5(a): return 1 if a%5 in (1,4) else -1
    def fixes_f(a): return c2(a)==c5(a)==1
    specs=[('Qsqrt3',12,2,lambda a:c3(a)==1),
           ('Qsqrt6',24,2,lambda a:c2(a)*c3(a)==1),
           ('Qsqrt2sqrt3',24,4,lambda a:c2(a)==c3(a)==1)]
    for m,degree in ((9,3),(11,5),(15,4),(16,4),(20,4)):
        specs.append((f'K{m}',m,degree,lambda a,m=m:a%m in (1,m-1)))
    records=[]
    for name,conductor,degree,fixes in specs:
        modulus=lcm(40,conductor)
        units=[a for a in range(modulus) if gcd(a,modulus)==1]
        field_kernel=[a for a in units if fixes(a)]
        assert len(units)==degree*len(field_kernel)
        assert len(units)==4*sum(fixes_f(a) for a in units)
        witness=next(a for a in units if fixes_f(a) and not fixes(a))
        other_modulus=lcm(7,conductor)
        other_units=[a for a in range(other_modulus) if gcd(a,other_modulus)==1]
        other_witness=next(a for a in other_units if fixes(a) and a%7 not in (1,6))
        records.append(dict(field=name,degree=degree,conductor=conductor,
                            not_subfield_F=dict(modulus=modulus,automorphism=witness,
                                                meaning='Fixes sqrt2 and sqrt5, moves the target field'),
                            does_not_contain_K7=dict(modulus=other_modulus,automorphism=other_witness,
                                                    meaning='Fixes the target field, moves K7')))
    return records


def connected(vertices,edges):
    seen={vertices[0]}
    while True:
        nxt=seen|{v for v in vertices if any(edge(v,w) in edges for w in seen if v!=w)}
        if nxt==seen: return len(seen)==len(vertices)
        seen=nxt


def forest(vertices,edges):
    # Independent union-find cycle test.
    parent={i:i for i in vertices}
    def find(i):
        while parent[i]!=i: i=parent[i]
        return i
    for a,b in combinations(vertices,2):
        if (a,b) not in edges: continue
        aa,bb=find(a),find(b)
        if aa==bb: return False
        parent[aa]=bb
    return True


def graph_checks():
    examined=0
    passing=set()
    # Four black nodes in cyclic order0,1,2,3; minimal opposite diagrams
    # contain their black pair and either no white node or one white node.
    shapes=[((0,2),(1,3)),((0,2),(1,3,4)),((0,2,4),(1,3)),
            ((0,2,4),(1,3,5)),((0,2,4),(1,3,4))]
    for left,right in shapes:
        vertices=sorted(set(left+right))
        pairs=list(combinations(vertices,2))
        whites=sorted(set(vertices)-{0,1,2,3})
        corners=[sorted(whites+[a,b]) for a,b in ((0,1),(1,2),(2,3),(3,0))]
        for mask in range(1<<len(pairs)):
            examined+=1
            edges={pair for i,pair in enumerate(pairs) if (mask>>i)&1}
            if not connected(left,edges) or not connected(right,edges): continue
            if not set(left)&set(right) and not any(edge(a,b) in edges for a in left for b in right): continue
            if any(edge(a,b) in edges or any(edge(a,w) in edges and edge(b,w) in edges for w in whites)
                   for a,b in ((0,1),(1,2),(2,3),(3,0))): continue
            if any(not forest(corner,edges) for corner in corners): continue
            candidates=[]
            for minimal in (left,right):
                if len(minimal)==2:
                    candidates.append([None])
                    continue
                high=[]
                for a,b in combinations(minimal,2):
                    if edge(a,b) not in edges: continue
                    if all(not {a,b}<=set(corner) or all(
                        sum(edge(v,w) in edges for w in corner if w!=v)==1 for v in (a,b))
                           for corner in corners): high.append(edge(a,b))
                candidates.append(high)
            for high_left,high_right in product(*candidates):
                assert left==(0,2,4) and right==(1,3,5)
                assert high_left==(0,2) and high_right==(1,3)
                assert (4,5) in edges
                assert all(edge(a,b) not in edges for a in (0,2) for b in (1,3,5))
                assert all(edge(a,b) not in edges for a in (1,3) for b in (0,2,4))
                assert edges&{(0,4),(2,4)} and edges&{(1,5),(3,5)}
                passing.add(tuple(sorted(edges)))
    assert examined==35904 and len(passing)==9
    return dict(examined=examined,surviving_supports=[[list(pair) for pair in e] for e in sorted(passing)],
                explanation='Both minimal diagrams are triangles on separate white nodes, distinguished bonds are opposite black pairs, and the only connection is the white-white bond.')


def scalar_checks():
    # z=2cos(2pi/9) lies in(1,2) and solves z^3-3z+1=0.
    polynomial=lambda z:z*z*z-3*z+1
    assert polynomial(F(3,2))==F(-1,8) and polynomial(F(2))==3
    assert 3*F(3,2)**2-3>0
    assert (2+F(3,2))/4==F(7,8)
    assert 1/(4*(1-F(7,8)))==2
    # sqrt3/2>3/4 proves dense-bond cost>2 for every m>=6.
    assert F(3,4)>F(3,4)**2 and 1/(2*(1-F(3,4)))==2
    # Inverse of [[1,-c],[-c,1]] has positive entries. Two contacts
    # each>=1/2 cost>=(1+c)/(2(1-c^2))=1/(2(1-c)).
    # Clearing denominators gives (1+c)(1-c)=1-c^2.
    expanded=[F(0)]*3
    for i,a in enumerate((F(1),F(1))):
        for j,b in enumerate((F(1),F(-1))): expanded[i+j]+=a*b
    assert expanded==[F(1),F(0),F(-1)]
    return dict(sparse=dict(cos_squared_lower='7/8',cost_strict_lower='2',min_label=9),
                dense=dict(cos_lower='3/4',cost_strict_lower='2',min_label=6,contacts=2),
                schur='Each residual white diagonal<-1; its off-diagonal is the original white pairing, with square<1. Therefore the residual2x2 block is negative definite.')


def forest_checks():
    records=[]
    for rank in range(1,31):
        families=['A']+(['D'] if rank>=4 else [])+(['E'] if rank in (6,7,8) else [])
        for family in families:
            if family=='A': edges={(i,i+1) for i in range(rank-1)}
            elif family=='D': edges={(i,i+1) for i in range(rank-2)}|{(rank-3,rank-1)}
            else: edges={(i,i+1) for i in range(rank-2)}|{(2,rank-1)}
            short={pair for pair in combinations(range(rank),2) if pair in edges or any(
                edge(pair[0],w) in edges and edge(pair[1],w) in edges for w in range(rank) if w not in pair)}
            expected=0 if rank==1 else 2*rank-3 if family=='A' else 2*rank-2
            assert len(short)==expected and expected<=2*rank-2
            records.append(dict(family=family,rank=rank,count=len(short)))
    assert 1+1==2*2-2 # one bonus on isolated I2(m), m>=6.
    # The formulas cover all ranks; B,F,H are paths and D,E have one fork.
    return dict(records=records,all_rank_identity='edges+(sum_v choose(degree(v),2))=2r-3 for paths,2r-2 for one three-valent fork',
                vertex_upper_bound='2n-2',high_I2_total=2)


def counting_checks():
    n=21
    coefficients=[]
    for i in range(11):
        coefficient=comb(n-i,2)+comb(i,2)-100
        assert coefficient==(10-i)*(11-i)
        coefficients.append(coefficient)
    assert coefficients[0]==110 and 5*50-comb(21,2)==40==2*21-2
    # Polynomial coefficient identities in u for all higher dimensions.
    def mul(a,b):
        out=[F(0)]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b): out[i+j]+=x*y
        return out
    for start,wanted in ((21,[0,5,F(1,2)]),(22,[2,F(11,2),F(1,2)])):
        a,b=([start-1,2],[start-5,2]) if start%2 else ([start,2],[start-6,2])
        difference=[v/8 for v in mul(a,b)]
        difference[0]-=2*start-2
        difference[1]-=4
        assert difference==wanted
    return dict(dimension21_h_coefficients=coefficients,strict_threshold=40,vertex_maximum=40,
                odd_difference=['0','5','1/2'],even_difference=['2','11/2','1/2'])


def run():
    primary=json.loads((ROOT/PRIMARY).read_text())
    for path,digest in primary['sources'].items(): assert sha(path)==digest,path
    fields=field_checks()
    graphs=graph_checks()
    assert graphs['surviving_supports']==primary['local']['surviving_supports']
    canonical={'Qsqrt3':'Q(sqrt3)','Qsqrt6':'Q(sqrt6)','Qsqrt2sqrt3':'Q(sqrt2,sqrt3)'}
    assert [canonical.get(row['field'],row['field']) for row in fields]==primary['fields']['new_fields']
    paths=[PRIMARY,f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_EIGHT_FIELD_HIGH_BOND_REVIEW',sources={p:sha(p) for p in paths},
                fields=fields,graphs=graphs,scalars=scalar_checks(),finite_forests=forest_checks(),
                counting=counting_checks(),
                scope='Independent finite arithmetic/graph checks and explicit mathematical proof review. Published Lanner classification, ground-field inheritance, compact face geometry and symmetric nonnegative h-vector remain theorem inputs. External validation and priority unresolved.')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=run()
    path=ROOT/'results'/f'{TAG}.json'
    if args.verify: assert result==json.loads(path.read_text())
    else: path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],len(result['fields']),result['graphs']['examined'])
