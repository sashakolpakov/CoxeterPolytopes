#!/usr/bin/env python3
"""Triangle payments with finite distance8 and ideal distance5/feature5."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_ideal6_triangle_payments_20261017 as T

ROOT=T.ROOT
TAG='fresh_ideal5_finite8_triangle_payments_20261017'

def inverse(matrix):
    n=len(matrix);a=[list(map(F,r))+[F(i==j)for j in range(n)]for i,r in enumerate(matrix)]
    for i in range(n):
        p=a[i][i];assert p>0;a[i]=[v/p for v in a[i]]
        for j in range(n):
            if j!=i:
                q=a[j][i];a[j]=[x-q*y for x,y in zip(a[j],a[i])]
    return [r[n:]for r in a]

def corner(adj,labels,marks,a,b):
    nodes=set(adj)-(set(marks)-{a,b});routes=T.paths(adj,nodes,a)
    p=routes.get(b)
    if p is None:return dict(pair=[a,b],path=None,length=None,feature_witness=None,elementary_payment=False)
    nodes=set(routes);features={i for i in nodes if len(adj[i]&nodes)>=3}
    for (i,j),m in labels.items():
        if m>3 and i in nodes and j in nodes:features.update((i,j))
    witness=None
    for start in(a,b):
        routes=T.paths(adj,nodes,start)
        for f in sorted(features):
            if len(routes[f])-1<=5:
                witness=dict(endpoint=start,feature=f,path=routes[f],
                    feature_neighbors=sorted(adj[f]&nodes),
                    incident_labels=[[j,labels[min(f,j),max(f,j)]]for j in sorted(adj[f]&nodes)])
                break
        if witness:break
    distance=len(p)-1
    return dict(pair=[a,b],path=p,length=distance,feature_witness=witness,
        elementary_payment=distance<=5 or(distance<=8 and witness is not None))

def exceptional(diagram):
    adj,labels=T.graph(diagram);assert set(labels.values())=={3}
    g=[[2 if i==j else -int(j in adj[i])for j in range(10)]for i in range(10)]
    affine=[r[:9]for r in g[:9]];kernel=T.positive_kernel(affine)
    pivots=T.positive_pivots([r[:8]for r in g[:8]])
    e8=[r[:8]for r in g[:8]];J=inverse(e8)
    assert all(x>=0 for row in J for x in row)
    diag=[J[i][i]for i in range(8)];assert diag==list(map(F,[30,8,14,4,20,12,6,2]))
    assert all(diag[i]>2 for i in range(7))
    rows=[]
    specs=[((0,1,9),(0,9),3,2,'affine_E8_kernel'),
           ((0,2,9),(0,9),1,0,'affine_E8_kernel'),
           ((1,2,8),(1,8),1,0,'finite_E8_inverse'),
           ((1,2,9),(1,9),1,0,'affine_E8_kernel')]
    for marks,pair,leaf,neighbor,method in specs:
        nodes=set(range(10))-(set(marks)-set(pair));path=T.paths(adj,nodes,pair[0])[pair[1]]
        assert len(path)-1<=8
        assert adj[leaf]&nodes=={neighbor}and len(adj[neighbor]&nodes)==2
        assert leaf in T.paths(adj,nodes,pair[0])and labels[min(leaf,neighbor),max(leaf,neighbor)]==3
        protected=list(range(9))if method=='affine_E8_kernel'else list(range(7))
        assert leaf in protected and neighbor in protected
        core_corner=sorted(nodes);cp=T.positive_pivots([[g[i][j]for j in core_corner]for i in core_corner])
        rows.append(dict(marks=list(marks),method=method,protected_core_nodes=protected,
            enforcing_corner=sorted(set(marks)-{9 if method=='affine_E8_kernel'else 8}),
            alternate_pair=list(pair),alternate_path=path,finite_corner_pivots=cp,
            fixed_simple_leaf=leaf,fixed_degree_two_neighbor=neighbor,
            neighbor_edges=[[i,labels[min(i,neighbor),max(i,neighbor)]]for i in sorted(adj[neighbor]&nodes)],
            maximum_exceptional_affine_order=9))
    return dict(rows=rows,affine_kernel=list(map(str,kernel)),affine_nodes=list(range(9)),
        positive_principal_pivots=pivots,E8_nodes=list(range(8)),E8_inverse=[[str(v)for v in r]for r in J],
        E8_inverse_diagonal=list(map(str,diag)),
        argument='Affine E8 forces all extra supporting roots orthogonal to nodes0..8. For the finite E8 case, nonnegative inverse entries and diagonals>2 force orthogonality to nodes0..6. The selected simple leaf and its degree2 neighbor persist at the alternate corner. That component is finite (distance<=8 pays) or exceptional affine of order<=9 (whole-factor pays); no classical affine graph has this fixed leaf pattern.')

def run():
    data=json.loads((ROOT/f'results/{T.CENSUS}.json').read_text());rows=[];exceptions=[]
    for diagram in data['diagrams']:
        if diagram['order']<7:continue
        adj,labels=T.graph(diagram)
        for marks in combinations(range(diagram['order']),3):
            cs=[corner(adj,labels,marks,a,b)for a,b in combinations(marks,2)]
            paid=sum(c['elementary_payment']for c in cs)
            rows.append(dict(diagram=diagram['index'],order=diagram['order'],marks=list(marks),
                corners=cs,elementary_payment_count=paid))
            if paid<2:exceptions.append([diagram['index'],list(marks)])
    assert len(rows)==1025
    assert Counter(r['order']for r in rows)=={7:105,8:224,9:336,10:360}
    assert exceptions==[[25,[0,1,9]],[25,[0,2,9]],[25,[1,2,8]],[25,[1,2,9]]]
    special=exceptional(next(d for d in data['diagrams']if d['index']==25))
    assert [r['marks']for r in special['rows']]==[r[1]for r in exceptions]
    for r in rows:
        if r['elementary_payment_count']<2:
            e=next(e for e in special['rows']if e['marks']==r['marks'])
            assert r['elementary_payment_count']==1
            assert not next(c for c in r['corners']if c['pair']==e['alternate_pair'])['elementary_payment']
    stems=[T.CENSUS,T.TAG,'fresh_cusp_factor_payment_framework_20261017']
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_FINITE8_IDEAL5_FEATURE5_TRIANGLE_PAYMENTS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        rule=dict(finite_distance_at_most=8,ideal_distance_at_most=5,
            ideal_feature_distance_at_most=5,ideal_feature_pair_distance_at_most=8,
            ideal_whole_factor_order_at_most=9),rows=rows,
        counts_by_order={str(k):v for k,v in Counter(r['order']for r in rows).items()},
        elementary_exceptions=exceptions,exceptional_resolution=special,
        lower_order_argument='In a connected core of order<=6, a Steiner tree through the three marks has at least two marked-pair paths avoiding the third mark; these have at most5 edges.',
        scope='Universal triangle lower payment2 with finite8 and ideal5/feature5/wholefactor<=9, including arbitrary supporting whites and nonsimple ideal vertices. No quadrilateral coverage, ideal vertex budget or dimension improvement asserted.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['rows']),'triples',len(d['elementary_exceptions']),'E8 resolutions',flush=True)
