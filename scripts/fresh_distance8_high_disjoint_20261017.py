#!/usr/bin/env python3
"""Exact all-label exclusion of disjoint high cores at distance greater than8.

This is a local lemma, not a new ambient dimension theorem.  All overlap
and supporting-extension obligations are outside this package's scope.
"""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path
import fresh_selective_distance9_dimension108_20261017_review as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_distance8_high_disjoint_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
F=R.F


def strict_obstruction(a):
    """Nonnegative, strictly negative-square witness in an indefinite Z block."""
    for part in R.pieces(a):
        b=R.extract(a,part);o=R.obstruction(b)
        if not o:continue
        v=list(map(F,o['vector']));q=R.square(b,v)
        if q==0:
            bv=[sum((x*y for x,y in zip(row,v)),F())for row in b]
            k=next((i for i,x in enumerate(bv)if x<0),None)
            if k is None:continue
            v[k]+=-bv[k]/b[k][k];q=R.square(b,v)
        assert q<0 and all(x>=0 for x in v)
        z=[F(0)]*len(a)
        for i,x in zip(part,v):z[i]=x
        assert R.square(a,z)==q
        return dict(vector=list(map(str,z)),square=str(q))
    return None


def corner_obstruction(native,black,u,v):
    white=[i for i in range(len(native))if i not in black]
    for p,q in product(u,v):
        ids=sorted(white+[p,q]);a,sc=R.native_to_rational(R.extract(native,ids))
        o=strict_obstruction(a)
        if o:return dict(corner=[p,q],nodes=ids,gauge_squared=sc,**o)
    return None


def run():
    ds={d['index']:d for d in json.loads((ROOT/f'results/{CENSUS}.json').read_text())['diagrams']}
    graphs={i:R.bond_graph(d)for i,d in ds.items()}
    contexts=R.all_contexts(list(ds.values()))
    assert len(contexts)==384 and max(max(r['distances'])for r in contexts)==7

    @lru_cache(None)
    def deletion_diagonals(di,black,node):
        out=[]
        for b in black:
            ids=[i for i in range(len(graphs[di]))if i!=b]
            a,sc=R.native_to_rational(R.extract(graphs[di],ids));R.psd(a,True)
            k=ids.index(node);out.append(sc[k]*R.inverse(a)[k][k])
        return tuple(out)

    @lru_cache(None)
    def white_diagonal(di,black,node):
        ids=[i for i in range(len(graphs[di]))if i not in black]
        a,sc=R.native_to_rational(R.extract(graphs[di],ids));R.psd(a,True)
        k=ids.index(node);return sc[k]*R.inverse(a)[k][k]

    @lru_cache(None)
    def zero_cofactor(di,node):
        assert node in ds[di]['ideal_deletions']
        ids=[i for i in range(len(graphs[di]))if i!=node]
        a,sc=R.native_to_rational(R.extract(graphs[di],ids));R.psd(a)
        assert R.nullspace(a)
        # The inherited minimal-indefinite classification, together with
        # a finite black deletion, makes the full core nonsingular index1.
        return dict(diagram=di,node=node,deletion_gram=R.encode(a),
                    gauge_squared=sc,deletion_signature=R.inertia(a),
                    full_determinant=ds[di]['determinant'],inverse_diagonal='0')

    ports=[]
    for c in contexts:
        for node,d in enumerate(c['distances']):
            if d<1:continue
            vals=deletion_diagonals(c['diagram'],c['black'],node)
            ports.append(dict(diagram=c['diagram'],black=list(c['black']),node=node,distance=d,
                inverse_diagonals=list(map(str,vals)),allowed=max(vals)<2))
    allowed=[p for p in ports if p['allowed']]
    assert len(ports)==2404 and len(allowed)==754
    rows=[];good={};types={};counts=Counter()
    for x,y in product(allowed,repeat=2):
        if x['distance']+y['distance']<8:continue
        di,dj=x['diagram'],y['diagram'];u,v0=x['black'],y['black'];i,j=x['node'],y['node']
        native,mp=R.unite(graphs[di],graphs[dj],());native[i][mp[j]]=native[mp[j]][i]=-1
        v=[mp[k]for k in v0];black=u+v;white=[k for k in range(len(native))if k not in black]
        wp=white_diagonal(di,tuple(u),i)*white_diagonal(dj,tuple(v0),j)
        row=dict(left=[di,u,i],right=[dj,v0,j],white_inverse_product=str(wp))
        if wp>=1:
            a,sc=R.native_to_rational(R.extract(native,white));o=R.obstruction(a);assert o
            row['white_rejection']=dict(nodes=white,gauge_squared=sc,**o)
            counts['white_rejected']+=1;rows.append(row);continue
        maxima=[max(map(F,p['inverse_diagonals']))for p in(x,y)]
        labels=[3]
        if max(maxima)>=1 or 2*wp>=1:
            row['higher_label_gate']=dict(first_excluded_label=4,endpoint_maxima=list(map(str,maxima)),twice_white_product=str(2*wp))
        else:
            # For label>=5 the squared bond is >=(3+sqrt(5))/2>5/2.
            assert F(5,2)*max(maxima)>=2
            labels.append(4);counts['joins_allowing4']+=1
            row['higher_label_gate']=dict(first_excluded_label=5,endpoint_maxima=list(map(str,maxima)),five_halves_maximum=str(F(5,2)*max(maxima)))
        key=(di,tuple(u),dj,tuple(v0));good.setdefault(key,[]).append((i,j,tuple(labels)))
        row['labels']=labels;row['label_rejections']=[]
        for label in labels:
            g=[r.copy()for r in native];g[i][mp[j]]=g[mp[j]][i]=-(label-2)
            typ=tuple(sorted(((di,i),(dj,j))))+(label,)
            if typ not in types:
                if i in ds[di]['ideal_deletions'] or j in ds[dj]['ideal_deletions']:
                    side=0 if i in ds[di]['ideal_deletions']else 1
                    # Rank-one cross-block Schur subtraction is zero exactly.
                    types[typ]=dict(ports=[list(z)for z in typ[:2]],label=label,
                        method='zero_inverse_rank_one_schur',zero_endpoint=[di,i]if side==0 else[dj,j],
                        zero_cofactor=zero_cofactor(di,i)if side==0 else zero_cofactor(dj,j),
                        block_signatures=[[ds[di]['order']-1,1,0],[ds[dj]['order']-1,1,0]],
                        zero_inverse_diagonal='0',negative_index=2)
                else:
                    a,sc=R.native_to_rational(g);sig=R.inertia(a)
                    types[typ]=dict(ports=[list(z)for z in typ[:2]],label=label,
                        method='rational_inertia',gram=R.encode(a),gauge_squared=sc,signature=sig,negative_index=sig[1])
            cert=types[typ]
            if cert['negative_index']==2:
                row['label_rejections'].append(dict(label=label,kind='two_negative_directions',type=[*map(list,typ[:2]),label]))
                counts['single_label_index2']+=1
            else:
                assert cert['negative_index']==1
                ob=corner_obstruction(g,black,u,v);assert ob,(row,label)
                row['label_rejections'].append(dict(label=label,kind='actual_corner_not_psd',**ob))
                counts['single_label_bad_corner']+=1
        rows.append(row)
    multiple=[]
    for (di,u,dj,v0),edges in sorted(good.items()):
        base,mp=R.unite(graphs[di],graphs[dj],());black=list(u)+[mp[i]for i in v0]
        for size in range(2,len(edges)+1):
            for chosen in combinations(edges,size):
                g=[r.copy()for r in base]
                for i,j,_ in chosen:g[i][mp[j]]=g[mp[j]][i]=-1
                row=dict(left=[di,list(u)],right=[dj,list(v0)],edges=[[i,j]for i,j,_ in chosen])
                ob=R.edge_bad_native(g,black)
                if ob:
                    row['all_labels_actual_edge_rejection']=ob;counts['multiple_monotone_z']+=1
                else:
                    row['label_cases']=[];counts['multiple_needing_labels']+=1
                    for labels in product(*(e[2]for e in chosen)):
                        h=[r.copy()for r in base]
                        for (i,j,_),m in zip(chosen,labels):h[i][mp[j]]=h[mp[j]][i]=-(m-2)
                        z=R.edge_bad_native(h,black)
                        if z:case=dict(labels=list(labels),kind='actual_edge_not_pd',**z);counts['multiple_labeled_edge']+=1
                        else:
                            a,sc=R.native_to_rational(h);sig=R.inertia(a);assert sig[1]==2
                            case=dict(labels=list(labels),kind='two_negative_directions',gram=R.encode(a),gauge_squared=sc,signature=sig)
                            counts['multiple_labeled_index2']+=1
                        row['label_cases'].append(case)
                multiple.append(row)
    assert len(rows)==12135 and counts['white_rejected']==3037
    assert counts['joins_allowing4']==178 and counts['single_label_bad_corner']==20,dict(counts)
    assert len(multiple)==10990 and counts['multiple_needing_labels']==3
    assert counts['multiple_labeled_edge']==3 and counts['multiple_labeled_index2']==3
    paths=[Path(__file__),Path(R.__file__),ROOT/'scripts/fresh_selective_distance10_dimension116_20261017_review.py']
    paths +=[ROOT/f'{folder}/{stem}.{ext}'for stem in(CENSUS,'fresh_distance8_low_core_20261017')for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    out=dict(status='PASS_DISTANCE8_ALL_LABEL_HIGH_DISJOINT_EXCLUSION',
        sources={str(p.relative_to(ROOT)):R.sha(p)for p in paths},marked_contexts=len(contexts),ports=ports,
        allowed_ports=len(allowed),single_joins=rows,single_types=list(types.values()),multiple_edge_subsets=multiple,
        counts=dict(counts),context_pair_multiplicities=dict(Counter(map(len,good.values()))),
        conclusion='Under actual quadrilateral white/edge/corner positivity and no opposite-pair white path of length<=8, the two minimal indefinite cores cannot be disjoint when both have orders6 through10.',
        scope='Complete high-order disjoint-core local lemma, including every cross-edge subset and every Coxeter label. No overlap-extension or ambient dimension conclusion is asserted.')
    return json.loads(json.dumps(out))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],flush=True)
