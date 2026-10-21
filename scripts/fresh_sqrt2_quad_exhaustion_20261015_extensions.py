#!/usr/bin/env python3
"""All one-extra-white profiles on the surviving finite-label union cores.

Every vertex of an actual external white support must pass these tests, but
this alone does not classify how several external whites can meet each other.
"""
import argparse,importlib.util,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE='fresh_sqrt2_quad_exhaustion_20261015';TAG=BASE+'_extensions'

def load(suffix):
    path=ROOT/'scripts'/f'{BASE}{suffix}.py';spec=importlib.util.spec_from_file_location('quad'+suffix,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

E=load('');A=load('_cores')

def profiles(row):
    n=row['nodes'];W=row['white']+[n]
    g=dict(zip(combinations(range(n),2),row['labels']));contacts=[2]*n;out=[];nodes=0
    def rec(k):
        nonlocal nodes
        nodes+=1
        if k==n:
            if any(v!=2 for v in contacts):out.append(list(contacts))
            return
        p=(k,n)
        g[p]=2;contacts[k]=2;rec(k+1)
        for label in (3,4,6,8):
            g[p]=label;contacts[k]=label
            if E.no_short(W,g) and E.corners_finite(W,g):rec(k+1)
        del g[p];contacts[k]=2
    rec(0)
    return out,nodes

def extension(row,contacts):
    n=row['nodes'];g=dict(zip(combinations(range(n),2),row['labels']))
    for i,v in enumerate(contacts):g[(i,n)]=v
    return dict(nodes=n+1,white=row['white']+[n],labels=[g[p] for p in combinations(range(n+1),2)])

def run():
    base=json.loads((ROOT/'results'/f'{BASE}.json').read_text())
    core=json.loads((ROOT/'results'/f'{BASE}_cores.json').read_text());rows=[];dotted=[];counts=Counter()
    for rec in core['records']:
        if rec['status']=='DOTTED_PARAMETER_UNFILTERED':
            row=base['labelled_union_supports'][rec['union_index']];domain,nodes=profiles(row);n=row['nodes']
            g=dict(zip(combinations(range(n),2),row['labels']));trials=[]
            for contacts in domain:
                touched=[i for i,v in enumerate(contacts) if v!=2]
                trial=dict(contacts=contacts)
                if set(touched)<=set((0,2)):
                    trial['status']='ONLY_DOTTED_ENDPOINT_CONTACTS'
                elif touched in ([1],[3]):
                    b=touched[0];d=4-b
                    old=next(w for w in row['white'] if g[E.pair(b,w)]==3)
                    support=[b,d,old,n]
                    ext=extension(row,contacts);gg=dict(zip(combinations(range(n+1),2),ext['labels']))
                    assert [gg[E.pair(b,v)] for v in (d,old,n)]==[4,3,3]
                    assert all(gg[E.pair(u,v)]==2 for u,v in combinations((d,old,n),2))
                    trial.update(status='REJECT_AFFINE_433_STAR',support=support,
                                 unit_gram_kernel=['1','sqrt2/2','1/2','1/2'])
                else:
                    assert touched==[1,3] and contacts[1]==contacts[3]==3
                    ext=extension(row,contacts);gg=dict(zip(combinations(range(n+1),2),ext['labels']))
                    assert all(gg[E.pair(u,v)]==2 for u in (0,2) for v in (1,3,n))
                    trial.update(status='REJECT_ORTHOGONAL_INDEFINITE_PIECES',
                                 dotted_pair=[0,2],triangle=[1,3,n],triangle_unit_determinant='-sqrt2/4')
                trials.append(trial)
            dotted.append(dict(class_index=rec['class_index'],union_index=rec['union_index'],
                               structural_nodes=nodes,structural_profiles=len(domain),trials=trials))
            continue
        if rec['status']!='PASS_PRINCIPAL_UNION_NECESSARY_TESTS':continue
        row=base['labelled_union_supports'][rec['union_index']];domain,nodes=profiles(row)
        C=A.matrix(row,True)
        # A positive-semidefinite principal core forces every new row to
        # annihilate its nullspace. Check this before the larger determinants.
        kernel=C.nullspace().to_list()
        accepted=[];trials=[]
        for contacts in domain:
            ext=extension(row,contacts);trial=dict(contacts=contacts)
            counts['finite_corner_profiles']+=1
            cyc=A.cyclic_field(ext)
            if not cyc['passes']:trial['status']='REJECT_CYCLE_FIELD'
            else:
                M=A.matrix(ext,True);v=[M.to_list()[i][-1] for i in range(row['nodes'])]
                products=[sum((a*b for a,b in zip(k,v)),A.K.zero) for k in kernel]
                trial['conjugate_kernel_products']=list(map(A.enc,products))
                if any(products):trial['status']='REJECT_CONJUGATE_KERNEL'
                else:
                    check=A.exact_conditions(ext);trial.update(check)
                    if not check['physical']['passes']:trial['status']='REJECT_TWO_NEGATIVE_DIRECTIONS'
                    elif not check['conjugate']['passes']:trial['status']='REJECT_CONJUGATE_NOT_PSD'
                    else:trial['status']='PASS_ONE_EXTRA_WHITE';accepted.append(contacts)
            counts[trial['status']]+=1;trials.append(trial)
        rows.append(dict(class_index=rec['class_index'],union_index=rec['union_index'],
                         structural_nodes=nodes,structural_profiles=len(domain),
                         conjugate_kernel=[list(map(A.enc,r)) for r in kernel],accepted=accepted,trials=trials))
        print(rec['class_index'],'domain',len(domain),'accepted',accepted,flush=True)
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',f'scripts/{BASE}_cores.py',f'results/{BASE}_cores.json']
    assert sum(r['structural_profiles'] for r in dotted)==18
    return dict(status='PASS_COMPLETE_ONE_EXTRA_WHITE_PROFILE_FILTER',counts=dict(counts),records=rows,dotted_records=dotted,
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                scope='All single extra-white contacts to these nine finite-label cores. Multiple extra whites, components unattached directly to the core, and interactions among them are not exhausted.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'])
