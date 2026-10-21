#!/usr/bin/env python3
"""Necessary cusp-edge coverage scout for the changed-mirror construction."""
import itertools
import json
from collections import Counter
import numpy as np
import fresh_h20_quartic_orbit_completion_20261014 as P

def run():
    P.mp.mp.dps=90
    z=json.loads((P.ROOT/P.RAW).read_text())['coxeter_input']
    Q=-np.array(z['basis_gram'],dtype=np.int64)
    W=np.array(z['wall_coordinates_as_columns_transposed'],dtype=np.int64)
    d=np.einsum('ij,ij->i',W@Q,W)
    local=json.loads((P.ROOT/P.LOCAL).read_text())
    walk=json.loads((P.ROOT/P.WALK).read_text())['corner_orbits']
    source=f'results/{P.TAG}.json'
    scout=json.loads((P.ROOT/source).read_text())
    bankfile=f'data/enumerations/{P.TAG}_bank.npz'
    bank=np.load(P.ROOT/bankfile)['roots']
    dn=np.einsum('ij,ij->i',bank@Q,bank)
    cache={};rows=[];counts=Counter()
    for spec in scout['sections']:
        index=spec['local_record'];row=local['directions'][index]
        if index not in cache:
            ab=list(map(P.mpr,row['root_interval']))
            poly=list(map(P.mpr,row['minimal_polynomial']))
            alpha=P.mp.findroot(lambda x:P.mp.polyval(poly,x),sum(ab)/2,tol=P.mp.mpf('1e-75'))
            ev=lambda cs:P.mp.polyval(list(map(P.mpr,cs)),alpha)
            e0=[ev(p)for p in row['normal_constant']]
            D=list(map(P.mpr,row['null_shift']))
            L=ev(row['parameter_lower_magnitude']);q=float(ev(row['normal_square']))
            cache[index]=(e0,D,L,q)
        e0,D,L,q=cache[index]
        e=np.array([float(a-P.mpr(spec['parameter_u'])*L*b)for a,b in zip(e0,D)])
        seed=spec['seed_walls'];ids=spec['compatible_cap_indices']+spec['near_angle_cap_indices']
        delta=np.array(list(map(float,walk[spec['cusp']]['primitive_ray'])))
        B=bank[ids].copy()
        B[(B@Q@delta)>0]*=-1
        se=W[seed]@Q@e;be=B@Q@e
        H=W[seed]@Q@W[seed].T-se[:,None]*se[None,:]/q
        b=B@Q@W[seed].T-be[:,None]*se[None,:]/q
        n=dn[ids]-be*be/q
        factors=[];costs=[]
        for factor in walk[spec['cusp']]['affine_components']:
            walls=[i for i in factor['walls'] if i in seed]
            ff=[seed.index(i)for i in walls]
            factors.append(walls);cc=[]
            for omit in ff:
                subset=[i for i in ff if i!=omit]
                A=H[np.ix_(subset,subset)]
                assert np.min(np.linalg.eigvalsh(A))>1e-8
                v=b[:,subset]
                cc.append(np.einsum('ij,ij->i',v@np.linalg.inv(A),v))
            costs.append(np.array(cc).T)
        assert sum(len(f)-1 for f in factors)==19
        edges=[];uncovered=[]
        for omissions in itertools.product(*[range(len(f))for f in factors]):
            remaining=n.copy()
            for k,o in enumerate(omissions):
                remaining-=costs[k][:,o]
            winners=np.flatnonzero(remaining>=-1e-7)
            record={'omitted_seed_walls':[f[o]for f,o in zip(factors,omissions)],
                    'candidate_cap_indices':[ids[i]for i in winners]}
            if len(winners)==0:
                uncovered.append(record['omitted_seed_walls'])
            else:
                edges.append(record)
        counts['tested_sections']+=1
        counts['seed_edges']+=len(edges)+len(uncovered)
        counts['covered_seed_edges']+=len(edges)
        if not uncovered:
            counts['sections_passing_necessary_edge_coverage']+=1
        rows.append({'local_record':index,'parameter_u':spec['parameter_u'],
                     'seed_walls':seed,'bank_caps_including_near_angle_remainders':len(ids),
                     'covered_edges':edges,'uncovered_edge_omissions':uncovered})
    sources=[P.RAW,P.LOCAL,P.WALK,source,bankfile,
             f'scripts/{P.TAG}.py',f'scripts/{P.TAG}_edges.py']
    out={'status':'COMPLETED_NUMERICAL_CUSP_EDGE_COVERAGE_SCOUT',
         'sources':{p:P.sha(p)for p in sources},'counts':dict(counts),'sections':rows,
         'scope':'Necessary Gram-Schur cusp-edge coverage using the finite one-reflection bank, including near-angle numerical remainders. Floating-point only; no exact nonexistence or finite-volume certification.'}
    (P.ROOT/'results'/f'{P.TAG}_edges.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print(out['counts'])
    print('best',sorted((len(r['uncovered_edge_omissions']),r['local_record'],r['parameter_u'])for r in rows)[:12])

if __name__=='__main__':run()
