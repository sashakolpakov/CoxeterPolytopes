#!/usr/bin/env python3
"""Bounded numerical event search around the second S6 exterior component."""
import argparse
import json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
from pathlib import Path

import numpy as np
import mpmath as mp
import fresh_nonarithmetic_a2_escape_partition_20261017 as P

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_s6_component_b_bridge_scout_20261017'
CONE='fresh_nonarithmetic_pi8_quad_cone_20261017'
CAPS=['fresh_nonarithmetic_pi8_tetrahedron_completion_20261017','fresh_nonarithmetic_pi8_quad_completion_20261017']
LABELS=list(range(2,17))+[20,24,32,64]
TRIPLES=[(8,44,52),(30,44,52),(30,47,52),(34,47,52)]


def read(tag):
    d=json.loads((ROOT/f'results/{tag}.json').read_text())
    for p,h in d['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return d


def mid(x):
    if hasattr(x,'l'):
        return float((x.l+x.h)/2)
    return float((Fraction(x[0])+Fraction(x[1]))/2)


def high(x):
    value=(x.l+x.h)/2 if hasattr(x,'l') else (Fraction(x[0])+Fraction(x[1]))/2
    return mp.mpf(value.numerator)/value.denominator


def run():
    _,q,old,_,_,_,_,_=P.physical_data()
    mp.mp.dps=60
    q=mp.matrix([[high(x)for x in row]for row in q])
    physical=[[high(x)for x in row]for row in old]
    for tag in CAPS:
        physical.extend([[high(x)for x in r['physical']]for r in read(tag)['new_normals']])
    assert len(physical)==60
    evals,evecs=mp.eigsy(q)
    assert sum(v<0 for v in evals)==1
    metric=np.array([float(mp.sign(v))for v in evals])
    transform=evecs*mp.diag([mp.sqrt(abs(v))for v in evals])
    temporary=mp.matrix(physical)*transform
    h=[sum(temporary[i,j]for i in range(60))for j in range(10)]
    square=sum(metric[i]*h[i]*h[i]for i in range(10));assert square<0
    h=[x/mp.sqrt(-square)*mp.sign(h[0])for x in h]
    boost=mp.eye(10);boost[0,0]=h[0]
    for i in range(1,10):
        boost[0,i]=boost[i,0]=-h[i]
        for j in range(1,10):boost[i,j]+=h[i]*h[j]/(h[0]+1)
    transform=transform*boost
    normals=np.array((mp.matrix(physical)*transform).tolist(),dtype=float)
    assert np.max(abs(np.sum(normals*normals*metric,axis=1)-2))<1e-7
    saved=read(CONE)
    faceids=[i for i,r in enumerate(saved['final_rays'])if set(P.S)<=set(r['active'])]
    rays=np.zeros((len(saved['final_rays']),10))
    for i in faceids:
        v=mp.matrix([[high(x)for x in saved['final_rays'][i]['normalized_coordinate_bounds']]])*transform
        rays[i]=[float(x)for x in v]
    finite=next(i for i in faceids if saved['final_rays'][i]['kind']=='finite_vertices')
    witness=rays[finite]
    angle_grid=np.array([0.]+[2*np.cos(np.pi/m)for m in range(3,257)]+[2.])
    a2=[(2,2)]+[(m,2)for m in(3,4,5)]+[(2,m)for m in(3,4,5)]
    assignments=np.array([a+b for a in a2 for b in product(LABELS,repeat=3)],dtype=int)
    magnitudes=2*np.cos(np.pi/assignments)
    magnitudes[assignments==2]=0
    groups=[]
    for sides in TRIPLES:
        frame=P.S+list(sides)
        graph=normals[frame]@np.diag(metric)@normals[frame].T
        inverse=np.linalg.inv(graph)
        contacts=np.zeros((len(assignments),9));contacts[:,4:]=-magnitudes
        keep=np.ones(len(contacts),dtype=bool)
        # The first side is shared by the two actual finite edge supports.
        edge_costs=[]
        for endpoint in(7,8):
            indices=list(range(7))+[endpoint]
            gi=np.linalg.inv(graph[np.ix_(indices,indices)])
            assert np.linalg.eigvalsh(graph[np.ix_(indices,indices)])[0]>0
            costs=np.einsum('bi,ij,bj->b',contacts[:,indices],gi,contacts[:,indices])
            keep &= costs<2-1e-8
            edge_costs.append(costs)
        numbers=np.where(keep)[0]
        c=contacts[numbers]
        coeff=c@inverse
        particular=coeff@normals[frame]
        residual=2-np.einsum('bi,bi->b',coeff,c)
        assert np.all(residual>0)
        rayid=next(i for i in faceids if set(saved['final_rays'][i]['active'])==set(frame))
        kernel=rays[rayid]
        kernel/=np.sqrt(np.dot(kernel*metric,kernel))
        bank=[];rejected=Counter();reject_witnesses=[]
        for sign in(-1,1):
            candidates=particular+sign*np.sqrt(residual)[:,None]*kernel
            values=candidates@np.diag(metric)@normals.T
            for pos,(normal,products) in enumerate(zip(candidates,values)):
                num=int(numbers[pos]);target=contacts[num]
                residual_gate=float(np.max(abs(products[frame]-target)))
                norm_error=float(abs(np.dot(normal*metric,normal)-2))
                if residual_gate>2e-5 or norm_error>2e-5:
                    rejected['numerical_residual_unresolved']+=1;continue
                if np.max(products)>2e-6:
                    rejected['obtuse']+=1
                    reject_witnesses.append(dict(labels=assignments[num].tolist(),kernel_sign=sign,
                        reason='obtuse',wall=int(np.argmax(products)),pair=float(np.max(products))))
                    continue
                if np.dot(normal*metric,witness)>=-1e-8:
                    rejected['no_common_persistent_vertex']+=1;continue
                other=[i for i in range(60)if i not in frame]
                errors=[float(np.min(abs(angle_grid+products[i])))if products[i]>-2+1e-6 else 0 for i in other]
                if max(errors)>2e-6:
                    rejected['non_Coxeter_old_contact']+=1
                    j=other[int(np.argmax(errors))]
                    reject_witnesses.append(dict(labels=assignments[num].tolist(),kernel_sign=sign,
                        reason='non_Coxeter_old_contact',wall=j,pair=float(products[j]),error=max(errors)))
                    continue
                record=dict(labels=assignments[num].tolist(),kernel_sign=sign,
                    normalized_lorentz_normal=normal.tolist(),old60_products=products.tolist(),
                    maximum_unprescribed_angle_error=max(errors),norm_error=norm_error,
                    prescribed_residual_error=residual_gate,
                    endpoint_projection_costs=[float(v[num])for v in edge_costs])
                bank.append(record)
        group=dict(sides=list(sides),full_support=frame,assignments=len(assignments),
            two_endpoint_cost_survivors=len(numbers),support_condition_number=float(np.linalg.cond(graph)),
            rejections=dict(rejected),rejection_witnesses=reject_witnesses,survivors=bank)
        groups.append(group)
        print('COMPONENT_B_EVENTS',sides,len(numbers),len(bank),dict(rejected),flush=True)
    pairs=[];paircounts=Counter()
    for gi,gj,shared in((0,1,(44,52)),(2,3,(47,52))):
        for ai,a in enumerate(groups[gi]['survivors']):
            va=np.array(a['normalized_lorentz_normal'])
            for bi,b in enumerate(groups[gj]['survivors']):
                vb=np.array(b['normalized_lorentz_normal']);ab=float(np.dot(va*metric,vb))
                paircounts['pairs_tested']+=1
                if ab>2e-6 or ab<-2+1e-6:
                    paircounts['no_finite_mutual_contact']+=1;continue
                target=int(np.argmin(abs(angle_grid+ab)));error=float(abs(angle_grid[target]+ab))
                if error>2e-6:
                    paircounts['mutual_finite_angle_gap']+=1;continue
                corner_eigenvalues=[]
                for wall in shared:
                    nn=np.vstack([normals[P.S+[wall]],va,vb]);gg=nn@np.diag(metric)@nn.T
                    corner_eigenvalues.append(float(np.linalg.eigvalsh(gg)[0]))
                if min(corner_eigenvalues)<=1e-8:
                    paircounts['bridge_corner_not_PD']+=1;continue
                label=2 if target==0 else target+2
                pairs.append(dict(groups=[gi,gj],indices=[ai,bi],shared_sides=list(shared),
                    mutual_label=label,mutual_product=ab,angle_error=error,
                    finite_bridge_corner_minimum_eigenvalues=corner_eigenvalues))
    paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',f'results/{CONE}.json']+[f'results/{t}.json'for t in CAPS]
    return dict(status='COMPLETED_BOUNDED_NUMERICAL_COMPONENT_B_BRIDGE_SCOUT',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        prescribed_labels=LABELS,unprescribed_angle_recognition_through=256,
        groups=groups,pair_counts=dict(paircounts),bridge_hits=pairs,
        metric=metric.tolist(),physical_to_normalized_lorentz_transform=np.array(transform.tolist(),dtype=float).tolist(),
        coordinate_conditioning='60-digit midpoint eigendecomposition and Lorentz boost along the timelike sum of normals, followed by double-precision event discovery.',
        scope='Bounded floating-point event discovery only. Isolated high-label tails beyond the displayed grid and unrecognized exact finite coincidences remain open. No numerical survivor is an accepted Coxeter wall; exact recipes, all gates and whole-face/cone tests are required. No all-label nonexistence conclusion follows from an empty hit list.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['pair_counts'],'hits',len(result['bridge_hits']),flush=True)
