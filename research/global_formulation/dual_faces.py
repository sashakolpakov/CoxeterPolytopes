"""Exact facial implications for the normalized certificate SDP.

Coordinate coefficient directions that are diagonal PSD expose forced dual
zero rows. Later steps are checked on the face already proved by earlier
steps. This preserves the certificate feasibility problem. It is not an
equivalent primal feasibility reduction in a weakly infeasible case.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def mask_digest(active):
    import numpy as np
    return hashlib.sha256(np.packbits(active).tobytes()).hexdigest()


def geometry(layout):
    import numpy as np
    scalar=layout['scalar_cones']
    diagonal=np.zeros(scalar+sum(n*(n+1)//2 for n in layout['PSD_sizes']),dtype=bool)
    diagonal[:scalar]=True
    blocks=[]
    start=scalar
    for n in layout['PSD_sizes']:
        i,j=np.triu_indices(n)
        diagonal[start+np.flatnonzero(i==j)]=True
        blocks.append((start,n,i,j))
        start+=len(i)
    return diagonal,blocks


def discover(data,layout,active):
    import numpy as np
    diagonal,_=geometry(layout)
    A=data['A'][active].tocsc()
    diagonal=diagonal[active]
    rays={}
    for col in range(A.shape[1]):
        start,end=A.indptr[col:col+2]
        if start==end or not np.all(diagonal[A.indices[start:end]]):
            continue
        values=A.data[start:end]
        if np.all(values<=0):
            rays[col]=1
        elif np.all(values>=0):
            rays[col]=-1
    return rays


def exact_step(problem,layout,active,rays):
    """Read all relevant exact source coefficients, not numerical signs."""
    import numpy as np
    names=problem['variables']
    if any(int(col)<0 or int(col)>=len(names) for col in rays):
        raise ValueError('Ray column outside the exact model')
    chosen={names[int(col)]:(int(col),int(sign)) for col,sign in rays.items()}
    if len(chosen)!=len(rays) or any(sign not in [-1,1] for _,sign in chosen.values()):
        raise ValueError('Invalid exact coordinate directions')
    _,blocks=geometry(layout)
    scalar=0
    full=0
    forced=set()
    counts={col:0 for col in map(int,rays)}
    for block in problem['blocks']:
        entries=block['entries']
        n=len(entries)
        if n==1:
            if active[scalar]:
                for name,c in entries[0][0].items():
                    if name not in chosen:
                        continue
                    col,sign=chosen[name]
                    value=sign*F(c)
                    if not value:
                        continue
                    if value<0:
                        raise ValueError('A source scalar coefficient is negative on the proposed face ray')
                    forced.add(scalar)
                    counts[col]+=1
            scalar+=1
            continue
        base,size,ii,jj=blocks[full]
        full+=1
        if size!=n:
            raise ValueError('Exact model and numeric block layout differ')
        position=base
        for i in range(n):
            for j in range(i,n):
                if active[position]:
                    for name,c in entries[i][j].items():
                        if name not in chosen:
                            continue
                        col,sign=chosen[name]
                        value=sign*F(c)
                        if not value:
                            continue
                        if i!=j or value<0:
                            raise ValueError('A source PSD coefficient is not diagonal PSD on the current face')
                        forced.add(position)
                        counts[col]+=1
                position+=1
    if scalar!=layout['scalar_cones'] or full!=len(blocks) or any(not n for n in counts.values()):
        raise ValueError('Incomplete or zero proposed source direction')
    after=active.copy()
    for row in forced:
        after[row]=False
    for base,n,i,j in blocks:
        diagonal_rows=base+np.flatnonzero(i==j)
        keep=after[diagonal_rows]
        after[base:base+len(i)]=keep[i]&keep[j]
    if np.array_equal(active,after):
        raise ValueError('A claimed facial step did not expose a proper face')
    return after,{'rays_0_based':{int(c):int(s) for c,s in rays.items()},
        'ray_variable_names':{int(c):names[int(c)] for c in rays},
        'positive_exact_coefficient_counts':counts,
        'positive_diagonal_coordinates_0_based':sorted(forced),
        'support_before_sha256':mask_digest(active),'support_after_sha256':mask_digest(after),
        'dual_coordinates_removed':int(np.count_nonzero(active)-np.count_nonzero(after)),
        'all_source_coefficients_checked_over_Q':True}


def reduced_certificate(data,layout,active):
    import numpy as np
    import scipy.sparse as sp
    _,blocks=geometry(layout)
    scalar=list(np.flatnonzero(active[:layout['scalar_cones']]))
    singletons=[]
    full=[]
    sizes=[]
    for base,n,i,j in blocks:
        keep=active[base+np.flatnonzero(i==j)]
        order=int(np.count_nonzero(keep))
        retained=base+np.flatnonzero(keep[i]&keep[j])
        if order==1:
            singletons.extend(retained.tolist())
        elif order>1:
            full.extend(retained.tolist())
            sizes.append(order)
    embedding=np.asarray(scalar+singletons+full,dtype=np.int64)
    if len(embedding)!=np.count_nonzero(active):
        raise ValueError('The cone embedding omits retained dual coordinates')
    projected=data['A'][embedding].tocsc()
    columns=np.flatnonzero(np.diff(projected.indptr))
    rhs=data['b'][embedding]
    equalities=len(columns)+1
    A=sp.vstack([projected[:,columns].T,sp.csc_matrix(rhs.reshape((1,-1))),
                 -sp.eye(len(embedding),format='csc')],format='csc')
    b=np.zeros(A.shape[0])
    b[equalities-1]=-1
    return {'A':A,'b':b,'c':np.zeros(len(embedding))}, {
        'zero_cones':equalities,'scalar_cones':len(scalar)+len(singletons),'PSD_sizes':sizes,
        'mode':'normalized_certificate_feasibility_on_verified_dual_face',
        'original_dual_coordinate_dimension':len(active),
        'original_ordinary_variable_dimension':data['A'].shape[1],
        'embedding':'Insert zeros at removed original dual coordinates',
        'no_claim_of_original_primal_feasibility_from_reduced_dual_infeasibility':True}, {
        'dual_coordinate_embedding':embedding,'retained_ordinary_equation_columns':columns}


def check_steps(problem,layout,steps):
    import numpy as np
    diagonal,_=geometry(layout)
    active=np.ones(len(diagonal),dtype=bool)
    for step in steps:
        rays={int(c):int(s) for c,s in step['rays_0_based'].items()}
        after,checked=exact_step(problem,layout,active,rays)
        # JSON changes integer dictionary keys to strings. Compare canonically.
        if json.dumps(checked,sort_keys=True)!=json.dumps(step,sort_keys=True):
            from . import conic
            if conic.canonical(checked)!=conic.canonical(step):
                raise ValueError('A replayed exact facial step differs from its proof')
        active=after
    return active
