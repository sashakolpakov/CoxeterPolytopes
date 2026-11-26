"""Exact small-coordinate parametrization of the balanced problem's dual.

Let E_j be the trace coefficients of source F_j, q=Tz solve E_j q=0,
and tau>0 be the checked common unit. For G(u)=tau F_0+sum F_j u_j/c_j,
where F_j has the unique scalar pivot c_j, y=(1,svec(G(u))) satisfies
all balanced-dual linear identities. Its cone membership remains to be checked.
"""
from fractions import Fraction as F
import gzip
import json


def build(problem,original_layout,embedding,face_layout,steps,tau):
    import numpy as np
    names=problem['variables']
    selected={int(s['source_row']):int(s['pivot']) for s in steps}
    if len(selected)!=len(steps) or len(set(selected.values()))!=len(steps):
        raise ValueError('Repeated primal coordinate or certificate pivot')
    columns=sorted(selected)
    if any(col>=len(names) for col in columns):
        raise ValueError('The homogeneous parametrization eliminated the normalization')
    units=[None]*len(columns)
    position={names[col]:i for i,col in enumerate(columns)}
    pivots=[selected[col] for col in columns]
    lookup={int(embedding[p]):i for i,p in enumerate(pivots)}
    size=original_layout['scalar_cones']+sum(n*(n+1)//2 for n in original_layout['PSD_sizes'])
    inverse=np.full(size,-1,dtype=np.int64)
    inverse[embedding]=np.arange(len(embedding))
    forms=[None]*len(embedding)
    descriptions=[]
    scalar=0
    base=original_layout['scalar_cones']
    for block in problem['blocks']:
        M=block['entries'];n=len(M)
        if n==1:
            locations=[(scalar,0,0)]
            scalar+=1
        else:
            locations=[]
            for i in range(n):
                for j in range(i,n):
                    locations.append((base,i,j));base+=1
        kept=[]
        for old,i,j in locations:
            col=int(inverse[old])
            if col<0:continue
            forms[col]=M[i][j]
            kept.append(col)
            if old in lookup:
                k=lookup[old]
                name=names[columns[k]]
                if i!=j or set(M[i][j])!={name} or not M[i][j][name]:
                    raise ValueError('A reconstruction pivot is not a unique pure scalar diagonal')
                units[k]=F(M[i][j][name])
        if kept:
            lo,hi=min(kept),max(kept)
            # Face compression moves full-block singletons into scalar order.
            order=1 if lo<face_layout['scalar_cones'] else next(
                k for k in range(2,n+1) if k*(k+1)//2==len(kept))
            if hi-lo+1!=len(kept):raise ValueError('A projected cone is not contiguous')
            descriptions.append(dict(name=block['name'],start=lo,size=order))
    if any(q is None for q in forms) or any(c is None for c in units):
        raise ValueError('The exact dual parametrization is incomplete')
    # Source references are no longer needed after each expression is replaced.
    for col,q in enumerate(forms):
        result={}
        for name,value in q.items():
            if not name:
                result[-1]=tau*F(value)
            else:
                if name not in position:
                    raise ValueError('A discarded ordinary coordinate remains on the proved dual face')
                k=position[name]
                result[k]=F(value)/units[k]
        forms[col]=result
    for k,pivot in enumerate(pivots):
        if forms[pivot]!={k:F(1)}:
            raise ValueError('The dual coordinate is not exactly its chosen cone diagonal')
    descriptions.sort(key=lambda item:item['start'])
    expected=0
    for block in descriptions:
        if block['start']!=expected:raise ValueError('A dual cone was duplicated or omitted')
        expected+=block['size']*(block['size']+1)//2
    if expected!=len(forms):raise ValueError('Incomplete dual cone partition')
    meta=dict(schema='balanced-dual-coordinate-model-v1',
        variables=[names[col] for col in columns],source_variable_indices=columns,
        certificate_scalar_pivots=pivots,source_pivot_coefficients=list(map(str,units)),
        positive_common_unit=str(tau),scalar_cones=face_layout['scalar_cones'],
        PSD_sizes=face_layout['PSD_sizes'],forms=len(forms),
        matrix_coordinates='Actual entries over Q; svec only in the numerical export',
        certificate_map='y_equal=1; y_cone=svec(G(u)); G(u)=tau F0+sum Fj*u_j/c_j',
        original_primal_coordinates='x_j=u_j/(tau*c_j); only the compressed face model is asserted')
    return forms,descriptions,meta


def save_forms(path,forms):
    with gzip.open(path,'xt') as stream:
        for row in forms:
            stream.write(json.dumps([[k,str(v)] for k,v in row.items()],separators=(',',':'))+'\n')


def read_forms(path,number=None):
    from flint import fmpq
    with gzip.open(path,'rt') as stream:
        forms=[{int(k):fmpq(v) for k,v in json.loads(line)} for line in stream]
    if number is not None and len(forms)!=number:raise ValueError('Missing dual matrix entries')
    return forms


def encode(forms,meta):
    import numpy as np
    import scipy.sparse as sp
    factors=np.ones(len(forms));base=meta['scalar_cones']
    for n in meta['PSD_sizes']:
        i,j=np.triu_indices(n)
        factors[base:base+len(i)][i!=j]=np.sqrt(2.)
        base+=len(i)
    if base!=len(forms):raise ValueError('Wrong dual cone dimension')
    rr=[];cc=[];vv=[];rhs=np.zeros(len(forms))
    for r,row in enumerate(forms):
        rhs[r]=factors[r]*float(row.get(-1,0))
        for col,value in row.items():
            if col>=0:
                rr.append(r);cc.append(col);vv.append(-factors[r]*float(value))
    return dict(A=sp.coo_matrix((vv,(rr,cc)),shape=(len(forms),len(meta['variables']))).tocsc(),
                b=rhs,c=np.zeros(len(meta['variables']))),dict(zero_cones=0,
        scalar_cones=meta['scalar_cones'],PSD_sizes=meta['PSD_sizes'],
        mode='normalized_dual_certificate_for_balanced_certificate_feasibility')
