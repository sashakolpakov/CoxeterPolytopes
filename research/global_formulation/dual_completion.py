"""Exact completion of numerical cone weights to the original Farkas dual.

Numerical eigenvectors only propose rational Gram factors. Nonnegative
coordinate rows repair the adjoint; a triangular equality history lifts
its potentials. The complete original dual is then checked over Q.
"""
from fractions import Fraction as F

from . import conic


def adjoint(problem,matrices):
    total={}
    for block,M in zip(problem["blocks"],matrices,strict=True):
        for i,row in enumerate(block["entries"]):
            for j in range(i,len(row)):
                weight=M[i][j]*(1 if i==j else 2)
                if not weight:continue
                for v,c in row[j].items():total[v]=total.get(v,F(0))+c*weight
    return {v:c for v,c in total.items() if c}


def equality_potentials(original,q):
    order={v:i for i,v in enumerate(original["variables"])}
    basis={};history=[]
    for index,source in enumerate(original["equations"]):
        row=dict(source);deps=[]
        for pivot in sorted(basis,key=order.get):
            if pivot in row:
                k=basis[pivot];c=row[pivot]
                row=conic.plus(row,conic.times(-c,history[k]["row"]))
                deps.append((k,c))
        variables=[v for v in row if v]
        if not variables:
            if row:raise ValueError("The original equalities are inconsistent")
            continue
        pivot=min(variables,key=order.get);scale=row[pivot]
        basis[pivot]=len(history)
        history.append(dict(row=conic.times(1/scale,row),scale=scale,source=index,deps=deps))
    weights=[F(0)]*len(history); remainder=dict(q)
    for pivot in sorted(basis,key=order.get):
        if pivot in remainder:
            k=basis[pivot];c=remainder[pivot];weights[k]-=c
            remainder=conic.plus(remainder,conic.times(-c,history[k]["row"]))
    if any(v for v in remainder):raise ValueError("The cone adjoint is not in the original equality row space")
    potentials=[F(0)]*len(original["equations"])
    for k in reversed(range(len(history))):
        item=history[k];w=weights[k]/item["scale"]
        if not w:continue
        potentials[item["source"]]+=w
        for j,c in item["deps"]:weights[j]-=w*c
    check=conic.plus(q,*(conic.times(y,row) for y,row in zip(potentials,original["equations"],strict=True) if y))
    if check!=remainder:raise ValueError("The equality-potential lift failed its exact identity")
    return potentials,remainder.get("",F(0))


def rational_grams(matrices,digits):
    import numpy as np
    grid=10**digits;answer=[];grams={}
    for b,M in enumerate(matrices):
        if len(M)==1:
            answer.append([[max(F(0),F(round(M[0][0]*grid),grid))]])
            continue
        A=np.asarray([[float(v) for v in row] for row in M]);A=(A+A.T)/2
        if not np.all(np.isfinite(A)):raise ValueError("Nonfinite proposed dual cone matrix")
        eig,U=np.linalg.eigh(A)
        root=U*np.sqrt(np.maximum(eig,0.))[None,:]
        factors=[[int(round(float(v)*grid)) for v in row] for row in root]
        keep=[j for j in range(len(M)) if any(row[j] for row in factors)]
        factors=[[row[j] for j in keep] for row in factors]
        Q=[[F(sum(a*c for a,c in zip(left,right,strict=True)),grid*grid) for right in factors] for left in factors]
        answer.append(Q);grams[b]=dict(denominator=grid,factors=factors)
    return answer,grams


def complete(original,compiled,matrices,digits):
    weights,grams=rational_grams(matrices,digits)
    pivots={}
    for b,block in enumerate(compiled["blocks"]):
        if len(block["entries"])!=1:continue
        row=block["entries"][0][0]
        if len(row)==1:
            v,c=next(iter(row.items()))
            if v and c>0:pivots.setdefault(v,(b,c))
    if set(pivots)!=set(compiled["variables"]):
        raise ValueError("Some free moment lacks a pure nonnegative coordinate row")
    for b,_ in pivots.values():weights[b]=[[F(0)]]
    residual=adjoint(compiled,weights)
    negative=[]
    for v,(b,c) in pivots.items():
        value=-residual.get(v,F(0))/c
        if value<0:negative.append((v,value))
        weights[b]=[[value]]
    if negative:return None,dict(reason="negative_completed_coordinate_weights",count=len(negative),
        examples=[dict(variable=v,value=c) for v,c in negative[:5]])
    repaired=adjoint(compiled,weights)
    if any(v for v in repaired):raise ValueError("The pure-coordinate adjoint repair was not exact")
    if repaired.get("",F(0))>=0:return None,dict(reason="nonnegative_certificate_constant")
    potentials,constant=equality_potentials(original,adjoint(original,weights))
    if constant!=repaired[""]:raise ValueError("The equality lift changed the reduced constant")
    normalization=-constant
    point={f"potential_{i}":v/normalization for i,v in enumerate(potentials)}
    for b,M in enumerate(weights):
        for i,row in enumerate(M):
            for j in range(i,len(row)):point[f"dual_b{b}_{i}_{j}"]=row[j]/normalization
    for item in grams.values():item["weight"]=1/normalization
    return (point,grams),dict(reason="exact_candidate_completed",coordinate_pivots=len(pivots),
        rational_Gram_digits=digits,positive_normalization=normalization)


def verify(original,dual,point,grams):
    if conic.canonical(dual)!=conic.canonical(conic.farkas(original)):
        raise ValueError("The retained original dual differs from the complete Farkas alternative")
    if set(point)!=set(dual["variables"]):raise ValueError("Wrong full original dual coordinate set")
    if any(conic.value(q,point) for q in dual["equations"]):raise ValueError("An original Farkas equation failed")
    for b,block in enumerate(dual["blocks"]):
        r=len(block["entries"])
        if r==1:
            if conic.value(block["entries"][0][0],point)<0:raise ValueError("A scalar dual weight is negative")
            continue
        item=grams[b];B=item["factors"];d=item["denominator"];w=item["weight"]
        if w<=0 or d<=0 or len(B)!=r:raise ValueError("Invalid rational Gram proof")
        for i,row in enumerate(block["entries"]):
            for j,entry in enumerate(row):
                expected=w*F(sum(a*c for a,c in zip(B[i],B[j],strict=True)),d*d)
                if expected!=conic.value(entry,point):raise ValueError("An original dual cone entry differs from its Gram proof")
    return dict(status="EXACT_ORIGINAL_FARKAS_CERTIFICATE_VERIFIED",equalities=len(dual["equations"]),
        cone_blocks=len(dual["blocks"]),method="exact original adjoint and rational positive Gram decompositions",
        point_sha256=conic.digest_bytes(conic.canonical(point)))
