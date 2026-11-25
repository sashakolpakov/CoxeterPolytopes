"""Exact changes of units for an affine SDP; no cone is relaxed."""
from fractions import Fraction as F
import math

from . import conic


def power10(exponent):
    return F(10**exponent) if exponent>=0 else F(1,10**(-exponent))


def unit(value,half=False):
    if not value:return F(1)
    exponent=-math.log10(abs(float(value)))/(2 if half else 1)
    return power10(max(-24,min(24,round(exponent))))


def equilibrate(problem,reference):
    if problem["equations"]:raise ValueError("Eliminate equations before equilibration")
    columns={v:1/unit(max(abs(reference.get(v,F(0))),F(1,10**12)))
        for v in problem["variables"]}
    change={v:conic.form(**{v:c}) for v,c in columns.items()}
    blocks=[];diagonals=[];scales=[]
    for block in problem["blocks"]:
        entries=[[{v:c*columns.get(v,F(1)) for v,c in q.items()} for q in row]
            for row in block["entries"]]
        r=len(entries)
        diagonal=[unit(max((abs(c) for c in entries[i][i].values()),default=F(0)),half=True)
            for i in range(r)] if r>1 else [F(1)]
        magnitude=max((abs(c*diagonal[i]*diagonal[j]) for i,row in enumerate(entries)
            for j,q in enumerate(row) for c in q.values()),default=F(0))
        scale=unit(magnitude)
        scaled=[[conic.times(scale*diagonal[i]*diagonal[j],q) for j,q in enumerate(row)]
            for i,row in enumerate(entries)]
        blocks.append(dict(block,entries=scaled));diagonals.append(diagonal);scales.append(scale)
    result=conic.model(problem["variables"],blocks,metadata=dict(
        original_model_sha256=conic.digest_bytes(conic.canonical(problem)),
        positive_invertible_column_units=True,positive_diagonal_congruences=True))
    # Certify the maps for every affine entry, including its constant term.
    for old,new,d,s in zip(problem["blocks"],result["blocks"],diagonals,scales,strict=True):
        if s<=0 or any(v<=0 for v in d):raise ValueError("Nonpositive cone unit")
        for i,row in enumerate(old["entries"]):
            for j,q in enumerate(row):
                inverse={v:c/(s*d[i]*d[j]*columns.get(v,F(1)))
                    for v,c in new["entries"][i][j].items()}
                if inverse!=q:raise ValueError("A cone equilibration inverse identity failed")
    report=dict(column_scales=columns,block_scales=scales,diagonal_congruences=diagonals,
        all_affine_inverse_identities_checked=True,original_feasible_set_preserved=True,
        reference_only_selects_units_and_adds_no_constraints=True)
    return result,change,scales,diagonals,report


def self_check():
    P=conic.model(["x","y"],[dict(name="mixed_units",entries=[
        [conic.form(x=10**8),conic.form(1)],
        [conic.form(1),conic.form(y=F(1,10**8))]])])
    old={"x":F(2,10**8),"y":F(2*10**8)}
    Q,change,scales,diagonals,report=equilibrate(P,old)
    new={v:old[v]/change[v][v] for v in old}
    conic.check_point(P,old);conic.check_point(Q,new)
    if {v:conic.value(q,new) for v,q in change.items()}!=old:raise ValueError("Wrong equilibrated lift")
    return dict(all_affine_inverse_identities=True,positive_congruences=True,exact_point_preserved=True)
