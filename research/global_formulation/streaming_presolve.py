"""Exact SDP preparation without materializing all quadratic cone tests.

Only equality elimination, PSD zero diagonals, opposite scalar coordinate
rows, constant cones and duplicate positive cone scalings are used. This
is equivalent presolve, not a relaxation. Remaining singular cones are
allowed in the strictly feasible common-margin auxiliary SDP.
"""
from fractions import Fraction as F
import time

from . import conic,feasibility


def prepare(original,progress,seconds=300):
    started=time.monotonic()
    rows,proofs=feasibility.original_binomial_faces(original,progress)
    current,substitutions=conic.eliminate(original,rows)
    progress("streaming_equalities_eliminated",dict(free_variables=len(current["variables"])))
    maps=[dict(original_block=b,indices=list(range(len(block["entries"]))),scale=F(1))
        for b,block in enumerate(original["blocks"])]
    trace=[dict(stage="original_scalar_faces",proofs=proofs)]
    for step in range(8):
        derived,proofs,units={ },[],{}
        for b,block in enumerate(current["blocks"]):
            if len(block["entries"])!=1: continue
            q=block["entries"][0][0]
            if len(q)==1 and "" not in q:
                v,c=next(iter(q.items()))
                units.setdefault(v,{})[1 if c>0 else -1]=(b,c)
        for v,signs in units.items():
            if len(signs)!=2: continue
            b,c=signs[1];d,a=signs[-1]
            if conic.plus(current["blocks"][b]["entries"][0][0],
                    conic.times(-c/a,current["blocks"][d]["entries"][0][0])):
                raise ValueError("Invalid opposite coordinate identity")
            row=conic.form(**{v:1});derived[conic.canonical(row)]=row
            proofs.append(dict(kind="opposite_scalar_coordinates",blocks=[b,d],weights=[F(1),-c/a],variable=v))
        blocks,new_maps,seen=[],[],set()
        removed_zero_rows=removed_blocks=0
        for b,(block,mapping) in enumerate(zip(current["blocks"],maps,strict=True)):
            M=block["entries"];zero=[]
            for i in range(len(M)):
                if M[i][i]: continue
                zero.append(i)
                for q in M[i]:
                    if q: derived[conic.canonical(q)]=q
                if any(M[i]): proofs.append(dict(kind="PSD_zero_diagonal",block=b,index=i,forced_rows=M[i]))
            keep=[i for i in range(len(M)) if i not in zero]
            removed_zero_rows+=len(zero)
            if not keep: removed_blocks+=1;continue
            A=[[M[i][j] for j in keep] for i in keep]
            if not any(v for row in A for q in row for v in q if v):
                conic.psd_pivots([[q.get("",F(0)) for q in row] for row in A])
                removed_blocks+=1;continue
            scale=max(abs(c) for row in A for q in row for c in q.values())
            A=[[conic.times(1/scale,q) for q in row] for row in A]
            key=conic.canonical(A)
            if key in seen: removed_blocks+=1;continue
            seen.add(key);blocks.append(dict(name=block["name"],entries=A))
            new_maps.append(dict(original_block=mapping["original_block"],
                indices=[mapping["indices"][i] for i in keep],scale=mapping["scale"]/scale))
        item=dict(step=step,input_variables=len(current["variables"]),input_blocks=len(current["blocks"]),
            remaining_blocks=len(blocks),derived_equalities=len(derived),removed_zero_rows=removed_zero_rows,
            removed_blocks=removed_blocks,proofs=proofs)
        trace.append(item);maps=new_maps
        if derived:
            current,update=conic.eliminate(conic.model(current["variables"],blocks,list(derived.values())))
            substitutions={v:feasibility.substitute(q,update) for v,q in substitutions.items()}
        else: current=conic.model(current["variables"],blocks)
        progress("streaming_exact_presolve",{k:v for k,v in item.items() if k!="proofs"})
        if not derived or time.monotonic()-started>=seconds: break
    return current,substitutions,maps,trace


def self_check():
    P=conic.model(["x","y","z"],[dict(name="PSD_zero",entries=[
        [conic.form(x=1),conic.form(y=1)],[conic.form(y=1),{}]]),
        dict(name="scalar_positive",entries=[[conic.form(z=1)]]),
        dict(name="scalar_negative",entries=[[conic.form(z=-1)]])],[conic.form(-1,x=1)])
    Q,substitutions,maps,trace=prepare(P,lambda *_:None)
    if Q["variables"]: raise ValueError("The zero-row fixture retained an unconstrained coordinate")
    point={v:conic.value(q,{}) for v,q in substitutions.items()}
    if point!={"x":F(1),"y":F(0),"z":F(0)}: raise ValueError("Wrong exact streaming reconstruction")
    conic.check_point(P,point)
    return dict(exact_PSD_zero_row_and_opposite_coordinate_fixture=True,
        unique_point_verified_in_original=True,all_pairs_quadratic_tests_materialized=False)
