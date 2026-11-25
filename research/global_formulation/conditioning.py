"""Invertible exact conditioning of narrow paired scalar cone bounds."""
from fractions import Fraction as F

from . import conic, feasibility


def thin_intervals(problem):
    if problem["equations"]:raise ValueError("Eliminate equalities before interval conditioning")
    names={b["name"]:i for i,b in enumerate(problem["blocks"])}
    pairs=[];chosen=set()
    for i,block in enumerate(problem["blocks"]):
        name=block["name"]
        marker=next((term for term in ("_inverse_response_lower","_root_inverse_lower") if term in name),None)
        if marker is None or len(block["entries"])!=1:continue
        other=name.replace(marker,marker.replace("lower","upper"))
        if other not in names:continue
        j=names[other];low=block["entries"][0][0];high=problem["blocks"][j]["entries"][0][0]
        candidates=[]
        for v in low.keys()&high.keys():
            if not v or v in chosen or low[v]*high[v]>=0:continue
            ratio=-low[v]/high[v]
            width=conic.plus(low,conic.times(ratio,high))
            delta=max((abs(c) for c in width.values()),default=F(0))
            magnitude=max(max(abs(c) for c in low.values()),max(abs(ratio*c) for c in high.values()))
            if 0<delta<F(1,10**10)*magnitude:
                candidates.append((delta/magnitude,v,ratio,delta,width))
        if not candidates:continue
        _,pivot,ratio,delta,width=min(candidates)
        chosen.add(pivot)
        pairs.append(dict(lower=i,upper=j,pivot=pivot,ratio=ratio,delta=delta,
            coordinate="conditioned_interval_"+str(len(pairs)),defining_row=low,width=width))
    if not pairs:return problem,{v:conic.form(**{v:1}) for v in problem["variables"]},[F(1)]*len(problem["blocks"]),dict(pairs=[])
    order=[p["pivot"] for p in pairs]+[v for v in problem["variables"] if v not in chosen]+[p["coordinate"] for p in pairs]
    equations=[conic.plus(p["defining_row"],conic.form(**{p["coordinate"]:-p["delta"]})) for p in pairs]
    augmented=conic.model(order,problem["blocks"],equations)
    transformed,all_subs=conic.eliminate(augmented)
    change={v:all_subs[v] for v in problem["variables"]}
    inverse={v:conic.form(**{v:1}) for v in transformed["variables"] if v in problem["variables"]}
    inverse.update({p["coordinate"]:conic.times(1/p["delta"],p["defining_row"])
        for p in pairs if p["coordinate"] in transformed["variables"]})
    # These are polynomial identities, not checks at sampled points.
    for v,q in change.items():
        if feasibility.substitute(q,inverse)!=conic.form(**{v:1}):
            raise ValueError("The interval coordinate change has no exact inverse")
    for v,q in inverse.items():
        if feasibility.substitute(q,change)!=conic.form(**{v:1}):
            raise ValueError("The inverse interval change does not compose to the identity")
    scales=[F(1)]*len(problem["blocks"])
    for p in pairs:
        scales[p["lower"]]=1/p["delta"]
        scales[p["upper"]]=p["ratio"]/p["delta"]
        if conic.plus(feasibility.substitute(p["defining_row"],change),
                conic.times(-p["delta"],all_subs[p["coordinate"]])):
            raise ValueError("A narrow-band slack identity failed")
    blocks=[dict(b,entries=[[conic.times(scales[k],q) for q in row] for row in b["entries"]])
        for k,b in enumerate(transformed["blocks"])]
    result=conic.model(transformed["variables"],blocks,metadata=dict(
        original_model_sha256=conic.digest_bytes(conic.canonical(problem)),
        exact_invertible_interval_coordinates=True,positive_cone_scalings=True))
    report=dict(pairs=pairs,old_variables=len(problem["variables"]),new_variables=len(result["variables"]),
        both_affine_inverse_identities_checked=True,original_feasible_set_preserved=True)
    return result,change,scales,report


def self_check():
    delta=F(1,10**15)
    P=conic.model(["x","y"],[
        dict(name="H3_inverse_response_lower_test",entries=[[conic.form(x=1,y=-1)]]),
        dict(name="H3_inverse_response_upper_test",entries=[[conic.form(delta,x=-1,y=1)]]),
        dict(name="nonnegative_y",entries=[[conic.form(y=1)]]),
        dict(name="bounded_y",entries=[[conic.form(1,y=-1)]])])
    # Match the same general naming convention used for paired bounds.
    for b in P["blocks"]:b["name"]="fixture_"+b["name"]
    Q,change,scales,report=thin_intervals(P)
    if len(report["pairs"])!=1:raise ValueError("The narrow interval fixture was not conditioned")
    old={"x":F(1,3)+delta/2,"y":F(1,3)}
    values={v:old[v] for v in Q["variables"] if v in old}
    values[report["pairs"][0]["coordinate"]]=F(1,2)
    conic.check_point(P,old);conic.check_point(Q,values)
    if {v:conic.value(q,values) for v,q in change.items()}!=old:raise ValueError("Wrong exact narrow-band lift")
    values[report["pairs"][0]["coordinate"]]=F(2)
    try:conic.check_point(Q,values)
    except ValueError:pass
    else:raise ValueError("The conditioned fixture admitted a point outside the original interval")
    return dict(exact_inverse_identities=True,positive_cone_scales=True,
        feasible_point_preserved=True,outside_interval_point_rejected=True)
