#!/usr/bin/env python3
"""Exact classification and compactness obstruction for one-coordinate H8 lifts.

Old wall projections are preserved up to positive scaling. Acuteness and
Coxeter angles leave two nontrivial seeds up to diagram symmetry; both
contain an explicit positive affine null-root combination.
"""
import argparse
from itertools import combinations
import json
from pathlib import Path

from scripts.fresh_h4_blocks_independent_verify import (
    Z, O, P, parse, enc, add, sub, scale, mul, div, sign, total, dot, mv, pair,
    matrix, rref, kernel_one, determinant,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/known_examples/bugaenko_compact_h8.json"
OUTPUT = ROOT / "data/enumerations/fresh_h8_lift_h9_20261013.json"


def encoded_matrix(rows):
    return [[enc(x) for x in row] for row in rows]


def build():
    data = json.loads(SOURCE.read_text())
    Q = matrix(data["lattice_gram"])
    roots = matrix(data["simple_roots"])
    H = matrix(data["wall_gram"])
    assert tuple(tuple(pair(r, Q, s) for s in roots) for r in roots) == H
    assert len(roots) == 11 and len(rref(roots)[1]) == 9
    assert Q[0][0] == scale(P,-2) and all(Q[0][i] == Q[i][0] == Z for i in range(1,9))
    assert all(sign(determinant(tuple(tuple(Q[i][j] for j in range(1,k+1)) for i in range(1,k+1))))>0
               for k in range(1,9))
    reversal = (8,7,6,5,4,3,2,1,0,10,9)
    assert all(H[i][j] == H[reversal[i]][reversal[j]] for i in range(11) for j in range(11))
    assert all(H[i][i] == scale(O, 2) for i in range(11))
    assert all(sign(H[i][j]) <= 0 for i in range(11) for j in range(i))
    negative_edges = {frozenset((i,j)) for i in range(11) for j in range(i) if sign(H[i][j]) < 0}
    assert not any(all(frozenset(e) in negative_edges for e in combinations(triple,2))
                   for triple in combinations(range(11),3))
    unseen, components = set(range(11)), []
    while unseen:
        component, pending = set(), [min(unseen)]
        while pending:
            i = pending.pop()
            if i in component:
                continue
            component.add(i)
            pending.extend(j for j in range(11) if j not in component and H[i][j] == scale(O,-1))
        unseen -= component
        components.append(sorted(component))
    assert sorted(components) == [[0], [1,2,3,4,5,6,7,9,10], [8]]
    for i,j in ((0,1),(8,7)):
        assert H[i][j] == scale(P,-1)
        assert all(H[i][k] == Z for k in range(11) if k not in (i,j))
    k3, k4 = scale(P,2), sub(P,O)
    both_endpoint_values = []
    for a,b in ((k3,k3),(k3,k4),(k4,k4)):
        value = div(scale(mul(a,b),4),mul(add(scale(O,2),a),add(scale(O,2),b)))
        assert sign(value)>0
        assert sign(sub(O,value))>0 or sign(sub(value,O))>0 and sign(sub(scale(O,2),value))>0
        both_endpoint_values.append(enc(value))
    records = []
    for target, extra, affine_ids, name in (
        (3,k3,[0,1,2,3,4,5,6,9],"affine_E7"),
        (4,k4,[0,1,2,3,4,9],"affine_B5"),
    ):
        host = tuple(tuple((*row,Z)) for row in Q) + ((*(Z for _ in range(9)),extra),)
        assert sign(extra)>0
        lifted = tuple(tuple((*r,O if i==0 else Z)) for i,r in enumerate(roots))
        lifted = (tuple(div(x,P) for x in lifted[0]),*lifted[1:])
        assert len(rref(lifted)[1]) == 10
        gram = tuple(tuple(pair(r,host,s) for s in lifted) for r in lifted)
        assert len(rref(gram)[1]) == 10
        assert gram[0][0] == (scale(O,2) if target==3 else O)
        assert gram[0][1] == scale(O,-1)
        assert div(scale(mul(gram[0][1],gram[0][1]),4),mul(gram[0][0],gram[1][1])) == scale(O,1 if target==3 else 2)
        assert all(sign(gram[i][j])<=0 for i in range(11) for j in range(i))
        center = (*tuple(map(parse,data["geometry"]["interior_point"])),Z)
        assert sign(pair(center,host,center))<0 and sign(center[0])>0
        assert all(sign(pair(center,host,r))<0 for r in lifted)
        block = tuple(tuple(gram[i][j] for j in affine_ids) for i in affine_ids)
        coefficients = kernel_one(block)
        assert coefficients is not None
        if sign(coefficients[0])<0:
            coefficients = tuple(scale(x,-1) for x in coefficients)
        assert all(sign(x)>0 for x in coefficients)
        coefficients = tuple(div(x,coefficients[0]) for x in coefficients)
        delta = tuple(total(mul(c,lifted[i][k]) for i,c in zip(affine_ids,coefficients)) for k in range(10))
        assert delta != (Z,)*10 and pair(delta,host,delta)==Z
        assert sign(delta[0])>0 and sign(pair(delta,host,center))<0
        values = tuple(pair(delta,host,r) for r in lifted)
        assert all(sign(x)<=0 for x in values)
        assert all(values[i]==Z for i in affine_ids)
        assert len(rref(block)[1])==len(affine_ids)-1
        records.append({"old_endpoint":0,"changed_edge":[0,1],"old_label":5,"new_label":target,
            "extra_coordinate_square":enc(extra),"normalization":"lift old root0 by extra coordinate1, then divide that whole root by phi",
            "host_gram":encoded_matrix(host),"wall_normals":encoded_matrix(lifted),"wall_gram":encoded_matrix(gram),
            "wall_gram_rank":10,"wall_gram_signature":[9,1,1],"strict_center":[enc(x) for x in center],
            "strict_center_square":enc(pair(center,host,center)),"affine_type":name,"affine_indices":affine_ids,
            "affine_rank":len(affine_ids)-1,"positive_null_coefficients":[enc(x) for x in coefficients],
            "future_null_vector":[enc(x) for x in delta],"null_vector_pairings_with_all_old_walls":[enc(x) for x in values]})
    return {"schema":"coxeter.fresh_h8_single_coordinate_lifts.v1","source":str(SOURCE.relative_to(ROOT)),
        "label_three_components":components,"negative_pairing_graph_has_no_triangle":True,
        "same_sign_lift_support_maximum":2,"total_nonzero_lift_support_maximum":4,
        "both_endpoints_lifted_4_cosine_squared_cases_33_34_44":both_endpoint_values,
        "seeds":records,"conclusion":"Every nontrivial single-extra-coordinate Coxeter lift preserving all eleven old wall projections is one of these two seeds, up to endpoint symmetry and reversal of the extra coordinate. Each contains a future null vector that survives every wall acute-compatible with the retained seed, so no such completion is compact.",
        "scope":"One extra spatial coordinate; all old projected wall directions retained. This does not exclude changing or removing those old walls, nor noncompact finite-volume completions, nor other H9 constructions."}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--verify",action="store_true")
    args=p.parse_args();result=build()
    if args.verify:
        assert result==json.loads(OUTPUT.read_text())
    else:
        OUTPUT.write_text(json.dumps(result,indent=2)+"\n")
    print("PASS exact single-coordinate lift classification: two rank10 seeds, affine E7/B5 null obstructions")
    for seed in result["seeds"]:
        print("5->",seed["new_label"],seed["affine_type"],seed["positive_null_coefficients"])


if __name__=="__main__":
    main()
