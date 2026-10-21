#!/usr/bin/env python3
"""Fresh determinant neighbors <-alpha>+E8+<1>, exact bounded search.

The recursion is adapted from the audited OG10 engine, with every coordinate
filtered by reflection integrality for its actual totally positive norm.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from itertools import groupby
import json
from pathlib import Path
import time
from typing import Sequence
from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.exact.matrix import inertia
from coxeter_polytopes.search.bugaenko_vinberg import RationalMatrix, Vector, _pair, _mat_vec
from coxeter_polytopes.search.og10_golden_h9_vinberg import (
    OG10GoldenH9HeightLevel, _reconstruct_root, _positive_definite_enumeration_data,
    h9_field_gram, h9_cone_roots,
)
from coxeter_polytopes.search.h11_vinberg import (
    RingInteger, _coordinate_alphabet, _ring_integer, _from_ring_integer,
    _ri_add, _ri_multiply, _ri_scale, _ri_sign, _ri_subtract, _integral_form,
    _standard_e8_gram,
)
from coxeter_polytopes.search.golden_path_host import bounded_golden_integers
from scripts.search_golden_path_h9_20261013 import complete_norm_pool

def host_gram(alpha):
    old=h9_field_gram()
    return tuple(tuple(-alpha if i==j==0 else ONE if i==j==9 else old[i][j] for j in range(10)) for i in range(10))

def enumerate_level(
    control_coordinate: Golden,
    previous_nonzero_roots: Sequence[Vector],
    inverse_e8: RationalMatrix,
    schur: Sequence[RationalMatrix],
    norm: Golden,
    alpha: Golden,
) -> OG10GoldenH9HeightLevel:
    x = control_coordinate
    target = 2 * norm + 2 * alpha * x * x
    if target.sign(0) < 0 or target.sign(1) < 0:
        raise ValueError("the control coordinate violates conjugate positivity")
    alphabet = _coordinate_alphabet(target)
    visited = [0] * 10
    schur_pruned = [0] * 10
    wall_pruned = [0] * 10
    leaf_norm_mismatches = 0
    solutions: list[Vector] = []

    target_integer = _ring_integer(target)
    alphabet_integers = tuple(_ring_integer(entry) for entry in alphabet.values if (2 * entry / norm).is_integral)
    wall_coefficients: list[tuple[RingInteger, ...]] = []
    wall_bounds: list[RingInteger] = []
    e8 = _standard_e8_gram()
    for wall in previous_nonzero_roots:
        wall_weights = _mat_vec(e8, wall[1:9])
        if (
            any(entry.sign(0) < 0 for entry in wall_weights)
            or any(entry.sign(0) < 0 for entry in wall[1:9])
            or wall[9].sign(0) < 0
        ):
            raise AssertionError(
                "dominant-wall pruning requires nonnegative coefficients"
            )
        # In w-coordinates, (r,wall)<=0 is
        #     2*w.wall_z + 4*u*wall_u <= 4*phi*x*wall_x.
        wall_coefficients.append(
            tuple(
                _ring_integer(entry)
                for entry in (
                    *(2 * entry for entry in wall[1:9]),
                    2 * wall[9],
                )
            )
        )
        wall_bounds.append(_ring_integer(2 * alpha * x * wall[0]))

    integral_schur = tuple(_integral_form(matrix) for matrix in schur)

    def extend(
        prefix: tuple[RingInteger, ...],
        minimum: tuple[RingInteger, int] | None,
        parent_partial_products: tuple[RingInteger, ...],
    ) -> None:
        nonlocal leaf_norm_mismatches
        depth = len(prefix)
        visited[depth] += 1
        if depth:
            assert minimum is not None
            minimum_numerator, minimum_denominator = minimum
            remainder = _ri_subtract(
                _ri_scale(target_integer, minimum_denominator), minimum_numerator
            )
            if _ri_sign(remainder, 0) < 0 or _ri_sign(remainder, 1) < 0:
                schur_pruned[depth] += 1
                return
            entry = prefix[-1]
            partial_products = tuple(
                _ri_add(
                    parent_partial_products[i],
                    _ri_multiply(entry, wall_coefficients[i][depth - 1]),
                )
                for i in range(len(wall_coefficients))
            )
            for partial, bound in zip(partial_products, wall_bounds, strict=True):
                if _ri_sign(_ri_subtract(bound, partial), 0) < 0:
                    wall_pruned[depth] += 1
                    return
        else:
            partial_products = parent_partial_products

        if depth == 9:
            assert minimum is not None
            minimum_numerator, minimum_denominator = minimum
            if minimum_numerator != _ri_scale(target_integer, minimum_denominator):
                leaf_norm_mismatches += 1
                return
            golden_prefix = tuple(_from_ring_integer(entry) for entry in prefix)
            solutions.append(_reconstruct_root(x, golden_prefix, inverse_e8))
            return

        denominator, next_form = integral_schur[depth]
        base: RingInteger = (0, 0)
        for i in range(depth):
            for j in range(depth):
                base = _ri_add(
                    base,
                    _ri_scale(
                        _ri_multiply(prefix[i], prefix[j]), next_form[i][j]
                    ),
                )
        linear: RingInteger = (0, 0)
        for i in range(depth):
            linear = _ri_add(
                linear,
                _ri_scale(prefix[i], 2 * next_form[i][depth]),
            )
        diagonal = next_form[depth][depth]
        for entry in alphabet_integers:
            child_minimum = _ri_add(
                _ri_add(base, _ri_multiply(linear, entry)),
                _ri_scale(_ri_multiply(entry, entry), diagonal),
            )
            extend((*prefix, entry), (child_minimum, denominator), partial_products)

    extend((), None, tuple((0, 0) for _ in wall_coefficients))
    gram = host_gram(alpha)
    if any(_pair(root, gram, root) != norm for root in solutions):
        raise AssertionError("height enumerator produced a vector of the wrong norm")
    if any(
        _pair(left, gram, right).sign(0) > 0
        for i, left in enumerate(solutions)
        for right in solutions[i + 1 :]
    ):
        raise AssertionError("same-height Vinberg roots do not meet acutely")
    return OG10GoldenH9HeightLevel(
        control_coordinate=x,
        height=alpha * x * x / norm,
        target=target,
        alphabet=alphabet,
        visited_by_depth=tuple(visited),
        schur_pruned_by_depth=tuple(schur_pruned),
        wall_pruned_by_depth=tuple(wall_pruned),
        leaf_norm_mismatches=leaf_norm_mismatches,
        roots=tuple(solutions),
    )


def run(alpha, height, output):
    start=time.monotonic()
    gram=host_gram(alpha)
    assert inertia(gram,embedding=0).positive==9
    assert inertia(gram,embedding=1).positive==10
    norms=complete_norm_pool(-alpha)
    controls=sorted((alpha*x*x/q,q,x) for q in norms
        for x in bounded_golden_integers((height*q/alpha,-q/alpha))
        if x.sign()>0 and (2*alpha*x/q).is_integral)
    roots=[]
    levels=[]
    inverse,schur=_positive_definite_enumeration_data(1)
    spatial=_mat_vec(inverse,(ONE,)*8)
    def save(completed,last_height):
        walls=(*h9_cone_roots(),*roots)
        epsilon=Golden(Fraction(1,2))
        while True:
            interior=(ONE,*(epsilon*v for v in spatial),epsilon)
            if _pair(interior,gram,interior).sign()<0 and all(_pair(interior,gram,w).sign()<0 for w in walls):
                break
            epsilon/=2
        payload={"alpha":alpha.to_json(),"dimension":9,"a1_norm":1,
            "norm_pool":[q.to_json() for q in norms],"requested_height":height.to_json(),
            "completed_height":(height if completed else last_height).to_json(),"completed":completed,
            "root_search_seconds":round(time.monotonic()-start,3),"levels":levels,
            "wall_count":len(walls),"gram":[[v.to_json() for v in row] for row in gram],
            "all_walls":[[v.to_json() for v in row] for row in walls],
            "strict_interior":[v.to_json() for v in interior],
            "scope":"Exact bounded root search with complete primitive norm candidate alphabet; initial E8+A1 follows from spatial unimodularity and trace-minimum argument. No global nonreflectivity conclusion."}
        output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    save(False,ZERO)
    for i,(h,group) in enumerate(groupby(controls,key=lambda row:row[0])):
        batch=[]
        for _,q,x in group:
            level=enumerate_level(x,roots,inverse,schur,q,alpha)
            batch.extend(level.roots)
            levels.append({"height":h.to_json(),"norm":q.to_json(),"control":x.to_json(),"roots":len(level.roots),"nodes":sum(level.visited_by_depth)})
        # Unit-associated candidates can coincide. Remove identical oriented
        # hyperplanes before checking same-height acuteness.
        unique={}
        for r in batch:
            unique.setdefault(tuple(v/r[0] for v in r),r)
        batch=list(unique.values())
        assert all(_pair(a,gram,b).sign()<=0 for j,a in enumerate(batch) for b in (*roots,*batch[:j]))
        roots.extend(batch)
        save(False,h)
        if batch:
            print(json.dumps({"height":str(h),"new_roots":len(batch),"walls":9+len(roots)}),flush=True)
    save(True,height)
    print(json.dumps({"alpha":str(alpha),"walls":9+len(roots),"norms":[str(q) for q in norms],"seconds":round(time.monotonic()-start,3)}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--alpha',required=True,help='a,b for alpha=a+b phi')
    p.add_argument('--height',type=int,default=30)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    run(Golden(*(int(v) for v in args.alpha.split(','))),Golden(args.height),args.output)
