#!/usr/bin/env python3
"""Independent integer conflict-graph and Fraction escape replay for H9.

The common-escape certificate already rejects every subset. This extra
audit checks the archived combinatorial search and its twelve examples.
"""
import argparse
import json
from pathlib import Path

from scripts.fresh_h4_blocks_independent_verify import mv
from scripts.fresh_multibasepoint_h9_targeted_verify import escape
from scripts.fresh_rank4_h9_selection_audit_20261013 import ROOT, OUTPUT as COMMON, read, setup

OUTPUT = ROOT / "results/fresh_rank4_h9_selection_audit_graph_20261013.json"
COMPLETION = ROOT / "results/fresh_rank4_h9_completion_20261013.json"


def sign(a, b):
    x, y = 2*a+b, b
    if y == 0:
        return (x > 0) - (x < 0)
    if x == 0 or (x > 0) == (y > 0):
        return (y > 0) - (y < 0)
    return ((x > 0) - (x < 0))*((x*x > 5*y*y) - (x*x < 5*y*y))


def coefficients(v):
    assert all(a.denominator == b.denominator == 1 for a, b in v)
    return tuple((int(a), int(b)) for a, b in v)


def ip(a, b):
    return (sum(x*u+y*v for (x,y), (u,v) in zip(a,b)),
            sum(x*v+y*u+y*v for (x,y), (u,v) in zip(a,b)))


def conflicts(p, denominator):
    a, b = p
    if sign(a, b) > 0:
        return True
    x, y = 4*(a*a+b*b), 4*(2*a*b+b*b)
    if sign(x-4*denominator, y) >= 0:
        return False
    return (x, y) not in {(0,0), (denominator,0), (2*denominator,0), (3*denominator,0),
                          (denominator,denominator), (2*denominator,denominator)}


def build():
    d, base, gram, old, walls, center, _, _, _ = setup(read(COMMON))
    roots = walls[len(old):]
    R = tuple(coefficients(r) for r in roots)
    N = tuple(coefficients(mv(gram, r)) for r in roots)
    norms = []
    for i in range(len(roots)):
        q = ip(R[i], N[i])
        assert q[1] == 0 and q[0] == d["norms"][i] and q[0] in (1,2)
        norms.append(q[0])
    edges, neighbors = set(), [set() for _ in roots]
    for i in range(len(roots)):
        for j in range(i):
            if conflicts(ip(R[i], N[j]), norms[i]*norms[j]):
                edges.add((j,i))
                neighbors[i].add(j)
                neighbors[j].add(i)
    assert edges == {tuple(edge) for edge in d["conflict_edges"]}
    assert len(edges) == len(d["conflict_edges"])
    completion = read(COMPLETION)
    reports = []
    for case in completion["cases"]:
        ids = case["selected_root_indices"]
        chosen = set(ids)
        assert sorted(chosen) == ids
        assert all(0 <= i < len(roots) for i in ids)
        assert all(not neighbors[i].intersection(chosen) for i in ids)
        assert all(neighbors[i].intersection(chosen) for i in range(len(roots)) if i not in chosen)
        candidate = (*old, *(roots[i] for i in ids))
        assert len(candidate) == case["wall_count"]
        assert case["status"] == "EXACT_SPACELIKE_ESCAPE"
        escape(case["hyperideal_witness"], candidate, gram)
        reports.append({"method": case["selection_method"], "walls": len(candidate),
                        "escape_norm": case["hyperideal_witness"]["norm"]})
    return {"status": "PASS_INDEPENDENT_COMPLETE_GRAPH_AND_ALL_RECORDED_SELECTIONS",
            "dictionary_roots": len(roots), "complete_graph_conflicts": len(edges),
            "selections": reports, "scope": "Replays the finite conflict graph and greedy maximal selections; the separate common strict escape excludes every subset of this dictionary."}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--verify", action="store_true")
    args = p.parse_args()
    result = build()
    if args.verify:
        assert result == read(OUTPUT)
    else:
        OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "selections"}, sort_keys=True))
    print("Selections checked:", len(result["selections"]))


if __name__ == "__main__":
    main()
