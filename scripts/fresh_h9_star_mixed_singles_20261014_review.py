#!/usr/bin/env python3
"""Independent core-domain and nonnegative-witness review.

Reuses only the prior independent Fraction field, never the producer's
SymPy/ANP arithmetic. New fallback matrices use determinant expansion and
Cramer's rule, rather than the producer's LDL recurrence. The old 7,840
affine determinant computations remain inherited through their review.
"""
import argparse
import json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
from fresh_h9_star_zero_20261014_own_three_review import (
    Z, O, COS, add, sub, neg, scale, mul, inverse, sign, total, dec_anp)

ROOT = Path(__file__).resolve().parents[1]
TAG = "fresh_h9_star_mixed_singles_20261014"
S = [3, 4, 6, 7, 9, 10, 12, 13]
CORE = sorted(S + [11, 14, 17])


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


@lru_cache(None)
def determinant(A):
    """Laplace expansion, memoizing submatrices; no pivot arithmetic."""
    if not A:
        return O
    if len(A) == 1:
        return A[0][0]
    value = Z
    for j, x in enumerate(A[0]):
        if x == Z:
            continue
        minor = tuple(tuple(v for k, v in enumerate(row) if k != j)
                      for row in A[1:])
        term = mul(x, determinant(minor))
        value = add(value, neg(term) if j % 2 else term)
    return value


def prefix(A, n):
    return tuple(tuple(row[:n]) for row in A[:n])


def enc(v):
    return list(map(str, v))


@lru_cache(None)
def nonnegative_failure(A):
    minors = [O]
    for k in range(1, len(A) + 1):
        q = determinant(prefix(A, k))
        minors.append(q)
        if sign(q) <= 0:
            break
    else:
        raise AssertionError("selected matrix is positive definite")
    # Earlier leading minors are positive, so the preceding block is PD.
    # Cramer's rule gives the exact Schur witness (x,1), x=-A0^-1 b.
    B = prefix(A, k)
    A0 = prefix(A, k - 1)
    x = []
    for j in range(k - 1):
        replaced = tuple(tuple(neg(B[i][-1]) if t == j else A0[i][t]
                               for t in range(k - 1)) for i in range(k - 1))
        x.append(mul(determinant(replaced), inverse(minors[-2])))
    v = tuple(x + [O])
    assert all(sign(t) >= 0 for t in v)
    value = total(mul(v[i], mul(B[i][j], v[j]))
                  for i in range(k) for j in range(k))
    assert sign(value) <= 0
    assert value == mul(minors[-1], inverse(minors[-2]))
    assert all(total(mul(B[i][j], v[j]) for j in range(k)) == Z
               for i in range(k - 1))
    return k, tuple(minors), v, value


def run():
    seen = {}
    def read(p):
        if p in seen:
            return seen[p]
        d = json.loads((ROOT / p).read_text())
        seen[p] = d
        for path, digest in d.get("sources", {}).items():
            assert sha(path) == digest, path
            if path.endswith(".json"):
                read(path)
        return d

    primary = f"results/{TAG}.json"
    d = read(primary)
    assert d["status"] == "PASS_CORE_SHARED_SINGLE_ALL_INTERNAL_LABELS"
    g = read("results/fresh_h9_tree_cube_wedges_20261014.json")["representatives"][0]["geometry"]
    V = [set(v["incident"]) for v in g["vertices"]]
    missing = {frozenset(m) for m in g["minimal_nonfaces"]}
    assert len(V) == 378 and d["rank_core"] == CORE
    assert any(set(S) <= v for v in V)
    assert all(any(set(pair) <= v for v in V) for pair in combinations(CORE, 2))
    faces = {frozenset(v & set(CORE)) for v in V}
    maximal = {f for f in faces if not any(f < h for h in faces)}
    assert len(maximal) == 9
    assert maximal == set(map(frozenset, d["complete_maximal_core_faces"]))

    own = {(3,14):3, (4,14):3, (11,12):2, (11,13):3,
           (12,17):2, (13,17):3}
    blocks = {tuple(S[2*i:2*i+2]):4 for i in range(4)}
    cross = {tuple(sorted((a,b))) for i,j in combinations(range(4),2)
             for a in S[2*i:2*i+2] for b in S[2*j:2*j+2]}
    got_cross, got_mutual, got_upper, kernels = set(), set(), set(), []
    for row in d["core_positive_kernel_guards"]:
        path = row["ordered_support"]
        assert len(path) == len(set(path)) and set(path) <= set(CORE)
        assert set(path) == set(row["support"]) <= V[row["actual_vertex"]]
        pair, trial = tuple(row["test_pair"]), row["test_label"]
        available = own | blocks | {pair:trial}
        bonds = {tuple(sorted((path[i],path[i+1]))):m
                 for i,m in enumerate(row["labels"])}
        if row["complete_triangle"]:
            bonds[tuple(sorted((path[0],path[-1])))] = 3
        # Every nonright matrix entry must follow the literal hypotheses or trial.
        assert all(available.get(p) == m for p,m in bonds.items())
        kind = row["kind"]
        if kind == "cross_block_right":
            assert pair in cross and trial == 3 and row["labels"] == [4,3,4]
            got_cross.add(pair)
        elif kind == "mutual_anchor_right":
            assert pair == (11,17) and trial == 3 and set(path) == {11,13,17}
            assert row["complete_triangle"]
            got_mutual.add(pair)
        else:
            assert kind == "third_contact_upper" and trial == 4
            assert path == ([3,14,11,13,17] if pair == (11,14)
                            else [3,14,17,13,11])
            assert pair in ((11,14),(14,17)) and row["actual_vertex"] == 295
            got_upper.add(pair)
        B = tuple(tuple(O if i == j else neg(COS[bonds[tuple(sorted((a,b)))]])
                        if tuple(sorted((a,b))) in bonds else Z
                        for j,b in enumerate(path)) for i,a in enumerate(path))
        assert B == tuple(tuple(dec_anp(x) for x in rr) for rr in row["gram"])
        k = tuple(map(dec_anp,row["positive_kernel"]))
        assert all(sign(x) > 0 for x in k)
        laplacian = tuple(tuple(mul(k[i],mul(B[i][j],k[j]))
                               for j in range(len(path))) for i in range(len(path)))
        assert all(total(rr) == Z for rr in laplacian)
        assert all(x[1:] == (0,0,0) for rr in laplacian for x in rr)
        assert all(sign(laplacian[i][i]) > 0 for i in range(len(path)))
        assert all(sign(laplacian[i][j]) <= 0 for i in range(len(path))
                   for j in range(len(path)) if i != j)
        kernels.append(dict(kind=kind,pair=list(pair),actual_vertex=row["actual_vertex"],
                            rational_laplacian=[[str(x[0]) for x in rr] for rr in laplacian]))
    assert got_cross == cross and got_mutual == {(11,17)}
    assert got_upper == {(11,14),(14,17)} and len(kernels) == 27

    expected_budgets = {11:(12,13),17:(12,13),14:(3,4)}
    for row in d["projection_budgets"]:
        leaf = row["leaf"]
        drop, keep = expected_budgets.pop(leaf)
        assert (drop,keep) == (row["removed_own_root"],row["retained_own_root"])
        assert own[tuple(sorted((leaf,keep)))] == 3
        support = (set(S)-{drop}) | {leaf}
        assert len(support) == 8 and support == set(row["support"]) <= V[row["actual_vertex"]]
    assert not expected_budgets and len(d["projection_budgets"]) == 3
    assert F(1,4) + 2*F(1,2) == F(5,4) > 1
    assert F(1,4) + F(1,2)*2 == F(5,4) > 1
    assert F(1,4) + 4*F(1,4) == F(5,4) > 1
    row = d["high_last_third_right"]
    assert set(row["support"]) == {12,13,14} <= V[row["actual_vertex"]]
    assert F(3,4)+F(1,4) == 1
    assert all(frozenset([12,13,x]) in missing for x in (11,17))
    # c>=cos(pi/6) gives D>=4, X=D(1+c)/2>=D/2>=2.
    assert 1/(1-F(3,4)) == 4 and F(4,2) == 2
    # Formal expansion in a=sqrt(F), b=sqrt(H).
    def pmul(p,q):
        out = {}
        for (i,j),x in p.items():
            for (k,l),y in q.items():
                out[i+k,j+l] = out.get((i+k,j+l),F(0)) + x*y
        return {k:v for k,v in out.items() if v}
    oneH, twoF, cross_square = {(0,0):F(1),(0,2):F(1)}, {(0,0):F(2),(2,0):F(1)}, {(0,0):F(1),(1,1):F(1)}
    lhs = pmul(oneH,twoF)
    for key,value in pmul(cross_square,cross_square).items():
        lhs[key] = lhs.get(key,F(0))-value
    assert {k:v for k,v in lhs.items() if v} == {(0,0):F(1),(2,0):F(1),(0,2):F(2),(1,1):F(-2)}

    oldpath = "results/fresh_h9_star_shared_high_20261014.json"
    old = read(oldpath)
    oldreview = read("results/fresh_h9_star_shared_high_20261014_review.json")
    assert oldreview["sources"][oldpath] == sha(oldpath)
    assert oldreview["complete_cases"] == 7840
    expected = {tuple(block)+(ia,ib,ic,ma,mb)
                for block in product((4,5),repeat=3)
                for ia,ib,ic in product([None]+list(range(6)),[None]+list(range(6)),[None,2,3,4,5])
                for ma,mb in product((2,3),repeat=2)}
    got, positive, fallback_indices = set(), 0, set()
    for i,row in enumerate(old["cases"]):
        key = tuple(row["block_labels"]+row["foreign_contacts"]+row["third_labels"])
        assert key not in got
        got.add(key)
        # Arithmetic identities are inherited, but their literal coverage and
        # encoded signs are checked afresh; this is not a determinant recensus.
        slope, intercept, bound = (dec_anp(row[k]) for k in ("slope","intercept","boundary_D4"))
        assert sign(slope) > 0 and bound == add(intercept,scale(slope,4))
        if row["status"] == "STRICT_POSITIVE_FOR_ALL_D_GREATER_THAN_FOUR":
            assert sign(bound) > 0
            positive += 1
        else:
            assert row["status"] == "NONPOSITIVE_ACTUAL_VERTEX_AT_RELAXED_LAST_INTERNAL_FOUR"
            fallback_indices.add(i)
    assert got == expected and len(got) == 7840 and positive == 7392
    assert len(fallback_indices) == 448

    records, used, type_ids, witnesses = [], set(), {}, []
    for row in d["new_core_fallbacks"]:
        index = row["old_case_index"]
        assert index in fallback_indices and index not in used
        used.add(index)
        oldrow = old["cases"][index]
        for key in ("block_labels","foreign_contacts","third_labels"):
            assert row[key] == oldrow[key]
        labels = own | blocks | {tuple(S[2*i:2*i+2]):m for i,m in enumerate(row["block_labels"])}
        labels |= {(11,14):row["third_labels"][0],(14,17):row["third_labels"][1]}
        for leaf,extra in zip((11,17,14),row["foreign_contacts"]):
            if extra is not None:
                labels[tuple(sorted((leaf,S[extra])))] = 3
        assert sorted([list(p)+[m] for p,m in labels.items()]) == row["relaxed_core_labels"]
        assert all(set(p) <= set(CORE) for p in labels)
        w = row["witness"]
        nodes = w["component_order"]
        assert set(nodes) <= set(w["containing_core_face"]) <= V[w["actual_vertex"]]
        assert frozenset(w["containing_core_face"]) in maximal
        expected_labels = [labels.get(p,2) for p in combinations(nodes,2)]
        assert expected_labels == w["component_pair_labels"]
        A = tuple(tuple(O if a == b else neg(COS[labels.get(tuple(sorted((a,b))),2)])
                        for b in nodes) for a in nodes)
        k,minors,v,value = nonnegative_failure(A)
        assert k == w["prefix_order"] and nodes[:k] == w["support"]
        pivots = [mul(minors[j],inverse(minors[j-1])) for j in range(1,k+1)]
        assert pivots == list(map(dec_anp,w["pivots"]))
        if A not in type_ids:
            type_ids[A] = len(type_ids)
            witnesses.append(dict(matrix_type=type_ids[A],matrix_size=len(A),failure_order=k,
                matrix=[[enc(x) for x in rr] for rr in A],
                leading_determinants=[enc(x) for x in minors[1:]],
                nonnegative_quadratic_vector=[enc(x) for x in v],quadratic_value=enc(value)))
        records.append(dict(old_case_index=index,matrix_type=type_ids[A],
                            actual_vertex=w["actual_vertex"],support=nodes[:k]))
    assert used == fallback_indices and len(records) == 448 and len(type_ids) == 39
    # Each explicit nonnegative vector remains nonpositive when unspecified
    # acute contacts are added or the relaxed internal labels increase.

    swap = list(range(18))
    swap[12],swap[13] = 13,12
    assert d["literal_last_root_swap"] == swap
    assert {frozenset(swap[i] for i in v) for v in V} == set(map(frozenset,V))
    assert {frozenset(swap[i] for i in m) for m in missing} == missing
    bounds = [S[2*i:2*i+2]+[4 if i<3 else 7] for i in range(4)]
    assert len(d["families"]) == 2
    for i,row in enumerate(d["families"]):
        required = own if i == 0 else {tuple(sorted((swap[a],swap[b]))):m for (a,b),m in own.items()}
        assert row["fixed"] == sorted([list(p)+[m] for p,m in required.items()])
        assert row["bounds"] == bounds
    paths = [primary,f"scripts/{TAG}.py",f"scripts/{TAG}_review.py",
             "scripts/fresh_h9_star_zero_20261014_own_three_review.py"]
    return dict(status="PASS_INDEPENDENT_MIXED_SINGLES_CORE_REVIEW",
        sources={p:sha(p) for p in paths},recursive_manifests=len(seen),
        arithmetic="Separate Fraction basis1,sqrt2,sqrt5,sqrt10. Laplace determinants and Cramer's-rule nonnegative witnesses; no producer arithmetic imports. Old affine determinants inherited through their independent review.",
        actual_vertices=len(V),rank_core=CORE,positive_kernel_guards=kernels,
        core_budgets=3,high_last_triangle=sorted(row for row in (12,13,14)),
        high_first_identity="1+(sqrt(F)-sqrt(H))^2+H>0",
        inherited_case_count=len(got),inherited_positive_cases=positive,
        fresh_fallback_count=len(records),selected_matrix_types=len(type_ids),
        matrix_witnesses=witnesses,fallback_placements=records,
        literal_swap_checked=True,
        scope="Star0 only: first3 finite internal labels>=4,last>=7;14 own3/3;11 and17 share one own3 and are right to the other last root. All noncore contacts free. Both last-root orientations. No whole-incidence exclusion or new label traversal.")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--verify",action="store_true")
    a = p.parse_args()
    answer = run()
    out = ROOT / "results" / f"{TAG}_review.json"
    if a.verify:
        assert answer == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(answer,indent=2)+"\n")
    print(answer["status"],answer["fresh_fallback_count"],answer["selected_matrix_types"])
