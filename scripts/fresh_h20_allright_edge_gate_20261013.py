#!/usr/bin/env python3
"""Exact incidence certificates for all sections crossing all-right edges.

Only standard-library integer/Fraction arithmetic is used. The continuous
argument is in the companion report; this checks its complete finite input,
including actual edge incidences and both all-right vertices.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = "fresh_h20_allright_edge_gate_20261013"
SOURCE = ROOT / "data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json"
OUT = ROOT / "results" / (TAG + ".json")


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def mv(A, v):
    return [dot(row, v) for row in A]


def run():
    raw = json.loads(SOURCE.read_text())["coxeter_input"]
    Q = [[-int(x) for x in row] for row in raw["basis_gram"]]
    W = [list(map(int, row)) for row in raw["wall_coordinates_as_columns_transposed"]]
    cov = [mv(Q, r) for r in W]
    G = [[dot(r, c) for c in cov] for r in W]
    h = [sum(r[j] for r in W) for j in range(22)]
    g = gcd(*h)
    h = [x // g for x in h]
    hc = mv(Q, h)
    assert dot(h, hc) == -14 and max(dot(r, hc) for r in W) < 0
    short = [i for i in range(210) if G[i][i] == 2]
    assert len(short) == 42
    color = {short[0]: 0}
    queue = [short[0]]
    while queue:
        i = queue.pop()
        for j in short:
            if i == j or not G[i][j]:
                continue
            assert G[i][j] == -1
            if j in color:
                assert color[j] != color[i]
            else:
                color[j] = 1 - color[i]
                queue.append(j)
    assert set(color) == set(short)
    sides = [sorted(i for i in short if color[i] == k) for k in range(2)]
    assert list(map(len, sides)) == [21, 21]
    side_rows = []
    count = Counter()
    for side, points in enumerate(sides):
        lines = sides[1 - side]
        assert all(G[i][j] == (2 if i == j else 0) for i in points for j in points)
        a = [sum(W[i][j] for i in points) for j in range(22)]
        v = [F(2, 7) * (h[j] + F(a[j], 2)) for j in range(22)]
        vc = mv(Q, v)
        vp = [dot(r, vc) for r in W]
        assert dot(v, vc) == -2 and dot(h, vc) < 0
        assert [i for i, p in enumerate(vp) if p == 0] == points
        assert max(vp) == 0 and all(vp[l] == -1 for l in lines)
        # Twenty-one orthogonal positive normals and v certify signature(21,1).
        delta = {p: [v[j] - W[p][j] for j in range(22)] for p in points}
        zeros = {}
        for p in points:
            dc = mv(Q, delta[p])
            values = [dot(r, dc) for r in W]
            assert dot(delta[p], dc) == 0 and dot(h, dc) < 0
            assert max(values) == 0
            zeros[p] = {i for i, z in enumerate(values) if z == 0}
            assert zeros[p] & set(points) == set(points) - {p}
            # Open segment v--delta[p] has exactly the other20 zero walls.
            count["allright_incident_edges"] += 1
        pair_rows = []
        for p, k in combinations(points, 2):
            incident = [l for l in lines if G[l][p] == G[l][k] == -1]
            assert len(incident) == 1
            l = incident[0]
            neighbors = [i for i in points if G[l][i] == -1]
            assert len(neighbors) == 5
            common = zeros[p] & zeros[k]
            assert common == (set(points) - {p, k}) | {l}
            # x=u*delta[p]+(1-u)*delta[k], x²=-4u(1-u).
            assert dot(delta[p], mv(Q, delta[k])) == -2
            leaves = sorted(set(neighbors) - {p, k})
            assert len(leaves) == 3
            assert all(G[l][i] == -1 for i in leaves)
            assert all(G[i][j] == (2 if i == j else 0) for i in leaves for j in leaves)
            pair_rows.append({"p": p, "k": k, "line": l, "d4_leaves": leaves,
                              "edge_zero_walls": sorted(common)})
            count["cusp_to_cusp_edges"] += 1
        avoidance = Counter()
        for p in points:
            pencil = [l for l in lines if G[p][l] == -1]
            assert len(pencil) == 5
            others = sorted(set(points) - {p})
            # A Coxeter cut at the firstedge has at most2 nonzero pairings.
            for size in range(3):
                for support in combinations(others, size):
                    allowed = [l for l in pencil if all(G[l][a] == 0 for a in support)]
                    assert len(allowed) >= 5 - size
                    avoidance[str(size)] += 1
                    count["support_avoidance_checks"] += 1
        side_rows.append({"point_walls": points, "line_walls": lines,
                          "finite_vertex": list(map(str, v)),
                          "cusp_rays": {str(p): list(map(str, d)) for p, d in delta.items()},
                          "second_edges": pair_rows,
                          "support_avoidance_counts": dict(avoidance)})
    out = {"status": "PASS_EXACT_CONTINUOUS_ALLRIGHT_EDGE_GATE",
           "source_sha256": sha256(SOURCE.read_bytes()).hexdigest(),
           "producer_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
           "counts": dict(count), "sides": side_rows,
           "scope": "Every real spacelike normal whose hyperplane crosses transversely an open edge incident to either displayed all-right finite vertex. Hyperplanes containing such an edge use the earlier all-right-vertex theorem. No claim for sections avoiding all42 open edges or for other ambient polytopes.",
           "continuous_argument": {
               "first_edge_point": "P=v-u*r_p, 0<u<1; v²=-2; r_p²=2; v.r_p=0",
               "sparse_pairings": "A PD Coxeter Gram projecting20 orthogonal roots has at most2 nonzero section-normal pairings: opposite signs give an obtuse pair and three equal signs give a non-PD negative triangle.",
               "normal": "e=c*r_p-c*u*v+sum_{a in S} b_a*r_a, |S|<=2, c!=0",
               "second_edge_point": "x=u*delta_p+(1-u)*delta_k=v-u*r_p-(1-u)*r_k",
               "line_pairing": "e.r_l=-c*(1-u)!=0",
               "projected_line_norm": "N=2-(e.r_l)^2/e²; 3/2<N<2 from the PD projected D4 block",
               "forbidden_actual_ridge": "cos²(theta)=1/(2*N), strictly between1/4 and1/3; no Coxeter label"}}
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"status": out["status"], "counts": out["counts"]}))


if __name__ == "__main__":
    run()
