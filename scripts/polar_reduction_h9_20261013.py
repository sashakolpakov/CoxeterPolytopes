#!/usr/bin/env python3
"""Rescreen the unchanged H9 templates after removing invisible polar cuts.

The retention theorem is proved in the companion report. For a simple
hyperideal vertex, its polar wall is visible in the hyperbolic slice iff
its incident old-wall diagram is connected (and therefore Lorentzian).
The template generation and numerical tolerances remain those of the
previous fresh_gram_polar_h9_20261013 search. No labels or heights are added.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import numpy as np

import fresh_gram_polar_h9_20261013 as original


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/polar_reduction_h9_20261013.json"
CONNECTIVITY = {}


def connectivity_table(gram):
    adjacent = abs(gram) > original.TOL
    np.fill_diagonal(adjacent, False)
    key = np.packbits(adjacent).tobytes()
    if key not in CONNECTIVITY:
        table = {}
        for omitted in combinations(range(11), 2):
            incident = set(range(11)) - set(omitted)
            seen = {min(incident)}
            todo = list(seen)
            while todo:
                i = todo.pop()
                for j in incident - seen:
                    if adjacent[i, j]:
                        seen.add(j)
                        todo.append(j)
            table[omitted] = seen == incident
        CONNECTIVITY[key] = table
    return CONNECTIVITY[key]


def analyze_reduced(gram, null_vector, counts):
    null_vector = null_vector / np.linalg.norm(null_vector)
    if min(abs(null_vector)) < original.TOL:
        counts["parked_zero_kernel_entry"] += 1
        return None
    positive = np.where(null_vector > 0)[0]
    negative = np.where(null_vector < 0)[0]
    if min(len(positive), len(negative)) < 2:
        counts["parked_redundant_facet_relation"] += 1
        return None
    projector = np.outer(null_vector, null_vector)
    inverse_on_image = np.linalg.inv(gram + projector) - projector
    pairings = []
    omitted_pairs = []
    for i in positive:
        for j in negative:
            row = np.zeros(11)
            row[i], row[j] = -abs(null_vector[j]), -abs(null_vector[i])
            pairings.append(row)
            omitted_pairs.append(tuple(sorted((int(i), int(j)))))
    pairings = np.array(pairings)
    vertex_gram = pairings @ inverse_on_image @ pairings.T
    squares = vertex_gram.diagonal()
    hyperideal = np.where(squares > original.TOL)[0]
    if any(abs(squares) <= original.TOL):
        counts["parked_cusp_bearing_skeleton"] += 1
        return None

    connected = connectivity_table(gram)
    retained = [int(i) for i in hyperideal if connected[omitted_pairs[i]]]
    removed = [int(i) for i in hyperideal if not connected[omitted_pairs[i]]]
    counts["retained_polar_walls"] += len(retained)
    counts["removed_invisible_polar_walls"] += len(removed)
    if removed:
        counts["skeletons_with_polar_reduction"] += 1
    if len(hyperideal) and not retained:
        counts["unexpected_no_visible_polar_wall"] += 1
        return None

    for i in retained:
        for j in omitted_pairs[i]:
            if not original.coxeter_angle(pairings[i, j], 2 * squares[i]):
                counts["retained_polar_old_angle_failure"] += 1
                return None
    counts["passes_all_retained_polar_old_angles"] += 1
    for i, j in combinations(retained, 2):
        if not original.coxeter_angle(vertex_gram[i, j], squares[i] * squares[j]):
            counts["retained_polar_polar_angle_failure"] += 1
            return None
    counts["all_reduced_angles_survive_requires_exact_audit"] += 1
    return {
        "hyperideal": len(hyperideal),
        "retained_polar_vertex_indices": retained,
        "removed_polar_vertex_indices": removed,
        "vertex_omitted_pairs": omitted_pairs,
        "vertex_squares": squares.tolist(),
        "vertex_pairing_vectors": pairings.tolist(),
    }


def main():
    # Reuse the archived generator without modifying its source or outputs.
    previous = original.analyze
    original.analyze = analyze_reduced
    try:
        result = original.search("complete_3456")
    finally:
        original.analyze = previous
    result.update({
        "schema": "coxeter.polar_reduction_h9.v1",
        "generator_source": "scripts/fresh_gram_polar_h9_20261013.py",
        "generator_sha256": sha256(Path(original.__file__).read_bytes()).hexdigest(),
        "retention_rule": "Keep a hyperideal-vertex polar cut iff its nine incident old walls have a connected Gram diagram.",
        "retention_proof": "results/polar_reduction_h9_20261013.md",
        "scope": "The same bounded 11-node templates as the previous scan, now with all invisible polar cuts removed before angle tests. Numerical Gram/ray decisions remain unproved. Singular/degenerate/cusp cases are parked; other diagrams, labels and nonpolar cuts remain open.",
    })
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"output": str(OUTPUT), "counts": result["counts"],
                      "survivors": len(result["hits_requiring_exact_verification"])}))


if __name__ == "__main__":
    main()
