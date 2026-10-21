#!/usr/bin/env python3
"""Fresh numerical H9 Gram deformations followed by polar truncation.

This is a bounded heuristic search, not an exact nonexistence certificate.
The original wall labels are explicitly bounded. Polar-wall labels are
tested against the nearest integer m>=2, without imposing a label cutoff.
Degenerate and cusp-bearing skeletons are parked, not ruled out globally.
All proposed polar walls are treated as retained facets; redundancy among
polar cuts is not removed. Original walls do remain facets (see report).
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
TOL = 1e-7


def magnitude(label):
    return {2: 0., 3: 1., 4: np.sqrt(2), 5: (1 + np.sqrt(5)) / 2,
            6: np.sqrt(3)}.get(label, 2 * np.cos(np.pi / label))


def coxeter_angle(pairing, norm_product):
    if pairing > TOL:
        return False
    squared_cosine = pairing * pairing / norm_product
    if squared_cosine >= 1 - TOL or squared_cosine <= TOL:
        return True
    proposed = max(2, int(round(np.pi / np.arccos(np.sqrt(squared_cosine)))))
    return abs(squared_cosine - np.cos(np.pi / proposed) ** 2) < TOL


def chain_choices(mode):
    if mode == "complete_3456":
        return list(product((3, 4, 5, 6), repeat=8)), (3, 4, 5, 6)
    original = (5, 3, 3, 3, 3, 3, 3, 5)
    alphabet = (2, 3, 4, 5, 6, 7, 8, 10)
    choices = {original}
    for count in (1, 2):
        for positions in combinations(range(8), count):
            for new_labels in product(alphabet, repeat=count):
                row = list(original)
                for position, label in zip(positions, new_labels):
                    row[position] = label
                choices.add(tuple(row))
    return sorted(choices), (3, 4, 5, 6, 7, 8, 10)


def analyze(gram, null_vector, counts):
    """Enumerate all cone rays through their two supported wall pairings."""
    null_vector = null_vector / np.linalg.norm(null_vector)
    if min(abs(null_vector)) < TOL:
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
            omitted_pairs.append((int(i), int(j)))
    pairings = np.array(pairings)
    vertex_gram = pairings @ inverse_on_image @ pairings.T
    squares = vertex_gram.diagonal()
    hyperideal = np.where(squares > TOL)[0]
    ideal_count = int(sum(abs(squares) <= TOL))
    counts[f"rays_hyperideal_{len(hyperideal)}_ideal_{ideal_count}"] += 1
    if ideal_count:
        counts["parked_cusp_bearing_skeleton"] += 1
        return None
    if not len(hyperideal):
        counts["numerical_compact_untruncated"] += 1
        return {"hyperideal": 0, "vertex_omitted_pairs": omitted_pairs,
                "vertex_squares": squares.tolist()}
    for i in hyperideal:
        for j in omitted_pairs[i]:
            if not coxeter_angle(pairings[i, j], 2 * squares[i]):
                counts["polar_old_angle_failure"] += 1
                return None
    counts["passes_all_polar_old_angles"] += 1
    for i, j in combinations(hyperideal, 2):
        if not coxeter_angle(vertex_gram[i, j], squares[i] * squares[j]):
            counts["polar_polar_angle_failure"] += 1
            return None
    counts["all_angle_survivor_requires_exact_and_volume_audit"] += 1
    return {"hyperideal": len(hyperideal), "vertex_omitted_pairs": omitted_pairs,
            "vertex_squares": squares.tolist(), "vertex_pairing_vectors": pairings.tolist()}


def search(mode):
    started = time.monotonic()
    chains, branch_labels = chain_choices(mode)
    counts = Counter()
    hits = []
    positions = list(combinations_with_replacement(range(9), 2))
    for chain_number, chain in enumerate(chains):
        base = 2 * np.eye(9)
        for i, label in enumerate(chain):
            base[i, i + 1] = base[i + 1, i] = -magnitude(label)
        eigenvalues = np.linalg.eigvalsh(base)
        negative_count = int(sum(eigenvalues < -TOL))
        if min(abs(eigenvalues)) < TOL:
            counts["parked_singular_base"] += 1
            continue
        if negative_count > 1:
            counts["base_multiple_negative_directions"] += 1
            continue
        counts["positive_base" if negative_count == 0 else "lorentz_base"] += 1
        inverse = np.linalg.inv(base)
        for left, right in positions:
            for left_label, right_label in product(branch_labels, repeat=2):
                a, b = magnitude(left_label), magnitude(right_label)
                s0 = 2 - a * a * inverse[left, left]
                s1 = 2 - b * b * inverse[right, right]
                cross = a * b * inverse[left, right]
                if negative_count == 1 and min(s0, s1) <= TOL:
                    continue
                if negative_count == 0 and max(s0, s1) >= -TOL:
                    continue
                for sign in (-1, 1):
                    dotted = -cross - sign * np.sqrt(s0 * s1)
                    if dotted <= 2 + TOL:
                        continue
                    counts["rank10_lorentz_skeleton"] += 1
                    last = np.array([np.sqrt(abs(s1)), -np.sign(s0) * sign * np.sqrt(abs(s0))])
                    coupling = np.zeros((9, 2))
                    coupling[left, 0], coupling[right, 1] = -a, -b
                    null_vector = np.r_[-inverse @ coupling @ last, last]
                    gram = np.zeros((11, 11))
                    gram[:9, :9] = base
                    gram[9, 9] = gram[10, 10] = 2
                    gram[left, 9] = gram[9, left] = -a
                    gram[right, 10] = gram[10, right] = -b
                    gram[9, 10] = gram[10, 9] = -dotted
                    survivor = analyze(gram, null_vector, counts)
                    if survivor is not None:
                        row = {"chain": list(chain), "attachment_positions": [left, right],
                               "branch_labels": [left_label, right_label], "schur_sign": sign,
                               "dotted_weight": float(dotted), "root_gram": gram.tolist(), **survivor}
                        hits.append(row)
                        print(json.dumps({"candidate": {k: v for k, v in row.items()
                                                       if k not in ("root_gram", "vertex_pairing_vectors")}}), flush=True)
        if chain_number % 4000 == 0:
            print(json.dumps({"chains_processed": chain_number,
                              "rank10_skeletons": counts["rank10_lorentz_skeleton"],
                              "seconds": round(time.monotonic() - started, 2)}), flush=True)
    return {
        "schema": "coxeter.fresh_gram_polar_h9.v1",
        "proof_status": "NUMERICAL_BOUNDED_SEARCH_ONLY",
        "mode": mode,
        "chain_count": len(chains),
        "branch_label_alphabet": list(branch_labels),
        "attachment_position_pairs": len(positions),
        "labelled_templates": len(chains) * len(branch_labels) ** 2 * len(positions),
        "polar_angle_test": "Nearest integer Coxeter label, with no imposed maximum label; floating tolerance 1e-7.",
        "counts": dict(sorted(counts.items())),
        "hits_requiring_exact_verification": hits,
        "seconds": round(time.monotonic() - started, 3),
        "scope": "Only the stated 11-node diagrams and wall lists retaining every proposed polar cut. Original-wall irredundancy is proved in the report, but redundancy among polar cuts is not removed; reduced Coxeter intersections may therefore be missed. Parked degeneracies, nonpolar cuts, other diagrams and other labels remain open.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("complete_3456", "two_twists"), default="complete_3456")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = search(args.mode)
    output = args.output or ROOT / f"data/enumerations/fresh_gram_polar_h9_{args.mode}_20261013.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"output": str(output), "survivors": len(result["hits_requiring_exact_verification"]),
                      "seconds": result["seconds"]}))
