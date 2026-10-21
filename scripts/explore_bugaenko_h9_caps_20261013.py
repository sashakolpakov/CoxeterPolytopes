#!/usr/bin/env python3
"""Heuristic H9 cap search; ONLY the separate exact verifier proves a result.

Dependencies for this exploratory script: NumPy and SciPy.  Numeric angle
tolerances, deduplication, Qhull, and finite dictionaries are deliberately
not treated as mathematical nonexistence certificates.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import time

import numpy as np
from scipy.spatial import HalfspaceIntersection


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/known_examples/bugaenko_compact_h8.json"
PHI = (1 + np.sqrt(5)) / 2


def build_stock(gram, labels, frame, base_inverse, useful_points=None):
    # Match the exact-looking standard constants without trigonometric drift.
    values = np.array([{2: 0., 3: -1., 4: -np.sqrt(2), 5: -PHI, 6: -np.sqrt(3)}[m]
                       for m in labels])
    controls = np.array(list(product(values, repeat=9)))
    found = {}
    basis_count = 0
    for omitted in combinations(range(11), 2):
        chosen = [i for i in range(11) if i not in omitted]
        block = gram[np.ix_(chosen, chosen)]
        if abs(np.linalg.det(block)) < 1e-7:
            continue
        basis_count += 1
        inverse = np.linalg.inv(block)
        extra = controls @ inverse @ gram[np.ix_(chosen, omitted)]
        square = np.einsum("ij,jk,ik->i", controls, inverse, controls)
        keep = square < -1e-8
        if useful_points is not None:
            keep &= square > -100
        for column in range(2):
            x = extra[:, column]
            keep &= np.any(abs(x[:, None] - values[None, :]) < 1e-7, axis=1) | (x <= -2 + 1e-7)
        indices = np.where(keep)[0]
        pairings = np.zeros((len(indices), 11))
        pairings[:, chosen] = controls[indices]
        pairings[:, list(omitted)] = extra[indices]
        if useful_points is not None:
            vectors = pairings[:, :9] @ base_inverse @ frame
            slack = (-vectors[:, 0, None] + vectors[:, 1:] @ useful_points[:, :-1].T
                     + np.sqrt(2 - square[indices])[:, None] * useful_points[:, -1])
            indices_in_batch = np.where(np.max(slack, axis=1) > 1e-6)[0]
        else:
            indices_in_batch = range(len(indices))
        for i in indices_in_batch:
            found[tuple(np.round(pairings[i], 8))] = square[indices[i]]
    stock = np.array(list(found))
    squares = np.array(list(found.values()))
    return stock, squares, {
        "labels": labels,
        "nonsingular_bases": basis_count,
        "assignments_tested": basis_count * len(controls),
        "retained_caps": len(stock),
        "additional_filter": None if useful_points is None else {
            "projected_norm_interval": [-100, 0],
            "must_cut_one_of_this_many_selected_horizon_points": len(useful_points),
        },
    }


class Search:
    def __init__(self, pairings, squares, gram, inverse, frame, interior, points):
        self.pairings = pairings
        self.inverse = inverse
        self.t = np.sqrt(2 - squares)
        self.vectors = pairings[:, :9] @ inverse @ frame
        self.forms = np.c_[self.vectors[:, 1:], self.t]
        self.roots = gram[:, :9] @ inverse @ frame
        self.interior = interior
        self.points = points
        self.hits = self.cut_lists(points)
        self.compatibility = {}
        self.coverage = {}
        self.nodes = 0
        self.deadline = time.monotonic() + 300

    def cut_lists(self, points):
        answer = []
        for begin in range(0, len(points), 100):
            slack = self.forms @ points[begin:begin + 100].T - self.vectors[:, 0, None]
            answer.extend(np.where(slack[:, i] > 1e-6)[0] for i in range(slack.shape[1]))
        return answer

    def compatible(self, index):
        if index not in self.compatibility:
            x = (self.pairings[:, :9] @ self.inverse @ self.pairings[index, :9]
                 + self.t * self.t[index])
            targets = np.array([-2 * np.cos(np.pi / m) for m in range(2, 21)])
            allowed = (x <= -2 + 1e-7) | np.any(abs(x[:, None] - targets) < 1e-6, axis=1)
            allowed[index] = False
            self.compatibility[index] = allowed
        return self.compatibility[index]

    def covered(self, index):
        if index not in self.coverage:
            self.coverage[index] = self.forms[index] @ self.points.T - self.vectors[index, 0] > 1e-6
        return self.coverage[index]

    def cover(self, chosen, live, covered):
        self.nodes += 1
        if time.monotonic() > self.deadline:
            raise TimeoutError
        if np.all(covered):
            return chosen
        best = None
        for point in np.where(~covered)[0]:
            options = self.hits[point][live[self.hits[point]]]
            if len(options) == 0:
                return None
            if best is None or len(options) < len(best):
                best = options
            if len(options) == 1:
                break
        scores = [int(np.sum(self.covered(i) & ~covered)) for i in best]
        for index in best[np.argsort(-np.array(scores))]:
            answer = self.cover(chosen + [int(index)], live & self.compatible(index),
                                covered | self.covered(index))
            if answer is not None:
                return answer
        return None

    def run(self):
        history = []
        status = "ITERATION_LIMIT"
        try:
            for iteration in range(50):
                self.coverage.clear()
                chosen = self.cover([], np.ones(len(self.t), bool), np.zeros(len(self.points), bool))
                if chosen is None:
                    status = "FINITE_NUMERICAL_DICTIONARY_INCOMPATIBLE"
                    break
                walls = np.concatenate([
                    np.c_[self.roots, np.zeros(11)],
                    np.array([[0.] * 9 + [-np.sqrt(2)]]),
                    np.c_[self.vectors[chosen], self.t[chosen]],
                ])
                halfspaces = np.c_[walls[:, 1:], -walls[:, 0]]
                rays = HalfspaceIntersection(halfspaces, self.interior).intersections
                excess = np.sum(rays * rays, axis=1) - 1
                bad = np.where(excess > 1e-6)[0]
                row = {"iteration": iteration, "sample_count": len(self.points),
                       "facets": len(walls), "rays": len(rays), "hyperideal": len(bad),
                       "maximum_excess": float(max(excess)), "chosen": chosen}
                history.append(row)
                print(json.dumps(row), flush=True)
                if len(bad) == 0:
                    status = "NUMERICAL_CANDIDATE_REQUIRES_EXACT_VERIFICATION"
                    break
                new_points = []
                for i in bad[np.argsort(-excess[bad])[:40]]:
                    direction = rays[i] - self.interior
                    a = direction @ direction
                    b = 2 * self.interior @ direction
                    c = self.interior @ self.interior - 1
                    parameter = (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)
                    new_points.append(self.interior + parameter * direction)
                new_points = np.array(new_points)
                self.hits.extend(self.cut_lists(new_points))
                self.points = np.vstack([self.points, new_points])
        except TimeoutError:
            status = "TIMEOUT_UNDECIDED"
        return {"status": status, "nodes": self.nodes, "history": history}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expand", action="store_true", help="also search labels 2,3,4,5,6 with the stated bounded useful-cap filter")
    parser.add_argument("--output", type=Path, default=ROOT / "data/enumerations/bugaenko_h9_caps_20261013_search.json")
    args = parser.parse_args()
    data = json.loads(SOURCE.read_text())
    gram = np.array([[float(x["a"]) + PHI * float(x["b"]) for x in row] for row in data["wall_gram"]])
    inverse = np.linalg.inv(gram[:9, :9])
    vertices = []
    for face in data["geometry"]["vertices"]:
        _, _, vh = np.linalg.svd(gram[np.ix_(face, range(9))] @ inverse)
        vector = vh[-1]
        if sum(vector) > 0:
            vector = -vector
        vertices.append(vector / np.sqrt(-(vector @ inverse @ vector)))
    vertices = np.array(vertices)
    eigenvalues, eigenvectors = np.linalg.eigh(gram[:9, :9])
    frame = eigenvectors @ np.diag(np.sqrt(abs(eigenvalues)))
    frame[:, 0] *= np.sign((vertices.mean(axis=0) @ inverse @ frame)[0])
    center = vertices.mean(axis=0) @ inverse @ frame
    interior = np.r_[center[1:] / center[0], .001]
    klein = vertices @ inverse @ frame
    klein = klein[:, 1:] / klein[:, :1]
    weights = np.random.default_rng(42).dirichlet(np.full(41, .08), 1000)
    points = np.vstack([klein, weights @ klein])
    points = np.c_[points, np.sqrt(1 - np.sum(points * points, axis=1))]
    stock, _, census = build_stock(gram, [2, 3, 5], frame, inverse)
    squares = np.einsum("ij,jk,ik->i", stock[:, :9], inverse, stock[:, :9])
    search = Search(stock, squares, gram, inverse, frame, interior, points)
    output = {"schema": "coxeter.bugaenko_h9_cap_heuristic_search.v1",
              "source_sha256": sha256(SOURCE.read_bytes()).hexdigest(),
              "proof_status": "NUMERICAL_EXPLORATION_ONLY",
              "scope": "A bounded finite cap dictionary. Incompatibility is not an H9 nonexistence theorem.",
              "initial_stock": census, "initial_search": search.run()}
    if args.expand:
        useful = np.vstack([search.points[:41], search.points[-80:]])
        extra, extra_squares, extra_census = build_stock(gram, [2, 3, 4, 5, 6], frame, inverse, useful)
        expanded = Search(np.vstack([stock, extra]), np.r_[squares, extra_squares],
                          gram, inverse, frame, interior, search.points)
        output.update(expanded_stock=extra_census, expanded_search=expanded.run())
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), "initial_status": output["initial_search"]["status"]}))


if __name__ == "__main__":
    main()
