#!/usr/bin/env python3
"""Fraction-only replay of an escape satisfying the entire 670-root bank.

This verifier uses no project Golden arithmetic, numerical optimization,
network, or external solver. Build mode adds a strictly feasible perturbation
to the archived exact extreme ray; --verify only checks saved data.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from scripts.fresh_h4_blocks_independent_verify import (
    Z, O, parse, enc, add, sub, scale, sign, dot, mv, pair, matrix, inverse,
)
from scripts.fresh_multibasepoint_h9_targeted_verify import escape

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/fresh_rank4_h9_selection_audit_20261013.json"
DICTIONARY = ROOT / "results/fresh_rank4_h9_dictionary_20261013.json"


def read(path):
    return json.loads(path.read_text())


def maximum(values):
    ans = values[0]
    for value in values[1:]:
        if sign(sub(value, ans)) > 0:
            ans = value
    return ans


def setup(data):
    d = read(DICTIONARY)
    base = read(ROOT / d["base_source"])
    gram = matrix(d["gram"])
    assert gram == matrix(base["gram"])
    old = matrix(base["cases"][d["base_case"]]["all_walls"])
    roots = matrix(d["roots"])
    assert len(old) == data["old_wall_count"] == 382
    assert len(roots) == data["dictionary_root_count"] == 670
    walls = (*old, *roots)
    assert len(walls) == data["wall_count"] == 1052
    center = tuple(map(parse, d["strict_interior"]))
    assert center == tuple(map(parse, data["strict_interior"]))
    assert sign(center[0]) > 0 and sign(pair(center, gram, center)) < 0
    covectors = tuple(mv(gram, r) for r in walls)
    center_values = tuple(dot(center, n) for n in covectors)
    assert all(sign(x) < 0 for x in center_values)
    v = escape(data["hyperideal_witness"], walls, gram)
    assert sign(pair(v, gram, center)) < 0
    return d, base, gram, old, walls, center, center_values, covectors, v


def build(data):
    d, base, gram, old, walls, center, center_values, covectors, v = setup(data)
    epsilon = F(1)
    while True:
        x = tuple(add(a, scale(b, epsilon)) for a, b in zip(v, center))
        if sign(pair(x, gram, x)) > 0:
            break
        epsilon /= 2
    values = tuple(dot(x, n) for n in covectors)
    assert all(sign(p) < 0 for p in values)
    data["dictionary_sha256"] = sha256(DICTIONARY.read_bytes()).hexdigest()
    data["strict_spacelike_witness"] = {"epsilon": str(epsilon), "identity": "x = v + epsilon * center",
        "vector": [enc(y) for y in x], "norm": enc(pair(x, gram, x)),
        "center_pairing": enc(pair(x, gram, center)), "maximum_wall_pairing": enc(maximum(values))}
    data["strict_center_norm"] = enc(pair(center, gram, center))
    data["maximum_center_wall_pairing"] = enc(maximum(center_values))
    data["independent_replay"] = "PASS_FRACTION_ONLY_ALL_1052_WALLS_AND_STRICT_SPACELIKE_PERTURBATION"
    data["scope"] = "The exact strict spacelike witness satisfies every old wall and every saved direct/swap root. It excludes finite-volume completion by ANY subset of this finite dictionary, irrespective of mutual compatibility. It does not exclude additional roots outside this bank, including roots in the same parameter box omitted by the old-target cutting filter."
    return data


def verify(data):
    d, base, gram, old, walls, center, center_values, covectors, v = setup(data)
    assert sha256(DICTIONARY.read_bytes()).hexdigest() == data["dictionary_sha256"]
    record = data["strict_spacelike_witness"]
    epsilon = F(record["epsilon"])
    assert epsilon > 0
    x = tuple(map(parse, record["vector"]))
    assert x == tuple(add(a, scale(b, epsilon)) for a, b in zip(v, center))
    assert sign(x[0]) > 0
    q, pc = pair(x, gram, x), pair(x, gram, center)
    assert q == parse(record["norm"]) and sign(q) > 0
    assert pc == parse(record["center_pairing"]) and sign(pc) < 0
    values = tuple(dot(x, n) for n in covectors)
    assert all(sign(p) < 0 for p in values)
    assert maximum(values) == parse(record["maximum_wall_pairing"])
    assert pair(center, gram, center) == parse(data["strict_center_norm"])
    assert maximum(center_values) == parse(data["maximum_center_wall_pairing"])
    print(json.dumps({"status": "PASS_INDEPENDENT_ALL_DICTIONARY_ESCAPE", "walls": len(walls),
                      "original_escape_norm": data["hyperideal_witness"]["norm"],
                      "strict_epsilon": str(epsilon), "strict_escape_norm": record["norm"]}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--certificate", type=Path, default=OUTPUT)
    args = parser.parse_args()
    data = read(args.certificate)
    if not args.verify:
        data = build(data)
        args.certificate.write_text(json.dumps(data, indent=2) + "\n")
    verify(data)


if __name__ == "__main__":
    main()
