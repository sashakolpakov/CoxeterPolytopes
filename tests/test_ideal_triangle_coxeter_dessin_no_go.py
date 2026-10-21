from __future__ import annotations

import json

from scripts.certify_ideal_triangle_coxeter_dessin_no_go import (
    DEFAULT_OUTPUT,
    boundary_roots,
    build,
    congruent,
    pairing,
    semicircle_root,
)


def test_frozen_archive_is_exact_replay() -> None:
    certificate = build(12)
    assert json.loads(DEFAULT_OUTPUT.read_text()) == certificate
    assert certificate["all_checks_pass"]


def test_general_norm_and_dotted_identities() -> None:
    n0 = semicircle_root(0)
    for j in range(80):
        nj = semicircle_root(j)
        assert pairing(nj, nj) == 1
        assert pairing(n0, nj) == 1 - 2 * j * j


def test_every_fixed_principal_level_eventually_collides() -> None:
    for modulus in range(2, 80):
        for j in range(5):
            assert congruent(
                semicircle_root(j), semicircle_root(j + modulus), modulus
            )


def test_every_boundary_gram_is_coxeter_acute() -> None:
    for N in range(1, 30):
        roots = boundary_roots(N)
        assert all(pairing(root, root) == 1 for root in roots)
        assert all(
            pairing(roots[i], roots[j]) <= 0
            for i in range(len(roots))
            for j in range(i)
        )

