# Golden AGCat four-node packets: independent hostile audit

**Verdict: PASS.** A fresh exact implementation reproduces the complete
connected four-node tranche with either no dotted edge or exactly one dotted
edge. It also reproduces the rank-23--30 terminal-dotted path family. The
audit does not import the owner module, its certifier, or the triangle packet
helpers.

## Independent reconstruction

The replay independently rebuilds the 6,875-value two-place edge alphabet,
including its 6 finite/orthogonal and 6,869 dotted values. It implements
arithmetic in `Z[phi]`, inert/split/ramified prime-ideal valuation parity,
the conjugate-sign unit squareclass, and the corrected exact square solver

`5(n^2)^2 - (4A+2B)n^2 + B^2 = 0`

for `A+B*phi=(m+n*phi)^2`. Every accepted cycle product is reconstructed as
an identity-positive golden root and multiplied back on every packet.

The finite tranche is canonicalized under a separately generated `S4`
action. The unique dotted edge is moved to `01`, then its complement is
canonicalized under the independently generated order-four setwise
stabilizer `S2 x S2`. This gives exactly 2,380 connected finite orbits and
2,080 dotted-complement orbits, the latter representing
`2,080*6,869=14,287,520` expanded inputs.

All present triangles and all three Hamilton cycles are generated from the
graph rather than copied from the owner. Their element-squareclass equations
are imposed simultaneously. These equations are sufficient at the declared
field-level symmetrizability stage: a tree gauge determines vertex scalings,
and each chord is consistent exactly when its fundamental-cycle product is a
square. This does not by itself certify an integral reflection order.

## Determinants, orientations, and inertia

Instead of using the owner's closed determinant formula, the replay expands
the `3x3` and `4x4` determinants over all permutations. An odd acute cycle
uses the negative identity-positive root; an even cycle uses the positive
root. Conjugation is then applied to that same algebraic root, so no
independent second-place sign choice is made.

At the conjugate place every `2x2` minor is positive. Exact positivity of the
leading triple and full determinant is therefore Sylvester's criterion. In
the finite selected tranche, positive leading and full determinants give a
positive-definite restriction; a negative full determinant, together with a
positive two-dimensional restriction and interlacing, gives exactly one
negative direction. In the dotted tranche the leading `2x2` determinant is
negative; negative leading-three and full determinants make the remaining
LDL pivots positive and hence give exact inertia `(3,1)`.

The original owner schema conflated two disjoint boundary events as “24
degenerate packets.” The audit found 21 packets with zero leading-triple
determinant and 3 different packets with zero full-four-node determinant, at
each real place. The repaired schema records this `21/3` split explicitly.
Both dotted zero counts are zero, and the retained sets were unaffected.

## Exact replay

| tranche | compatible | conjugate positive / retained | selected split |
|---|---:|---:|---:|
| finite/orthogonal only | 327 | 119 / 119 | 5 positive, 114 Lorentzian |
| exactly one dotted | 532,734 | 180,402 / 180,402 | 180,402 Lorentzian |

The replay matches all 48 exact topology/sign strata and every topology
aggregate. In particular, retained diamonds and complete graphs are nonzero,
so overlapping cycle constraints were genuinely exercised. It also matches
the four owner streams byte-for-byte through their SHA-256 digests:

```text
finite compatible: 687762835339ac1f56fed84bcc41975e51da471abba4d97a7ab4cafc64590e7a
finite retained:   a8a04a852ce054ebfa1668156200998263939d7417d7c437662130b2320c50e8
dotted compatible: 0c988e10b93d804457dcfb9f238a832c546fd749b0a66751e67c875462153cc5
dotted retained:   7de03d51631e744ba34c248a23a933da855e841256f7f27a33c1cf09b8a0b0a1
```

## Global path-family check and scope

For the path with terminal product `c` and all later products one, the replay
independently verifies `D_j=2j-(j-1)c`, exact selected inertia one for `c>4`,
and conjugate positive definiteness iff `c'<2r/(r-1)`. For ranks 23 through
30 it reproduces the counts

`3590, 3584, 3577, 3571, 3566, 3561, 3556, 3552`

and all eight admissible-value stream hashes. The integral Cartan gauge has
`A01=-1`, `A10=-c`, later adjacent entries `-1`, and symmetrizer
`diag(c,1,...,1)`. Conjugate definiteness makes this `K`-form anisotropic, so
the full arithmetic orthogonal host is cocompact. The reflection subgroup is
not known to have finite index, and the displayed path cone is not a
finite-volume Coxeter simplex.

This PASS is exhaustive only for the stated local at-most-one-dotted
four-node tranche, plus the explicit path family. It does not cover local
packets with two or more dotted edges, complete rank-23--30 graph packets,
integral orders/conductors, or Vinberg completion. It proves neither a new
compact reflection chamber nor golden nonexistence in dimensions 22--29;
that question remains **OPEN**.

Reproduce with:

```sh
python3 scripts/independent_hostile_audit_golden_agcat_four_node_at_most_one_dotted.py --verify
PYTHONPATH=src python3 scripts/certify_golden_agcat_four_node_at_most_one_dotted.py --verify
PYTHONPATH=src pytest -q tests/test_golden_agcat_four_node_at_most_one_dotted.py tests/test_golden_agcat_four_node_at_most_one_dotted_hostile_audit.py
```
