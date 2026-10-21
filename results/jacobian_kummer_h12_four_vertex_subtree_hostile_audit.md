# TOOT: hostile audit of the Jacobian–Kummer H12 four-vertex subtree census

Status: **PASS/FREEZE; all exact checks pass** (2026-08-08).

## Positive theorem under audit

Let `P` be the pinned 172-facet finite-volume noncompact
Jacobian–Kummer Coxeter polytope in `H^12`.  Its walls `W61,...,W172`
form a 112-wall stock of pairwise disjoint doubling walls.  Four seed
chambers joined by a colored `P4` or `K1,3` subtree therefore give a
fourfold Allcock gluing.

Up to the exact folded seed action and tree automorphisms, the complete
four-vertex census contains

```text
P4       1,964
K1,3       680
total    2,644.
```

All 2,644 resulting ordinary Coxeter matrices are pairwise nonisomorphic.
Their facet histogram is

```text
607 :   23
616 :  170
625 :  844
634 : 1607.
```

Each is a finite-volume noncompact hyperbolic Coxeter polytope in `H^12`,
has reflection-subgroup index and volume multiplier four, and retains a
weak full-Picard divisor-reflection carrier in

```text
U(2) + D5(-1)^2 + <-4>.
```

This is not a full-nef realization and makes no surface-automorphism claim.
Worldwide priority remains open.

## Independence of the hostile reconstruction

The verifier imports neither the production census nor any project module.
It starts again from the pinned Dolgachev–Shimada source and reconstructs:

- the saturated rank-13 fixed lattice of signature `(1,12)` and determinant
  `256`;
- all 172 oriented seed roots and a strict timelike interior;
- the order-23,040 ambient action, its order-768 involution centralizer, and
  the order-384 folded facet action;
- an independent six-generator action on the 112-wall stock.

Only after finishing the independent census does it compare the frozen
production rows.

## Complete coloring denominator

Generator BFS partitions all `1,379,952` admissible ordered `P4` words.
After reversal there are `696,192` unrooted words and 1,964 orbits.  A
separate Burnside calculation gives the same orbit count from fixed sums

```text
identity part   1,446,528
reversal part      61,824.
```

For `K1,3`, generator BFS partitions all `C(112,3)=227,920` unordered
triples into 680 orbits.  The independent cycle-index Burnside fixed sum is
`261,120`, again giving 680 after division by 384.  Thus all `924,112`
unrooted four-tree colorings are accounted for.

The raw, pre-quotient facet histogram is

```text
607 :   2,480
616 :  47,616
625 : 293,376
634 : 580,640.
```

## Independent boundaries and ordinary Coxeter matrices

For every representative, the hostile verifier transports all 172 seed
walls into the four chambers.  It identifies internal mirrors by cancellation
of opposite local orientations, rather than by the production
transformation-adjacency rule.  Every row has exactly three internal mirrors,
a strict timelike summed interior, and primitive integral reflective boundary
roots.  The direct boundary hashes, facet counts, ordinary relation
inventories, and degree-profile hashes agree with all 2,644 frozen production
rows.

The fixed walls `W1,...,W60` occur in every boundary.  Their root matrix has
Smith diagonal `1^13`; hence every boundary spans the full fixed Picard
lattice, not a proper sublattice.

## Structurally different nonisomorphism certificate

Sorted ordinary degree profiles already yield 2,535 types, leaving 102
collision groups containing 211 matrices; the largest group has size four.
The hostile verifier does not use the production color-refinement word on
these residues.

Instead it assigns weights

```text
m2,m3,m4,m6,infinity = 1,2,3,5,7
```

to form a symmetric zero-diagonal matrix `M`.  For every edge-labelled
Coxeter matrix it takes the sorted multiset

```text
( ordinary_label(i,j), (M^2)(i,j) ),  i<j.
```

This multiset is invariant under simultaneous row/column permutation.
Different multisets therefore rigorously imply nonisomorphism.  It separates
every member of every residual group.  All 211 matrices also pass an explicit
vertex-reversal permutation test.  The largest observed entry of `M^2` is
`30,089`, below the pre-multiplication bound `634*7^2=31,066`.

Combining the 2,535 degree classes with this independent residue invariant
gives exactly 2,644 distinct ordinary Coxeter types.

## Finite volume, noncompactness, and carrier scope

Finite volume is not inferred from the Picard lattice.  The seed is pinned
by both a literal Vinberg edge/endpoint/parabolic-neighborhood certificate
and an independent exact GMP cone enumeration, each classifying it as
finite-volume noncompact.  The finite four-vertex Allcock gluing contains
four interiors-disjoint seed chambers, proving index four, volume multiplier
four, and finite volume.  It contains a closed noncompact seed chamber, so
the glued polytope is noncompact.

The full-rank Smith gate and the pinned seven-dimensional K3 carrier family
give the weak algebro-geometric statement: every new reflection arrangement
lives in the same full Néron–Severi lattice.  Proving that a particular glued
chamber is an actual nef cone would require additional effectivity and
no-extra-wall arguments; none is asserted here.

## Production hardening found during audit

The hostile review caught two auditability gaps in the initial production
freeze candidate, without finding a changed mathematical row:

- the claim of a pre-bound before every `int64` multiplication omitted the
  possible four-reflection product and the transformed-interior operation;
- canonical numeric refinement IDs needed an explicit induction statement
  and permutation tests.

The production lane repaired both before this hostile freeze.  The final
script supplies eight named preflight entry bounds and tests one discrete and
one nondiscrete refinement case under reverse and even/odd vertex
permutations.  Its rebuilt archive and six focused tests pass.  The hostile
verifier pins that coherent repaired script/JSON/test/report quartet.

## Reproduction

```bash
python3 scripts/audit_jacobian_kummer_h12_four_vertex_subtree_hostile.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_jacobian_kummer_h12_four_vertex_subtree_hostile_audit.py
```

Machine-readable archive:
[`jacobian_kummer_h12_four_vertex_subtree_hostile_audit.json`](../data/enumerations/jacobian_kummer_h12_four_vertex_subtree_hostile_audit.json).

Verifier:
[`audit_jacobian_kummer_h12_four_vertex_subtree_hostile.py`](../scripts/audit_jacobian_kummer_h12_four_vertex_subtree_hostile.py).

Current hostile SHA-256 values:

```text
script  a6c92d37fdcf34243570ad01a062029974b386b46fd832840c5185ee6e629aec
JSON    83efbd60b5ded8886d1cca201a323f6a3ae4eacdd66fe110f41adb028796f7d8
test    df178dc4b6c35ad17c6e82127ab68c907cd3351f72ca0affb8b17df96b96d527
```

Fresh exact archive replay: `archive_check=PASS`, payload SHA-256
`83efbd60b5ded8886d1cca201a323f6a3ae4eacdd66fe110f41adb028796f7d8`.

Fresh focused hostile test result: `4 passed in 368.99s`.
