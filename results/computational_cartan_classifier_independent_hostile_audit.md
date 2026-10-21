# Independent hostile audit: computational Cartan-cycle classifier

**Date:** 2026-08-12  
**Verdict:** **PASS for the explicitly stated finite-box proposition.**

I audited the classifier frozen in the architecture report as if it were an
external implementation, including its exact scope:

* ground field \(\mathbf Q\);
* connected acute symmetrizable Coxeter simplices;
* rank \(3\le m\le8\), subject to the explicit assignment cap;
* a user-supplied finite set \(S\) of recognized rational Coxeter
  length-two products; and
* quotient by diagonal root gauge and all facet permutations.

The pilot censuses use only

\[
                  S=\{1,2,3,4\},
\]

corresponding to labels \(3,4,6,\infty\), with zero corresponding to
orthogonality/label (2).  They do not include dotted products (>4),
unrecognized finite angles, arbitrary number fields, non-simplex full wall
sets, or reflection-maximality.

## Clean-room census replay

I independently wrote a second enumerator without importing the classifier or
the chamber verifier.  It directly enumerated edge-color assignments,
connectivity, graph simple cycles, rational-square/integral cycle conditions,
the full symmetric-group orbit minimum, a separately constructed symmetric
Gram matrix, exact congruence inertia, and every deleted principal link.

| stage | rank 3 | rank 4 |
|---|---:|---:|
| all labelled assignments | 125 | 15,625 |
| connected | 112 | 15,104 |
| integral rational cycle packets | 68 | 2,976 |
| permutation orbits | 18 | 207 |
| Lorentzian | 13 | 200 |
| finite volume | 13 | 18 |
| compact | 4 | 1 |
| noncompact | 9 | 17 |

These counts agree exactly with both frozen JSON files.  The repository's
separate chamber engine then independently replays all 31 accepted chambers,
including interior point, relations, extreme rays, vertices, and cusps.

## Boundary-condition audit

1. **Recognized Coxeter labels.**  The census function accepts an arbitrary
   positive rational set, but certificate emission correctly rejects any
   value outside \(\{0,1,2,3,4\}\).  Proposition 5.1 explicitly assumes that
   (S) consists of recognized Coxeter relations, so this is an interface
   precondition rather than a completeness gap.  It would be safer for a
   production API to reject the set before enumeration.
2. **Dotted products.**  A dotted/ultraparallel pair has product (>4) and
   cannot occur between two facets of a finite-volume Coxeter simplex: every
   pair of simplex facets meets in its finite or ideal projective closure.
   Thus their omission is correct for the simplex proposition.  Such products
   are essential for non-simplex bounded frames and are expressly deferred to
   the host/Vinberg-completion backend.
3. **Singular forms.**  The full Gram must have inertia
   \((m-1,1,0)\), so globally singular Euclidean packets are excluded.
   Deleted principal matrices are allowed precisely as positive semidefinite
   corank-one affine links.  More singular or indefinite links are rejected.
4. **Finite-volume criterion.**  For a Lorentzian Coxeter simplex, all
   deleted principal subdiagrams spherical or affine is exactly the
   finite-volume criterion.  The independent chamber verifier provides a
   stronger geometric replay and detected no sheet, hyperideal, or cusp-rank
   defect in the 31 survivors.
5. **Integral gauge.**  At each prime, the inequalities
   \(x_j-x_i\ge-v_p(a_{ij})\) are difference constraints.  Their cycle
   obstruction is exactly negative valuation of a directed closed product.
   Requiring every simple cycle, including length-two backtracks, integral is
   sufficient because every closed walk decomposes into these cycles.  The
   relaxation constructs an integral gauge rather than assuming individual
   Cartan entries integral.
6. **Gauge/orbit completeness.**  A complete cycle packet determines the
   connected Cartan representation up to diagonal gauge.  Exhausting every
   vertex permutation gives the exact unmarked quotient for ranks at most
   eight; no graph-isomorphism heuristic is involved.

## Scope verdict

The implemented proposition is correct.  It is not a global arithmetic
Coxeter classification theorem and does not claim to find new examples.  Its
outputs are exhaustive only inside the displayed rational simplex/edge box.
The architecture correctly identifies the missing extensions: explicit field
manifests, prime-ideal gauges, non-simplex frame-to-host completion,
order/conductor lattice enumeration, certified Vinberg completion, and
reflection-maximality.

No mathematical defect was found within the stated scope.
