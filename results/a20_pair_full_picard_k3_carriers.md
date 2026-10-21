# TOOT: both high `A20` chambers now have very-general full-Picard K3 carriers

## Positive theorem

The known finite-volume noncompact chambers `A20-C2a-H12-19` and
`A20-C2b-H11-20` each occur as a full-rank divisor-reflection subsystem on a
very general algebraic K3 surface.  More precisely, there are respectively
seven-dimensional and eight-dimensional lattice-polarized K3 period families
whose very general members have

```text
                         H12/19                    H11/20
rank / signature        13 / (1,12)               12 / (1,11)
determinant              256                       -1024
Smith invariants         1^6 2^6 4                 1^2 2^10
full-Picard moduli       7                         8
primitive wall span      index 1                   index 1
```

Every archived facet normal is primitive in this full Néron--Severi lattice,
every corresponding reflection is integral, and the walls span the whole
lattice.  This upgrades two special-period realizations to positive-dimensional
full-Picard algebraic families.  It is not a new abstract-polytope discovery.

## From one special K3 to two generic carrier families

The source construction puts both chambers on the same singular rank-20 K3:
the `H12/19` section is fixed by a non-symplectic order-four automorphism and
the `H11/20` section by a non-symplectic involution.  The new certificate uses
that surface only as an exact primitive-embedding witness.

For each saturated fixed kernel `S`, it reconstructs the common marked K3
lattice, writes an integral primitive embedding

\[
 S\hookrightarrow \Lambda_{K3}=U^3+E_8(-1)^2,
\]

and computes the saturated orthogonal complement.  The two complements have
signatures `(2,7)` and `(2,8)`, and the gluing indices are exactly 256 and
1024, the discriminant orders of their respective `S`.  Period surjectivity,
together with avoidance of the countable collection of extra algebraic
hyperplanes, then gives very general projective K3 surfaces with `NS(X)=S`.

No source automorphism is asserted to deform to those very general surfaces.

## The carrier denominator is forced

For each complete wall system the certifier exhausts all integral
reflection-stable isotropic subgroups of the discriminant module.  In both
cases the only possible subgroup orders are

```text
integral stable orders    1, 2
even stable orders        1
```

The unique proper stable integral overclosure has index two, but its Gram
lattice is odd.  Its determinant is 64 in the `H12/19` case and -256 in the
`H11/20` case.  Consequently, the original saturated kernel is the **unique
even reflection-stable carrier** and hence the only one of these stable
denominators that can be the Néron--Severi lattice of a K3 surface.

This is useful beyond these two rows: the continuous surface-existence step
can be attached only after the discrete discriminant-form denominator has
been fixed.  Here the wall group itself fixes it.

## Exact wall and volume gates

The primitive mirror types are

```text
H12/19    7*(-2,div1) + 12*(-4,div2)
H11/20    3*(-2,div1) + 6*(-2,div2) + 11*(-4,div2)
```

Explicit primitive positive classes have squares 32,068 and 5,386 and pair
strictly positively with every oriented wall.  The finite-volume conclusion
is pinned in two independent exact lanes: a literal Vinberg all-edge and
full affine-neighbourhood audit, and rational GMP double description.

```text
                         H12/19                    H11/20
finite vertices          527                       311
cusps                     25                        29
hyperideal rays           0                         0
elliptic codim-1 edges    3,846                     2,229
bad endpoint tests        0                         0
bad cusp tests            0                         0
```

Thus both chambers are finite-volume and noncompact; no hyperideal ray is
being mistaken for a cusp.

## Exact boundary

These are period-level existence theorems, not equation-level families.  The
mixed square-minus-four mirrors are divisor reflections, so neither chamber
is claimed to be the full nef cone.  Picard-lattice reflections are not
asserted to extend to Hodge isometries or surface automorphisms.  The special
order-four automorphism and involution remain valid on the source singular
K3 only.  The chambers themselves were already certified, so the gain is two
new project-level full-Picard carriers and zero new geometric slots.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_a20_pair_full_picard_k3_carriers.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_a20_pair_full_picard_k3_carriers.py
```

Frozen payload SHA-256:

`2d6c68b8aed6fbe7c39b527d94a23a1b7e466504d3a6097b13536f9edc96cbfe`
