# Exact lineage of the 206-wall characteristic-three H11 chamber

## Result

`F3-C2-H11-206` is a genuine finite-volume noncompact hyperbolic Coxeter
polytope, but its ambient integral lattice and maximal reflection group are
known.  More precisely, its saturated fixed lattice is

\[
U\oplus E_8(-2)\oplus A_1(-2)^2,
\qquad (r,a,\delta)=(12,10,1),
\]

and all 206 facet reflections form a finite-index reflection subgroup of the
14-wall full reflection group published by Alexeev for the even
2-elementary lattice `(12,10,1)`.

This closes the lattice and maximal-group lineage.  It does **not** close the
priority of the particular 206-wall finite-index subgroup chamber or of its
Fermat-quartic divisor construction.

The same actual fixed-divisor lattice now also carries the dual-Alexeev
29-wall chamber.  The unique stable odd index-eight closure on the `A16_6`
side has doubled dual equal to this lattice by an explicit determinant-one
isometry; all 29 oriented walls transport and generate it.  This is a second
subgroup chamber on the known lattice, not a claim that the 29- and 206-wall
chambers coincide.  See the
[`H11/29` Fermat carrier](a16_h11_29_fermat_char3_scaled_dual_carrier.md).

## Integral lattice identity

The Smith-kernel computation gives an even saturated lattice of signature
`(1,11)`, determinant `1024`, and Smith diagonal `1^2 2^10`.  The lineage
certificate goes beyond that fingerprint: it archives a determinant-one
matrix `P` and checks entry by entry that

```text
P^T G_fixed P = Gram(U + E8(-2) + A1(-2)^2).
```

Thus the Nikulin invariants are exactly `(12,10,1)`.  An independent rational
Hasse-invariant calculation also identifies the underlying rational quadratic
space with the classical odd-unimodular `I_(1,11)` space; this is not a new
arithmetic commensurability class.

## Exact match to Alexeev's maximal chamber

The pinned primary source is Valery Alexeev,
[*Reflexive hyperbolic 2-elementary lattices, K3 surfaces and hyperkahler
manifolds*](https://arxiv.org/abs/2209.09110).  Its distributed Sage data give
`coxeter_matrix[(12,10,1)]`, a 14-by-14 primitive-root Gram matrix.  The paper
states that this is the full reflective diagram `Gamma_r`, has 14 roots and
five maximal parabolic subdiagrams, and has trivial diagram automorphism
group.

The exact source freezes are:

| object | SHA-256 |
|---|---|
| arXiv source archive | `10e4002f6518f8bf6d669946775bbe42731f69277a52847c5602fa60ea33c2a6` |
| distributed matrix file | `4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457` |
| canonical 14-root matrix | `a12618b59ad911429dc4dfe1a18a1d7b6b64a7f5dae53fcf0c677f28c43eb9e5` |

The first twelve Alexeev roots have Gram determinant `-1024`, so they can be
an integral basis rather than merely a finite-index root sublattice.  The
certificate supplies an explicit determinant-minus-one basis matrix in the
Fermat fixed lattice.  It reconstructs the final two roots from their Gram
pairings and verifies the complete 14-by-14 matrix entrywise.  Their primitive
types are

```text
11 x (square -4, divisibility 2)
 2 x (square -2, divisibility 1)
 1 x (square -2, divisibility 2).
```

All fourteen reflections are integral involutive isometries.  Therefore the
match is an explicit integral root realization, not a collision of rank,
determinant, and facet count.

## The 206 induced mirrors

Every induced facet is reconstructed from the complete 5,944-wall
Kondo--Shimada stock in the saturated fixed basis.  The exact primitive-root
histogram is

```text
164 x (square -2, divisibility 2)
 30 x (square -4, divisibility 2)
 12 x (square -2, divisibility 1).
```

All 206 reflection matrices are integral, square to the identity, preserve
the fixed Gram, and negate their defining root.  Hence `W_206` is literally a
reflection subgroup of Alexeev's full 14-wall group.  Since both chambers
have finite covolume, the subgroup index is finite.  Its exact numerical
index is not obtainable from the even-dimensional Solomon-Euler ratio used
for the `H10` examples: in hyperbolic dimension 11 both finite-covolume groups
have rational Euler characteristic zero.  It instead requires an explicit
coset/tiling computation or an odd-dimensional volume comparison.

The independent literal Vinberg certificate for the induced chamber remains:
327,566 finite vertices, 18,354 cusps, 2,688,957 elliptic rank-ten edges, and
zero bad or unrecognized endpoint neighborhoods.  Thus none of the lineage
identification weakens the positive hyperbolic result.

## Bounded subgroup-chamber audit

The official pinned CoxIter snapshot at commit
`b800be48d8240fbabd2d1c9542e61054584ad7fe` has 1,214 ordinary records.  Its
entire `H11` facet histogram is `12:14, 13:6, 15:2`; the global maximum in the
snapshot is 50 facets.  The separate 2,959-file reflection-doubling archive
has `H11` histogram `15:2, 16:2, 18:2, 22:2, 30:2, 46:2`.  Consequently
neither pinned corpus contains a literal 206-facet record.

This is a bounded exclusion, not a worldwide priority theorem.  The precise
status is:

| identity level | status |
|---|---|
| hyperbolic Coxeter construction | proved finite-volume noncompact |
| integral lattice | known: `U + E8(-2) + A1(-2)^2` |
| maximal reflection chamber/group | known: Alexeev `(12,10,1)`, 14 walls |
| arithmetic commensurability class | known odd-unimodular `H11` rational space |
| literal 206-wall subgroup chamber | absent from two pinned corpora; global priority open |
| divisor-lattice realization | positive and potentially new; priority open |

## Reproducibility

- certificate:
  [`ss_fermat_quartic_char3_c2_h11_lineage.json`](../data/enumerations/ss_fermat_quartic_char3_c2_h11_lineage.json)
- generator:
  [`certify_ss_fermat_char3_c2_h11_lineage.py`](../scripts/certify_ss_fermat_char3_c2_h11_lineage.py)
- tests:
  [`test_ss_fermat_char3_c2_h11_lineage.py`](../tests/test_ss_fermat_char3_c2_h11_lineage.py)

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_ss_fermat_char3_c2_h11_lineage.py --check-archive

PYTHONPATH=src:. python3 \
  scripts/certify_ss_fermat_char3_c2_h11_lineage.py \
  --check-archive --verify-external

PYTHONPATH=src:. pytest -q \
  tests/test_ss_fermat_char3_c2_h11_lineage.py
```
