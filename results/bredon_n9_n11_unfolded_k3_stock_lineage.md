# Exact K3-stock lineage of the Bredon `N9` and `N11` triality rows

**Status (2026-08-09): one exact positive return and one exact stock
nonmatch.**  The two Lorentzian unfoldings singled out by the ADEG triality
precompiler have now been compared with the complete finite-volume K3
root-chamber classification, not merely with diagrams sharing the same
dimension and number of walls.

- `ADEG-N11-01-branched-path`, orientation 1, is exactly the `C3` triality
  fold of Roulleau's complete 15-wall
  `U+E8(-1)+D4(-1)` K3 nef chamber.  The weighted Grams, integral wall
  lattices, triality actions, and folded fixed lattices agree through explicit
  unimodular changes of basis.  The existing order-six Weierstrass family
  therefore gives a named

  ```text
  Bredon N11 polytope -> projective K3 -> same Bredon N11 polytope
  ```

  realization.

- `ADEG-N9-01-branched-path`, orientation 0, is **not** the complete
  13-wall `U+A2(-1)+E8(-1)` stock.  It fails already as a weighted graph and
  also as an integral wall lattice.  This rejects the proposed stock reuse;
  it does not prove that the Bredon `N9` polytope has no carrier of another
  kind.

The `N11` result is a new named polytope-to-surface identification in this
project, not a new abstract Coxeter polytope.  Its abstract chamber was
already archived as `R14_UE8D4_C3_S3`, Scharlau--Kirschmer `List11` entry 4,
and Tumarkin's `11-tum04_12_05` chamber.

## 1. Exhaustive stock denominator

The high-rank K3 root classification contains nineteen complete
finite-volume chamber types for Picard ranks 11 through 19.  Rank and wall
count leave exactly one possible archived stock for each unfolding:

| Bredon row | unfolded rank / roots | unique classified stock with that pair |
|---|---:|---|
| `N9`, orientation 0 | `12 / 13` | `R12_UA2E8`, `U+A2(-1)+E8(-1)` |
| `N11`, orientation 1 | `14 / 15` | `R14_UE8D4`, `U+E8(-1)+D4(-1)` |

The exact source matrices are not inferred from this numerical filter.  They
are loaded from their original frozen wall archives and cross-checked against
the independent literal Vinberg archive and the Roulleau stock inventory.
Thus the first column selects the comparison denominator; entrywise Gram
isomorphism decides the result.

## 2. `N9`: exact nonmatch

The `N9` orientation with only vertex 0 long unfolds its 11 Coxeter facets to
13 square-two roots.  Its graph and the complete `R12_UA2E8` source have the
following exact invariants:

| invariant | Bredon `N9` unfolding | `R12_UA2E8` complete stock |
|---|---:|---:|
| vertices | 13 | 13 |
| edges | 12 | 13 |
| cycle rank | 0 | 1 |
| triangles | 0 | 1 |
| degree histogram | `1^5 2^6 3 4` | `1^2 2^9 3^2` |
| weighted-graph automorphism order | 12 | 2 |
| wall-lattice determinant | `-12` | `-3` |

Consequently there are zero weighted Gram isomorphisms.  More concretely,
the Bredon unfolding is a tree with a degree-four vertex: three of its arms
are the triality roots over the unique long wall.  Roulleau's complete stock
has maximum degree three and contains a triangle.  The determinant mismatch
also excludes an integral unimodular change of wall-lattice basis.

There is no compatible triality action on the proposed stock: the target
graph contains a visible `S3` permuting the long-wall triple, whereas the
complete `R12_UA2E8` weighted graph has automorphism group `C2`.  This is a
failure at the source graph, lattice, and action levels simultaneously.

The conclusion is deliberately narrow:

```text
N9 does not reuse the classified complete R12_UA2E8 nef-wall stock.
```

The `N9` same-polytope carrier question remains open.  A different divisor
lattice, a non-nef wall subsystem, or a genuinely new completion would need
a separate search and a complete-root/no-hidden-wall proof.

## 3. `N11`: exact unfolded wall identification

The successful Bredon orientation has vertices `0,...,11` short and vertex
`12` long.  Its long wall unfolds to target roots `13,14,15` in one-based
notation.  There are exactly six weighted Gram isomorphisms to Roulleau's
15-wall source, differing only by permutation of the three outer `D4` roots.
The lexicographically first target-to-source wall permutation is

\[
 (2,3,5,6,7,4,8,10,9,1,11,12,13,14,15).             \tag{3.1}
\]

It sends the target triple `(13,14,15)` exactly to source walls
`(13,14,15)`.  Those are the three outer roots attached to source wall 12 in
the displayed `D4` diagram.  Hence (3.1) respects the structure that produced
the Bredon square-six wall; this is not an accidental uncoloured graph match.

Let `G_B` be the positive-normal Bredon unfolded Gram and `G_R` the negative
intersection Gram of the Roulleau curves.  If `pi` denotes (3.1), the
certifier checks all 225 entries of

\[
                         (G_B)_{ij}=-(G_R)_{\pi(i)\pi(j)}.           \tag{3.2}
\]

Both systems use their first fourteen roots as integral wall-lattice bases.
Writing the permuted Roulleau wall coordinates in its archived
`A1,...,A14` basis gives a `14 x 14` matrix `C` with

\[
       \det C=-1,\qquad C G_{R,\mathrm{basis}}C^t=-G_{B,\mathrm{basis}}. \tag{3.3}
\]

Every dependent-wall coordinate also maps entrywise.  Thus (3.2) extends to
an integral lattice isometry

\[
       L_B\simeq U+E8(-1)+D4(-1),\qquad \det L_B=-4,                \tag{3.4}
\]

not just a Coxeter-graph isomorphism.

## 4. Triality and folded-Gram compatibility

On target unfolded roots, let

\[
                  c_B=(13\ 14\ 15)
\]

and fix roots 1 through 12.  The Roulleau `C3` generator is the same cycle on
source walls.  The permutation (3.1) satisfies

\[
                         \pi c_B=c_R\pi.                            \tag{4.1}
\]

The stronger lattice statement is also literal.  If `A_B,A_R` are the two
rank-fourteen integral action matrices in the bases of (3.3), then

\[
                         A_BC=CA_R,
       \qquad A_B^3=A_R^3=I.                                      \tag{4.2}
\]

Orbit-summing the three outer roots gives the long folded root of square six.
The induced target-to-source folded-facet permutation is

\[
                  (2,3,5,6,7,4,8,10,9,1,11,12,13).                \tag{4.3}
\]

Under (4.3), all 169 entries of the folded integral Gram agree, as do all
entries of the normalized Coxeter relation matrix.  In particular the sole
label-six edge is carried to the source pair consisting of wall 12 and the
orbit `{13,14,15}`.

The folded walls generate the full fixed lattice.  A second explicit change
of basis has determinant `+1` and identifies the two rank-twelve lattices:

\[
                 L_B^{C3}\simeq U+E8(-1)+A2(-1),\qquad \det=-3.    \tag{4.4}
\]

Thus graph, ambient divisor lattice, group action, folded Gram, and fixed
lattice all agree.

## 5. Surface and same-polytope return

The matched source is Roulleau's complete 15-wall nef chamber.  The existing
equation-level carrier is the very general minimal resolution of

\[
 y^2=x^3+t^3(t-1)(t-a)(t-b)(t-c),                                  \tag{5.1}
\]

with pairwise distinct nonzero `1,a,b,c`.  It has

\[
                         NS=U+E8(-1)+D4(-1),                        \tag{5.2}
\]

and the order-six automorphism

\[
                  (x,y,t)\longmapsto(\omega x,-y,t)                \tag{5.3}
\]

acts on `NS` by the order-three outer triality used in (4.1)--(4.2).

The archived fixed-nef-section audit proves that its 13 orbit walls are
irredundant and satisfy Vinberg's finite-volume criterion:

- dimension 11 and 13 facets;
- 27 finite vertices and one cusp;
- 162 elliptic edges, each with exactly two finite or ideal endpoints;
- both full affine-neighbourhood audits pass; and
- zero bad edges or hyperideal vertices.

Because (4.3) identifies the entire folded Gram, this is the same marked
Bredon `N11` polytope, not merely a chamber with the same header.  Hence the
new lineage statement is

\[
 \boxed{\text{Bredon `N11` is the K3 `D4`-triality chamber
 `R14_UE8D4_C3_S3`.}}                                               \tag{5.4}
\]

Its maximal reflection group has index one.  The same chamber is already
recorded as Scharlau--Kirschmer `List11` entry 4 and by the exact CoxIter
file `graphs/Tumarkin pyramids n+2/11-tum04_12_05.coxiter`.  Therefore
(5.4) is a positive polytope-to-surface interpretation of a known abstract
polytope, not a new hyperbolic example.

## 6. Reproduction

- Certifier:
  `scripts/certify_bredon_n9_n11_unfolded_k3_stock_lineage.py`.
- Tests:
  `tests/test_bredon_n9_n11_unfolded_k3_stock_lineage.py`.
- Triality input:
  `data/enumerations/bredon_adeg24_triality_census.json`.
- Complete source stocks:
  `data/enumerations/roulleau_non2elementary_involution_folds.json` and
  `data/enumerations/roulleau_ue8d4_all_subgroup_folds.json`.
- Equation-level K3 action:
  `data/enumerations/ue8d4_order6_weierstrass_family.json`.

Run

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_bredon_n9_n11_unfolded_k3_stock_lineage.py
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_bredon_n9_n11_unfolded_k3_stock_lineage.py
```

The tests freeze the unique stock denominator, both graph invariants, the
six `N11` wall isomorphisms, the long-root/triple correspondence, both
unimodular changes of basis, integral `C3` equivariance, folded Gram equality,
and the equation-level same-polytope handoff.
