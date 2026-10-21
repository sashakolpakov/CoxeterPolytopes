# Exact lineage of the 17-wall H12 triality chamber

## Result

The 17-wall finite-volume chamber obtained from the K3 divisor-lattice
triality fold is **not a new abstract hyperbolic Coxeter polytope**.  It is
the full reflection chamber of the determinant-six lattice already present
as entry 6 of Scharlau--Kirschmer's `List12` census.

This is stronger than a numerical match:

1. an explicit unimodular matrix proves equality of the integral lattices;
2. every face number and the cusp count agree;
3. an exhaustive all-reflective-root Vinberg replay proves that the 17 K3
   roots generate the full lattice reflection group.

The distributed census archive contains no chamber roots or weighted
Coxeter Gram, so it does not permit an explicit archive-facet-to-K3-facet
permutation.  The chamber identity is nevertheless proved abstractly, up to
conjugacy by the full reflection group.  No surface automorphism enters the
argument.

## Frozen primary source

The source is the official Scharlau--Kirschmer
[`Reflective Lorentzian lattices`](https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html)
page and its
[`Lorentz.zip`](https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/Lorentz.zip)
archive, retrieved on 2026-08-06.

| object | SHA-256 |
|---|---|
| `Lorentz.zip` | `74a4dd3e534f1e793cba45765ff33575dbd331d5cdce6f34cdd64340b3f0d4a5` |
| `List12` | `5a1a4972fffc83ec072e38573bb3848a84b4b42d1747abbc3c708eaa6260ab88` |
| `Gram12` | `acc8f92bf41c533ea041440eccbdc333dfd17eb68c8448db6d4b694215877c0d` |

Entry 6 of `List12` is

```text
<6, "II_{12,1}([2^{1}]_{1} 3^{-1})",
 [17,130,596,1835,4017,6444,7680,6807,4433,2054,623,93], 3>
```

The source documentation says that `faces[k]` counts faces of dimension
`12-k`; the final `93` therefore counts finite vertices, while the last
entry `3` records cusps separately.  The matching entry of `Gram12` is a
13-by-13 integral Gram matrix of determinant `-6` and signature `(12,1)`.
The source also states that each displayed genus has class number one.

## Exact lattice identity

The K3 fixed lattice is

\[
S=U\oplus E_8(-1)\oplus A_2(-1)\oplus A_1(-1),
\qquad \operatorname{sig}(S)=(1,12),\quad |\det S|=6.
\]

The archived matrix has the opposite sign convention and visibly splits as

\[
E_8(+1)\oplus U\oplus A_1(+1)\oplus A_2(+1).
\]

The audit does not infer equality merely from signature and determinant.  It
constructs a unimodular 13-by-13 matrix `P`, with determinant `1`, and checks
entry by entry that

\[
P^{\mathsf T}(-G_{K3})P=G_{\mathrm{SK}}.
\]

Its matrix hash is
`b5892112c08e28409e0223c391dab586c8d3a3a5a336beffaaaa7111d9f91322`;
the full matrix and two checked orthogonal decompositions are in the JSON
certificate.  Thus the same integral lattice, and hence the same genus and
class, is proved.

## Complete face and cusp match

Restoring the three ideal vertices to the census's 93 finite vertices gives
the full face vector

\[
(f_0,\ldots,f_{12})=
(96,623,2054,4433,6807,7680,6444,4017,1835,596,130,17,1).
\]

This is exactly the independently enumerated K3 chamber vector.  Equivalently,
the comparison in the source's ordering is

```text
facets through finite vertices:
17,130,596,1835,4017,6444,7680,6807,4433,2054,623,93
cusps: 3
```

So the match is not a collision of determinant, facet count, and cusp count:
all twelve proper-face numbers agree.

## Exhaustive reflective-root proof

Let a primitive reflecting root satisfy `r^2=-2k` and have divisibility `d`.
Reflection integrality gives `k | d`, while `r/d` has exact order `d` in the
discriminant group.  Here

\[
A_S\cong A_{A_2(-1)}\oplus A_{A_1(-1)}\cong \mathbf Z/3\oplus\mathbf Z/2,
\]

whose class-square data are

| exact order | square modulo `2Z` | multiplicity |
|---:|---:|---:|
| 1 | `0` | 1 |
| 2 | `3/2` | 1 |
| 3 | `4/3` | 2 |
| 6 | `5/6` | 2 |

Testing the finite possibilities `k | d | 6` leaves exactly

\[
(k,d)=(1,1),(1,2),(3,3).
\]

Consequently every primitive reflective root has norm `-2` and divisibility
1 or 2, or norm `-6` and divisibility 3.  There are no other root norms to
search.

Choose the finite control vector

\[
p=(0,2,4,2,4,4,4,4,4,4,2,2,1),\qquad p^2=2.
\]

Its saturated orthogonal complement is negative definite of determinant 12.
An exact PARI/GP `qfminim` enumeration through norm 6 finds 91 pairs of
norm-2 reflective roots and one pair of norm-6 reflective roots.  Relative
to the twelve height-zero candidate roots, every one is either a nonnegative
or a nonpositive linear combination.  Their simple Gram is exactly

\[
D_{10}\oplus A_1\oplus A_1,
\]

where the last isolated root has norm 6.  Thus the finite stabilizer chamber
at `p` is complete.

For positive height, every integral root satisfying the preceding Vinberg
inequalities was searched inside rigorous positive-definite coordinate
bounds.  Each “extra root” query was exactly `UNSAT`:

| `-r^2` | `(p,r)` | height `(p,r)^2/(-r^2)` | recovered facets | extra root |
|---:|---:|---:|---|---|
| 2 | 1 | `1/2` | 12 | none |
| 6 | 3 | `3/2` | 14 | none |
| 2 | 2 | `2` | 1 | none |
| 2 | 3 | `9/2` | none | none |
| 6 | 6 | `6` | none | none |
| 2 | 4 | `8` | 16, 17 | none |

The norm-6 pairings need only be checked in multiples of 3 because their
divisibility is 3.  The final 17-wall polyhedron has already passed the
literal finite-volume Vinberg--CoxIter criterion: 93 finite vertices, three
full parabolic vertices, no hyperideal rays, and every elliptic edge has
exactly two endpoints.  Vinberg's finite-volume stopping criterion therefore
terminates the root algorithm at these 17 roots.

Hence these roots are the complete simple-root system of the full reflection
group `W(S)`, rather than merely generators of an unknown finite-index
reflection subgroup.

## Identity hierarchy and novelty status

| comparison level | status |
|---|---|
| same genus | proved |
| same integral lattice | proved by explicit unimodular isometry |
| same complete face vector and cusps | proved |
| explicit census-facet permutation / weighted-root Gram equality | unavailable: the archive distributes no roots |
| full reflection-group identity | proved independently by the all-root Vinberg replay |

The strongest defensible conclusion is therefore:

> The K3 triality construction recovers the full reflection chamber of the
> exact determinant-six Scharlau--Kirschmer lattice.  It is the same abstract
> finite-volume H12 Coxeter chamber, up to `W(S)`-conjugacy.

Accordingly, this is a valid positive hyperbolic example but not a new
abstract polytope, lattice, or reflection chamber.  Whether the particular
divisor-lattice/triality derivation is new as an algebro-geometric realization
is a separate literature question and is not claimed here.

## Reproducibility

- exact certificate:
  [`k3_rank15_order3_h12_scharlau_kirschmer_audit.json`](../data/enumerations/k3_rank15_order3_h12_scharlau_kirschmer_audit.json)
- audit script:
  [`audit_k3_rank15_order3_h12_scharlau_kirschmer.py`](../scripts/audit_k3_rank15_order3_h12_scharlau_kirschmer.py)
- independent K3 chamber certificate:
  [`k3_rho15_h12_cyclic_fold_audit.json`](../data/candidates/k3_rho15_h12_cyclic_fold_audit.json)
- tests:
  [`test_k3_rank15_order3_h12_scharlau_kirschmer.py`](../tests/test_k3_rank15_order3_h12_scharlau_kirschmer.py)

```bash
PYTHONPATH=src:. python3 \
  scripts/audit_k3_rank15_order3_h12_scharlau_kirschmer.py \
  --verify-source-archive /path/to/Lorentz.zip

PYTHONPATH=src:. python3 \
  scripts/audit_k3_rank15_order3_h12_scharlau_kirschmer.py \
  --verify-certificate

PYTHONPATH=src:. pytest -q \
  tests/test_k3_rank15_order3_h12_scharlau_kirschmer.py
```
