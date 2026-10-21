# Maximal-reflection lineage of the 84-wall Jacobian--Kummer H10 chamber

## Result

The 84-wall finite-volume chamber is not the fundamental chamber of the full
integral reflection group of its intrinsic lattice

\[
L=U\oplus E_8(-2)\oplus\langle-4\rangle.
\]

An exhaustive all-root Vinberg calculation gives a 13-wall maximal
integral-reflection chamber with 29 finite vertices and three cusps.  Its
normalized Coxeter relation matrix is entrywise identical, after an explicit
facet permutation, to the classical chamber for

\[
f_2^{10}=-2x_0^2+x_1^2+\cdots+x_{10}^2.
\]

Thus the maximal abstract chamber and reflection group are known.  The 84
facet reflections generate a finite-index reflection subgroup of that group.
This does **not** yet prove that the particular 84-wall subgroup chamber was
previously printed; that finer subgroup-lineage question remains open.

## Exhaustive integral-root search

Reflection integrality and the discriminant form reduce the primitive root
norms to `-2`, `-4`, and `-8`; all three types are searched.  With control
vector `(1,1,0,...,0)`, the height-zero stabilizer is
`A1 + E8 + A1`.  The positive-height search accepts exactly three further
roots:

| control pairing | square | accepted roots |
|---:|---:|---:|
| 2 | -4 | 2 |
| 4 | -4 | 1 |

The resulting 13 primitive simple roots have types

| square and divisibility | roots |
|---|---:|
| `(-2,1)` | 1 |
| `(-4,2)` | 11 |
| `(-4,4)` | 1 |

The exact cone stopping calculation finds 29 timelike rays, three lightlike
rays, and no spacelike/hyperideal or opposite-sheet rays.  Hence Vinberg's
finite-volume stopping criterion fires and proves that no later root can add
a facet.

## Exact classical identification

The pinned CoxIter snapshot contains four 13-facet records in `H10`.  Exact
edge-labelled graph isomorphism gives a unique match:

```text
graphs/group units qf/10-2,1,1,1,1,1,1,1,1,1,1.coxiter
```

The certificate records the permutation from the Vinberg root order to the
zero-based corpus order:

```text
9,1,2,3,4,5,6,7,10,11,12,8,0
```

Every relation is then checked entry by entry.  This is equality of the
normalized Coxeter Gram matrices, not a comparison by dimension, determinant,
or face counts.

The historical source is
[Vinberg's 1972 paper](https://www.mathnet.ru/eng/sm3032), which treats the
groups of units of these quadratic forms.  Later summaries record that
`f_2^n` is reflective exactly through `n=14`; in particular its `H10`
chamber long predates this project.

## What is and is not left open

All 84 mirrors were independently proved primitive and reflective in `L`, so
their group `W_84` lies in the full integral reflection group `W(L)`.  Both
its chamber and the 13-wall chamber have finite hyperbolic volume.  Therefore

\[
[W(L):W_{84}]<\infty.
\]

Solomon's exact Euler sum over every finite special subgroup gives

\[
\chi(W_{84})=-\frac{527}{11520},\qquad
\chi(W(f_2^{10}))=-\frac{31}{11147673600}.
\]

Therefore

\[
[W(f_2^{10}):W_{84}]
=\frac{\chi(W_{84})}{\chi(W(f_2^{10}))}
=16{,}450{,}560
=2^{10}3^3\cdot5\cdot7\cdot17.
\]

The same calculation gives

\[
\operatorname{vol}(P_{84})=\pi^5\frac{527}{340200}.
\]

This was independently replayed in CoxIter 1.1.8.  CoxIter returned the full
face vector

```text
(16560,107233,263403,347634,281244,146599,49643,10650,1350,84,1)
```

with 862 ideal vertices, the same Euler characteristic, and the same
covolume.  Thus the remaining priority question is strictly finer: whether
this particular index-16,450,560 reflection subgroup, or its 84-wall chamber,
occurs in a previous publication or database.

The bounded CoxIter audit now includes both the ordinary 1,214-file corpus
and its separate `testing_index2.tar.gz` reflection-doubling archive.  The
latter contains 2,959 Coxeter files.  Its `H10` chambers have 12, 13, 14, 15,
16, 19, 20, 28, or 44 facets; none has 84.  Exact web searches for the full
face vector, Euler characteristic, and index also returned no relevant prior
record.  This strengthens the bounded absence evidence but is not a global
priority proof.  Searching older Sterk/Enriques and unpublished
reflection-subgroup computations remains the next lineage gate.

## Reproducibility

- exact certificate:
  [`jacobian_kummer_h10_maximal_reflection_lineage.json`](../data/enumerations/jacobian_kummer_h10_maximal_reflection_lineage.json)
- exact Euler/index certificate:
  [`jacobian_kummer_h10_reflection_index.json`](../data/enumerations/jacobian_kummer_h10_reflection_index.json)
- Vinberg implementation:
  [`jk_h10_maximal_vinberg.py`](../src/coxeter_polytopes/search/jk_h10_maximal_vinberg.py)
- certificate generator:
  [`certify_jacobian_kummer_h10_maximal_reflection_lineage.py`](../scripts/certify_jacobian_kummer_h10_maximal_reflection_lineage.py)
- index generator:
  [`certify_jacobian_kummer_h10_reflection_index.py`](../scripts/certify_jacobian_kummer_h10_reflection_index.py)
- tests:
  [`test_jacobian_kummer_h10_maximal_reflection_lineage.py`](../tests/test_jacobian_kummer_h10_maximal_reflection_lineage.py)

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_jacobian_kummer_h10_maximal_reflection_lineage.py \
  --check-archive

PYTHONPATH=src:. pytest -q \
  tests/test_jacobian_kummer_h10_maximal_reflection_lineage.py

PYTHONPATH=src:. python3 \
  scripts/certify_jacobian_kummer_h10_reflection_index.py \
  --check-archive
```
