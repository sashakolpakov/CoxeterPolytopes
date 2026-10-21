# Hostile audit of the Type-II downstairs anti-Gram compiler

**Verdict:** `PASS/FREEZE`.

The hostile archive is
[`enriques_type_ii_downstairs_anti_gram_census_hostile_audit.json`](../data/enumerations/enriques_type_ii_downstairs_anti_gram_census_hostile_audit.json),
rebuilt by
[`audit_enriques_type_ii_downstairs_anti_gram_census_hostile.py`](../scripts/audit_enriques_type_ii_downstairs_anti_gram_census_hostile.py).
It pins the production script, archive, report, tests, Shimada datum, both
primary-source PDFs and provenance notes, and Martin's source archive.  No
production file or research ledger was changed by this audit.

## Independent `32,768 -> 16` replay

The audit does not call the project's GAP-record parser.  A record-local
balanced-bracket parser extracts the three literal arrays needed from
Shimada's `12B` record and verifies

\[
  W G_{L_{10}} W^t=C_{12B}.
\]

It independently obtains twelve walls, fifteen intersection-one edges, a
connected support graph of cycle rank four, rational wall rank ten, and
mod-two wall rank nine.  It then generates all

\[
  2^{15}=32,768
\]

anti-Grams from scratch.

The production code makes sorted breadth-first-tree edges negative.  The
hostile replay instead makes sorted depth-first-tree edges positive.  This
different gauge again gives exactly sixteen classes.  For every class, all
`2^11=2,048` diagonal-sign conjugates are generated directly, and every one
returns to the same hostile canonical key.  The complete hostile key set is
the production key set after regauging.

The audit also avoids the project's inertia routine.  It factors each exact
characteristic polynomial and uses exact Sturm root counts.  The result is:

| inertia `(positive,negative,zero)` | raw assignments | switch classes |
|---|---:|---:|
| `(0,9,3)` | 4,096 | 2 |
| `(1,9,2)` | 2,048 | 1 |
| `(1,10,1)` | 12,288 | 6 |
| `(2,10,0)` | 14,336 | 7 |

Thus the signature gate independently leaves exactly two negative-semidefinite
rank-nine classes.

## Independent full half-sum lattices

The mod-two left kernel of the `12 x 10` wall matrix is exhaustively computed,
not imported.  It has eight elements and dimension three.  For each surviving
anti-Gram `D`, the audit chooses a reverse-lexicographic rank-nine delta basis,
different from the production basis, and factors

\[
  D=C G C^t.
\]

The congruence lattice image is then built in one step from the twelve images
of `2e_i` and the images of all eight kernel classes:

\[
  A=\left\{\frac12\sum_i z_i\delta_i:
  z\in\mathbb Z^{12},\ zW=0\pmod 2\right\}.
\]

A single exact HNF gives the image basis.  This is independent of the
production route, which first takes an HNF of the congruence lattice and then
a second HNF after projection.  Direct fractional-coset counting checks the
HNF covolume: one class modulo the reverse-lex delta span for the survivor and
two for the rejected orbit.

Exact rational LDL ellipsoid recursion, independent of PARI/GP, gives:

| residue | determinant | Smith diagonal | short vectors | class |
|---|---:|---|---:|---|
| survivor | `2048=2^11` | `2^8,8` | no `(-2)` vectors; 144 signed `(-4)` vectors | `D9(-2)` |
| rejected | `512=2^9` | `2^9` | 18 signed `(-2)` vectors | `<-2>^9` |

For the survivor, an independently chosen positive root system has nine
simple roots, unimodular coordinate determinant, and Dynkin arm lengths
`1,1,6`.  For the rejected class, the nine root pairs themselves form a
unimodular basis whose Gram after division by `-2` is `I9`.  These are exact
isometry proofs, not classifications from determinant alone.

## The 36 witnesses and sheet invariance

The ternary search is independently exhaustive.  For a vector in
`{-1,0,1}^12`, its support is exactly its reduction modulo two, so only the
eight kernel supports can occur.  The audit enumerates every signing of those
supports, fixing the first nonzero sign to quotient global negation.  It finds

```text
surviving orbit:  0 witnesses
rejected orbit:  36 witnesses, all of support 9.
```

It separately reruns the production-gauge search and reproduces the archived
coefficient hash
`6c7a35b3bafc46c5569742bcbe83c915681e7ef641175e0112a350654e357415`.
Transporting all 36 coefficients by the exact sheet switch from the production
gauge gives the complete hostile-gauge set, not merely the three displayed
examples.

The 36 coefficient presentations map into the one-stage HNF lattice as exactly
nine root pairs, four presentations per root pair.  Finally, for every one of
the 2,048 sheet switches and every one of the 36 witnesses, the audit directly
checks both congruence parity and square `-2`: 73,728 exact transport checks.

One useful warning emerged.  A `0/1`-only half-sum search is **not** invariant
under sheet relabelling.  Both production negative-BFS representatives have
zero binary witnesses, but three of the rejected roots become binary in the
positive-DFS hostile gauge.  This is not a disagreement in the theorem: the
full congruence lattice and the signed transport are invariant.  It does mean
that “no binary witness” must never be used as a gauge-free rootlessness
criterion.

## Primary-source calibration

The pinned Brandhorst--Shimada PDF was independently text-extracted and checked
at the cited locations:

- Definition 1.1 and Propositions 1.2--1.3 define the half-sum lattice `M_R`
  and identify the anti complement as `R-tilde(2)`;
- Table 1.1 row 184 has `R=D9`, primitive closure and anti ADE type again
  `D9`, existence flag `1`, and chamber type `12B`;
- Definition 1.4(ii) identifies the K3 Neron--Severi lattice of a generic
  carrier with `M_R`; and
- Section 7.4 explicitly states that a very general Type-II surface is
  `(D9,D9)`-generic Number 184 and has the `12B` nef chamber.

Martin's pinned TeX source is checked in memory from the hashed arXiv archive.
It contains the Type-II K3-cover equation, the involution `s -> -s-1` followed
by translation by `(0,0)`, and the fixed-point-free range
`beta != 0,-64`.  This is existence/equation lineage only; no lift Gram is
imported from Martin.

## Frozen conclusion and scope

The production compiler survives the hostile gate:

> The `12B` two-sheet problem has 32,768 labelled anti-Grams and sixteen sheet
> switching classes.  Exact signature leaves two rank-nine classes.  Their
> full half-sum lattices are respectively `D9(-2)` and `<-2>^9`; the latter is
> eliminated by genuine, sheet-invariant `(-2)` roots.  Hence exactly one
> wall-lift orbit survives, and it is realized by the known very general
> Type-II carrier.

This is a positive surface-to-wall compiler result, but it is not a new
surface or a new abstract hyperbolic Coxeter polytope.  It also does not yet
reconstruct the period or ample cone from the downstairs datum alone.
