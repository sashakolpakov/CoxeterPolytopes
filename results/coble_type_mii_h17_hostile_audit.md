# Hostile audit: the Type-MII `H17/48` claim survives

**Decision (2026-08-09): green, with a stronger literal volume
certificate.**  I rebuilt the Type-MII chamber without importing its project
certificate or its frozen wall dictionary.  The full `48`-wall cone has

```text
196,052 extreme rays
= 190,808 finite + 5,244 ideal + 0 hyperideal.
```

Every ray is on the time sheet selected by a strict timelike chamber point.
This is a direct exact double-description proof that the projectivized cone is
a finite-volume, noncompact polytope in `H17`.  In particular, the result is
not an inference from finite sidedness.  The independently obtained `5,244`
ideal vertices agree exactly with Alexeev's `5,244` maximal parabolic
subdiagrams.

## What was independently reconstructed

Starting with the standard basis of
`Bl_16(P1 x P1)`, I recovered

```text
F1.F2=1, E_ij^2=-1,
B_i=F1-sum_j E_ij,
B'_j=F2-sum_i E_ij,
C_sigma=F1+F2-sum_i E_i,sigma(i).
```

The canonical class is

```text
K=-2F1-2F2+sum E_ij,   K^2=-8,
sum(B_i+B'_i)=-2K.
```

The resulting curves have the asserted split

```text
8*(-4) boundaries + 16*(-1) exceptionals + 24*(-2) graph curves.
```

They span the odd unimodular `NS=I_(1,17)` with index one.  Their Gram is
acute, has only `m=2,3,4,infinity`, and every inequality has an exact timelike
relative-interior facet witness.

Adjoining the eight boundary half-classes has index `2^7`.  Doubling the form
gives an even lattice of determinant `-16`, Smith invariants `1^14 2^4`, and
Nikulin invariant `(18,4,0)`.  The primitive arithmetic walls have split

```text
24*(-2, divisibility 1) + 24*(-4, divisibility 2),
```

span the lattice integrally, and all `48` reflections are integral.

## The hidden `(-1)`-curve gate

Kondō's pinned primary source proves that the displayed `40` Coble--Mukai
roots are exactly the effective irreducible roots.  This exhausts all
`(-2)`-curves and all `(-1)`-curves meeting two different boundaries, but by
itself it does **not** exclude a `(-1)`-curve meeting one boundary twice.

For such a curve `E`, adjunction gives total boundary incidence `B.E=2`, and

```text
f = E + B_i/2
```

is a primitive isotropic CM class.  Moreover `2f=2E+B_i` is nef: it has zero
intersection with `E` and `B_i`, and nonnegative intersection with every
other irreducible curve.

The independent CM reconstruction gives `CM(S)=II_(1,9)`.  Exact double
description of its `40`-root chamber gives

```text
4,954 rays = 4,752 finite + 202 ideal + 0 hyperideal.
```

All `4,954` rays lie on one causal sheet.  This last check is logically
important: it proves that every nef isotropic class in the chamber is one of
the `202` null extreme rays, rather than merely assuming it.  All
`202*8=1,616` candidates `f-B_i/2` were tested in the full Picard lattice;
none is integral.  Their affine-type census is

```text
72  A7~ + 2A1~
96  A5~ + A2~ + 2A1~
18  A3~ + A3~ + 2A1~ + 2A1~
16  A2~ + A2~ + A2~ + A2~.
```

Therefore no hidden `(-1)`-curve exists.  Together with Kondō's root theorem,
the `48` displayed curves are exactly all negative curves.  Since the full
48-wall cone is causal, any effective curve of nonnegative square is on the
same closed positive sheet and cannot add another nef inequality.  Thus this
really is the full nef cone.

## Raw Alexeev match and the volume logic

The audit parses the pinned raw `(18,4,0)` ancillary Gram, rather than loading
the project's Alexeev reconstruction archive.  A fresh weighted-graph search
finds exactly `1,152` Gram isomorphisms.  Choosing the lexicographically first
one and using Alexeev's first `18` roots as a basis produces an integral
determinant-`±1` marking which carries all `48` roots entry by entry and
identifies the lattice forms exactly.

The published claims have distinct logical roles:

- “Vinberg's algorithm terminates after two steps” proves completeness of the
  finite simple-root list; termination alone does not prove finite volume.
- The maximal-parabolic list proves noncompactness once the chamber is known
  complete; its existence alone is not a volume proof.
- Alexeev's theorem
  `O+(L)=Aut(K4,4) semidirect W_r`, with group order `1,152`, proves that the
  reflection group has finite index in the arithmetic lattice and hence
  finite covolume.
- Independently of that published index theorem, the literal `196,052`-ray
  census above proves finite volume directly and reproduces all `5,244`
  cusps.

The conclusion is therefore robust: this is a genuine finite-volume
noncompact full-nef `H17/48` realization.  Its abstract Coxeter polytope and
arithmetic lattice were already known; no new abstract-polytope priority is
claimed.

## Reproduction

The independent script is
[audit_coble_type_mii_h17_hostile.py](../scripts/audit_coble_type_mii_h17_hostile.py),
and its frozen output is
[coble_type_mii_h17_hostile_audit.json](../data/enumerations/coble_type_mii_h17_hostile_audit.json).

```bash
PYTHONPATH=src:. python3 scripts/audit_coble_type_mii_h17_hostile.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_type_mii_h17_hostile_audit.py
```

The replay requires `cddexec_gmp` and the pinned local primary-source
extractions used elsewhere in the project.
