# Hostile audit: exact AGCat/finiteness synthesis

**Audited synthesis bytes:**
`97abdea54e201bebda64917299eb905318c43b73009a7a061f0c320678d57887`

**Verdict:** **PASS.**  The synthesis accurately separates the globally
proved exact functor, the two scoped finiteness pushbacks, the conditional
K3/Mukai capacity theorem, and the still-open all-dimensional boundedness
step.  It does not substitute classical reflection-finiteness for the
requested AG pushback.

## Dependency and scope checks

- **Global exactness:** a marked smooth projective quadric determines the
  form up to scalar, polarity returns every root line, and the finite real
  framing selects inward rays and the hyperbolic component.  Normalized
  inner products and signs recover angles, parallelism, ultraparallel
  distances, and the same half-space intersection.  Hence the claimed return
  functor is exact on the image.  This is an algebraic arrangement carrier
  with finite real framing, not a claim that the carrier image is bounded.

- **Vinberg translation:** the repaired text correctly says that the
  nonselected ambient form is positive definite while the redundant
  (f\times f) facet Gram is positive semidefinite of rank (n+1).

- **Dimension two:** the synthesis identifies the exact real-DM orbifold
  curve and the Long--Maclachlan--Reid bridge, but explicitly labels the
  latter arithmetic--spectral.  It discloses congruence, the automorphic
  spectral gap, Zograf's inequality, and bounded-area arithmetic finiteness,
  and does not call this a pure KSBA/BAB consequence or import
  ABSW/Nikulin reflection finiteness.

- **Lambert:** the repaired text includes the essential primitivity step
  (k=\mathbf Q(A,B)=\mathbf Q(u)), so the integer norm has exponent
  ([k:\mathbf Q]-1).  The interval contraction and exact cutoff
  ([k:\mathbf Q]\le92) agree with the owner theorem.  It is stated only for
  the exact Lambert stable-pair stratum.

- **Higher-dimensional curve obstruction:** the compact arithmetic
  (\mathbf H^6) quotient has virtual cohomological dimension six, whereas a
  Fuchsian orbifold group has vcd at most two.  This correctly obstructs an
  exact quotient-curve replacement retaining the orbifold fundamental group.
  The synthesis expressly allows an unrelated decorated curve to encode the
  Gram and explains why that loses Fuchsian finiteness leverage.

- **`AUTO-TRACE`:** the inequality
  ([k:\mathbf Q](n+1)\le23) is stated only after a rational K3/Mukai
  cohomological action preserves a common standard Coxeter slice.  The text
  lists the unresolved integral gluing, period/crystal, lift, hidden-wall,
  characteristic, finite-type, and kernel gates.  It calls the inequality a
  necessary carrier-capacity test, not an unconditional Coxeter bound.

- **`MAX-ARR`:** both the finite-type envelope and zero-dimensional
  arithmetic exact-Gram closure are visibly labelled an open theorem.  The
  implication from `MAX-ARR` to finiteness is valid: a zero-dimensional
  closed substack of a finite-type stack has finitely many geometric points,
  and exact return separates the marked normalized Grams.  The alternative
  bounded degree/height formulation is only sufficient and is not presented
  as proved.

The synthesis therefore answers the original program without inflation:
the exact functor exists globally; an AG-mediated arithmetic--spectral
finiteness theorem exists in dimension two; the Lambert surface stratum is
finite by a new integral-point contraction; and global AG-derived finiteness
remains precisely open at `MAX-ARR`.
