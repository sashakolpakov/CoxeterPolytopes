# Gamma16 is reflection-maximal without a volume formula

**Status (2026-08-10): direct proof closed.**  Let

\[
 S=U\oplus E_8(-1)\oplus E_7(-1)
\]

and let `W=W_2(S)` be the reflection group of Tumarkin's unique
finite-volume noncompact Coxeter polytope `P` in `H^16` with 19 facets.
Then `W` is maximal among **all discrete hyperbolic reflection groups**.
In particular it is maximal among arithmetic reflection groups.  No
Prasad volume formula, minimal-covolume theorem, parahoric classification,
or lattice-integrality restriction is needed.

## Direct proof

Suppose that `W < G`, where `G` is an arbitrary discrete reflection group
in `Isom(H^16)`.  Since `W` has finite covolume, every discrete overgroup
has finite index, so \([G:W]<\infty\).  A finite cover preserves
compactness, hence a chamber `Q` of `G` is noncompact and has finite volume.

Felikson--Tumarkin prove that if a finite-volume reflection group with a
`k`-facet chamber contains a finite-index reflection subgroup, then a
chamber of the subgroup has at least `k` facets.  Apply this with ambient
group `G` and subgroup `W`.  The 19 facet reflections of `P` generate `W`,
so

\[
  17\leq \#\operatorname{facets}(Q)\leq19.
\]

The three cases are already classified.

1. `17` facets would make `Q` a finite-volume Coxeter simplex in `H^16`.
   Such simplices exist only through dimension 9.
2. `18` facets would make `Q` an `n+2`-facet polytope.  The complete
   Kaplinskaja--Esselmann--Tumarkin classification has no such polytope in
   `H^16`: prisms stop by dimension 5; Tumarkin's noncompact product case
   is in `H^4`; and his pyramid tables occur in dimensions
   `3,...,13,17`, with no rows in dimensions `14,15,16`.
3. With `19` facets, Tumarkin's `n+3` theorem says that the noncompact
   finite-volume `H^16` polytope is unique.  Therefore `Q` is isometric to
   `P`.  The covolume identity

   \[
      \operatorname{vol}(P)=[G:W]\operatorname{vol}(Q)
   \]

   now gives `[G:W]=1`.

Thus no proper `G` exists.

## Hostile boundary checks

- Proving that every reflection in `O(S)` has norm two rules out only
  same-lattice mixed-root enlargement.  By itself it says nothing about a
  reflection overgroup obtained after a commensurable lattice change.  The
  facet-count proof really does cover those groups.
- Maximality of an integral lattice or of a parahoric normalizer would not,
  by itself, exclude a non-normal finite-index reflection overgroup of
  `W`.  The finite-index facet theorem is the missing geometric input.
- The diagram involution and its `9+8` eigenspace split are not needed for
  reflection maximality.  They are relevant only when identifying the full
  normalizer and explaining why a chamber symmetry is not a hyperplane
  reflection.
- The direction of the Felikson--Tumarkin inequality matters: `W` is the
  finite-index **subgroup**, so its 19-facet chamber bounds the number of
  facets of the overgroup chamber from above.
- Tumarkin's `H^16` uniqueness statement is for noncompact finite-volume
  polytopes.  That is sufficient here because a finite quotient of the
  noncompact orbifold `H^16/W` cannot be compact: if `H^16/G` were
  compact, its finite orbifold cover `H^16/W` would be compact as well.

## Primary sources

- A. Felikson and P. Tumarkin,
  [*On finite index reflection subgroups of discrete reflection groups*](https://arxiv.org/abs/math/0305093),
  Funct. Anal. Appl. 38 (2004), 313--314.  The abstract states precisely
  the finite-volume facet inequality used above.
- P. Tumarkin,
  [*Hyperbolic Coxeter n-polytopes with n+2 facets*](https://arxiv.org/abs/math/0301133),
  Math. Notes 75 (2004), 848--854.  Together with Kaplinskaja and
  Esselmann it completes the `n+2` classification; the source tables have
  no `H^16` row.
- P. Tumarkin,
  [*Hyperbolic Coxeter n-polytopes with n+3 facets*](https://arxiv.org/abs/math/0311272),
  Russian Math. Surveys 58 (2003), 805--806.  It proves existence and
  uniqueness in `H^16`.
