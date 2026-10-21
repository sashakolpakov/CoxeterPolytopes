# Hostile audit of the provisional record-12 residual census

## Verdict

The archived `344`-orbit computation is reproducible **inside its declared
one-branch blocker filter**, but it is not yet an exhaustive census of the
whole record-12 graph shadow.

The problematic step rejects a one-sided branch whenever the complete
blocker of that branch contains a set whose order is not four or five.  The
Lanner order bound applies to the minimal nonfaces of the **coupled global
complex**.  No proof currently identifies the blocker of either one-sided
branch with the minimal-nonface clutter of the global complex or of a
specific deletion or link.  Consequently the one-sided test is not known to
be necessary and cannot be used in an exhaustive theorem.

The `344` objects remain valid, fully audited incidence candidates.  What is
withdrawn is only the claim that they exhaust all candidates over record 12.

## What was independently reproduced

The replay did not import
`coxeter_polytopes.combinatorics.h9_14_record12_residual`.

1. A fresh `geng -q -t 10` stream has `12,172` triangle-free graph
   representatives and the archived raw SHA-256
   `162cd507b93e74307ed9e88052da7bfb6f68e6cfeae6ccd2254c949c84ce0214`.
2. Recompiling and running the corrected branch generator over every
   normalized switching cut gives `56,310` switched-triangle-free
   `(X,p)` states and `294` survivors of the provisional one-sided blocker
   gate.  Its output SHA-256 is
   `21694c23cd1b8ef29b621afc7010268ce7396dc22647fc5928f39ed0c0672ee2`.
3. For every one of those `294` rows, a separate direct dualization of the
   *full* branch triple clutter reproduced the archived blocker histogram.
   This checks the corrected nonzero-`y` update, including all three
   extensions
   \[
        b\mathbin\cup Y,\qquad b\mathbin\cup\{s\},\qquad
        b\mathbin\cup\{t\}.
   \]
   There were no discrepancies.
4. Applying the local coupled-ridge feasibility equations leaves `39` raw
   states.  An independent coloured-graph canonicalization gives exactly
   the archived `28` branch orbits, including their raw multiplicities.
5. Testing both a target cut and its complement gives exactly `108` ordered
   branch-type pairs: `100` balanced and `8` trivial.  There are no
   nonzero-`y` pairs and no pair admitted only by complementary-cut
   normalization **within the provisionally filtered branch list**.
6. The balanced relative-embedding enumeration has `2,418` fixed-orientation
   candidates and raw SHA-256
   `0dd2f10bfa56619081dd29a49b120c34e648fd26fb8ebf2152c7bd8f392c3edb`.
   The full row/column, transpose, independent layer-complement, and
   layer-swap quotient gives `344` representatives and raw SHA-256
   `479832e0e28c86527b9642dc3022ba784e0d74f91f796f6c0f74e2d6c0ed44b1`.
   The resulting `(X,Z,multiplicity)` triples equal the archived set and
   their multiplicities sum to `2,418`.
7. The `14` trivial-cut double cosets all fail the global incidence gate.

Thus the implementation of the restricted census, the corrected
nonzero-`y` formula, the cut/complement normalization, and the finite group
quotient are internally sound.

## The completeness gap

Write `P` and `Q` for the two clutters of selected five-set complements and
`B(-)` for the blocker.  The global minimal-nonface clutter is

\[
 B(P\mathbin\cup Q)
   = \min_{\subseteq}\{p\mathbin\cup q:p\in B(P),\ q\in B(Q)\}.
\]

It does **not** follow that every member of `B(P)` or `B(Q)` is a global
minimal nonface.  A large `p in B(P)` can disappear after coupling because
all unions containing it are dominated by a union formed from another
one-sided transversal.  The elementary example

\[
 P=\{\{a,x_1\},\{a,x_2\}\},\qquad Q=\{\{a\}\}
\]

has `B(P)={{a},{x_1,x_2}}` but
`B(P union Q)={{a}}`: the large one-sided blocker member vanishes.  This
does not itself construct a record-12 residual, but it proves that blocker
duality alone cannot justify the discarded branches.

For record 12, the one-sided triples describe facets avoiding one of the two
size-two covers of the missing `K_{2,2}`.  The current proof does not show
that these facets are the complete facets of a deletion or link to which the
Lanner bound can independently be applied.  Until such a containment theorem
is proved, the one-sided blocker filter must be removed.

## Safe route to a complete census

An exhaustive replacement has a clear finite order of operations:

1. enumerate all triangle-free `X` and all normalized cuts `p` for which
   `X^p` is triangle-free;
2. enumerate every `y` satisfying the per-pair coupled-ridge feasibility
   equations, without imposing a one-sided blocker-order restriction;
3. quotient these coloured branches and pair both sides, checking both cut
   normalizations;
4. impose every coupled six-set ridge equation;
5. form the complete fourteen-vertex facet-complement clutter `P union Q`;
6. compute `B(P union Q)` directly and only then enforce the exact missing
   edges, allowed orders four and five, the Dehn--Sommerville schedule,
   ridge degrees, connectedness, local `g`/Macaulay constraints, and all
   required link audits; and
7. quotient the survivors by the full `D_8 times S_10` action, freeze every
   raw stream hash, and replay the final representatives independently.

Until that unfiltered coupled computation is frozen, theorem-facing text
must describe `344` as a restricted subclass, not as the whole record-12
completion census.
