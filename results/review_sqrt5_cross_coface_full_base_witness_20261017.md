# Independent exact review: one cap across nine orthogonal cofaces

The saved restricted dual violates a genuine cross-coface incidence constraint by a factor of **9**. This review independently reconstructs the first witness in `Q(sqrt(5))`; it does not recompute the broader cap family or repair the dual.

Generation and a fresh saved replay pass in [the review script](../scripts/review_sqrt5_cross_coface_full_base_witness_20261017.py); [the exact certificate](review_sqrt5_cross_coface_full_base_witness_20261017.json) binds its input hashes. Arithmetic uses only standard-library `Fraction`, with no producer arithmetic imports.

## Literal identification

The containing maximal simplex has base type 54999, row 80143, and dimension 4. Its full fourteen-root star consists of the five-root Lorentz core 1 and nine orthogonal `A1` components. The core has norm-2 diagonal and bonds `(0,1)=3`, `(0,4)=5`, `(1,2)=3`, `(2,3)=3`. All five one-root deletions are positive definite. Exhausting all 120 permutations proves that this weighted core has trivial automorphism group.

Actual prism provenance 440, row 26259, uses canonical prism core 19 and raw source 23. Its two base records have types 54998 and 54999, so exactly one contributes to the latter type. For that base the deleted cap normal is raw root 2, the opposite base facet is raw root 0, and the remaining core identification is `1→4, 3→3, 4→1, 5→2, 6→0`. Root 0 is an isolated `A1` in the containing simplex star. The cap pairs by −1 with canonical core root 3, by a dotted product with root 0, and by zero with every other old root. Its squared dotted product is

\[
\tau=(7+\sqrt5)/2.
\]

These identifications are checked against the frozen [raw local cores](fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017.json), [base incidence data](fresh_sqrt5_prism_simplex_base_incidence_20261017.json), [pointed profiles](fresh_sqrt5_pointed_simplex_base_capacity_20261017.json), and [actual base map](fresh_sqrt5_labelled_coface_base_packing_20261017.json). Local realization and the maximal-prism/full-star interpretation are inherited from those certificates.

## Exact geometric incompatibility

Let `G` be the core Gram and let `p=−e_3`. Direct independent inversion gives

\[
\det G=4-2\sqrt5,\qquad v=G^{-1}p,\qquad
v^2=p^TG^{-1}p=(1-\sqrt5)/4.
\]

The positive definite four-root deletions and negative determinant give signature `(4,1)`. Thus the full core-plus-nine-`A1` star is nondegenerate and spans the ambient fourteen-dimensional Lorentz space. Contacts with that star determine a cap normal uniquely. For its possible occurrence at the i-th isolated norm-2 root `q_i`, the normal must be

\[
r_i=v-\frac{\sqrt\tau}{2}q_i,\qquad r_i^2=v^2+\tau/2=2.
\]

Two distinct isolated roots are orthogonal, so for different copies

\[
(r_i,r_j)=v^2=(1-\sqrt5)/4\in(-1,0).
\]

This is forbidden for distinct norm-2 walls of an acute Coxeter chamber: their product is zero, at most −1 for a nonright finite angle, or at most −2 for a nonintersecting pair. Repeated use of the same isolated root gives the same cap normal, hence cannot give another facet. The trivial core automorphism eliminates alternative core placements; permutations of the nine identical `A1` components are exactly the alternatives tested above. Consequently **at most one cap of this exact profile can occur across all nine cofaces**.

Counting the unique type-54999 base of each prism of this provenance yields the valid selected-type inequality

\[
f_{26259}\le f_{80143}.
\]

Omitting other prism types only weakens the positive side. No independent choice of provenance has been made: this is the actual expanded row and its literal saved base.

## Saved dual violation and scope

The unchanged masses in the [five-port diagnostic](fresh_sqrt5_joint_five_port_diagnostic_20261017.json) are

\[
\alpha_{26259}=7551/10^{10},\qquad
\alpha_{80143}=839/10^{10}.
\]

The new zero-budget signed column has load `alpha[26259]−alpha[80143]=839/1250000000>0`, exactly nine times the allowed simplex mass. This invalidates that fixed dual witness once the new incidence constraint is included. It establishes neither the outcome of a repaired dual nor a stronger dimension bound.

This is an independent arithmetic and mapping check of one witness, together with the geometric uniqueness/product argument. It is not a duplicate census of all cross-coface profiles, an external human review, or a formal proof-assistant verification.
