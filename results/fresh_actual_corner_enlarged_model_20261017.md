# Exact enlarged actual-corner model: optimum still above12

One bounded HiGHS interior-point solve terminated successfully in53.52 seconds. Exact rational primal and dual certificates give

\[
\frac{338427277204734900}{27027022505771209}
\ \le T_*\le\
\frac{3130455146279}{250000000000},
\]

or

\[
12.52181135130475\ldots\le T_*\le12.521820585116.
\]

Thus this explicitly enlarged **finite partial model** still cannot attain the target12. The earlier actual-corner family reduced one fixed witness to8.80486; reoptimization moves to a different obstruction. Neither calculation is a compact arithmetic dimension theorem.

The [source](../scripts/fresh_actual_corner_enlarged_model_20261017.py) and [saved certificates](fresh_actual_corner_enlarged_model_20261017.json) preserve every coefficient of the inherited [cross-coface full model](fresh_sqrt5_cross_coface_cap_full_model_20261017.md), then add the fully exported subsequent constraints. No frozen input is modified and no duplicate matrix is written to disk.

## Preserved domain and exact additions

All82,806 inherited face obligations, 6,558 vertex-forest budgets, 428,768 old variables, and94,722 zero-budget signed transfers remain unchanged. The common injective map from1,858 actual Q provenances to face rows is checked before adding geometric capacities.

The successor adds:

- 1,454 exact cap columns not already present. The six frozen cap packages contain2,071 records and1,695 distinct complete coefficient tuples;241 were already inherited. Deduplication retains the simplex row, capacity, and every positive face coefficient. It does not merge different capacities or choose a new provenance.
- All52 [core23/25 forced-quadrilateral incidence constraints](fresh_sqrt5_core23_25_spectator_incidence_20261017.md).
- All294 [actual raw23 corner-clique inequalities](fresh_dominant_cap_compatibility_family_20261017.md), including the first dominant-vertex inequality as one member.
- All98 saved critical joint-port cuts from the [profile-clique scan](fresh_sqrt5_profile_clique_joint_ports_20261017.md).

The final finite matrix has430,666 variables, of which164,840 have a positive face entry, and807,707 face-coefficient entries. The optimizer's combined face/budget matrix has89,364 rows,164,841 columns including the objective variable, and882,929 nonzeros. Every rational positive-cost column is scaled uniformly by a positive common denominator in both its face coefficients and its vertex budget.

## Signed columns with a vertex cost

The52 forced-incidence inequalities have the form

\[
\sum_Q a_Q N_Q-\sum_S b_S N_S\le cN_Y,
\qquad c>0.
\]

Their variables have negative S entries **and a positive single-vertex cost**. They are deliberately excluded from the zero-budget transfer list. The dual condition is the ordinary signed load bound

\[
\sum_Q a_Q\alpha_Q-\sum_S b_S\alpha_S\le c\mu_Y.
\]

The final verifier explicitly recounts all52 such signed loads. Treating them as zero-budget constraints would impose the wrong inequality.

The96,176 genuine zero-budget columns retain the earlier Q-positive/S-negative structure. Exact dual recovery first raises S masses as necessary for those columns, then recomputes every positive-cost vertex load, including the52 signed ones. Increasing S masses therefore causes no unchecked vertex-cost constraint: all such loads are evaluated after the repair. The robust inherited verifier checks nonnegative masses, indices, the exact objective and normalizer, and the target flag.

Exact primal recovery rounds the numerical variables to rational values and raises only inherited pair/port repair variables. Every such coefficient is asserted strictly positive wherever present, and every face has a positive repair contribution. None of the new signed columns is in this repair set. The verifier then checks every face inequality and every vertex budget exactly, including all negative entries. Variables absent from all positive face entries may safely be fixed to zero.

## Scope of the new certificates

The new dual has4,932 positive face masses and1,724 positive vertex allowances. Of its Q rows,397 lie in the explicitly split actual-provenance domain. A separate support calculation identifies the remaining merged Q rows and reconstructs an exact restriction before assigning any common geometric features.

The complete local simplex-packing covers and the millions of joint-family inequalities checked for earlier fixed masses are **not** a complete exported constraint matrix here. This solve retains their saved finite cuts, together with all the explicit constraints listed above. Their old fixed-alpha certificates cannot be reused at the new masses. Any subsequent full local or joint-family claim needs a new feature reconstruction and separation; it is outside this LP certificate.

Likewise, the face/forest domain and ambient-field coverage remain partial. An exact target12 primal would still not establish a dimension bound without resolving those omissions. The present exact lower certificate instead shows that the exported finite constraint system itself remains above12.

## Validation

Prepared-domain generation45118, the single solve1008, and fresh saved replay89108 all terminated successfully. The fresh replay rebuilds the entire preserved matrix, matches its saved digest, and checks both rational certificates without invoking the optimizer. There was exactly one optimization call, using `highs-ipm` with a300-second limit; it returned optimal after15 iterations.

The model digest is `5ad462457e20ed3955d05caa77a5ab501619cab13ffc2ec19abc11c87cde480b`. The source SHA256 is `16c72b6d97399361a31315e5e3b7c65d06d07db2306968d2cca87f11cb9c5641`. All direct prepared-source hashes and source syntax checks pass. The saved JSON binds the exact added cap sources, positive denominator scalings, and all52 vertex-cost signed columns.

Construction-audit's narrow independent source read passes cap remapping and exact deduplication, old-matrix preservation, denominator scaling, placement of the52 signed vertex-cost columns, their explicit load recount, and the strictly positive inherited primal-repair indices. It does not duplicate the matrix build or optimizer, and the general inherited solver engine remains a separately checked dependency. These are internal computer-assisted checks, not external refereeing or proof-assistant formalization.
