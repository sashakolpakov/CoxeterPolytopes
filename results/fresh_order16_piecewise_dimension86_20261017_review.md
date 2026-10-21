# Independent scalar and composition review of the dimension86 argument

**FINAL/FROZEN: independent hash-bound generation and fresh saved replay both PASS.** The [independent checker](../scripts/fresh_order16_piecewise_dimension86_20261017_review.py) uses only standard-library integer/Fraction arithmetic. It imports neither the producer nor any project arithmetic module. A complete written/source composition read of the [dimension86 argument](fresh_order16_piecewise_dimension86_20261017.md) finds no identified gap, with its inherited classification and local payment proofs retained as explicit inputs.

This review independently reconstructs the piecewise scalar comparison. It does not claim another260-case induced-embedding census or another reconstruction of the nine exceptional glued cores.

## Independent enumeration and exact curvature

The envelope under review is

    b(r)=r/2 for1<=r<=14; b(15)=7; b(16)=127/17;
    b(r)=217/29 forr>=17.

For each link dimension d=86,...,90, the checker chooses the optional high part in{0,18,...,d}, the optional low part in{0,...,13}, and whether a15 is present. It then chooses the multiplicities of17 and16; the remaining rank must be a nonnegative multiple of14, forcing its multiplicity. Zero parts are omitted and sorted shapes are deduplicated. This is a different enumeration from the producer's loops over the three repeated-rank counts, followed by low/middle choices and a forced high part. Both give exactly the same shapes.

Sliding-window convolution reconstructs

    H(t)=product_i(1+t+...+t^r_i)=sum_j h_j t^j,
    V=product_i(r_i+1),  G=sum_j h_j min(j,d−j).

The checker verifies polynomial length, symmetry and H(1)=V. It independently expands product_i(1−t^(r_i+1))=sum_l a_l t^l, combining repeated shifts exactly. For s factors and q=floor(d/2), a cumulative-tail/hockey-stick identity gives

    G=qV−2 sum_l a_l binom(q−l+s,s+1).

The coefficient and cumulative computations agree in all3,355 profiles and the boundary profile. It also checks the absolute-deviation form

    kappa/V=(5G−dV)/(2V)
           =3d/4−(5/2)E|sum_i U_i−d/2|.

All polynomial coefficients, G, curvature, budget, margin, dimension minima and boundary values match the primary saved JSON exactly.

| Ambient dimension | Profiles | Minimum curvature minus cost | Finite slack beta_n−9n |
|---:|---:|---:|---:|
|87|621|146417/911250|197/2|
|88|646|203099/607500|110|
|89|671|136309/253125|123|
|90|696|34528/50625|135|
|91|721|606679/708750|297/2|

The minimum is146417/911250 at ranks(2,14,14,14,14,14,14). The ambient86 boundary profile(1,14,14,14,14,14,14) has margin−112409/1518750. This formal link shape is a limitation of the particular coefficient comparison, not a realized finite-volume chamber or a sharpness assertion.

## Why the reduced shapes cover every partition

The inherited [convex-order theorem](fresh_distance14_unrestricted_dimension156_20261017.md) applies to centered independent uniform factors: for1<=a<=b, transfer(a,b)->(a−1,b+1) increases dispersion and can only decrease normalized curvature, with the other factors fixed.

Only three types of transfer are used. Two ranks below14 are concentrated until one vanishes or the other becomes14; their linear budget is preserved. Two ranks15 are replaced by14 and16; their budget increases by8/17. Two ranks greater than17 are concentrated until one becomes17; their constant budget is preserved. In all three cases curvature minus budget cannot increase. In particular, the second step requires no false assumption of cost preservation: the increase in cost strengthens the desired reduction.

Every nontrivial transfer increases the integer sum of squared ranks at fixed total d, so the process terminates. It leaves arbitrary copies of14,16,17 and at most one part in each exceptional class:1,...,13;15; and18 or larger. Transfers do not create additional troublesome classes: low-rank moves stop at14, high-rank moves stop at17, and the doubled15 move produces only14 and16. Thus every integer partition reduces to a checked shape. Positivity at the reduced shape implies positivity at the original one.

## Local payment composition and budgets

The [order16 payment theorem](fresh_order16_finite_cross_payments_20261017.md) binds the [complete260-domain preferred-component embeddings](fresh_order16_preferred_embeddings_20261017.md) and the [universal antipodal-cycle half-coupon theorem](fresh_cycle16_half_coupon_transfer_20261017.md). Its finite coupons refer to whole components. They do not infer an A1/A15 or A1/D15 payment merely from those subdiagrams occurring within larger components.

The no-short9-path reduction remains valid: finite corners pay a short path directly; at ideal corners the preceding replacement theorem handles larger factors, and the new cycle lemma handles the sole unpaid antipodal order16 A pair. The complete inherited preferred-component classification can therefore be applied with unchanged hypotheses. All induced images in order16 factors are included, and the nine dropped D cases are treated separately.

For those D cases, the actual affine kernel forces every extra supporting white orthogonal to the full16-root affine factor. Contacts to both remaining black roots3,13 give a paid path2. A single nonzero contact has squared Lorentz projection at least4 by the inverse block4[[1,−1],[−1,1]], contradicting norm2 because the orthogonal complement is positive definite. Thus an otherwise unpaid star is sealed, preserving its whole finite components and the required coupons. The cycle lemma uses the arbitrary-spectator determinant equality to obtain the analogous sealing for A1/A15-center coupons. These are universal support arguments, rather than observations about one small example.

For finite cost, let a,b count A1,A2, let t count D15,D16,A15, and v count D13,D14. The source's upper bound is

    9n+a(3t/2−9)+b(3v−17)−24t−35v,
    a+2b+15t+13v<=n.

For n<=108, t<=7 and v<=8. If t<=6 and v<=5, all surplus terms are nonpositive. If t=7, then v=0 and a<=3, so its possible9/2 surplus is dominated by168. For v=6,7,8 the respective bounds b<=15,8,2 give positive contributions at most15,32,14, dominated by210,245,280; the a term is nonpositive. This exhausts the domain. The A15 coupon was safely overestimated in the displayed expression. Other finite components retain their inherited9-rank budget.

The order16 factor costs112,107,103,111 give normalized bound7. Order17 retains127/17. For larger orders the inherited all-rank spine/fork formulas have finite breakpoints and a final slope7 branch. Their maximal normalized cost is217/29, attained by D at order29. Smaller factors use all pairs and give r/2. This justifies the piecewise envelope for every rank; the finite BFS sample is not the tail proof.

At nonsimple cusps the [actual factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md) multiplies each factor's pair count by V/(r+1), giving total ideal cost V sum_i b(r_i). No minimum-rank or factor-count restriction has entered.

## Dimension conclusion and verification scope

The inherited [shape inequality](fresh_cusp_link_shape_curvature_20261017.md) gives

    D5 >= beta_n*v_f + sum_cusps kappa(L) + 5(n−1).

Actual payments give the corresponding upper bound, with finite cost9n and the stated ideal envelope. The independent margins and finite slacks are strictly positive for ambient87,...,91, so the bounds contradict one another. The completed [dimension91 theorem](fresh_ideal7_cross_dimension91_20261017.md) excludes every larger dimension. The resulting unrestricted ceiling is86; the separate compact-facet ambient ceiling79 is unchanged.

Independent reconstruction and complete comparison with all3,355 primary rows pass. Final hash-bound generation and a separate fresh saved replay both terminated exit0/PASS, including all input hashes and every primary comparison. A separate frontier mathematical/source review of the full primary composition also passed; it does not duplicate this scalar reconstruction. The root independently read this checker and found no issue in its enumeration or identities. This written/source review covers the full composition, concentration, finite budget, ideal envelope, actual cusp multiplicity and inherited tail. The260-domain embedding census, local core computations and complete historical classification remain explicit inputs rather than duplicated computations. External refereeing, publication priority and proof-assistant formalization are not asserted.
