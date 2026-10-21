# Every residual E4 cutter must be dotted to wall 144

**Any real spacelike wall normal that is right to the retained sixteen-root
support, acute to the 145-wall chamber, and has nonnegative E4 pairing is
strictly dotted to wall 144.** This covers all finite, parallel and dotted
choices of contacts to the other three residual caps, with no angle cutoff.
It does not exclude a cutter that is dotted to wall 144, or one that leaves
the sixteen-root orthogonal complement.

The [exact verifier](../scripts/fresh_h20_residual_cutter_domain_20261016.py)
and [certificate](fresh_h20_residual_cutter_domain_20261016.json) use the
[independently reconstructed 145-wall geometry](fresh_h20_residual_h4_corner_20261016_review.md).
Only three inequalities and three rational multipliers are needed. A
numerical polytope scout suggested the certificate but is not a proof premise.

Write n_0,...,n_3 for walls 141,...,144, let q_i=n_i^2>0, and let C be their
common finite corner with d=-C^2>0. These five vectors form an orthogonal
basis of the Lorentz five-space perpendicular to the support. Acuteness to
the n_i makes the four coefficients nonpositive. Acuteness to old wall 133,
which has negative product with each basis vector, forces the coefficient
of C to be positive for a spacelike candidate. After a positive rescaling,

    r = C - sum u_i*n_i,    u_i >= 0,
    z_i = u_i*sqrt(q_i/d) >= 0.

Old walls 133 and 140 give A_133.z<=1 and A_140.z<=1, where

    A_j,i = (n_i.n_j)/(C.n_j) * sqrt(d/q_i).

Nonnegative E4 pairing gives B.z>=1, with the same formula replacing n_j by
E4. All denominators are strictly negative. Take the positive multipliers

    529/100, 268/25, 759/50.

The independent rational intervals prove that

    (529/100)*A_133 + (268/25)*A_140 - (759/50)*B
        > (1,1,1,0) coordinatewise.

The approximate coefficients, for readability only, are
(1.00236676,2.32599748,1.00156003,0.00713549). The exact right side of the
combined inequality is 83/100. Therefore

    z_0+z_1+z_2 <= 83/100,
    z_0^2+z_1^2+z_2^2 <= 6889/10000 < 1.

Projection of r into n_3-perpendicular consequently has square

    -d + sum_{i<3} q_i*u_i^2 <= -(3111/10000)*d < 0.

Since r itself is spacelike and r.n_3<=0, this implies
r.n_3 < -sqrt(r^2*n_3^2): its contact with wall 144 is strictly dotted.
The proof only needs the four cap inequalities, old walls 133/140, and the
E4 halfspace. It imposes no field, norm, height, or retained-center condition.

This excludes the entire four-finite-contact cutter route, rather than just
a finite grid of labels. It identifies the required dotted contact for the
next construction search; it is not a finite-volume nonexistence theorem.

```sh
python3 scripts/fresh_h20_residual_cutter_domain_20261016.py --verify
```
