# Complete parabolic endpoint profiles for the three E8 partner edges

The exact deletion census has **20 profiles: 6, 8, and 6** on the three rank-19 spherical supports. Exactly **five** have a positive normal square after projection perpendicular to the fixed rank-17 support T, and hence can meet the target hyperbolic 3-space. This is a contact-profile classification, not a chamber completion or a finite-volume claim.

[Generator](../scripts/fresh_h20_parabolic_endpoint_profiles_20261017.py) and [saved certificate](fresh_h20_parabolic_endpoint_profiles_20261017.json) passed generation and a fresh saved replay. The inputs are the frozen [finite endpoint supports](fresh_h20_e8_partner_endpoints_20261017.json), [exact 162-wall cluster](fresh_h20_e8_exact_cluster_20261017.md), and [168-wall cusp completion](fresh_h20_e8_f4_cusps_20261017.md). The prior [finite-endpoint exclusion](fresh_h20_all_endpoint_separation_20261017.md) deliberately did not cover these profiles.

Let W be one of the three spherical rank-19 supports, with normalized Gram G and inverse J. Every old root and the new root have square 2. Write d_i = -(r,w_i). A one-root parabolic extension has d >= 0 and d^T J d = 2. Its Gram has rank 19, with kernel (Jd,1). The connected component containing the new root must be an irreducible affine Coxeter diagram; all other components remain spherical. Deleting the new node therefore leaves a union of whole connected components of W. Conversely, each label-preserving deletion match supplies precisely such a profile. The generator exhausts the standard affine Coxeter diagrams A, B, C, D, E, F, G up to rank 19, including the parallel affine A1. It retains every labelled graph isomorphism, deduplicates only identical literal contact lists, and checks the cost, kernel, and rank exactly. This standard affine classification is the completeness input, rather than a cutoff on angle labels. No affine noncrystallographic family is omitted.

The literal old supports have types A12+A4+A2+A1, B12+A4+A2+A1, and A12+A4+A2+A1. Thus they have no noncrystallographic component obstructing a possible ideal endpoint. The five profiles that can meet T are:

| Old partner pair | New nonright contacts | Affine component | Projected normal square |
| --- | --- | --- | --- |
| 159,160 | 159 and 160, both pi/3 | affine A12 | 2 |
| 159,161 | 10 at pi/3 and 151 at pi/4 | affine C13 | 1/11 |
| 159,161 | 10 at pi/4 | affine C12 | 2/11 |
| 159,161 | 9 at pi/3 | affine B12 | 4/11 |
| 160,161 | 1 and 161, both pi/3 | affine A12 | 12/11 |

Each support also has the same five profiles of projected square zero: close the A4 using roots 12 and 15 at pi/3; add a root parallel to 151; close the A2 using 61 and 87 at pi/3; or extend the A2 to affine G2 by a pi/6 contact to either 61 or 87. These null projected normals do not define facets meeting the hyperbolic T-space. They can still be relevant to completion of an ambient cusp.

For each profile the normal has the form r = -WJd + h nu, where nu is a null vector in W-perp. The classification alone neither selects the feasible future ray nor proves a permitted value of h. Those require all old ambient inequalities and a common strict future interior. The affine relation then reads r + WJd = h nu; its positive coefficients constrain the orientation. The height, mutual new-wall Coxeter relations, and bounded full-rank cusp link are separate obligations. In particular, the five positive projected squares do not by themselves establish five compatible new walls or close an old edge.

The numerical height diagnostic following this census is outside this frozen package. No global H20 obstruction, complete ambient vertex census, or finite-volume result is asserted.
