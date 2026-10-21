# A single opposite prism cap is excluded by the old21-root core

No coherent acute Coxeter wall can simultaneously make all three old partner-edge endpoints finite or ideal while retaining the old core T+{158,159,160,161}. This result has **no dotted-other hypothesis and no angle-label cutoff**. It therefore tests a topology outside the previous endpoint-family no-go. It does not exclude replacing a blocking old wall or using several intermediate facets.

The [source](../scripts/fresh_h20_opposite_prism_gate_20261017.py) and [certificate](fresh_h20_opposite_prism_gate_20261017.json) reconstruct the complete intersection of the three [finite endpoint catalogues](fresh_h20_e8_partner_endpoints_20261017.md) and the three [parabolic endpoint catalogues](fresh_h20_parabolic_endpoint_profiles_20261017.md). The intersection has40 fixed finite profiles and5 parabolic profiles, plus the shared unbounded singleton151 family. Matching preserves literal roots and every bond label; no geometrically favorable contact is discarded. The actual core comes from the [exact162 cluster](fresh_h20_e8_exact_cluster_20261017.md).

Let B be the20 old roots T+{159,160,161}, and let u be the component of158 perpendicular to B. The exact core Gram gives

\[
u^2=1/30.
\]

For any prescribed products with B, all norm-two completions have the form r=p+lambda u, where p lies in span(B). The old three finite vertices are perpendicular to u. For44 of the45 fixed/affine intersection profiles, the source also finds p perpendicular to all three old vertices. The same identity holds symbolically for **every** singleton151 magnitude d: its particular vector is -d n151/2. Thus these completions vanish on old finite vertices and cannot define coherent acute added walls.

Here is the strictness argument. If a new acute normal vanishes at an old finite vertex, it belongs to the positive-definite span of that corner's20 wall normals. Its products with those normals are nonpositive, so its expansion coefficients are nonpositive because the spherical inverse Gram is entrywise nonnegative. It consequently pairs nonnegatively with every strict old interior point. This contradicts a common strict interior for the new halfspace. The certificate checks the three actual spherical corner Grams by exact positive LDL pivots and entrywise nonnegative inverses. No chosen fixed center is used to reject a normal.

Exactly one profile remains:

\[
r\perp T,159,161,\qquad (r,160)=-1.
\]

Its particular vector and two norm branches satisfy

\[
p^2=\frac{-5+\sqrt2}{23},\qquad (p,158)=0,
\]

\[
(r,158)=\pm\sqrt{\frac{51-\sqrt2}{690}}.
\]

Since1<sqrt2<2, the quantity under the last square root lies strictly between49/690 and50/690, hence between0 and1. One branch is obtuse to158. The other has a nonzero product in(-1,0), strictly between a right contact and the smallest nonright norm-two Coxeter magnitude1. Neither is Coxeter-compatible. This exact test uses only wall158; the numerical old96 obstruction is unnecessary.

The five common parabolic profiles have contacts(12,15) at pi/3,151 parallel,(61,87) at pi/3,61 at pi/6, or87 at pi/6. They are included in the44 old-vertex-zero cases. There is no additional mixed finite/ideal intersection profile. In particular, allowing ideal endpoints does not rescue the one-wall opposite-cap topology.

The calculation identifies a separate replacement experiment rather than another label search. The two surviving algebraic branches are mutually dotted:

\[
(r_+,r_-)=2p^2-2<-2.
\]

If their blocking old walls are explicitly omitted, the pair has the correct three spherical endpoint profiles to bound opposite triangular faces. Compatibility with every retained ambient wall and the whole target-face cone must then be checked anew. That replacement branch is outside this theorem and outside the frozen168-wall control.

Generation and a fresh saved replay completed successfully. Source syntax, direct source manifests, and report links were checked. This package is final/frozen.
