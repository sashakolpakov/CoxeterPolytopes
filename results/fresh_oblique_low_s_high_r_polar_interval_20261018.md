# Remaining high R/M and R/S polar obstructions

No system in the complete frozen low-S domain can have a high finite
R/M or R/S bond on either of the two fixed parent chains, for either
high-N family and any finite N/X label at least6. This closes the
remaining high-R cases left by the
[physical low-S/low-R proof](fresh_oblique_low_s_continuous_n_interval_20261018.md).
It gives a necessary obstruction for this parent/incidence family, not
a new chamber or a finite-volume example.

Files: [producer](../scripts/fresh_oblique_low_s_high_r_polar_interval_20261018.py),
[certificate](fresh_oblique_low_s_high_r_polar_interval_20261018.json),
[standalone checker](../scripts/check_oblique_low_s_high_r_polar_interval_saved_20261018.py).

## Fixed physical domain and finite duties

The parents, their four selected origins and the exact physical definitions
of J67,L68,X69,M71 are unchanged. Ancestral R70 is removed after the five
parent extensions are reconstructed. Old63/64 remain recipe witnesses;
they are not retained facets. N73 has profile0 and frame10/J/X labels
(2,2,m) or(3,2,m), m>=6. Write n=2cos(pi/m). All proofs cover the larger
closed interval [sqrt3,2], including its upper endpoint.

S72 has one of the seven inherited white profiles and frame52/L/M labels
in2–5, with actual finite duties52/L/S and L/M/S. No52/M corner is
assumed. The preceding physical proof covers all84 recipes and168 norm
sheets on the two parents. It strictly excludes92 sheets using retained
native walls, leaving76 overdomain candidates. Those earlier exclusions
are independent of R and N, so they remain valid here. The exact
finite-duty catalogue and complete seven-profile premise are inherited;
no floating positive-definiteness test selects the present domain.

The finite X/N/R duty forces RX=RN=right. To see the general high-bond
constraint directly, take a norm2 bond of magnitude h>=sqrt3. Adding
any nonright native contact of magnitude t>=1 gives a three-root
principal determinant at most

\[
8-2(h^2+t^2)\le0.
\]

The third mutual product is nonpositive, and any additional negative
product only decreases this determinant. Strict positive definiteness
therefore forces every other root in that actual finite duty to be
orthogonal to both endpoints of the high bond. This includes the white
roots, not just the three displayed black walls.

Consequently:

* High R/M in X/M/R forces R's profile0 and RX=right. The M/R/S duty
  also forces SM=SR=right. S's profile need not be zero.
* High R/S in M/R/S forces both R and S profiles0 and RM=SM=right.
  The high N/X duty separately supplies RX=RN=right.

Exactly24 of the76 low-S candidates have SM=right, twelve on each
parent. Exactly14 of these have profile0, seven on each parent. The
producer and checker select these complete lists directly from the
frozen physical certificate. There is no new label bank or numerical
candidate selection. Simultaneously high RM and RS is already forbidden
by the same finite three-root determinant.

## Continuous physical N reconstruction

The saved refined original Gram and five exact parent extensions are
restored from the physical predecessor. No rounded new product is
declared native. Each selected S uses its saved physical extension,
including its actual six-white projection and norm. These are enclosing
physical intervals, whose provenance and common-witness sheet
completeness were established by the predecessor's independent checker.

For frame10/J/X, N has covector(-a,0,-n), a=0 or1, and norm2. The
old64 residual has strictly positive square; the N norm remainder is
strictly positive in every box. The retained N sheet is reconstructed
from the resulting square root. The opposite N sheet has strictly
positive product with retained old8 in every box and is excluded here
again. No singular or zero-remainder case is silently omitted.

Sixteen directed rational boxes cover all four continuous parent/family
domains, with exactly matching endpoints. All calculations use outward
rounding at scale10^70. The first screen used numerical evaluations only
to choose these polar tests; none of those samples is a proof input.

## High R/M: the full three-frame adjugate

Let G be the projected Gram of X,N,S, and v the projected covector of M.
Every box has a positive diagonal and det G<0, hence signature(2,1).
The orthogonal complement in the four-dimensional projected Lorentz
space is a positive line. R has profile0 and norm2, so its two possible
normals have the same squared M product:

\[
(R\cdot M)^2=2\bigl(M^2-v^TG^{-1}v\bigr).
\]

Here M^2=2 and v=(0,M\cdot N,0), because XM and SM are defining right
duties. S may have a nonzero white profile; its projected norm and all
Gram entries retain that contribution. R's zero profile means its
projected products with every wall are its full physical products.

The certificate evaluates

\[
H=(3-2M^2)(-\det G)-2v^T\operatorname{adj}(G)v
  =(-\det G)\bigl(3-(R\cdot M)^2\bigr).
\]

For the middle covector, only the middle cofactor is needed. All192
box/candidate proofs have H>11.924 and det G<0. Thus both R signs have
(R\cdot M)^2<3, excluding every finite label at least6. The argument
uses the full three-frame determinant, with no division by4−n^2, so it
also covers n=2. The separately saved interval for RM squared is
diagnostic; the strict adjugate margin is the rejection certificate.

## High R/S: both signs of a different polar line

Now S and R have profile0, and R is orthogonal to X,M,N. The Gram of
X,M,N is nonsingular Lorentzian on every box. This polar line is
independent of the selected S. Project old8 orthogonally to it, writing
its positive squared residual as q. The two norm2 R roots have old8
products plus or minus sqrt(2q). The positive sign is excluded by a
strict old8 product greater than3.0874.

For the other sign, let w and s be the frame covectors of old8 and S.
The exact physical product is

\[
R\cdot S=-\sqrt{2/q}\,(8\cdot S-w^TG^{-1}s).
\]

All112 box/candidate proofs give either R\cdot S>0.0558 or
R\cdot S<−2.0298. The first violates the retained-wall nonpositive
condition. The second is dotted and cannot occur at the required finite
M/R/S corner. In particular it cannot equal a finite high-label product
in(−2,−sqrt3]. No native label cutoff or inverse-cosine rounding is used.

| Complete domain | Count |
|---|---:|
| Inherited low-S recipes / sheets |84 /168|
| Inherited strict exclusions / overdomain survivors |92 /76|
| Right-SM candidates / profile0 subset |24 /14|
| High-RM / high-RS continuous families |48 /28|
| Covering N boxes / attempted subdivisions |16 /28|
| High-RM / high-RS box proofs |192 /112|
| Unresolved / accepted systems / new cones |0 /0 /0|

## Consequence, debt and verification scope

Together with the frozen physical low-S/low-R proof and the
[all-high-S continuous-N proof](fresh_oblique_continuous_high_n_high_s_interval_20261018.md),
this excludes every finite S/R label allowed by the inherited finite-duty
and seven-profile classification on these two fixed parents, for both
specified N families and all finite high N/X labels. It does not exclude
other parents, N families or disk incidences. A next construction must
change one of those fixed geometric conditions; raising S/R label cutoffs
within this family cannot help.

The original deletion set is still30/44/63/64/65. Its exact debt remains
19 lost finite supports,12 touched compact two-faces and one compact
three-face, including three lost supports on8, from all199 original
finite supports. This packet inherits that ledger unchanged. The disk
remains V/E/F=10/21/12 with boundary8–10–N–R–S–52–8. No compact face
is claimed restored. All systems are excluded by necessary native or
finite-corner conditions before any common-center/cone/whole-face
acceptance stage, so no full ambient H9 volume conclusion is made.

Source-bound generation91185/807930 and standalone saved check69152/3cca4b
passed. The standard-library checker imports no project modules. It
independently implements directed rational arithmetic, reconstructs the
saved parent/S Gram enclosures, verifies the complete selected S lists,
all continuous N boxes, both polar arguments and every saved strict
bound. It inherits the earlier physical enclosures and92 native
exclusions; it does not repeat their original-scaffold reconstruction,
finite-duty/seven-profile proof, normal-bank enumeration or cone checks.
Separate fresh producer replay90190/d9523d also passed. A peer completed
a full producer/report/standalone-checker mathematical read, including
the immediate inherited physical restoration, extension and event-domain
interfaces. The read checked the high-bond profile/right consequences,
complete24/14 S selections, both N sheets, the three-frame Lorentz and
adjugate signs, nonzero S profiles, n=2, both R signs and the gap-free
16-box coverage. It did not execute the code, replay interval arithmetic,
reconstruct the original scaffold or normal/cone banks, or independently
prove the inherited profile classifier. The exact displayed decimal
bounds were checked against the saved rational endpoints (44db6c).
The packet is final and frozen after successful generation, fresh replay,
standalone check and this internal mathematical review.
