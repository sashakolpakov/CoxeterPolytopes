# One finite M6 edge and its actual chamber reversal

20 October 2026. Complete finite reduction and checker contract;
execution pending. No new example or global compactness claim.

The exact whole H3 face of the new dimension-nine lattice R_5 is the
literal `M6=<-phi>+I6`, proved in
`results/golden_dyadic_d4_fourth_vertex_and_sealed_h3_face_20261020.md`,
SHA256 `d6b148954b09e259d0e0598cc4f754904ba1977b373a3f491196ef87f93f9501`.
Use its norm-`-phi` time basis V0, spatial basis `e1=-c,e2,...,e6`
(the last five are the original I5), and actual simple roots

    beta_i=-e_i+e_(i+1) (i<6), beta6=-e6.

The old actual a1 is the intrinsic M6 root

    u=(phi-1;phi,0,0,0,0,0), u^2=2.

Its actual provenance is imported from the four-edge proof, not
inferred merely from a finite Gram. Norm1/2 roots of the unimodular
M6 extend to R_5, and the H3 seal excludes other transverse contacts.

## Finite upper and all nine profiles

Release beta6 at V0, retaining the whole A5 on beta1,...,beta5 and
the whole ambient H3. The ray in M6 coordinates is `(1;t,t,t,t,t,t)`.
The primitive norm1 integral root

    h=(2phi^3;2phi^2,phi^3,phi^3,phi^3,phi^3,phi^3)

has contact -1 only with beta1, so its retained completion is B6.
It is an upper, not yet an actual facet. Put

    aH=2phi^3, ZH=12phi+7,
    ZH^2-6phi*aH^2=1.

It meets the ray at `tH=phi*aH/ZH` before light. Every first new
simple root has time a>0, norm q1 or2, and one of these profiles:

|q|A5 contact product1 at node k|finite type|R|
|---|---|---|---|
|1|none|A5+A1|6|
|2|none|A5+A1|12|
|2|1 or5|A6|7|
|2|2 or4|D6|4|
|2|3|E6|3|
|1|1 or5|B6|1|

There are nine literal profiles. A5 has no multiple-contact finite
extension; q2 product phi would give a nonexistent H6, and q1 can
have only the endpoint4-bond. Sealed H3 contacts are impossible.

Let Z be the sum of the six spatial coordinates. For k0 they are
all Z/6; for k>0 they are `(y+1)^k,y^(6-k)`, with `y=(Z-k)/6`.
Thus the exact equation and glue are

    Z^2-6phi*a^2=R, Z>0, all six coordinates in O.       (1)

Comparison with the upper and the other embedding gives

    0<a, a^2<=4phi^6*R, a'^2<=phi*R/6.                (2)

These are complete first-wall bounds, not a chosen height cutoff.
For every R in the table, write `a=A+Bphi`. Since R<=12, (2) gives
`a<30`, `|a'|<9/5`, and consequently `0<=B<=13`, `-1<=A<=9`.
For B>=14, `a>14sqrt(5)-9/5>59/2`; the sharper physical bound is
`a<=2phi^3 sqrt(12)<59/2`, which is already below that expression.
The lower radical comparison is exact since `5*140^2-313^2=31>0`.
For B<=-1, positivity and the conjugate bound are incompatible.
The A bounds follow from `a'=A+Bphi'`, `|a'|<9/5` and phi<81/50.
The checker applies the EXACT bounds(2) inside this conservative
154-pair box, separately for all six R values.

For each field square in(1), both square roots are tested, followed
by positive Z, full O-integral reconstruction, norm and all old
beta_i/u acuteness. Unknown accepted candidates abort. The exact
accepted interface consists of TWO vectors at the SAME upper:

    h in the q1,k1 B6 row,
    g=2h+beta1
     =(4phi^3;(4phi+3)^2,(4phi+2)^4),
       in the q2,k2 D6 row.

Here `g^2=2`, its time is2aH and its spatial sum is2ZH, so it
cannot close earlier. More importantly, g is not a simple wall:
the full chamber lies in the oriented h halfspace since h is an
integral reflector, even before its facet status is known. Together
with the actual beta1 halfspace, `g=2h+beta1` gives g<=0, with
equality only where BOTH independent walls h,beta1 vanish. Thus g
can only meet a codimension-two stratum, not be an opposite facet.
The complete first-wall inventory consequently leaves h as the
unique new SIMPLE wall. This is a full-germ saturation step, not
a deletion of a failed candidate by convention.

The spatial unit coordinate phi^3 proves h primitive. Its six-root
B6 Gram has determinant1, which proves saturation; its endpoint is

    H=(3phi-1;2,2,2,2,2,2), H^2=-phi^-7.

The initial exact generation failed closed at the incorrect expectation
of one accepted vector; a bounded diagnostic found precisely h and
this nonsimple g, with no earlier root. It wrote no certificate.
The corrected source explicitly expects and verifies both vectors,
their tie and their positive-root identity. Firstness is promoted only
after corrected generation and a separate replay pass.

## Integral reversal, conditional on that finite result

Define additionally

    w=(phi^2;phi,phi,phi,phi,phi,0), w^2=2.

There is an exact Gram-preserving map A of M6 with

    beta_i -> beta_(6-i) (i1,...,5), beta6 -> h,
    h -> beta6, u <-> w.

The seven vectors `(u,beta1,...,beta6)` form an O-basis: its
determinant is the unit phi-1. Their proposed images are integral
with the same full Gram, so A is integral, and the fixed identities
give A^2=Id. In particular

    A(V0)=phi^4 H=(12phi+7;(6phi+4)^6), tr A=1.

After the actual finite firstness proof, A sends the full positive
B6+H3 germ at V0 to the actual full germ at H, fixing H3 pointwise.
It therefore preserves the same full ambient chamber. This proves
w=A(u) an ACTUAL facet and also proves the whole native list below
actual, without importing a full-arrangement claim from a finite Gram.
The map is an involution; no infinite period follows.

The nine roots `beta1,...,beta6,u,w,h` have a compact Coxeter chamber
as a separate known-example input. Their reflection double across h
matches exactly the pinned P6_10 matrix in
`results/p6_10_bugaenko_double_lineage_and_classification_audit.md`,
SHA256 `651176b56e7d0d0a5e7596b1c889e187f0f5a9703b8336c456aa2e0b445d56c3`.
In its numbering the doubled roots are
`(w,beta5,beta4,beta3,beta2,beta1,s_h(beta6),beta6,s_h(beta1),u)`.
For example the new dotted weights are `sqrt(2)*phi^3` and
`2phi^6-1`, exactly the recorded values. This input proves compactness
of that finite chamber; no disputed uniqueness theorem is imported.
Once all nine roots are proved actual, every vertex of this compact
finite chamber is feasible by the actual positive acute projection.
Its convex hull then proves that it is the FULL M6 chamber.

## Checker scope and execution

The new standalone worker uses only Fraction arithmetic and the pinned
pure quadratic-field helper. It checks the complete box, all profiles,
field-square roots, integral reconstruction, old-wall contacts, the
two-vector accepted set and unique simple primitive h, saturated
endpoint, full reversal matrix/Gram/integrality,
and both fixed doubling identities. Sources are streamed-hashed before
and after. It owns no lock; root runs it with the existing parent-lock
59-second/1-GiB wrapper, a45-second child alarm, and exclusive new
certificate and ledger paths. Generation and a separate replay are
required. All geometry, finite-class completeness and published compact
P6 lineage remain explicit written inputs, not consequences of hashes.
