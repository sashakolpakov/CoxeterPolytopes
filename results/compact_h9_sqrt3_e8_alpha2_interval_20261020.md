# A finite complete interval on the E8 alpha-two edge

20 October 2026. Exact paper algebra and checker source only. No census,
proof job, Lean module or build was run by this author. This is one edge
of a reflection chamber, not a compactness or finite-volume certificate.

Use the integral splitting and numbering in
`compact_h9_sqrt3_e8_integral_splitting_20261020.md` and
`compact_h9_sqrt3_e8_first_vinberg_batch_20261020.md`. Put

    s=sqrt(3), epsilon=2+s, eta=1+s, lambda=s-1,
    M=E8(Z[s]) orthogonal-sum Z[s] v,       v^2=-lambda.

The E8 simple roots alpha_1,...,alpha_7 form a path, with beta attached
at alpha_3. All have square 2, and omega_i are the integral dual weights;
omega_beta=omega_8. Release alpha_2. Its seven retained walls form
A1+A6, with isolated alpha_1 and ordered path

    (beta, alpha_3, alpha_4, alpha_5, alpha_6, alpha_7).       (1)

The outgoing edge is [v-tau omega_2], with 0<=tau<sqrt(lambda/14).
The earlier twelve-wall escape and the complete eighteen-profile
inventory are in
`compact_h9_sqrt3_twelve_wall_escape_and_completion_20261020.md`.
That note proves existence of large finite upper mirrors. The small
mirror below gives a direct bounded interval, without a unit search.

## An explicit finite upper mirror

Set

    a*=23+13s, h*=4+3s,
    r*=a* v-h* omega_2-omega_1-omega_beta.

The weight products needed here are

    omega_2^2=14, omega_2 omega_1=7, omega_2 omega_beta=10,
    omega_1^2=4, omega_beta^2=8, omega_1 omega_beta=5.

Consequently

    a*^2=1036+598s,       lambda a*^2=758+438s,
    (h* omega_2+omega_1+omega_beta)^2=760+438s,
    r*^2=2.

Its E8 wall products are -h* at alpha_2, -1 at alpha_1 and beta,
and zero elsewhere. It is primitive because one of these products is
-1. Its reflection preserves M, since 2/r*^2=1. The seven retained
walls together with r* have spherical diagram A8.

Write Z*=14h*+17=73+42s. Directly,

    Z*^2=10621+6132s,       Z*^2-14lambda a*^2=9.

Thus its collision parameter

    tau*=lambda a*/Z*=(16+10s)/(73+42s)

is positive and strictly before the null point. This is an upper mirror,
not a claim that r* is the first wall of the full arrangement.

For completeness it is compatible with all four displayed caps

    r1=eta epsilon v-epsilon omega_6,
    r2=eta epsilon v-(3+2s)omega_7,
    r3=2epsilon v-eta omega_1-omega_7,
    r4=(9+5s)v-epsilon omega_1-(3+2s)omega_7.

Their products with r*, in order, are

    -14-8s,       -35-20s,       -3-2s,       -22-13s.        (2)

All are ultraparallel: their squared products exceed the corresponding
products of root norms. This compatibility is checked separately; it
will not be used as a filter in the complete first-edge census.

## Complete profiles and their equations

Every primitive reflector in M can be unit-normalized to square q=2
or 2epsilon; its v coefficient a belongs to eta Z[s]. These are inputs
already proved in the scalar-facet and integral-splitting notes. Write

    r=a v-h omega_2-U,        a>0, h>=0,
    kappa=omega_2 U, C=U^2, Z=14h+kappa,
    R=kappa^2+14(q-C).

Then

    Z^2-14lambda a^2=R,       tau=lambda a/Z.                 (3)

In the following table a path weight means the actual E8 weight at the
specified position of (1). Every position is a separate profile.

| U | q | Positions | R |
|---|---|---|---|
| 0 | 2, 2epsilon | none | 28, 28epsilon |
| omega_1 | 2 | none | 21 |
| s omega_1 | 2 | none | 7 |
| eta omega_1 | 2epsilon | none | 14epsilon |
| epsilon omega_1 | 2epsilon | none | 7 |
| path weight | 2 | 1,2,3,4,5,6 | 16,8,4,4,8,16 |
| eta times path weight | 2epsilon | 1,6 | 4epsilon,4epsilon |
| omega_1+path weight | 2 | 1,2,5,6 | 9,1,1,9 |

Finite spherical classification supplies completeness: a new node has
at most one contact with the A6 path, since a second creates a cycle.
A simple attachment gives A7, D7 or E7; a 4-bond is possible only at a
path end. If it also contacts the isolated A1, the only finite cases
are the four listed simple bridges, giving A8 or E8. A 3-bond between
the two norm-unit classes would require a raw product outside Q(s).
The isolated A1 permits precisely the listed 3/6 or 4/12 contacts.

In particular every R and its conjugate R' are positive, with

    0<R<105,       0<R'<=28.                               (4)

All eighteen profiles will be enumerated, including the six excluded
at every height by the independent mod-seven argument in the completion
note. That exclusion is not needed to shrink the finite census.

## Certified coefficient boxes

An earlier or simultaneous collision means a/Z<=a*/Z*. By (3), this
is equivalent to either of the exact inequalities

    9a^2<=R a*^2,       9Z^2<=R Z*^2.                       (5)

Using 17/10<s<26/15 gives a*<46 and Z*<146. Therefore

    0<a<158,       0<Z<500.

At the conjugate embedding (3) becomes

    Z'^2+14eta a'^2=R'.                                   (6)

It follows that |Z'|<6 and |a'|<1: indeed a'^2<=2/eta=lambda<1.
Writing Z=z0+z1 s and a=a0+a1 s, the following integer boxes therefore
contain every earlier or simultaneous finite-endpoint root:

    -2<=z0<=252,       -1<=z1<=148,
     0<=a0<=79,         0<=a1<=46.                         (7)

For example z0=(Z+Z')/2 lies strictly between -3 and 253;
z1=(Z-Z')/(2s) lies strictly between -30/17 and 2530/17.
Similarly a0 lies between -1/2 and 159/2, and a1 between
-5/17 and 795/17. The additional condition a in eta Z[s] is exactly
a0 congruent a1 modulo 2.

The source checker builds a dictionary of a^2 from its two-dimensional
box and loops over Z in its separate two-dimensional box. For each
profile it matches a^2=(Z^2-R)/(14lambda), then checks the equation,
both embeddings, (5), a/Z<=a*/Z*, and coefficientwise integrality of
h=(Z-kappa)/14 with h>=0. All signs use exact integer-pair arithmetic.
There is no floating-point cutoff or unproved height truncation.

## First-wall and simultaneous-root interpretation

No sign test against r3 or r4 is permitted to discard a candidate:
their status as facets of this prefix does not prove that each is a
facet of the full arrangement. All cap products are reported afterwards.
The minimum positive a/Z among the complete candidates is the first
endpoint parameter. The full-arrangement first wall must lie in this
finite list, because r* already cuts the edge strictly before the null
point; the fixed whites and the first wall have a finite spherical
extension. Conversely each accepted vector is an integral reflector,
so it cannot precede the first intersection with the full arrangement.

For simultaneous roots l,r, form

    d=a_l r-a_r l.

Its v coefficient and its alpha_2 simple-root coefficient vanish. In
weight coordinates its simple-root coefficients are C_E8^{-1}d.
If its other seven coefficients are nonnegative, l and the retained
white inequalities imply r. The checker reports the roots satisfying
this positive-cone test against every simultaneous candidate. It claims
a selected simple wall only when they constitute one positive projective
ray; otherwise it records UNKNOWN rather than guessing a chamber wall.

The homogeneous endpoint associated to a selected (Z,a) is

    p=Z v-lambda a omega_2,       p^2=-lambda R<0.           (8)

The checker verifies its products with every accepted root, the eight
old walls and the four prefix caps. Complete reflection-group and
spherical-inventory arguments remain written inputs; generation and
replay of one implementation are not independent geometric proofs.

Even a successful finite endpoint proves only this single edge. A later
compactness claim must prove a pointed cone with a bounded future affine
section and control ALL its extreme rays. No such global completion is
asserted here.
