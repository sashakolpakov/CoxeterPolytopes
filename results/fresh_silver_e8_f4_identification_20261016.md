# The silver E8 host duplicates the silver F4-squared host

There is an explicit integral isometry

\[
 \langle-(1+\sqrt2)\rangle\perp F_4\perp F_4\perp\langle1\rangle
 \ \cong_{\mathbf Z[\sqrt2]}\
 \langle-(1+\sqrt2)\rangle\perp E_8\perp\langle1\rangle.
\]

Here each F4 block has diagonal 2 and path products
\((-1,-\sqrt2,-1)\); this is the literal equal-norm block used in the old
search. The E8 host is therefore a duplicate construction host. Its full
integral reflection group is conjugate to the group already proved
nonreflective in the [older full-chamber certificate](fresh_silver_full_chamber_20261015.md)
and [independent review](fresh_silver_full_chamber_20261015_review.md).
The new E8 symmetry proof remains valid, but does not exclude a new integral
isometry class.

The [standalone checker](../scripts/fresh_silver_e8_f4_identification_20261016.py)
uses only pairs of rational Fractions, with no producer or project algebra
imports. Its [certificate](fresh_silver_e8_f4_identification_20261016.json)
stores both literal ten-dimensional Grams, the complete basis matrix P and
its integral inverse. Generation and saved replay pass. Exactly,

\[
 P^T G_{E_8}P=G_{F_4^2},\qquad \det P=1,
 \qquad P,P^{-1}\in M_{10}(\mathbf Z[\sqrt2]).
\]

The first column of P, with \(s=\sqrt2\), is

\[
 (17+12s,\ 48+34s,\ 96+68s,\ 144+102s,\ 116+82s,
 88+62s,\ 60+42s,\ 30+21s,\ 72+51s,\ 0)^T.
\]

It has squared norm \(-(1+s)\), and its time coordinate is positive.
Thus P preserves the future sheet. It mixes the original time direction
with spatial coordinates; comparing the positive direct summands alone
would miss this identification.

The matrix was found using vertex 51 of the
[complete height-100 E8 finite-corner census](fresh_silver_e8_h9_continuation_20261016_h100_corners.json).
Its complete outward sector has the literal F4+F4+A1 Gram. P sends the
ordered initial F4-host normals \(-e_1,\ldots,-e_9\) to E8-prefix normals

\[
 (N_1,N_0,N_{12},N_7,N_4,N_3,N_{41},N_6,N_8).
\]

It sends the initial basepoint to \((17+12s)V_{51}\), where \(V_{51}\) has
time coordinate 1. The checker verifies these identities and all 48
inequalities at that image point. The integral Gram identity itself proves
the lattice identification without requiring completeness of the corner
search or any further Vinberg enumeration.

An integral isometry conjugates every lattice-preserving reflection to a
lattice-preserving reflection in the other host. Hence it conjugates the
full reflection groups, and the older infinite-covolume result transfers
to the E8 host. This is stronger than an unspecified rational or
commensurability comparison: a literal determinant-one integral matrix is
provided. No equivalence or nonreflectivity claim for any further lattice
is implied.

```sh
python3 scripts/fresh_silver_e8_f4_identification_20261016.py
python3 scripts/fresh_silver_e8_f4_identification_20261016.py --verify
```

Both runs pass, as do recursive source manifests and report links. All prior
obstruction packages remain unchanged. This separate identification should
be consulted before treating either Gram decomposition as a new host.
