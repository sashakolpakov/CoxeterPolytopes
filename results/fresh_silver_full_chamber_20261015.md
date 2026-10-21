# Full-group nonreflectivity of the older silver F4² host

An exact integral isometry with trace **8+3√2>10** stabilizes the full
reflection chamber of

\[
L=\langle-1-\sqrt2\rangle\perp F_4\perp F_4\perp\langle1\rangle,
\qquad O=\mathbb Z[\sqrt2].
\]

Consequently this literal rank-ten lattice is **nonreflective**: its full
integral reflection chamber in H9 has infinite volume. This closes the
full-group question left open in the
[older construction report](fresh_silver_f4_h9_report_20261013.md).
It supplies no new polytope or dimension bound. Independent review of the
new symmetry certificate is pending.

The [checker](../scripts/fresh_silver_full_chamber_20261015.py) and
[certificate](fresh_silver_full_chamber_20261015.json) reconstruct the witness
from the old92-wall list. Every generation or saved replay also reruns the
complete height30 Vinberg producer in a temporary directory. All mathematical
fields agree exactly with the old archive; only the two elapsed-time fields
are excluded from comparison. No earlier file is changed.

## Complete reflector scope

Each displayed F4 block has diagonal entries2, successive off-diagonal entries
−1,−√2,−1, and determinant1. The full lattice determinant is −1−√2, a unit.
It has physical signature(9,1), and its other real embedding is positive
definite. Both assertions are checked from the literal Gram.

The earlier norm-alphabet proof applies to every primitive integral reflector.
If its squared length is q, reflection integrality and the adjugate identity
give q dividing2det(L). The ring is norm-Euclidean, and its totally positive
units are squares. The only possible primitive norm classes, modulo square
units, are therefore **1,2+√2,2**. The complete zero-height enumeration has
sector F4²+A1, and the positive-height enumeration includes all three classes,
with parallel duplicates identified. Thus this is a full-reflection-group
argument, not a conclusion about the older norm-two subgroup alone.
These coverage facts and the earlier independent review are bound as explicit
dependencies; the new full producer replay checks the complete bounded run.

The prefix contains92 actual facets through height30, over17 control levels.
The new checker verifies reflection integrality for all normals, every pair
of walls, a strict timelike center, and a relative-interior timelike point on
each of the92 facets. Its pair counts are756 right angles,148 anglesπ/3,
104 anglesπ/4,72 anglesπ/8 and3,106 dotted pairs. A finite prefix alone would
not prove nonreflectivity; the next step concerns its persistent corners.

## Two persistent complete corners

Using zero-based indices in the old wall list, take the ordered supports

\[
S=(0,1,2,5,6,10,12,56,75),\qquad
D=(0,1,2,57,17,5,35,6,26).
\]

Their literal norm-colored9×9 Gram matrices are equal and positive definite,
of type A1+B4+I2(4)+I2(4). Their inverse has nonnegative entries. The checker
solves their perpendicular lines to obtain future timelike points U,V,
verifies all92 inequalities, and verifies that exactly the specified nine
walls meet each point. The saved choices satisfy

\[
U^2=24-17\sqrt2<0,
\qquad V^2=(60-43\sqrt2)/49<0.
\]

Let H=(1,0,…,0) be the original Vinberg basepoint. For each of U,V the
certificate gives the exact decomposition

\[
H=aU+\sum_{i\in S}b_i n_i,
\qquad a>0,\quad b_i\le0,
\]

or its D,V counterpart. Here the coefficients b are the inverse spherical
Gram applied to the old wall pairings with H. Every later simple wall r has
(H,r)<0 and (n_i,r)≤0. Pairing the decomposition with r yields

\[
0>(H,r)=a(U,r)+\sum b_i(n_i,r)\ge a(U,r),
\]

so (U,r)<0, and likewise (V,r)<0. Thus the complete outward corner sectors
at U and V survive in the full reflection chamber. This argument permits
arbitrarily many later walls; it is not a height extrapolation.

## An infinite-order integral chamber symmetry

Put F=[n_S,U] and B=[n_D,(3−√2)V], using columns in the displayed order.
Their Grams agree, since U²=(3−√2)²V². The checker forms

\[
T=BF^{-1},
\]

and verifies the full literal matrix and inverse lie in M10(O),
TᵀGT=G, det(T)=1, Tn_S=n_D, and TU=(3−√2)V with positive scale. Thus T
preserves the lattice and the chosen future sheet. It normalizes the full
integral mirror arrangement and maps one complete outward corner sector to
the other, so it preserves their full reflection chamber.

Finally tr(T)=8+3√2>10. A finite-order real10×10 matrix has all eigenvalues
on the unit circle and trace of absolute value at most10, so T has infinite
order. More strongly, this trace bound for a Lorentz isometry forces a
nonunit expanding eigenvalue. Its powers give infinitely many disjoint
translates of a sufficiently small interior ball, all within the preserved
chamber. That chamber therefore has infinite volume.

## Consequences and limits

The [earlier exact carrier identification](fresh_silver_f4_identifications_20261013.json)
proves that the time-index and scalar-index sublattices have the same full
integral orthogonal group under their displayed embeddings. The older
silver-K3 carrier is integrally identified with the scalar-index lattice.
The same nonreflectivity conclusion consequently transfers along those
specific established identifications. It is not asserted for every other
commensurable lattice or every form over Q(√2).

The exploratory corner search stopped after finding this witness. It does
not supply a complete finite-vertex, ideal-vertex, chamber-automorphism or
facet census beyond the stated bounded Vinberg prefix. There is no reason
to continue searching this literal full group for a finite-volume chamber.
Other lattices and nonintegral geometric mirror constructions remain open.

```sh
python3 scripts/fresh_silver_full_chamber_20261015.py
python3 scripts/fresh_silver_full_chamber_20261015.py --verify
```
