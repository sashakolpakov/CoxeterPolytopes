# A stronger third-corner bound for the heavy rank transfer

22 October 2026. Complete new manual candidate by root, pending full
outside-author review and acceptance. No mathematical worker, enumeration,
symbolic import, test, Lean or TeX job ran. The existing simultaneous
payment rule and all frozen manuscripts are unchanged.

## 1. Actual statement and inherited finite classification

Use EXACTLY the actual quadrilateral and root conventions of

`results/finite57_heavy_rank_transfer_quadrilateral_local_discharge_20261021.md`,

SHA-256
`8abe292196a6792ed350c85c0c1e66ce397f312af2a025c2502a287a601f3f53`.

Its full argument and reviews were accepted in the7ceefb66 batch. The
supporting-white forest is `W=A13+2A1+D17+D23`, of rank55. Its
black sides in cyclic order are a,b,c,d. The first corner has whole
forest2D17+D23, with a the long endpoint and b the centre of the
displayed active D17. The next wall c attaches by simple3 only to
the inert D23 long leaf. Its corner with b is D16+D17+D24.
Both incident corners of the actual fourth wall d are finite.
All actual normals have squared norm2. Write the A13 path as q1..q13
and its two old isolated leaves as f1,f2. Thus a.q1=-1 and
b.q13=b.f1=b.f2=-1.

**Proposition.** Under these same hypotheses, the third corner obeys

    C0(c,d)<=22633/64=353.640625.                  (1)

This strengthens the accepted357.125 bound. It does not give a
uniform vertex cap, a global discharging map or a dimension56 theorem.
The additional geometric input below is the actual opposite-side
indefiniteness lemma in `manuscripts/finite_volume57.tex`, source
SHA-256
`c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`.
It concerns original faces; no projected metric face is assumed Coxeter.

Here are the precise inherited finite obligations used from Sections3--4
of the accepted rank-transfer source. They are its complete finite
one-node classification and exact charge formulas, not an unexecuted
enumeration:

* If d joins either large D17 or D24 block at(c,d), including a
  possible contact with the new long leaf c itself, its six listed
  finite forests have charge at most353.640625.
* Otherwise D17+D24 remains, with charge263.5. The active
  `A13+2A1+d` has exactly the nine listed finite forms.
* The three heavy forms D16,A15+A1,B15+A1 are impossible by that
  source's full-white Schur/native-gap proof.
* If d does not contact A13, the remaining active forms are
  A13+3A1,A13+T2+A1,A13+T3. Their largest exact active charge
  is85.859375, hence total349.359375.

It therefore remains ONLY to exclude the active forms A14+2A1,
B14+2A1 and D14+2A1 when d joins neither large block. In
particular c.d=0 in these cases. A contact of d with c is NOT
discarded: it belongs to the first, already bounded alternative.

## 2. Full-white rank equation and every orientation

Put g=-a.c>=0 and h=-b.d>=0. These are products of distinct
ACTUAL norm-two walls, so

    g^2 in {0,1} union [2,infinity),
    h=0 or h>=1.                                 (2)

The55 white roots are independent and positive definite. Together
with four black roots they lie in the actual58-dimensional Lorentz
space of H57, so the resulting4-by-4 Schur matrix has determinant0.
After eliminating W its common entries are

    a^2=15/14, b^2=1/14, c^2=1,
    a.b=-1/14, a.c=-g, b.c=c.d=0.

If d contacts A13 by an endpoint3 or4, its endpoint is q1 orq13.
For next-to-end3 it is q2 orq12. No exceptional finite extension
contains A13 at this rank. In all six orientations a and d are
already connected through A13, so a nonzero a.d would close a
cycle in the actual finite corner(a,d). Hence a.d=0 BEFORE white
elimination, for every orientation. The near endpoint4 at q1 is impossible
at(a,d): q1 with leaves a,d,q2 has squared star weights
1/4,1/2,1/4, an actual affine principal. The near next-to-end3 at
q2 would make the full finite component at(a,d) a tree with arms
2,1,11. It is neither A,B,D nor any finite exceptional type.
Thus those two near orientations cannot occur.

The near endpoint3 at q1 IS a finite D15 completion at(a,d) and
must be retained; Section4 treats it explicitly. The orthogonality
used throughout follows from actual finite diagrams, not an imposed
raw product convention.

## 3. The three far orientations violate the native gap

Let B0 be the projected b,d matrix and let u=a.b, v=a.d be the
projected couplings. Eliminating the positive c coordinate leaves
a diagonal `15/14-g^2`. Since B0 below is nonsingular, determinant0
is exactly

    g^2=15/14-(u,v) B0^(-1) (u,v)^T.              (3)

All entries in the following calculations are after eliminating the
FULL W, including both inert D blocks.

**Endpoint3 at q13: A14+2A1.** Its entries are

    d^2=15/14, v=-1/14,
    B0=[[1/14,-h-13/14],[-h-13/14,15/14]],
    -det B0=Q=h^2+(13/7)h+11/14>0.

Equation(3) gives

    g^2=15/14+(3+2h)/(196Q).

The identity

    (42/11)Q-(3+2h)=(42/11)h^2+(56/11)h>=0

therefore gives `15/14<g^2<=12/11<2`, contrary to(2).

**Endpoint4 at q13: B14+2A1.** Its entries are

    d^2=1/7, v=-sqrt2/14,
    B0=[[1/14,-h-13sqrt2/14],[-h-13sqrt2/14,1/7]],
    -det B0=Q=h^2+(13sqrt2/7)h+12/7>0.

Now(3) gives

    g^2=15/14+(2+sqrt2*h)/(98Q).

Since

    (7/6)Q-(2+sqrt2*h)=(7/6)h^2+(7sqrt2/6)h>=0,

one has `15/14<g^2<=13/12<2`, again forbidden.

**Next-to-end3 at q12: D14+2A1.** The entries are

    d^2=2/7, v=-1/7,
    B0=[[1/14,-h-6/7],[-h-6/7,2/7]],
    -det B0=Q=h^2+(12/7)h+5/7>0.

Equation(3) is

    g^2=15/14+(h+1)/(49Q).

Here

    (7/5)Q-(h+1)=(7/5)h^2+(7/5)h>=0,

so `15/14<g^2<=11/10<2`, also contradicting(2).

All three exclusions cover EVERY opposite product h>=0, whether
finite, parallel or dotted. They are exact interval arguments.

## 4. The near endpoint3 is excluded separately

For d.q1=-1 the projected data instead are

    d^2=15/14, v=-13/14,
    B0=[[1/14,-h-1/14],[-h-1/14,15/14]],
    det B0=(1-2h-14h^2)/14.

If h0, B0 is positive definite. The original W+b+d would then
be positive definite and give an actual common finite point for
the OPPOSITE sides b,d, contrary to the inherited actual-face
lemma. Thus h is nonzero, and(2) forces h>=1. In particular B0
is nonsingular. Equation(3) now simplifies to

    g^2=h(15h+4)/(14h^2+2h-1).

The denominator is positive. Subtracting it from the numerator
gives `(h+1)^2>0`; subtracting the numerator from twice the
denominator gives `13h^2-2>0`. Therefore

    1<g^2<2,

which is impossible by(2). This argument does not mistake the
near orientation for the far one and does not omit its initially
allowed finite D15 corner.

## 5. Charge consequence and remaining global work

The complete active list now has no A13 attachment left except
the already bounded alternatives joining a large D block. The
unjoined active alternatives have total charge at most349.359375.
The six large-D alternatives have total charge at most353.640625.
Their maximum proves(1), with the unchanged held C0 rule and its
absence-of-order11-increments convention inherited from the source.

Numerically, one charge372 corner paired with one such third corner
has mean at most

    (372+22633/64)/2=46441/128=362.8203125<364.

Here364 is the dimension57 finite curvature coefficient
`(57-1)(57-5)/8`, not the even-rank formula. This explains why the
stronger local bound could be useful in a global allocation. It
does NOT construct such a pairing. The source still requires the
specified quadrilateral, c-role and finite endpoints; other edge
completions, ideal endpoints, higher-sided faces and receiving
multiplicities remain unresolved. No global372 improvement, dimension56
ceiling, actual example or runtime/formal verification is claimed.
