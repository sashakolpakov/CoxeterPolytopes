# The complete first Vinberg batch at the E8 vertex

20 October 2026. Paper proof only; no scientific job, root enumeration,
import, Lean, or build. This uses the exact integral splitting in
`compact_h9_sqrt3_e8_integral_splitting_20261020.md`, SHA256
`0e8f62f9354fb3be65738db2bb6b8ffcbfdd864d535a70c10292f71d4210bacb`.

Write s=√3, ε=2+s, η=1+s, λ=s−1. Use

    M=E8(O)⊥O v,       v²=−λ,       O=Z[s].

Choose the E8 simple roots α₁,…,α₇ along an A7 path and β attached at
α₃. The simple roots have square 2. Let ωᵢ denote their dual fundamental
weights (and ω₈ the weight for β). The weights are integral because E8
is unimodular.

**Result.** The first positive-height batch consists of exactly TWO
W(E8)-dominant walls, not one. In normalized norm 2ε they are

    P=ηε v−εω₆,
    Q=ηε v−(3+2s)ω₇.                                 (1)

Their common squared height is 5+3s. Together with the eight E8 walls
they give a genuine ten-wall initial Vinberg prefix. This does not prove
that the entire chamber has ten walls or finite volume.

## 1. Complete exclusion below the displayed height

Every primitive M reflector can be unit-normalized to square q=2 or 2ε.
Write it u+a v, with a>0. The preceding splitting note proves a∈ηO.
At the conjugate positive-definite embedding,

    |a′|²≤q′/η.                                      (2)

Its squared hyperbolic height above [v] is λa²/q. For (1) this is
λ(ηε)²/(2ε)=5+3s. Thus a candidate at or below this height satisfies
a²/q≤ε².

For q=2 set a=ηc. These bounds give

    0<c≤√ε<2,       |c′|≤√(η/2)<5/4.

Writing c=A+Bs gives |B|<13/(8s)<1, hence B=0 and c=1. But then
u²=2+λη²=2ε, which is not represented in E8(O): for u=X+sY the
constant-coefficient equation X²+3Y²=4 would force Y=0, contrary
to the required nonzero s coefficient. Thus there is no norm-2 root
in this entire height interval.

For q=2ε the bounds instead give

    0<c≤ε,       |c′|≤1/√η<1.

Now |B|<5/(2s)<2. If B=0 the second inequality is impossible for a
positive integer A. If B=−1, c>0 requires A≥2 and then c′=A+s>1.
If B=1, the two bounds leave A=1 or 2; A=1 is excluded because
η(s−1)²=2(s−1)>1. Therefore c=2+s=ε is the only possibility.
It attains the upper height exactly, and gives a=ηε and

    u²=42+24s.                                       (3)

## 2. Both E8 coefficient orbits are necessary

Write u=X+sY, with X,Y∈E8(Z). Equation (3) says

    X²+3Y²=42,       X·Y=12.

Set Y²=2k. Cauchy's inequality gives k(7−k)≥12, so k=3 or 4,
and equality holds in either case. Consequently exactly two forms occur:

    Y²=6, X=2Y,          u=εY;
    Y=2θ, X=3θ, θ²=2,   u=(3+2s)θ.                  (4)

The second line is essential. Integrality of θ follows from 2X=3Y
and coprimality of 2 and 3 in the free Z-module E8.

For completeness, the E8 inverse Cartan diagonal in our numbering is

    (4,14,30,20,12,6,2,8),

and (C⁻¹)₁₇=2, (C⁻¹)₆₇=3. These follow by adjoining β to A7:
the A7 inverse column at node 3 is
(5,10,15,12,9,6,3)/8, and the Schur complement is 1/8.
All entries of C⁻¹ are positive. Hence the only dominant integral
weight of square 2 is ω₇, and the only one of square 6 is ω₆.
Indeed the only diagonal entries ≤6 occur at nodes 1,6,7;
2ω₇ already has square 8, and ω₁+ω₇ has square 10. All other
nonnegative combinations have still larger square. Weyl reduction
of (4) therefore gives exactly the two roots in (1).

Both are primitive integral reflectors. Their squares are 2ε and
their reflection coefficients 2/q=ε⁻¹ are integral units. A common
nonunit factor could only have prime η and would leave a forbidden
represented unit norm, so cannot occur.

Their mutual product and E8 contacts are

    P·Q=−ε,
    P·α₆=−ε,             P·other E8 simple roots=0,
    Q·α₇=−(3+2s),        Q·other E8 simple roots=0.

Thus P--Q has label 3, P--α₆ has label 12, and Q--α₇ is dotted.
All first-batch walls are mutually acute. There is no omitted third
dominant orbit at this height.

The standard first-batch conclusion also has a direct metric check.
Every nonzero-height root hyperplane has distance at least the stated
height from [v]. The perpendicular foot on either (1) therefore lies
on the allowed side of every other such mirror, strictly unless it is
the same mirror. It also satisfies the E8 halfspaces by the nonpositive
contacts above. The seven E8 walls orthogonal to P form E6+A1;
those orthogonal to Q form E7. Adjoining P or Q gives a positive
definite rank-eight corner in either case. Local finiteness consequently
makes both feet actual finite vertices and both displayed mirrors actual
walls. This establishes the prefix without assuming termination of an
algorithm or completeness of any finite-height numerical list.

## 3. A smaller full-arrangement stopping test

The pair (α₆,P) has positive definite Gram matrix

    [[2,−ε],[−ε,2ε]],       determinant 1.             (5)

It therefore splits integrally from M. Its orthogonal complement N
has signature (6,1) and determinant class 1−s, and every primitive
reflection norm in N is again 2 or 2ε. The unit norm classes are
excluded already as represented norms in M, and the determinant
divisibility supplies the remaining necessary bound. Every N reflection
extends to M, because 2/q is an O-unit.

The associated codimension-two face is an actual nonempty face: the
inverse of (5) is entrywise positive, so orthogonally projecting a strict
interior point of the chamber to the pair preserves strict inequalities
at every other actual wall and preserves timelikeness. Any wall meeting
this face is orthogonal to the pair. Indeed a positive definite rank-three
finite Coxeter configuration containing I2(12) is I2(12)+A1; no finite
irreducible group of rank at least three has a 12-bond. For arbitrary
mirrors, the same statement follows by looking at their finite local
reflection group. Mirrors within the dihedral span contain the whole
face; they do not cut its relative interior. Thus the full N arrangement
is the intrinsic face arrangement. In particular compact reflectivity
of M would imply compact reflectivity of this exact N.

Here is a small integral generating model, retaining the glue explicitly.
Put

    x=α₅+2εα₆+εP,          y=α₇+2εα₆+εP,
    w=v+2ε²α₆+4εP.

The five roots α₁,…,α₄,β form D5 and lie in N. Furthermore

    x²=y²=−2η,       x·y=−2ε,
    y⊥D5,            x·α₄=−1,
    x·other D5 simple roots=0,
    ηw=4α₁+8α₂+12α₃+10α₄+8x+3y+6β.                  (6)

These eight generators and the one relation generate the ENTIRE N:
they are the projections of a full M basis to the unimodular complement
of (5). The span of D5,x,y has index ideal (η), not index one. In
particular discarding w would give the wrong integral lattice. Its
determinant is −38−22s; dividing by η² changes this to ε²(1−s), as
expected up to a unit square.

No nonreflectivity contradiction has yet been established for N. In
particular the orthogonal D5 root e and timelike y already allow the
relative-Pell family B e+A y of square 2 whenever B²−ηA²=1.
Finite-prefix negative-cone witnesses are not whole-lattice obstructions.
The useful new reduction is a definite pair of first walls followed by
an exact six-dimensional full-arrangement face, not an asserted global
completion or exclusion of the original H⁹ host.
