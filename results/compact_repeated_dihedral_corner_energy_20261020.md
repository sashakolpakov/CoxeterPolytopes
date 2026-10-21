# Repeated dihedral vertices and arbitrary mixed labels at least four

20 October 2026. Direct paper proof by root. Independent review of the
initial repeated-label theorem passed; the stronger mixed-label revision
below awaits its final source-bound review. No execution, census, TeX or
Lean job.

**Theorem.** At a farthest vertex of a compact hyperbolic Coxeter polytope:

1. The whole diagram cannot be I2(m)^r for ANY m>=2 if n=2r>=14.
2. More generally, it cannot be a sum of r arbitrary I2(m_j) components
   with EVERY m_j>=4 if n=2r>=14.

This is a pointwise vertex restriction. Other vertices are unrestricted.
It is not a uniform dimension ceiling. Mixed families containing A1 or
A2, and connected factors of rank at least three, are outside theorem2.

## 1. Actual block formulas and the energy budget

Use the covering energy of
`compact_farthest_vertex_energy_obstruction_20261020.md`. It is below one
at the chosen farthest vertex and adds over orthogonal components. A
block containing an actual edge with cosh L>11/8 has energy>3/19, because

    tanh^2(L/2)=(cosh L-1)/(cosh L+1).

For an old block put theta=pi/m, c=cos theta, d=csc^2 theta. Its unit
Gram and inverse are B=[[1,-c],[-c,1]] and C=d[[1,c],[c,1]]. For an
actual outside unit wall, let b,g be the contacts at the released and
retained node of this block, a>0 its height, and E its full spatial cost
in all OTHER blocks. Then

    a^2=d*(b^2+2*c*b*g+g^2)+E-1,
    delta=1-E-g^2>0,
    cosh L=sqrt(d)*(b+c*g)/sqrt(delta).                (1)

Each nonzero contact is a native cos(pi/k), k>=3, or exceeds one.
Contact one is forbidden by compactness, since its positive affine null
combination gives a feasible ideal point. Native contacts begin1/2 and
1/sqrt2. If both b,g are positive, (1) gives, for every m>=4,

    cosh L>=sqrt(d)*(1+c)/sqrt3
            =cot(theta/2)/sqrt3
            >=(1+sqrt2)/sqrt3>11/8.                   (2)

The last inequality is128sqrt2>171; its square is32768>29241.
Thus only single contacts in the distinguished block need attention.

A wall's positive definite replacement at another direction gives its
ACTUAL endpoint. Indeed the inverse of the retained finite Gram is
entrywise nonnegative, and projection of the old vertex along those
normals satisfies every other actual halfspace by acuteness. Equivalently
one can use the reviewed reciprocal-shadow lemma. This prevents treating
formal corners as actual, or counting one wall as two distinct walls.

When g=0 and b is finite, its mate replacement is positive definite
precisely if 1-E-b^2>0; in that case the SAME wall closes the mate and

    cosh L_mate=sqrt(d)*c*b/sqrt(1-E-b^2).             (3)

## 2. Complete off-block cost inventory for all m_j>=4

A complete spectator I2(l) block has inverse as above. Its minimum single
contact cost is csc^2(pi/l)/4. Two nonzero contacts cost at least
csc^2(pi/l)*(1+cos(pi/l))/2>1, so are forbidden in the positive definite
white set. For l=4 the ONLY allowed contact is a single3, cost1/2; the
next single contact already costs one. For l=5 the ONLY allowed contact
is a single3, cost

    k=(5+sqrt5)/10,

and the next costs>1. For l>=6 even its minimum cost is at least one.
Since each allowed cost is at least1/2, no wall can touch two spectator
blocks. Therefore the COMPLETE inventory is

    E=0, 1/2, or k.                                  (4)

The case b=g=0 cannot be an outside spacelike wall: all three values of
E are below one, so its squared height would be negative. For b,g>0,
(2) already supplies the required energy in every case.

## 3. No spectator contact: E=0

A single dotted contact has a direct edge with cosh L>sqrt(d)>=sqrt2,
which is greater than11/8. For a single FINITE contact b, both block
replacements are positive definite, so both endpoints are actual. The
mate has cosh L=cot(theta)*b/sqrt(1-b^2).

* m=4: positive height excludes labels3 and4. Hence b>=cos(pi/5),
  and cosh L_mate>=cot(pi/5).
* m=5: label3 has negative height. Hence b>=1/sqrt2, and
  cosh L_mate>=cot(pi/5).
* m=6: label3 has zero height. Hence b>=1/sqrt2, and
  cosh L_mate>=sqrt3.
* m>=8: b>=1/2 gives cosh L_mate>=cot(pi/8)/sqrt3
  =(1+sqrt2)/sqrt3>11/8.

For the first two cases,

    cot^2(pi/5)=1+2/sqrt5>(11/8)^2,

since128>57sqrt5, whose square is16384>16245. Thus each listed block
has energy>3/19.

The only exception is m=7 with a single3 contact. If its finite contact
is at least4 instead, the mate has cosh L>=cot(pi/7)>2. Its double or
dotted cases also have cosh L>2. An exceptional wall is supported ONLY
in its7 block, and has squared height

    h=csc^2(pi/7)/4-1 in (0,1/2).

There is at most ONE exceptional block: walls from two such blocks would
be distinct with disjoint spatial supports, giving product -h in(-1/2,0),
which violates the Coxeter gap. This counts blocks even though the same
exceptional wall closes both directions of its own block.

Its actual mate has cosh L=cot(pi/7)/sqrt3>109/91, hence energy>9/100.
For completeness put q=cos(2pi/7). The elementary bound
623/1000<q<5/8 follows by the increasing polynomial8q^3+4q^2-4q-1 on
[1/2,1], as in the reviewed A2 capacity proof. Then

    cot^2(pi/7)/3>(1623/1131)>(109/91)^2;

the second cross-product difference is2652. The upper bound q<5/8
also gives h<1/2, and the lower gives cot(pi/7)>2.

## 4. One spectator contact: E=1/2 or k

Keep only a single distinguished-block contact. If b=0, white positivity
forces g=1/2, because 1-E<=1/2 and the next native square is1/2. For
E=1/2 this has

    a^2=d/4-1/2,       cosh L=cot(theta).

At m=4 its height is zero and is forbidden. For m>=5 this is at least
cot(pi/5)>11/8. For E=k its direct edge instead satisfies

    cosh^2 L=(d-1)/(4*(3/4-k))>5> (11/8)^2.

Indeed d>=2, k>7/10 and 3/4-k<1/20 give the displayed bound.
The denominator is positive, since k<3/4.

Now take g=0. If E=1/2 and m=4, contact3 gives zero height and the next
native contact4 gives cosh L>=sqrt2. At m=5, contact3 gives a positive
definite mate with delta=1/4 and cosh L_mate=cot(pi/5)>11/8; any contact
at least4 gives direct cosh L>=sqrt(d)>11/8. At m>=6 even contact3 gives
cosh L>=sqrt(d/2)>=sqrt2.

For E=k and m>=5 the direct edge has

    cosh L>=sqrt(d/(4*(1-k)))
            >=sqrt(k/(1-k))=phi>11/8.

For m=4, contact3 has a positive definite mate with delta=3/4-k, whose
cosh-square is1/(4*(3/4-k))>5. Any contact at least4 has direct
cosh L>=1/sqrt(1-k)>11/8. This exhausts the inventory(4); every block
with E>0 is strong, with energy>3/19.

## 5. Sum of component energies and repeated small labels

For r>=7, if there is no exceptional7 block the energy is>3r/19>1.
If there is one, it is

    >3*(r-1)/19+9/100>=18/19+9/100=1971/1900>1.

This contradicts the farthest-vertex energy below one and proves theorem2.
The repeated cases m>=4 are included. Repeated m=2 and m=3 follow from
the reviewed A1 and A2 notes, with farthest-vertex thresholds12 and14,
respectively. Their shared linear tests are valid local tests at this
same vertex. This proves theorem1 for all m>=2.

The argument is independent of arithmetic. It excludes an infinite set
of whole vertex diagrams, including arbitrary mixtures of B2 and higher
dihedral labels. It does not exclude all finite diagrams. The general
compact29 and arithmetic compact20 ceilings are unchanged.
