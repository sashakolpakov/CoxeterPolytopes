# Twelve walls, an exact escape, and finite endpoint completion

20 October 2026. Written algebra only; no scientific job, import, root
census, Lean, or build. The exact lattice and first positive-height batch
are proved in `compact_h9_sqrt3_e8_integral_splitting_20261020.md`
(SHA256 `0e8f62f9354fb3be65738db2bb6b8ffcbfdd864d535a70c10292f71d4210bacb`)
and `compact_h9_sqrt3_e8_first_vinberg_batch_20261020.md`
(SHA256 `9d71d939945ca043851ffe53b7eea64c7f7abd53067b94ece6f600c537c1b36b`).

Use s=√3, ε=2+s, η=1+s, λ=s−1 and
M=E8(O)⊥Ov, v²=−λ. The E8 numbering is the A7 path α₁,…,α₇ with
β attached at α₃. Its fundamental weights are ω₁,…,ω₈.

## Two more compatible roots and a genuine prefix escape

In addition to the first walls

    r₁=ηεv−εω₆,                 r₂=ηεv−(3+2s)ω₇,

take

    r₃=2εv−ηω₁−ω₇,             r₄=(9+5s)v−εω₁−(3+2s)ω₇.

Their squares are respectively 2ε,2ε,2,2ε. The nonzero products with
E8 simple roots are

    r₁·α₆=−ε,
    r₂·α₇=−(3+2s),
    r₃·α₁=−η,         r₃·α₇=−1,
    r₄·α₁=−ε,         r₄·α₇=−(3+2s).

All six pair products among the four new roots are

    r₁r₂=−ε, r₁r₃=−ε, r₁r₄=−ε,
    r₂r₃=−2ε, r₂r₄=0, r₃r₄=0.

Thus every pair is orthogonal, has label 3 or 12, or is ultraparallel.
All reflections preserve M, since their norms are 2 times units. The
timelike E8 vertex has strict inequalities at these four walls; a small
inward perturbation gives a strict interior point for the twelve-wall
polytope. The acute projection lemma then makes every listed wall an
actual facet. No completeness beyond the proved first batch is assumed.

Nevertheless this twelve-wall Coxeter polytope has infinite volume.
The integral vector

    X=(5+2s)v−2ω₂

has square 33−17s>0. Its E8 products are −2 at α₂ and zero elsewhere;
its products with r₁,r₂,r₃,r₄ are, in order,

    −2s,       −8−2s,       0,       −2η.

Thus it lies in the twelve-wall cone. Adding a sufficiently small positive
multiple of a strict future timelike interior vector makes all inequalities
strict while retaining spacelikeness. The resulting segment toward a
timelike interior point meets the future null cone in a strict point,
giving an open ideal patch and hence infinite hyperbolic volume.
This is a prefix escape, not an escape from the FULL reflection chamber.

## The complete finite endpoint possibilities on this edge

The escaped edge releases α₂ from the E8 vertex. Its retained white
diagram is A1+A6: the isolated A1 is α₁, and the ordered A6 nodes are

    β, α₃, α₄, α₅, α₆, α₇.                          (1)

Let β₀≥0 be the raw contact with the isolated A1. If there is an A6
contact, let γ≥0 be its raw value at position j in (1). Finite spherical
classification and the norm pool {2,2ε} leave exactly eighteen necessary
profiles, counting the two choices of endpoint position separately:

| Contact pattern | q | β₀ | γ and positions j | Number |
|---|---:|---:|---|---:|
| None |2 or2ε|0|0|2|
| A1 only |2|1 or s|0|2|
| A1 only |2ε|η or ε|0|2|
| A6 only |2|0|1, any j=1,…,6|6|
| A6 only |2ε|0|η, j=1 or6|2|
| Bridge |2|1|1, j=1,2,5,6|4|

An A6 single attachment gives A7, D7 or E7 at the respective path
positions; a 4-bond gives B7 only at an end. A bridge gives A8 at j=1,6
and E8 at j=2,5. Positions 3,4 are not finite. A putative B8 bridge
would require a 3-bond between different norm-unit classes, whose raw
product is not in K. Multiple A6 contacts create a cycle and cannot be
finite. The isolated A1 permits labels 3,6 for q=2 and 4,12 for q=2ε.
This inventory is necessary at every height, not merely in a search box.

Write a candidate root

    r=a v+Tω₂−β₀ω₁−γω_j,       a=ηA.

Here ω_j in this formula means the E8 fundamental weight of the actual
node at position j in (1). Omit its term when γ=0. Set

    R=7β₀+R_jγ,       (R₁,…,R₆)=(10,20,16,12,8,4),
    Z=14T−R.

Completing the square gives the exact necessary equation

    Z²−28ηA²=K,
    K=14q−7β₀²−2j(7−j)γ²,       Z≡−R mod14O.          (2)

The term involving j is zero if γ=0. For A6-only short profiles the
values of K are 16,8,4 at distances 1,2,3 from a chain end; the long
end profiles give 4ε. The bridges give 9 or1 at an end or next-to-end.
All are positive. In this sign convention an outgoing finite endpoint
has a>0 and Z<0: its edge parameter is −λa/Z, and (2) puts it strictly
before the null point. A convention using r=a v−hω₂−white instead
uses the opposite completed square Z=14h+R>0.

## A uniform prime-seven exclusion, but not a full obstruction

All SIX profiles without an A6 contact are impossible at every height.
For γ=0, equation (2) reduces directly to

    7(2T−β₀)²−4ηA²=2q−β₀².                          (3)

The ring O/(7) is F49 because 3 is nonsquare modulo7. Reducing (3)
forces −(2q−β₀²)/(4η) to be a square in F49. An element of F49 is
a square exactly when its norm to F7 is a square. For

    (q,β₀)=(2,0),(2ε,0),(2,1),(2,s),(2ε,η),(2ε,ε)

the six norms modulo7 are respectively

    3,3,3,5,6,5.

All are nonsquares in F7. The facts used here are N(η)=−2,
N(ε)=1, and N(4)=16. Thus twelve A6-touching profiles remain.
It would be incorrect to regard this reduction as a nonreflectivity proof.

## In fact finite endpoint mirrors exist at unbounded height

The bridge β₀=γ=1, j=5 has R=15 and K=1. In the positive completed
square convention put

    r=a v−hω₂−ω₁−ω₆,
    Z=14h+15,          Z²−14λa²=1.                   (4)

There are arbitrarily large integral solutions suitable for (4). This
uses the standard Dirichlet unit theorem, explicitly, as an existence
input: the quartic field F=K(√(14λ)) has signature (2,1), so its unit
rank is 2, while K has unit rank 1. The relative norm-one subgroup
therefore has positive rank. After taking a power it contains a
nontrivial unit of the order O[√(14λ)]. Further powering puts it in
the identity of the finite unit quotient modulo 14η. Choose its
physical embedding greater than1 and, if necessary, invert it so that
its √(14λ) coefficient is positive. Its powers yield

    Z>0, a>0, Z→∞, Z≡1 mod14η, a≡0 mod14η.

Then h=(Z−15)/14 is integral and eventually positive. Equation (4)
gives r²=2. It is a primitive integral reflector; its retained white
extension is E8 and its released α₂ contact is −h≤0.

These mirrors eventually coexist with ALL twelve prefix walls. Indeed
h/a→√(λ/14), so

    r/a → Y=v−√(λ/14)ω₂.

The vector Y is the future null point of the outgoing edge. It has
strictly negative products with every nonwhite prefix wall: all four
extra walls meet that edge only after Y, as witnessed by the strict
prefix escape above. The retained white products of r are fixed at
0 or−1. Thus every remaining product is strictly negative for large
enough powers. The same limit preserves any prescribed finite collection
of future timelike points, including previously certified finite corners.

This is an existence theorem, not a small-height answer or a compactness
proof. It rules out the proposed strategy of excluding every finite
endpoint on this particular edge. A constructive search should seek the
smallest compatible solution or a different geometric mechanism, rather
than an impossible all-height exclusion of all twelve remaining profiles.
