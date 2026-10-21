# No farthest vertex with all finite components of rank at most four in n >= 15

20 October 2026. Complete paper candidate; independent review pending.
No scientific computation, census, optimizer, build or formal job.
Arithmeticity is not assumed.

**Theorem.** Let P be a compact hyperbolic Coxeter n-polytope, n >= 15,
and z any strict interior point. A vertex p maximizing d(z,p) cannot
have every irreducible component of its finite diagram of rank at most
four. Other vertices are unrestricted.

This is a restriction on the actual farthest vertex, not a new universal
dimension ceiling. In particular, large connected A, B, D and E factors
remain outside the theorem.

The argument combines the following fixed written component proofs:

|Source|SHA256|
|---|---|
|`compact_universal_small_a_component_energy_20261020.md`|`801581c4b855289aca7982b2d0f4bb6f5849f151a78f40b86c9d403e84f415e6`|
|`compact_universal_b3_b4_d4_f4_component_energy_20261020.md`|`2eb274288127f1bd2b7df36499bec1da0fedc607d3c7ac8bd7734ef8f80c44db`|
|`compact_universal_h3_h4_component_energy_20261020.md`|`3fc546347f0488200e36dbb6407d3a66ebe4e35bc353603936799841747ff0bb`|

All paths are under `results/`. The supporting universal geometric
mechanism is the no-common-opposite theorem at
`compact_farthest_vertex_orthogonal_component_opposite_obstruction_20261020.md`,
SHA256 `4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`.
The sources credit the joint work of root, batch_guard and
resource_cleanup. The new parts here are the universal dihedral floor
with its one named defect, the short rank-four completeness argument,
and the combined n >= 15 bookkeeping. No unreviewed polytope census is
an input.

## 1. Framework and the imported component estimates

At p let B be the old unit finite Gram, C=B^-1, d_i=C_ii, and
u_i=omega_i/sqrt(d_i). Writing z=E(p-x), the farthest-vertex comparisons
give

    ||x||^2<1,       x.u_i>=tau_i=tanh(L_i/2).

The orthogonal whole components have additive energies. The fixed
proofs establish, with arbitrary other finite components,

|Component|Energy lower bound|
|---|---|
|A2|greater than 319/2000, except a possible isolated A2_7 block|
|A3, B3, H3|greater than 1/4|
|A4|greater than 2/7|
|B4, D4, F4|greater than 1/3|
|H4|greater than 11/36|

Their actual-wall tools also give these universal facts.

* Distinct opposite walls force tanh(L_i)tanh(L_j) >= u_i.u_j.
  The pair/log-concavity estimate in the fixed sources makes the two
  directions of any I2(m), m >= 4, share one actual opposite if their
  block energy is at most 1/7. It suffices to use the stronger common-
  wall threshold 1/3 for correlations at least 1/sqrt(2).
* An opposite wall cannot close directions in two orthogonal whole
  components. If it closes one edge, its positive-definite white
  endpoint forces at most one contact in each other whole component.
  Every positive spectator cost is at least 1/4, and is strictly LESS
  than its squared height h; equality would be an ideal endpoint.
* If h<1/4, it has no spectator contact at all. This is the universal
  small-height isolation argument, independent of their types.
* At a farthest vertex the ONLY short A1 walls, for
  L < L_* = asinh(1/sqrt(2)), are an isolated dotted wall I, or a dotted
  wall with a single A1_3 spectator J. There is at most ONE short A1
  direction. A J pair funds its own two A1 components above their
  combined 2/10 baseline. Otherwise nonshort A1 energy is greater than
  1/10 per component. An I wall has 0<h<1/2; if h>1/5, the total A1
  energy is greater than a/10-2/35, where a is the number of A1 factors.

These interfaces are derived in the pinned small-A proof. The two
possible 7 defects will be kept distinct below.

## 2. A universal dihedral floor and its sole exception

We prove that an I2(m) block with m >= 4 pays more than 1/7, except
possibly when m=7 and an ISOLATED single 3 wall closes the block.

Put k=cos(pi/m), d=1/(1-k^2). Its unit inverse is

    C=d[[1,k],[k,1]].

Suppose its energy is at most 1/7. Section1 gives a common actual
opposite wall c for both directions. Let h be its squared height and Q
the total cost in other components. Its two contacts are b,g >= 0.
If both are positive, both are at least 1/2. At the actual release of
the b direction the positive white residual is delta=1-Q-g^2, and

    cosh L=sqrt(d)(b+kg)/sqrt(delta)
        >=sqrt(d)(1+k)/sqrt(3)
         =cot(pi/(2m))/sqrt(3)
        >=(1+sqrt(2))/sqrt(3)>11/8.

The last comparison follows from 128sqrt(2)>171. Its individual energy
already exceeds 3/19>1/7. Thus c has exactly one contact beta.

Its own cost q=d beta^2 has smallest derivative square k^2 q. Energy
at most 1/7 bounds every actual tangent square by 7/16; hence

    q-1 <= h <= (7/16)k^2 q,
    q(1-7k^2/16)<=1.                                  (1)

Since beta^2>=1/4, the left side is at least (9d+7)/64. This exceeds
one for every m>=8, since d>=4+2sqrt(2)>19/3. We are left with just
the four OLD labels m=4,5,6,7; this is a derived finite interval, not
an assumed bound on native contacts or a height search.

For m=4, (1) gives beta^2<=16/25<cos^2(pi/5). Only labels3/4 remain,
with q=1/2 or1 and h<=7/32<1/4. Small-height isolation gives Q=0 and
the impossible h=q-1<=0.

For m=5, d=2+2/sqrt(5) lies between25/9 and3. A contact of label4 or
larger makes the left side of(1) at least (9d+7)/32>1. Thus beta=1/2,
q=d/4<1 and h<=7(d-1)/64<7/32<1/4. Isolation again contradicts h>0.

For m=6, only label3 passes(1). It has q=1, so the old component plus
c is positive-null and forbidden by compactness.

For m=7, d>5 excludes label4 and larger using the same (9d+7)/32
bound. Therefore beta=1/2. Equation(1) also gives

    Q=1+h-d/4 <= (57-9d)/64 < 3/16.

The minimum positive spectator cost is1/4, so Q=0. This is exactly the
named exception: OLD I2(7), NEW single 3 contact, supported only in
that one block. Dotted contacts were included among the larger beta
values and have not been dropped.

The exception still pays more than9/100. Its actual mate has

    cosh L_mate=cot(pi/7)/sqrt(3)>109/91,

so tanh^2(L_mate/2)>9/100. For an exact check, put q_7=cos(2pi/7).
The elementary bounds623/1000<q_7<5/8 follow from the increasing cubic
8q_7^3+4q_7^2-4q_7-1 on[1/2,1]. Then
cot^2(pi/7)/3>1623/1131>(109/91)^2; the final cross-product difference
is2652. Its squared height is

    h_D=(1/4)csc^2(pi/7)-1 in (0,1/3).                 (2)

There is at most one exceptional BLOCK. Two in different blocks have
isolated spatial supports, so their wall product is -h_D in(-1/2,0),
forbidden by the Coxeter gap.

## 3. Exhaustion of the rank-at-most-four finite types

For completeness the small-rank finite list can be derived directly.
A positive-definite Coxeter graph is a forest: on any cycle all bond
cosines are at least1/2, so its all-ones vector has nonpositive square.
Ranks1 and2 therefore give A1 and I2(m), m>=3.

A connected rank3 tree is a path. Positivity requires the sum of its
two bond cosine squares to be less than one. Thus its bonds are
(3,3),(4,3) or(5,3), giving A3,B3,H3.

A rank4 tree is a star or a path. For the star positivity requires
the sum of its three bond squares to be less than one; only three3
bonds qualify, giving D4. For a path with bond cosines a,b,c its
determinant is (1-a^2)(1-c^2)-b^2. Positive rank3 subpaths already
force each bond to be at most5 and forbid adjacent bonds >=4.
If the middle bond is4 or5, both end bonds must be3; middle5 has
determinant9/16-cos^2(pi/5)<0, leaving F4. If the middle is3, both
ends >=4 give determinant at most zero, so at most one end is4 or5.
This gives A4,B4,H4. The entire finite list used here is consequently

    A1; I2(m),m>=3; A3,B3,H3; A4,B4,D4,F4,H4.

All nonexceptional blocks in this list have energy strictly greater
than their rank/14, apart from the aggregated A1 convention already
specified in Section1. In particular no large-rank finite classification
or list of hyperbolic chambers has been imported.

## 4. Two different isolated 7 defects and their interactions

The defect in Section2 is OLD I2(7), NEW contact3, with height h_D in(2).
The one in the pinned small-A proof is OLD A2=I2(3), NEW contact7:

    h_A=(4/3)cos^2(pi/7)-1 in (0,1/12).

It also pays more than9/100. Its actual mate has
tanh^2 L>164/541>3600/11881, which gives half-tangent greater than3/10;
the strict rational cross-product difference is884.

Both defects are isolated to their own old block, but their labels,
heights and propagation rules are DIFFERENT. For walls isolated in
different blocks the spatial product is zero, so their full product
is minus the product of their heights. All h_A,h_D and short-A1 heights
are less than1/2. Thus an h_D defect cannot coexist with another h_D
block, an h_A block, or ANY short A1 wall I or J: every such pair would
have product in(-1/2,0). There is no hidden second deficit in its case.

## 5. The common propagation lower rate

The small-A proof establishes two UNIVERSAL facts, valid with arbitrary
other finite component types:

* If the OLD A2/NEW7 defect occurs, EVERY direction in the other blocks
  has squared half-tangent greater than1/5. Its proof retains the whole
  old A2 as white and allows every positive-definite other-wall contact
  there; it does not assume the other blocks are A types.
* If an isolated dotted A1 wall I has h<=1/5, every direction in the
  other blocks has squared half-tangent at least1/5. Again the actual
  wall-product and derivative proof does not restrict the other types.

Under this direction floor, EVERY block of rank at most4 has energy
greater than (11/144) times its rank. Here is the complete accounting:

* A1 has energy at least1/5. For any I2(m), its direction Gram is
  [[1,k],[k,1]], and the positive covering vector gives energy at least
  (2/5)/(1+k)>1/5, more than the required2*11/144.
* A3,B3,H3 already pay more than1/4>3*11/144 intrinsically.
* For A4 the uniform direction-cover coefficient is
  T=(28-8sqrt(6))/5>5/3. Its positive K^-1*1 gives energy at least
  T/5>1/3>4*11/144. This is an ACTUAL cover argument at the same vertex.
* B4,D4,F4 intrinsically pay more than1/3; H4 pays more than11/36,
  precisely4*11/144.

Using the maximum of the intrinsic and covering estimates is essential:
no falsely large uniform covering coefficient for H4 or D4 is claimed.

## 6. The n >= 15 contradiction

If an OLD I2(7)/NEW3 defect occurs, Section4 excludes every short A1
direction and the OLD A2/NEW7 defect. All other blocks therefore pay
at least rank/14, with the strict component estimates. The total is

    ||x||^2 > 9/100+(n-2)/14
             >=9/100+13/14=713/700>1.

Suppose no such defect occurs. If the OLD A2/NEW7 defect occurs,
Section5 bounds its other n-2 ranks at rate11/144. Including its own
energy gives

    ||x||^2 > 9/100+11(n-2)/144
             >=9/100+143/144>1.

Neither defect can consequently occur. If there is an isolated dotted
A1-short wall I of height h<=1/5, Section5 bounds the other n-1 ranks
alone by more than11(n-1)/144>=154/144>1. Thus any remaining I has
h>1/5 and the full A1 energy exceeds a/10-2/35. Summing the baseline
component bounds now gives

    ||x||^2 > n/14-2/35 >=15/14-2/35=71/70>1.

If there is no I, nonshort A1 or the funded J case has no A1 debit,
and the simpler bound ||x||^2>n/14>1 applies. These cases exhaust the
whole rank-at-most-four finite forest and contradict ||x||^2<1.

This completes the pointwise theorem. It neither excludes arbitrary
farthest vertices with larger connected factors nor promotes the
overall compact Coxeter dimension ceiling. No attempt is made to close
the OLD I2(7) mixed n=14 branch by further local casework.
