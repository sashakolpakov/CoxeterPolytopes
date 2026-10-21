# Norm239: a unique first transverse edge with E7+I2(5) germ

21 October2026. **New complete manual supporting candidate, pending root
and outside-author review.** No worker, scientific/symbolic import,
enumeration process, test, Lean or TeX job ran. No frozen publication,
status guide or accepted proof was edited. All finite rows below follow
from displayed coefficient intervals and hand arithmetic.

## 1. Actual starting germ and claimed scope

Use O=Z[phi], phi²=phi+1, and the literal original module

    L=E8(O) perpendicular-sum [[-2phi,1],[1,8]],
    D=1+16phi, b=phi D=16+17phi, N(D)=-239.

The binary basis is X,Y. Its actual first edge has endpoint

    q=(21+33phi)X-(17+25phi)Y, q²=-2D,
    m=(3+7phi)X-(3+5phi)Y, m²=2, q.m=0,

and the full actual germ is E8+A1(m). Its native E8 walls are -e_j,
j0..7, with edges02,13,23,34,45,56,67 and native basis Gram2/-1.
The following source inputs are used at their declared boundaries:

| Input | SHA-256 | Interface |
|---|---|---|
| results/compact_h9_prime239_full_base_and_first_binary_edge_candidate_20261021.md | 879a3feba630072a4b127e68a3729c4da3fca4f3c44630d9aae623ddd899de02 | Literal q,m, original module, dyadic glue and actual full starting germ |
| results/golden_odd_glue_full_base_and_discriminant_reduction_20261021.md | b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce | Complete primitive norm alphabet2/2b up to unit rescaling; golden ring facts and full arrangement discreteness/local finiteness |

Finite spherical Coxeter classification and ordinary finite-germ
recognition remain explicit geometric inputs. No norm79 identity,
isometry, firstness result or invariant module is transported here.
No compactness or finite covolume of the full reflection chamber is
assumed. The author contributed to earlier norm79 work; it is unused.

Release native wall -e7 and retain E7 on0..6 and m. **The result proved
here is a unique first new simple wall**, namely

    R=phi(q-m)/2+phi³ Omega7,

with an actual finite endpoint germ E7+I2(5). This is one more actual
edge of the same literal full chamber, not a completion or exclusion.

## 2. The actual outgoing ray and an original-integral upper mirror

The E8 dual columns needed here are

    Omega7=(2,3,4,6,5,4,3,2),
    Omega6=(4,6,8,12,10,8,6,3).

Multiplication by the native Gram gives the respective coordinate
unit covectors at7 and6. Thus

    Omega7²=2, Omega6²=6, Omega7.Omega6=3.

The outgoing ray is p(tau)=q+tau Omega7, 0<tau<sqrt(D), of square
-2D+2tau². Its product with -e7 is -tau<0 and with every retained
old wall is0. Local finiteness and the ACTUAL starting germ give an
initial segment in the same chamber.

The proposed upper mirror R has ORIGINAL binary coordinates

    (13+22phi,-10-17phi)

and E8 part phi³ Omega7, so it belongs to L. Direct expansion gives

    R²=2, R.m=-phi, R.(-e7)=-phi³,
    R.(-e_j)=0 for j0..6,
    q.R=-D phi, Omega7.R=2phi³.

Norm2 and evenness make its reflection integral on the original module
and make R primitive. Its positive collision is

    tau*=D/(2phi²),
    1-tau*²/D=(7-4phi)/(4phi⁴)>0.

Here 7-4phi>0 follows from phi<7/4. Thus R is a genuine FINITE
mirror upper, established without assuming it is first or a facet.

## 3. Complete first-simple-wall list and original dyadic parity

Because the upper lies at a finite hyperbolic point, the compact closed
segment to it meets finitely many arrangement mirrors. A first new
simple chamber wall therefore occurs by tau*. Its finite local Coxeter
germ retains the eight old simple walls E7+m. Finite classification
gives only these possibilities: E7+I2(k), or E8+A1(m), the latter
attaching to native node6 by simple3. A connected rank9 finite type
cannot contain E7, and no other finite one-node completion of E7 exists.

Normalize a primitive first wall to square2 or2b. In square2 the
E7+I2(k) contact with m is -kappa, where kappa=0,1,phi only.
For a nonzero native contact kappa>=1, kappa<2, and conjugate
definiteness gives |kappa'|<2. Writing kappa=A+Bphi gives B=0 or1;
the physical strip then gives kappa=1 orphi. This includes every
native finite label without a root-scale assumption.

A square2b wall cannot have a nonzero finite contact p with ANY
retained norm2 wall. Reflection integrality gives p=b eta, eta in O.
Its angular invariant p²/b=b eta² has positive integer norm at least239
if nonzero. Physical finite intersection and conjugate definiteness
instead give both embeddings strictly below4, hence norm<16. This
contradiction proves orthogonality, including the E7 attachment case.

Consequently a first simple wall has exactly one of FIVE profiles:

    square2: epsilon=0, kappa=0,1,phi;
    square2: epsilon=1, kappa=0;
    square2b: epsilon=0, kappa=0.

Write its vector as

    r=(T q-kappa m)/2+h Omega7+epsilon Omega6.       (1)

The old released wall remains an actual chamber wall, so acuteness
gives h=-r.(-e7)>=0. If epsilon0, a positive collision requires h>0.
The E8 coefficients force h in O. Original binary coordinates force
T in O and, EXACTLY,

    T=kappa modulo2O.                              (2)

Indeed q and m both reduce to (1+phi)(1,1) modulo2 in the ORIGINAL
X,Y coordinates, and1+phi is a unit modulo2. Conversely (2) makes
the original binary coordinates integral. This uses the actual
dyadic quotient, not a fictitious unimodular (q,m) basis.

For square2 the norm and collision equations are

    D T²=4h²+12epsilon h+12epsilon²+kappa²-4,
    tau=D T/(2h+3epsilon), T>0.                    (3)

In particular T is even for kappa0. These necessary profiles may
overlist walls not acute to other old facets, which is harmless:
every first actual simple wall is included. The argument does NOT
assert that all mirrors meeting the transverse ray have finite native
products with the retained walls; it classifies the first SIMPLE wall.

## 4. The three epsilon0 norm2 profiles: complete finite rows

For epsilon0, put U=T/2 when kappa0. Equations(3), tau<=tau* and
the negative conjugate D' give

    h² <= (4-kappa²)phi⁴/(7-4phi),
    |h'|² <= (4-kappa'²)/4.                        (4)

The identity `(7-4phi)(1+phi)=3-phi` is useful here. With
8/5<phi<81/50, the physical upper for kappa0 is<7.3, that for
kappa1 is<6.3, and for kappaphi it is EXACTLY phi³=1+2phi.
For the first two bounds, compare respectively
`4phi⁶/(3-phi)<7.3²` and `3phi⁶/(3-phi)<6.3²` by the same
positive rational phi bounds.

Write h=A+Bphi, B=(h-h')/sqrt5 and A=B(phi-1)+h'. For kappa0,
|h'|<=1 and sqrt5>2.2 give B0..3, with exactly the following
nonzero possibilities. The equation is D U²=h²-1.

| h | h²-1 | Rational norm |
|---|---|---:|
|1|0|0, forcing U0 and no collision|
|phi|phi|-1|
|1+phi|1+3phi|-5|
|1+2phi|4+8phi|-16|
|2+2phi|7+12phi|-11|
|1+3phi|9+15phi|-9|
|2+3phi|12+21phi|-45|

Every nonzero norm fails divisibility by239, necessary for D U².
Thus this entire profile contributes no positive cut before the upper.

For kappa1, |h'|<.87, h<6.3 give B0..3. The B0 row is empty;
B1 has A0,1; B2 has A1,2; B3 has A1,2, but2+3phi>6.3 excludes
its last row. These are exactly the five cases below. The norm
equation is D T²=4h²-3.

| h | 4h²-3 | Rational norm |
|---|---|---:|
|phi|1+4phi|-11|
|1+phi|5+12phi|-59|
|1+2phi|17+32phi|-191|
|2+2phi|29+48phi|-71|
|1+3phi|37+60phi|-11|

Again no norm is divisible by239. No parity test could restore a row.

For kappaphi, |h'|<1 and h<=1+2phi give B0..2. B0 has no positive
row, B1 has A0,1, and B2 has A1 (A2 exceeds the physical upper).
The complete remaining table is

| h | 4h²+phi-3 | Rational norm |
|---|---|---:|
|phi|1+5phi|-19|
|1+phi|5+13phi|-79|
|1+2phi|17+33phi=D phi²|-239|

Only the last can be D times a square. It gives T²=phi², so T>0
forces T=phi; (2) holds. This is exactly the integral upper R.
It has collision tau*, not an earlier collision.

## 5. E8-completion and square2b profiles are absent in the interval

For epsilon1,kappa0 write T=2U. Equations(3) give

    D U²=(h+1)(h+2),
    tau²=D*[1-1/(2h+3)²].

The upper and the kappa0 bound above yield0<=h<2.15. Conjugate
positivity gives -2<h'<-1. Thus B=(h-h')/sqrt5 is exactly1,
and A=B(phi-1)+h' is exactly-1. The sole row is h=phi-1, where
(h+1)(h+2)=phi³ has norm-1, not divisible by239. No such wall exists.

For square2b, reflection integrality and E8 unimodularity force
h=b j, j in O. Equations(1)–(2) have T even, so r=U q+b j Omega7
with U in O. Its norm and collision give

    D U²=h²-b,
    j² <=4phi⁴/[b(7-4phi)]<1.1²,
    |j'|²<=1/b'<1/5.                              (5)

The middle denominator is b(7-4phi)=44-13phi; clearing with
phi<81/50 proves the displayed strict rational upper. Also
b'=33-17phi>5. A positive collision has j>0. Hence
`0<j<1.1`, `|j'|<.45` force B=(j-j')/sqrt5=0. Then j is a
positive integer with |j|<.45, impossible. Thus this fifth profile
has no cut at or before tau*.

These are necessary conditions on every original-integral reflector;
no missing square-root membership test or nonprimitive norm class is
being declared absent merely because it was not enumerated.

## 6. Actual firstness, endpoint and remaining construction scope

The complete first-simple-wall list and all its finite coefficient
intervals leave only R. Since R itself is an original-integral
mirror with finite positive upper collision, some first simple wall
exists by then; therefore it is R and the collision is exactly tau*.
This proves firstness in the FULL arrangement, with no assumption
that the chamber is compact and no use of the binary-only mixed-E8
exclusion on this transverse ray.

An explicit future endpoint ray is

    Q=2phi² q+D Omega7,
    Q²=-2D(7-4phi)<0.

It is a positive multiple of p(tau*), is orthogonal to E7,m,R, and
has Q.(-e7)=-D<0. It is feasible against every original chamber
wall: the initial actual edge stays in the chamber up to its first
collision. The retained eight walls and the unique new simple wall
have positive definite Gram E7+I2(5). Ordinary finite spherical germ
recognition thus gives the FULL actual simple facet germ at Q.
Other mirrors in the generated finite local arrangement do not
constitute extra chamber facets or a competing first-simple-wall tie.

This is a complete manual next-edge theorem for the literal norm239
host, subject to the two named inherited inputs and ordinary geometric
foundations. Transverse continuation from E7+I2(5), whole-chamber finite
volume and a chamber period remain unresolved. No new example,
nonreflectivity, commensurability-class exclusion, dimension bound,
runtime PASS or Lean statement is claimed. The original dyadic module
must still be respected in every later continuation.
