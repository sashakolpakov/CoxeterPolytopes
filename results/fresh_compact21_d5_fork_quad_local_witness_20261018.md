# A D5 fork quadrilateral: exact local witness and complete arithmetic right-corner branch

18 October 2026. **New exact local configuration and a complete unbounded
arithmetic classification of one finite-corner pattern.** There is a compact
quadrilateral face in an acute seven-wall hyperbolic five-dimensional local
chamber, with corner diagrams D5, A4+A1, A3+2A1, A4+A1. Its only basic
distance-two payment is the original D5 fork-leaf unit. The configuration
passes the full local arithmetic tests over Q(sqrt5).

This is **not a complete compact or finite-volume chamber**, either in H5 or
after adjoining orthogonal E8 spectators in H21. No global compact extension
is proved. It is a local proof-system obstruction to deleting the D5 fork
unit without checking its consumers. Other bonus families could change its
payment, and a genuinely global obstruction could still exclude its use.

More positively, every arithmetic quadrilateral in the same finite-corner
pattern has the same two opposite squared products, up to exchange. It
therefore forces a designated D5 triangle, and cannot coexist at a
D5+2E8 vertex with either E8 nonadjacent triangle type. This is a new
conditional incidence restriction, not yet a complete charge-saving rule.

## Exact witness

Use the D5 edges 01,12,23,24, with old roots0,...,4 of square two.
Put s=sqrt5, phi=(1+s)/2 and

  X=phi^2=(3+s)/2,   Y=(23+9s)/2.

Let c=root5 have square two and the single old contact c dot root4=-phi.
Let d=root6 be scaled to square2Y, with the single old contact
d dot root3=-Y. Let c dot d=0. All other unspecified products are zero.
The normalized opposite product for root3,d is -sqrtY; scaling d keeps
the entire symmetric Gram matrix field-valued:

```
       0   1   2   3   4    c     d
 0     2  -1   0   0   0    0     0
 1    -1   2  -1   0   0    0     0
 2     0  -1   2  -1  -1    0     0
 3     0   0  -1   2   0    0    -Y
 4     0   0  -1   0   2   -phi   0
 c     0   0   0   0  -phi  2     0
 d     0   0   0  -Y   0    0    2Y
```

The [producer](../scripts/fresh_compact21_d5_fork_quad_local_witness_20261018.py)
and [exact output](fresh_compact21_d5_fork_quad_local_witness_20261018.json)
check physical signature(5,1,1), nonidentity signature(6,0,1), and integrality
of all49 entries of the generalized Cartan matrix 2Gij/Gjj in Z[phi].
The nonzero off-diagonal support is a tree. All cyclic products therefore
belong to Z[phi]; the phi^2 cycle generates precisely Q(sqrt5).
The full span can be realized over that field using the first six roots
as a Lorentz basis. Under the nonidentity embedding their Gram is positive
definite. These are local admissibility and integrality statements; no
finite-covolume conclusion is inferred from an infinite-index reflection
subgroup passing them.

Set the white support to W={0,1,2}, and use cyclic black order4,3,c,d.
The four adjacent corners are positive definite. Eliminating W gives a
four-normal rank-three Lorentz Gram. The script constructs all four
vertex rays explicitly in a three-root Lorentz basis, checks their
timelike norms, all active equalities, all strict inactive inequalities,
and negative mutual products putting them on one future sheet. Their
four boundary edges form a closed convex quadrilateral in W-perpendicular.

Its vertex-ray sum has strictly negative product with each black wall.
A sufficiently small negative white-dual perturbation makes all three
white inequalities strict while preserving timelikeness and the strict
black inequalities. Thus these actual normals bound a nonempty acute H5
chamber, and the quadrilateral is its actual compact two-face, not merely
four unrelated positive matrices.

All127 nonempty principal subsets are covered:95 are checked exactly
and32 containing the strictly dotted pair3,d have a negative two-root
minor. None is affine. The four basic short-path payments are

|Corner diagram|Marked sides|Basic payment|
|---|---|---:|
|D5|3,4|1|
|A4+A1|3,c|0|
|A3+2A1|c,d|0|
|A4+A1|d,4|0|

The quadrilateral is intrinsically a Lambert quadrilateral with three
right angles and the remaining angle having cosine3/5. A compact face
does not imply that the whole seven-wall chamber is compact. In particular,
other root subsets are indefinite; no finite-volume chamber is produced.
The other three corners can receive nonzero payments from additional
families in a larger global weight rule, which have not been recomputed here.

## Complete classification of this finite-corner pattern

The following statement has no opposite-distance cutoff. Consider an actual
compact arithmetic chamber with ground field Q(sqrt5), a D5 fork face
marked3,4 and the same finite-corner pattern: each new side is orthogonal
to all of W, to its adjacent old black side, and to the other new side.
Any additional white spectators are assumed orthogonal to these four
black sides, as the two whole E8 factors force at a D5+2E8 vertex.
For clarity the norm-two normalized new roots are used in this proof.
Write their opposite products as -sqrtX and -sqrtY.

Here X and Y are the invariant length-two Cartan products
4(r dot t)^2/((r dot r)(t dot t)); their definitions do not require
the norm-two roots themselves to be field-valued. All normalizations
are by positive scalars at the physical embedding. Acuteness fixes
the negative signs of the opposite products. The zero products are
excluded below by opposite-corner indefiniteness. At the other real
embedding, a positive-definite arithmetic root form has positive root
squares, so its Cartan products are nonnegative. Equivalently, a real
normalized conjugate Gram can be used, up to diagonal sign congruence.
No negative conjugate square root or nonreal root normalization is
silently allowed.

Eliminating W gives old black diagonal5/4 and mutual product -3/4.
Eliminating the original positive D5 corner then gives the two-by-two
Schur matrix with diagonals2-5X/4 and2-5Y/4 and off-diagonal
-3sqrt(XY)/4. The complete quad star has corank one, so

  XY-(5/2)(X+Y)+4=0,
  (2X-5)(2Y-5)=9.                                      (1)

An opposite white-augmented pair is indefinite: otherwise positive
definiteness would give an actual common corner inside the two-face,
or singular positive semidefiniteness would violate compactness.
Here is the acute projection argument for the first implication. If B
is the positive-definite Gram matrix of the selected white and opposite
walls, its off-diagonal entries are nonpositive, and B^{-1} is entrywise
nonnegative. To see this, write B=D(I-C)D with D positive diagonal,
C entrywise nonnegative and spectral radius less than one; then use
(I-C)^{-1}=sum_{k>=0} C^k. For a strict interior timelike vector q,
write b=B^{-1}(q dot n_i). Every b_i is nonpositive. The projection
p=q-sum_i b_i n_i is timelike, since its squared norm is
q^2-(q dot n_i)^T B^{-1}(q dot n_i)<0. It is perpendicular to all
selected walls. Every other wall m satisfies
p dot m=q dot m-sum_i b_i(n_i dot m)<=q dot m<0.
Thus p is an actual point of precisely that face of the chamber.
For the selected opposite walls this would place their intersection
in the quadrilateral, contrary to their being opposite sides.
Its two-by-two residual determinant is5/2-X or5/2-Y. Hence

  X>5/2, Y>5/2.                                       (2)

Thus the physical Schur diagonals are strictly negative and its
off-diagonal -3sqrt(XY)/4 is strictly negative. Its rank-one branch is
the negative outer-product branch forced by the two outward walls at
the common old vertex. The other sign cannot supply an additional
quad branch; it would contradict acuteness of an opposite product or
the fixed Schur off-diagonal. Singular opposite corners are forbidden
by compactness, so the equality cases X=5/2 or Y=5/2 do not enter.

The arithmetic length-two cycles give X,Y in O=Z[phi]. The conjugate
positive-semidefinite Gram has nonnegative squared products X',Y',
and its old-D5 Schur diagonals give

  0<=X',Y'<=8/5.                                      (3)

Set U=2X-5, V=2Y-5. Equations(1)--(3) imply

  U,V in O, UV=9, U,V>0,
  -5<=U',V'<=-9/5.                                   (4)

These constraints can be solved exactly using only elementary quadratic
integer arithmetic, as follows.

First,3 is inert in O: the polynomial t^2-t-1 has no root modulo3.
If an algebraic integer has norm divisible by3, its residue in the field
O/(3) has zero norm and is zero, so the integer is divisible by3 in O.
Since N(U)N(V)=81 and both norms are negative, repeated division shows

  U=3^a epsilon, a in{0,1,2}, N(epsilon)=-1.

Every positive unit of O is a power of phi. For a short proof, multiply
such a unit by a suitable power of phi so that1<=beta<phi. If its norm
is+1 then beta-beta' lies in[0,1), so its sqrt5 coefficient is zero and
beta=1. If its norm is-1, then2<=beta-beta'=beta+1/beta<sqrt5, whereas
beta-beta' is an integer multiple of sqrt5; this is impossible. Thus
the original unit is phi^r. Negative norm makes r odd.

The condition X=(5+U)/2 in O implies U=1 modulo2. The image of phi has
order three in the field O/(2), since t^2-t-1 becomes t^2+t+1. Therefore
3 divides r, and altogether r=3 modulo6. The bound in(4) becomes

  9/5 <= 3^a phi^(-r) <= 5.                            (5)

There are exactly two possibilities:

|a|r|U|
|---:|---:|---|
|0|-3|phi^(-3)=sqrt5-2|
|2|3|9phi^3=18+9sqrt5|

Indeed phi^3=2+sqrt5 lies between4 and5, while phi^(-3)=sqrt5-2 lies
between1/5 and1/4. For a=1 the nearest allowed exponents already give
3phi^(-3)<3/4 or3phi^3>12, outside(5). For a=0 or2, moving r by another
six multiplies the value by phi^6>16 or its inverse, so the neighboring
listed exponents exclude every further r. This is an all-exponent
argument, not a scan of units up to a chosen height.

Substituting back into X=(5+U)/2,Y=(5+V)/2 yields precisely

  {X,Y}={ (3+sqrt5)/2, (23+9sqrt5)/2 }.                 (6)

The explicit witness above proves that this necessary parameter pair
really supplies the local compact quadrilateral. We have therefore
classified this single complete finite-corner pattern up to exchange.
No assertion that all D5 fork quadrilaterals have this pattern is made.

## Coupling to the new triangle theorems

The finite opposite product in(6) is -phi. Its new side is exactly the
previously classified D5 root L3 or L4, depending on exchange. In the
displayed witness it is L4, so the positive corners for marks(1,4) force
an actual triangle as well. The root's squared height is(-1+5sqrt5)/8.

At a D5+2E8 vertex the two possible E8 third-wall types have heights
4+2sqrt5 and1+sqrt5. Their squared mutual products with L4 are respectively

  (23+9sqrt5)/4,  3+sqrt5/2,

neither an algebraic integer in Z[phi]. Arithmetic cyclic integrality
therefore excludes both E8 types. Thus every quadrilateral in the
classified right-corner pattern is incompatible with every designated
E8 nonadjacent triangle at that same vertex.

Consequently there is an exact local incidence budget. At a D5+2E8
vertex let z be the indicator that its fork-leaf face has this classified
quad pattern, and let t_E be the total number of the eight possible
nonadjacent E8 triangles. The previous complete triangle domain gives
0<=t_E<=8, and the present theorem gives z=1 implies t_E=0. Therefore

  t_E + 8z <= 8,   equivalently z + t_E/8 <= 1.

This is a valid local budget assigning a unit to this specific quad
consumer and one eighth to each such triangle incidence. It does not
provide all the payments required by the global face theorem: other
D5 fork quad patterns are not covered, and E8 nontriangular corners
can still belong to quadrilaterals requiring payment. It is therefore
not labelled a proved reduction of any existing vertex allowance.

This is a useful possible priority condition for a new charging rule.
Using it to delete a unit still requires classification of the other
finite-corner patterns and preservation of every existing consumer, or
an independent replacement payment. New sides at a quad can meet a whole
untouched E7 through an E8 extension: the triangle's two-old-corner
rank22 isolation proof must not be reused for quadrilaterals without
its missing old-corner duty. Consequently no general A1+D5+E7+E8 quad
classification follows from the present orthogonal pattern.

## Reproduction and exact status

```sh
python3 -B scripts/fresh_compact21_d5_fork_quad_local_witness_20261018.py --verify
```

Generation passed in chunk `160297`, exit0, in0.058 seconds, with measured
peak RSS28,033,024 bytes. Separate-process replay passed in chunk
`c33aea`, exit0; the replay plus hash reporting took0.159 seconds, with
peak RSS28,000,256 bytes. The one-file job held the shared global proof
flock, used a59-second alarm, and checked peak RSS below one GiB. Its
state consists only of seven-root exact matrices and127 principal
subsets. No optimizer, NumPy, Lean, parallel worker, parameter search,
temporary proof output or bytecode cache was used.

The exact script checks the witness, all physical/conjugate matrix claims,
all vertex inequalities, the local Cartan integrality and the two E8
cyclic exclusions. The unrestricted unit classification in(1)--(6) is
the written mathematical proof above, not a conclusion extrapolated
from finite samples. Independent review and full Lean verification of
that proof are pending. No new global charge saving or dimension bound
is claimed.

Producer SHA-256:
`8cd609e79cf06f1bbea59faf57f9338fe35f06407784ce612efa4bc1e3eb59e0`.
Certificate SHA-256:
`559c3e7caa0fa73a531dc73e6371c5cb672a07b248306c85c12893b3a888469c`.
