# An explicit good lift package for the Type-I Enriques Coxeter chamber

**Status:** exact positive reconstruction; all finite `L1`--`L4` gates pass.

## 1. Result

Let `R` be Shimada's twelve-wall Type-I (`12A`) system in

\[
L_{10}=U\oplus E_8(-1).
\]

The archived rank-19 K3-cover Picard lattice in OSCAR's
`generic_enriques_surface(172)` contains an explicit good lift package for
`R` in the sense of the variable Enriques lift theorem.  More precisely:

1. every downstairs wall has exactly two integral square-minus-two lifts;
2. the two lifts are exchanged by one explicit integral involution and are
   disjoint;
3. all twenty-four lifted roots have nonnegative pairwise intersection;
4. their anti-invariant differences generate the full saturated lattice

   \[
   A\cong E_8(-2)\oplus\langle-4\rangle,
   \]

   which is negative definite, rank nine, and rootless;
5. the invariant ample class

   \[
   H=(38,23,46,70,57,45,34,24,15,7)
   \]

   has downstairs square `70`, K3 pullback square `140`, and wall degrees
   `1^10,8^2`; and
6. an exhaustive bounded-root search proves the irreducibility gate `L4`
   for all twenty-four lifts.

Together with the already frozen positive-sheet, complete-cone, and literal
Vinberg certificates for `12A`, this is a complete explicit positive example
of the surface-free lift-package theorem.  It reconstructs a known
finite-volume noncompact `H9/12` polytope on a known Type-I Enriques surface;
it is not a new abstract polytope or a new surface.  The project-level gain is
the first exact worked noncanonical simultaneous lift package and a concrete
specification for what the future downstairs compiler must recover.

The complete coordinate certificate is
[`enriques_type_i_explicit_lift_package.json`](../data/enumerations/enriques_type_i_explicit_lift_package.json),
recomputed by
[`certify_enriques_type_i_explicit_lift_package.py`](../scripts/certify_enriques_type_i_explicit_lift_package.py).

## 2. Pinned geometric datum and lineage

The input is OSCAR 1.4.0's serialized datum
`TauTaubarGenericEnriquesSurfaceNo172.mrdi`.  OSCAR's documentation identifies
Number 172 as a very general Type-I Enriques surface: its K3 cover has Picard
rank `19`, Picard determinant `4`, and root invariant `A1+E8`; its Enriques
nef cone has twelve walls, all defined by rational curves.  See the official
[OSCAR Enriques-surface documentation](https://docs.oscar-system.org/stable/AlgebraicGeometry/Surfaces/EnriquesSurfaces/)
and Brandhorst--Shimada, *Automorphism Groups of Certain Enriques Surfaces*,
FoCM 22 (2022), DOI `10.1007/s10208-021-09530-y`.

The exact upstream byte hash and GPL provenance are recorded beside the
archived datum in
[`data/external/oscar_enriques_type_i/README.md`](../data/external/oscar_enriques_type_i/README.md).
The OSCAR presentation is a lattice chain

\[
S_Y(2)\subset S_X\subset II_{1,25};
\]

the last lattice is the Borcherds computational ambient lattice and is **not**
the K3 cohomology lattice.  The certifier keeps these roles separate.

The first ten OSCAR basis vectors of `S_X` are exactly its `S_Y(2)` basis, and
their Gram is twice Shimada's pinned `L10` Gram.  Thus no numerical Gram
matching or floating isometry is being guessed: the two inputs are joined by
one exact integral marking.

## 3. The global Picard involution and its glue

Write an `S_X` coordinate row as

\[
z=(a,c),\qquad a\in\mathbb Z^{10},\quad c\in\mathbb Z^9.
\]

Exact orthogonal projection to the invariant lattice yields the integral
matrix

\[
D=\begin{pmatrix}
0&0&1&0&0&0&0&0&0&0\\
0&0&0&0&0&1&0&0&0&0\\
0&0&0&0&1&0&0&0&0&0\\
0&0&0&1&0&0&0&0&0&0\\
4&2&4&6&5&4&3&2&1&0\\
0&0&0&0&0&0&0&0&0&1\\
0&0&0&0&0&0&0&0&1&0\\
0&0&0&0&0&0&0&1&0&0\\
4&1&4&7&6&5&4&3&2&1
\end{pmatrix}.
\]

The downstairs class and the involution are therefore

\[
r(z)=2a+cD,
\qquad
\tau(a,c)=(a+cD,-c).
\]

The certifier proves integrally that `tau^2=1` and
`tau Gram(S_X) tau^t=Gram(S_X)`.  The invariant lattice is
`U(2)+E8(-2)`.  A basis for the integral anti-invariant lattice is

\[
Z=(-D\mid 2I_9).
\]

Its Gram is a signed presentation of `-2 C_(E8+A1)`: the non-isolated
Dynkin tree has arms of lengths `1,2,4`, and all signs are removed by vertex
sign switches.  Hence

\[
A=S_X^{-\tau}\cong E_8(-2)\oplus\langle-4\rangle,
\quad |\det A|=1024.
\]

The Smith diagonal of `Z` in `S_X` is `1^9`, so `A` is saturated.  Its Gram
has diagonal divisible by four and off-diagonal entries divisible by two;
therefore every norm is divisible by four and `A` contains no
square-minus-two root.  Finally,

\[
[S_X:S_X^\tau\oplus S_X^{-\tau}]=512,
\qquad
\frac{2^{10}2^{10}}{512^2}=4=|\det S_X|.
\]

This is the exact simultaneous half-sum glue hidden by the downstairs wall
Gram.

## 4. All twenty-four curve lifts

For a fixed downstairs root `r`, a lift is

\[
b=\left(\frac{r-cD}{2},c\right).
\]

The parity gate is `r-cD=0 mod 2`.  Orthogonality of the invariant and
anti-invariant parts reduces `b^2=-2` to

\[
c^t\left(-\frac{1}{2}\operatorname{Gram}(A)\right)c=2.
\]

The positive form in parentheses is `E8+A1`.  Exact PARI `qfminim`
enumeration gives all `242=240+2` roots.  Applying the parity gate to each of
the twelve downstairs roots leaves **exactly two** roots in every case.  The
two candidates sum to `(r,0)`, are exchanged by `tau`, and meet zero.

The two cases missed by a naive coefficient box `[-1,1]^9` illustrate why an
exact lattice enumeration matters.  Their `c`-vectors, up to sign, are

```text
wall 6:  (-2,-2, 3,-4,-3,-2,-1, 0,0)
wall 10: (-3,-2, 4,-6,-5,-4,-3,-2,0).
```

Among all unordered pairs of the twenty-four distinct lifts the exact
intersection histogram is

```text
intersection 0: 246
intersection 1:  28
intersection 2:   2
```

with no negative entry.  Consequently every one of the `2^12=4096` choices
of a preferred sheet passes the signed nonnegativity gate.  For a
deterministic package we choose the lexicographically first lift of every
wall.  The resulting twelve anti-differences have rank nine and Smith
diagonal `1^9` in the `Z` basis, hence generate the full saturated
anti-invariant lattice `A`.

This is the precise positive counterpart to the canonical same-sheet
obstruction: simultaneous lifts exist, but their anti-invariant data are
noncanonical and cannot be read off as twice the downstairs Gram.

## 5. Exact ample and irreducibility gates

Projecting OSCAR's serialized Borcherds Weyl vector to `S_Y(2)` gives the
small integral class

\[
H=(38,23,46,70,57,45,34,24,15,7).
\]

It has positive pairing with every wall.  The first ten wall pairs are one;
the last two are eight.  The degree-one lifted roots are automatically
irreducible once Riemann--Roch chooses the positive sign.  It remains to
exclude a decomposition of any of the four degree-eight lifts.

Let `G=Gram(S_X)` and let `p=G H^t`.  Hodge index makes

\[
Q_H=-140G+2pp^t
\]

positive definite.  If `e^2=-2` and `0<H.e<8`, then

\[
Q_H(e)=280+2(H.e)^2\le 378.
\]

Thus one complete positive-definite short-vector enumeration decides `L4`.
PARI `qfminim` finds 472 signed vectors of `Q_H`-norm at most 378.  Independent
integer replay shows that every one is a square-minus-two root with absolute
`H`-degree between one and seven.  The positive-degree histogram is

```text
degree: 1  2  3  4  5  6  7
count: 20 20 24 28 36 44 64
```

There is no degree-zero root, proving that `H` is root-regular.  None of the
236 positive smaller-degree roots pairs negatively with any degree-eight
lift.  Therefore the exact finite criterion `L4` passes for all twenty-four
lifts, including the isotropic-residual cases that a mere two-root
decomposition test would miss.

## 6. Polytope conclusion and exact scope

The pinned independent downstairs certificate gives

```text
dimension/facets: H9/12
extreme rays:      14 finite + 9 ideal + 0 hyperideal
Vinberg audits:    114 elliptic edges + 12 affine neighborhoods
decision:          finite volume, noncompact
```

Hence the Type-I chamber now passes, in one joined certificate:

```text
L1  exact integral glue and nonnegative signed intersections
L2  saturated negative-definite rootless anti-lattice
L3  exact rank-19 Picard lattice and root-regular ample orientation
L4  exhaustive bounded-root irreducibility
     plus complete cone and literal Vinberg finite-volume gates.
```

This closes the existence side for one explicit noncanonical lift package.
It does **not** yet solve the projection problem: the current reconstruction
starts with OSCAR's `S_X`.  The next compiler must start only with the
downstairs Coxeter datum, enumerate the finite anti-Gram/glue/embedding
orbits, and recover this package among its outputs.  The present example is
therefore both a positive theorem and an exact regression target for the
surface-softening program.
