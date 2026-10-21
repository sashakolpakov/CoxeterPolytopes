# A C4 x C2 K3 lift of the C18 H10/21 chamber

## Theorem

There is a special complex projective K3 surface X of Picard rank 20 with a
genuine automorphism group

\[
G\cong C_4\times C_2
\]

whose restriction to (N) is the archived (V_4) class `A18_4_0_C18`, where
(N) is Alexeev's even 2-elementary Picard stock `(18,4,0)`.  The common fixed
Néron--Severi lattice is integrally equal to the original rank-11 `C18`
fixed lattice.  Consequently the common fixed section of the preserved
48-wall reflection chamber is the already certified finite-volume,
noncompact 21-wall Coxeter polytope in `H10`.

The actual K3 group is important: it is not claimed to be a literal (V_4)
action on cohomology.  Two generators have order four and the same central
square; after restriction to (N), that square disappears and the image is
exactly (V_4).

This is a new project-level algebraic realization of the existing chamber.
It is not a new abstract polytope, arithmetic class, dimension, or general
finite K3-group type.

## Simultaneous transcendental lift

Set

\[
T=U(2)\oplus U(2)
\]

and take the two commuting isometries

\[
B_L=\begin{pmatrix}
0&0&0&-1\\
0&0&-1&0\\
0&1&0&0\\
1&0&0&0
\end{pmatrix},\qquad
B_R=\begin{pmatrix}
0&0&1&0\\
0&0&0&1\\
-1&0&0&0\\
0&-1&0&0
\end{pmatrix}.
\]

They satisfy

\[
B_L^2=B_R^2=-1,\qquad B_LB_R=B_RB_L,\qquad (B_LB_R)^2=1.
\]

In the exact ordering of the three archived nonidentity `C18` Picard
actions, use

\[
(B_LB_R,\ B_L,\ B_R).
\]

One common discriminant anti-isometry works simultaneously:

\[
P=\begin{pmatrix}
0&0&1&1\\
0&0&0&1\\
1&1&0&0\\
0&1&0&0
\end{pmatrix}.
\]

The certificate checks the quadratic form on all 16 discriminant classes
and all three equivariance equations.  Graph gluing with (P) produces the
even unimodular K3 lattice, and all three block actions become integral
isometries.

The two order-four glued actions commute, have the same nontrivial square,
and generate exactly eight matrices.  Their presentation is

\[
\langle x,y\mid x^4=y^4=1,\ x^2=y^2,\ [x,y]=1\rangle
\cong C_4\times C_2.
\]

The central square acts as `+1` on (N) and `-1` on (T).  Hence restriction
to (N) has kernel of order two and image (V_4), exactly the `C18` action.

## Common period and special Picard lattice

Take

\[
\omega=(-i,-i,1,1).
\]

For the three displayed transcendental actions its eigenvalues are

\[
(1,i,-i).
\]

Moreover

\[
\omega^2=0,\qquad \omega\bar\omega=8.
\]

For a rational vector ((a,b,c,d)\in T\otimes\mathbb Q), orthogonality is

\[
a+b=0,\qquad c+d=0.
\]

The saturated period-orthogonal kernel in the glued K3 lattice has rank 20,
signature `(1,19)`, determinant `-16`, and Smith diagonal `1^18 4^2`.
Thus X has Picard rank 20.

The induced action of (G) on this Picard lattice has common fixed rank 11.
An exact unimodular row-lattice comparison proves

\[
NS(X)^G=N^{C18}
\]

integrally.  Its determinant is 128 and its Smith diagonal is
`1^6 2^4 8`, exactly as in the original fixed-section certificate.

## Rootless coinvariant and the Kähler gate

The rank-9 algebraic coinvariant, the saturated orthogonal complement of
the rank-11 fixed lattice inside `NS(X)`, is negative definite.  It
has determinant `-512` and Smith diagonal

```text
1,1,2,2,2,2,2,2,8.
```

The certificate performs an exhaustive exact short-vector enumeration on
its positive negative.  Coordinate bounds are

```text
10,11,13,6,8,4,1,6,1,
```

and Schur-complement pruning leaves only nine complete leaves.  There are
zero vectors of square (-2).

This closes the group Kähler obstruction.  If a K3 root were orthogonal to
the entire common fixed positive cone, it would lie in that rank-9
coinvariant, contradicting the exact root census.  Every root therefore cuts
the nonempty strict fixed section in a proper hyperplane.  A generic common
fixed point in the source-section interior avoids all of them and determines
a (G)-invariant Kähler chamber.

The three glued maps are Hodge isometries preserving that chamber.  Strong
global Torelli realizes the full (C_4\times C_2) group by genuine
automorphisms of X.

## Hyperbolic and Vinberg conclusion

All three Picard actions permute the 48 source walls, and the archived common
fixed interior has positive pairing with every wall.  Because the common
fixed `NS` is exactly the archived `C18` lattice, its reflection-chamber
section is literally the same `H10/21` polytope.

The separate Vinberg certificate gives:

- 2,088 elliptic edges;
- 329 finite vertices;
- 22 ideal vertices;
- zero bad edges; and
- zero bad affine zero-neighborhoods.

Every elliptic edge has exactly two finite-or-ideal endpoints and every
affine full zero-neighborhood has the required cusp rank.  The polytope has
finite volume, and the 22 cusps prove that it is noncompact.

## Exact scope

Alexeev's 48-wall source is a complete mixed `(-2)`/`(-4)` reflection
chamber.  It is not the literal nef cone of this special K3 surface, and the
extra algebraic classes can further cut the actual nef cone.  What is proved
is an automorphism-invariant divisor/reflection chamber and its exact common
fixed Coxeter section; no quotient-nef claim is made.

The `H10/21` abstract chamber, its arithmetic lattice, and its containing
full reflection group were already certified.  Only the explicit K3-group
realization is new here.  Worldwide priority of the literal 21-wall chamber
remains open.

## Reproduction

The permanent files are:

- [exact C4 x C2 certificate](../data/enumerations/a18_c18_special_k3_c4xc2_lift.json);
- [certificate generator](../scripts/certify_a18_c18_special_k3_c4xc2_lift.py);
- [regression tests](../tests/test_a18_c18_special_k3_c4xc2_lift.py);
- [source construction certificate](../data/enumerations/alexeev_remaining_high_fixed_section_construction.json); and
- [Vinberg finite-volume certificate](../data/enumerations/alexeev_remaining_high_fixed_section_vinberg.json).

Run from the repository root:

```bash
PYTHONPATH=src:. python3 scripts/certify_a18_c18_special_k3_c4xc2_lift.py
PYTHONPATH=src:. pytest -q tests/test_a18_c18_special_k3_c4xc2_lift.py
```

The certificate payload SHA-256 is
`8b36468302ee6bc3b9cd3bd199578201cea8e26091ef0416cdd8937fe7da7cc3`.
The frozen JSON file SHA-256 is
`91dcfc1fec68826c5fb43b8f887cc1034d42cf24d2e054263fdb61331c4e3a17`.
