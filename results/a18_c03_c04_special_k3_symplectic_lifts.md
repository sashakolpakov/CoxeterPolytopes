# Special K3 lifts of the 28-wall H11 and H10 chambers

## Theorem

Two more fixed sections of Alexeev's complete 48-wall ((18,4,0))
reflection chamber have genuine K3-automorphism realizations:

| section | K3 Picard rank | automorphism | fixed `NS` | hyperbolic chamber |
|---|---:|---|---:|---|
| `A18_4_0_C03` | 20 | symplectic Nikulin involution | rank 12, exactly the archived (N^{C03}) | finite-volume noncompact `H11/28` |
| `A18_4_0_C04` | 19 | symplectic Nikulin involution | rank 11, exactly the archived (N^{C04}) | finite-volume noncompact `H10/28` |

Here (N) is the even 2-elementary lattice with invariants ((18,4,0)).
The `C03` surface has Picard rank 20, while the `C04` surface has Picard rank 19.
The K3 surfaces occur at special periods: their Picard lattices strictly
contain (N).  Crucially, every added algebraic class is anti-invariant.
Consequently the fixed Néron--Severi lattice does not grow, and the original
hyperbolic dimensions 11 and 10 are unchanged.

This is a positive algebraic-realization result for two already certified
Coxeter polytopes.  They are not new abstract polytopes.  Symplectic
involutions on K3 surfaces are classical, and
neither their general automorphism type nor the abstract polytopes are
claimed new.  For a modern lattice treatment of prime-order symplectic K3
automorphisms, see Garbagnati--Sarti,
[*Symplectic automorphisms of prime order on K3
surfaces*](https://arxiv.org/abs/math/0603742).

## Common K3-lattice construction

The discriminant group of (N) is ((\mathbb Z/2)^4).  For both cases take

\[
 T=U(2)\oplus U(2),\qquad
 G_T=\begin{pmatrix}
 0&2&0&0\\2&0&0&0\\0&0&0&2\\0&0&2&0
 \end{pmatrix}.
\]

For each Picard involution (A), the certificate gives an involution (B)
of (T) and a discriminant anti-isometry (P:A_N\to A_T) satisfying

\[
 AP=P B\pmod 2.
\]

All 16 discriminant classes are checked exactly.  Gluing (N\oplus T) along
the graph of (P) gives an even unimodular lattice of signature ((3,19)),
hence the K3 lattice.  Equivariance makes (A\oplus B) integral on that
glue.  In both cases its characteristic polynomial is

\[
 (x-1)^{14}(x+1)^8.
\]

The saturated anti-invariant lattice has rank 8 and determinant (2^8).
Dividing its Gram matrix by two gives an even negative-definite unimodular
rank-8 lattice.  Therefore

\[
 H^2(X,\mathbb Z)^-=E_8(-2).
\]

In particular it contains no class of square (-2), exactly as required for
a symplectic Nikulin involution.

## The rank-20 `C03` period

Use

\[
 B_{03}=\begin{pmatrix}
 0&1&0&0\\1&0&0&0\\0&0&0&1\\0&0&1&0
 \end{pmatrix},\qquad
 \omega_{03}=(1,1,i,i).
\]

Then

\[
 \omega_{03}B_{03}=\omega_{03},\qquad
 \omega_{03}^2=0,\qquad
 \omega_{03}\overline{\omega}_{03}=8.
\]

For ((a,b,c,d)\in T\otimes\mathbb Q), orthogonality to the period is

\[
 a+b=0,\qquad c+d=0.
\]

Thus the entire rank-two anti-invariant part of (T) becomes algebraic.  An
exact saturated Smith-kernel computation inside the glued K3 lattice gives

\[
 \rho(X_{03})=20,\qquad \det NS(X_{03})=-16,
\]

with Smith diagonal (1^{18},4^2).  The action on `NS` has fixed rank 12 and
anti-invariant rank 8.  A unimodular integral row-lattice comparison proves

\[
 NS(X_{03})^{g}=N^{A_{03}}
\]

exactly, not merely over (\mathbb Q).

The discriminant identification used here is

\[
 P_{03}=\begin{pmatrix}
 0&0&1&1\\0&0&0&1\\1&1&0&0\\0&1&0&0
 \end{pmatrix}.
\]

## The rank-19 `C04` period

Use

\[
 B_{04}=\begin{pmatrix}
 1&0&0&0\\0&1&0&0\\0&0&0&1\\0&0&1&0
 \end{pmatrix}
\]

and, writing (u=2^{1/4}),

\[
 \omega_{04}=(1,\sqrt2,iu,iu).
\]

Again

\[
 \omega_{04}B_{04}=\omega_{04},\qquad
 \omega_{04}^2=0,\qquad
 \omega_{04}\overline{\omega}_{04}=8\sqrt2.
\]

The irrational coefficients are purposeful.  For a rational vector
((a,b,c,d)), orthogonality gives

\[
 a\sqrt2+b=0,\qquad c+d=0,
\]

and hence (a=b=0) and (c+d=0).  Exactly one anti-invariant rational
direction in (T) becomes algebraic; no fixed rational direction does.  The
saturated Picard lattice therefore has

\[
 \rho(X_{04})=19,\qquad \det NS(X_{04})=16,
\]

and Smith diagonal (1^{16},2^2,4).  Its fixed lattice has rank 11 and is
integrally equal to the archived (N^{A_{04}}).

Here the discriminant identification is

\[
 P_{04}=\begin{pmatrix}
 0&1&1&1\\1&0&0&0\\0&0&1&1\\1&0&0&1
 \end{pmatrix}.
\]

## Kähler chamber and Torelli

Both Picard actions permute all 48 walls of Alexeev's source chamber, and
their fixed sections contain an exact strict positive interior point.  It
remains to ensure that the special extra algebraic classes do not force every
invariant positive class onto a K3 root wall.

That obstruction is absent for a structural reason.  If a `(-2)` class had
zero projection to the invariant real subspace, it would lie in the
anti-invariant lattice.  But that lattice is (E_8(-2)), whose norms are
divisible by four.  Hence every K3 root cuts the invariant positive cone in a
proper hyperplane.  A generic invariant point in the nonempty strict source
section avoids their countable locally finite union.  The Weyl chamber
containing that point is invariant.

Choose the marking so this is the Kähler chamber.  The two glued involutions
are effective Hodge isometries, fix their periods, and preserve a Kähler
class.  Strong global Torelli realizes them as genuine symplectic
automorphisms of order two.  Their (E_8(-2)) anti-lattices identify them as
Nikulin involutions.

## Vinberg handoff

Finite volume is inherited from the separate literal Vinberg audits, not
from the existence of 28 walls:

| section | elliptic edges | finite vertices | cusps | bad edges | bad affine neighborhoods |
|---|---:|---:|---:|---:|---:|
| `C03`, `H11/28` | 8,851 | 1,252 | 64 | 0 | 0 |
| `C04`, `H10/28` | 6,956 | 1,087 | 74 | 0 | 0 |

Every elliptic edge has exactly two finite-or-ideal endpoints, and every
affine full zero-neighborhood has the required cusp rank.  Therefore both
polytopes have finite volume.  Their positive cusp counts prove that they are
noncompact.

## Exact scope

The geometric object certified here is the intersection of an
automorphism-invariant divisor/reflection chamber with the fixed
Néron--Severi space.  Alexeev's source chamber is a mixed `(-2)`/`(-4)`
reflection chamber; it is **not the literal nef cone** of either special K3
surface, and extra algebraic classes may further cut the actual nef cone.  No
quotient-nef statement is made.

The `H11/28` and `H10/28` chambers, lattices, and maximal arithmetic groups
were already certified.  The new result is their genuine special-K3
symplectic-involution realization.  Worldwide priority of the literal
28-wall abstract chambers remains open.

## Reproduction

The permanent files are:

- [exact two-case certificate](../data/enumerations/a18_c03_c04_special_k3_symplectic_lifts.json);
- [certificate generator](../scripts/certify_a18_c03_c04_special_k3_symplectic_lifts.py);
- [regression tests](../tests/test_a18_c03_c04_special_k3_symplectic_lifts.py);
- [source construction certificate](../data/enumerations/alexeev_remaining_high_fixed_section_construction.json); and
- [Vinberg finite-volume certificate](../data/enumerations/alexeev_remaining_high_fixed_section_vinberg.json).

Run from the repository root:

```bash
PYTHONPATH=src:. python3 scripts/certify_a18_c03_c04_special_k3_symplectic_lifts.py
PYTHONPATH=src:. pytest -q tests/test_a18_c03_c04_special_k3_symplectic_lifts.py
```

The certificate payload SHA-256 is
`44fbce0d83222fa8f1e868b6e03516f726a1b3bcf9592effa07b01c3183fac0f`.
The frozen JSON file SHA-256 is
`2a31a6ca346d6387fb9454638cf2bcdb7e73c901499527d0e17983cd28accd12`.
