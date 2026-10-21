# A fourth actual vertex and a sealed diagonal face

20 October 2026. Paper geometry and integral-lattice argument only.
No search, worker, build or new compactness claim.

Import the first/second/third actual edges of `R_n=<-phi>+I_n+D4(O)`,
n>=5, from the following fixed sources:

* `golden_dyadic_d4_extension_first_edge_20261020.md`,
  SHA256 `3fc21952a3d119db4370ca46e39aa607666eaea7eaa54d835c3163e7a281035b`;
* `golden_dyadic_d4_extension_second_edge_20261020.md`,
  SHA256 `a894a30bb3b17e63a434049195931e8d6bd51a160b44d85141a820a08dd2e72e`;
* `golden_dyadic_d4_extension_third_edge_20261020.md`,
  SHA256 `11ccfde06a5ed9db05fafa803eb7a5fa8b399cf44234b8398f60fe487f350c3b`.

All paths are under `results/`. The conclusion below is a whole
actual-face identity, not identification by determinant alone.

## 1. A known actual facet with a finite completion is already opposite

The positive acute projection in
`compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`,
SHA256 `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`,
has no compactness requirement for the following use. If selected
actual facet normals M have positive definite Gram G and p is a
feasible future point, put `h=-M^T p>=0`. Then

    p0=p+M G^-1 h

is future timelike, orthogonal to M and feasible for every other
actual wall: the inverse of a positive definite Coxeter Gram is
entrywise nonnegative, and all other wall products with M are
nonpositive. At full selected rank this is an actual finite vertex.

Consequently an already actual facet which makes a positive definite
full corner with the retained walls closes that outgoing edge. The
segment between its projected vertex and the old vertex is feasible
and has those retained walls throughout. This argument applies also
to a locally finite chamber whose global volume is unknown.
Equivalently, it is the positive branch of the actual reciprocal
shadow lemma, `compact_actual_wall_reciprocal_shadow_20261020.md`,
SHA256 `c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`.

Thus the long second-edge census was unnecessary for firstness once
the actual old wall a4 and the positive definite corner
`B_n+A2+I2(5)` were recognized. Its saved complete census remains
valid independent evidence. This shortcut does not apply to an
arbitrary integral upper mirror which is not yet an actual facet;
in particular it does not turn the third-edge nonsimple root c
into an actual wall.

## 2. The fourth actual vertex

At the third endpoint T the full diagram is

    B_(n+1) + A1(a2) + I2(5)(r,s), s=a4.

Write `c=phi*(1;0^n;-1,0,0,-1)` and `v=c+e1`, so the B_(n+1)
component is the old B_n extended by v. Release a2. The retained
subspace has time coordinate1 and

    z1=-1/phi, z4=-1, z2+z3=-1/phi^2.

The old actual facet a3=2e3 is orthogonal to B_(n+1) and r, and has
product -2 with s. It completes the retained roots to the positive
definite whole corner `B_(n+1)+H3`. The preceding actual-wall lemma
therefore proves its firstness without another profile census.
The new actual vertex is the ray of the integral vector

    V=(1;0^n;-1/phi,-1/phi^2,0,-1),
    V^2=8-5phi=-phi^-5<0.                         (1)

The D4 sum is -2. The H3 simples are a3,s,r, all norm4, with bonds
3 and5. The norm1 vector c and the n old coordinate unit vectors
are orthogonal to this whole H3.

## 3. A primitive norm4 H4 cannot occur in this lattice

Suppose four primitive norm4 reflectors in R_n form a simple H4.
Their Gram is `2C`, where C has diagonal2 and the successive
off-diagonal entries `-phi,-1,-1`. Its determinant is

    det C=5-3phi=phi^-4,

a unit. Each norm4 reflector pairs every lattice vector into2O.
For any x in R_n, writing its four pairings as2k, its orthogonal
projection onto their root span N has coefficient vector

    (2C)^-1*(2k)=C^-1*k in O^4.

Therefore N splits integrally and orthogonally from R_n. But its
determinant ideal is `(2)^4`, whereas `det R_n` has ideal `(2)^2`.
An integral orthogonal complement cannot have a fractional
determinant ideal. This contradiction excludes the proposed H4.

This also applies to any finite H4 root group containing the present
H3, even if that H3 is not parabolic in it. All H4 reflections are
conjugate, so their primitive normal lengths have the same squareclass
as the norm4 H3 roots; an integral conjugacy preserves primitivity.
The norm4 unit normalization and the argument above then apply to
its simple roots.

A finite irreducible rank4 group containing H3 can only be H4.
Thus every mirror crossing the orthogonal hyperbolic subspace of
this H3 is orthogonal to H3. A nonorthogonal crossing would give
such a finite rank4 group at its timelike intersection. This proves
the seal for the FULL arrangement, not just for chosen native
continuations along one edge.

## 4. The whole orthogonal face is a literal diagonal chamber

Set `V0=phi^3 V`, so `V0^2=-phi`. The n original I_n unit vectors,
c and V0 are pairwise orthogonal, integral and span the H3-perpendicular
K-space. Their Gram is unimodular and exactly

    <-phi> + I_(n+1).                            (2)

For any vector of R_n in this K-space, pairings with these unit
basis vectors recover integral coefficients. Thus(2) is the exact
perpendicular lattice, not just a sublattice of the same rank.
The same projection argument gives an integral orthogonal split
`R_n=M_(n+1) + P3`, with P3 positive definite and determinant ideal
`(2)^2`.

Every intrinsic primitive norm1 or norm2 reflector of the unimodular
diagonal lattice extends to R_n through this integral split. Conversely,
an orthogonal primitive norm4 reflector would pair the unimodular
summand into2O and hence belong to2M_(n+1), making it imprimitive
in R_n. The imported norm pool leaves no further class. Together
with the full seal this proves that the actual H3 face has precisely
the FULL `M_(n+1)=<-phi>+I_(n+1)` chamber arrangement.

For the target n5 this is the known compact diagonal dimension6
face. It does not prove the dimension9 extension compact or
noncompact. A published literal diagonal cutoff, if separately
accepted with its precise scope, can be imported only for the
corresponding larger n; it is not part of this own-proof identity.
No period, global determinant census or new example is claimed.
