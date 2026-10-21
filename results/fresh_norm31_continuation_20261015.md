# Discriminant-31 continuation: an exact nonreflectivity witness

The new lattice
\[
L=E_8\perp\begin{pmatrix}-2\sqrt2&-1\\-1&2\end{pmatrix}
\quad\text{over }\mathbb Z[\sqrt2]
\]
is **nonreflective**. A literal integral isometry preserves the full reflection
chamber and has trace \(9+5\sqrt2>10\), hence infinite order. Consequently no
finite-volume Coxeter polytope can have all its wall reflections preserving this
literal lattice. This is a host-specific obstruction, not a dimension bound or
an exclusion of all commensurable lattices or arbitrary geometric walls.

This conclusion is proved by a complete bounded root enumeration and one exact
corner map. It is not inferred from a finite search failing to close a chamber.
Primary generation and saved replay pass; independent review is pending.

## Complete height-100 prefix

The [new producer](../scripts/fresh_norm31_continuation_20261015.py) uses the
frozen [height-20 pilot and its norm-alphabet proof](fresh_frontier_construction_20261015.md)
without modifying them. With \(s=\sqrt2\), \(\delta=1+4s\), and
\(\beta=9+5s\), every primitive integral reflector has a representative of norm
2 or \(2\beta\). The full zero-height sector is \(E_8+A_1\).
The basepoint is \(H=(2,0^8,1)\), of norm \(-2\delta\), and root height is
\(\delta x^2/(2q)\).

The exact enumeration through height 100 has **32 control levels, 11,882 recursion
nodes, and 28 walls**. It reproduces all 14 height-20 levels and all their accepted
roots literally before continuing. At this bound the universal control estimates
are \(x^2<1100\), \(|x'|<2\); the coordinate estimates are
\(z^2<6666\), \(|z'|<3\). Explicit integer boxes contain these two-embedding
regions, and exact filtering, parity, long-norm divisibility, and rational Schur
pruning make the enumeration complete. Every accepted root passes literal
reflection integrality and all prior-wall inequalities. The full Gram, strict
center, and 28 relative interior facet witnesses are saved in the
[certificate](fresh_norm31_continuation_20261015.json).

## Complete finite-corner census

The [finite-corner checker](../scripts/fresh_norm31_continuation_20261015_corners.py)
constructs the pairwise finite-neighbor graph directly from the exact Gram and
enumerates every nine-clique with a separate fixed-size recursion. There are
only 35 candidate supports. Full exact LDL tests leave **11 finite vertices**;
each is checked physically against every wall. Their shared-eight-wall pairs
give ten finite edges and 79 remaining finite-vertex germs.

The [corner certificate](fresh_norm31_continuation_20261015_corners.json) is a
complete finite-vertex census. It is not a complete ideal/hyperideal cone census,
and no such census is needed for the nonreflectivity argument.

## One infinite-order chamber symmetry

The source corner is the initial \(E_8+A_1\) corner \(H\), with ordered walls
\([0,1,2,3,4,5,6,7,8]\). The target has the same literal norm-colored Gram under
the ordered correspondence
\[
[0,1,2,3,4,5,6,7,8]
\longmapsto [5,1,4,3,2,0,21,20,7].
\]
Its exact future timelike point, using the nullspace normalization in the
certificate, is
\[
V=\left(
\frac{16+14s}{17},\frac{50+14s}{17},
\frac{75+21s}{17},\frac{100+28s}{17},
\frac{150+42s}{17},\frac{125+35s}{17},
\frac{100+28s}{17},\frac{75+21s}{17},
\frac{75+21s}{34},-1\right).
\]
Its norm is \((829-590s)/578<0\). Its active set is precisely
\([0,1,2,3,4,5,7,20,21]\), and every other wall inequality is strict.

The [symmetry producer](../scripts/fresh_norm31_continuation_20261015_symmetry.py)
reconstructs the matrix from these two root-plus-apex frames. The positive
perpendicular scale is \(18+14s\), so
\(T H=(18+14s)V\). Direct arithmetic checks give
\[
T^TGT=G,\qquad T,T^{-1}\in M_{10}(\mathbb Z[s]),\qquad
\det T=-1,\qquad \operatorname{tr}T=9+5s>10.
\]
The literal matrix and inverse, all root images, the two complete wall sign
lists, the spherical inverse, and both basepoint decompositions are stored in
the [exact symmetry certificate](fresh_norm31_continuation_20261015_symmetry.json).
The conjugate form's positive LDL pivots are checked separately.

The use of finite corners is global, despite the finite height cutoff. For any
finite corner \(P\) of a genuine acute Vinberg prefix, write
\[
H=aP+\sum_i b_i n_i,
\quad a=\frac{(H,P)}{P^2}>0,
\quad b=G_S^{-1}\big((H,n_i)\big)\le0.
\]
The last inequality follows from the nonnegative inverse of the positive
definite spherical Gram and the basepoint's nonpositive wall pairings. Any later
simple wall \(r\) has \((H,r)<0\) and \((n_i,r)\le0\), hence
\((H,r)\ge a(P,r)\), forcing \((P,r)<0\). Thus no later wall cuts or touches
either corner. The two exact decompositions are verified in the checker.

An integral future Lorentz isometry permutes the full lattice reflection
arrangement. The checked map sends the complete outward local sector at the
source corner to the complete outward sector at the target corner. Both sectors
belong to the same full reflection chamber, so \(T\) stabilizes that chamber.
Agreement of an abstract diagram alone would not have sufficed; here the literal
outward normal images, future direction, and lattice inverse are checked.

A finite-order ten-dimensional matrix has trace of absolute value at most 10.
The displayed trace proves infinite order. The admissible integral isometry
group is discrete: bounded sets in the physical group together with the compact
conjugate orthogonal group bound all coefficients in both real embeddings, and
\(\mathbb Z[s]\) is a lattice in their product. Infinitely many powers of this
chamber stabilizer therefore give disjoint translates of a sufficiently small
interior hyperbolic ball inside the chamber. Its volume is infinite.

Finally, a cofinite reflection group whose walls all preserve \(L\) would be a
subgroup of the full lattice reflection group. Its existence would force the
latter group to have finite covolume, contradicting the chamber symmetry. This
gives the stated host-wide consequence without claiming a result for other
lattices or nonintegral reflections.

## Reproduction and limits

```sh
python3 scripts/fresh_norm31_continuation_20261015.py --verify
python3 scripts/fresh_norm31_continuation_20261015_corners.py --verify
python3 scripts/fresh_norm31_continuation_20261015_symmetry.py --verify
```

All three generation and saved replays pass. No height beyond 100 was searched.
The automorphism search need not be exhaustive: the one exact infinite-order
witness is sufficient. All earlier pilot files remain unchanged. This closes
the proposed discriminant-31 host and calls for a different lattice or a
geometric construction whose reflections need not preserve this lattice.
