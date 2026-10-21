The all-high branch of the fixed compact H8 raw attachment is impossible. This holds for every real all-dotted cap direction whose twelve normalized projected squared norms are at least \(B=\sqrt5-1\), provided the twelve old-wall balls cover the placed source. There is no finite angle cutoff or numerical optimization in the proof.

The [producer](../scripts/fresh_h9_h8_all_high_20261014.py) and [exact certificate](fresh_h9_h8_all_high_20261014.json) pass generation and saved replay. They use the frozen independent quartic Fraction arithmetic and outward 128-bit dyadic intervals. The source vertices, old pairings, positive-support Grams and orthogonality relations are reconstructed from their bound inputs.

The [independent review](fresh_h9_h8_all_high_20261014_review.md) separately reconstructs the source vertices, all66 old pairings and all three norm obstructions with SymPy field arithmetic and outward 192-bit intervals.

Normalize the cap \(r\) and all twelve old normals \(n_i\) to squared norm2. Write \(a_i=-(r,n_i)>0\) and \(x_i=\kappa_i^2=a_i^2/2-2\). Thus the all-high hypothesis gives

\[
a_i^2\ge L^2=2(B+2)=2(1+\sqrt5).
\]

The two actual source vertices32 and33 have distance \(D\) with \(\cosh D=17+8\sqrt5\). A high ball has radius \(\rho(x)=\operatorname{arsinh}\sqrt{2/x}\), so its diameter is too small to contain both. For two old walls, put \(C=-(n_i,n_j)\). The projected center distance satisfies

\[
\cosh d(x,y;C)=\frac{\sqrt{(x+2)(y+2)}+C}{\sqrt{xy}}.
\]

Both this distance and the radii decrease as either positive variable increases. Any points in the two balls are at distance at most \(d+\rho(x)+\rho(y)\). For every nondotted old pair, \(0\le C\le\phi\); at \(x=y=B,C=\phi\), that bound is already strictly below \(D\). Consequently the balls covering32 and33 must come from one of the three literal old dotted pairs

\[
(4,11),\qquad(5,6),\qquad(8,11).
\]

For each pair the certificate evaluates the full ball-union diameter at \(x=A^2/2-2,y=B\). It is strictly below \(D\). Therefore coverage forces **both** anchor products \(a_i,a_j<A\). These are the three exact cases; decimals below are only readable summaries of outward rational interval bounds.

| Pair | Contact bound \(A\) | Orthogonal positive support \(S\) | \(\mathbf1^TG_S^{-1}\mathbf1\) | Lower bound on \(r^2\) |
|---|---:|---|---|---:|
|4,11|4|0,2,9|3/2|>5.0096|
|5,6|3|0,1,2,3|\(166+72\sqrt5\)|>1737.9235|
|8,11|3|0,2,5|3/2|>3.2160|

Here every support root has squared norm2, is orthogonal to both anchor roots, and its Gram is positive definite with entrywise nonnegative inverse. Consequently the projection of \(r\) onto \(S\) has square at least \(L^2\mathbf1^TG_S^{-1}\mathbf1\). The anchor pair has Lorentz Gram \(\left(\begin{smallmatrix}2&-C\\-C&2\end{smallmatrix}\right)\), where \(C>2\). Its projection contributes

\[
-\frac{2(a_i^2+a_j^2+C a_i a_j)}{C^2-4}
> -\frac{2A^2}{C-2}.
\]

The orthogonal complement of this nondegenerate Lorentz block and the positive support is positive definite, because the ambient form has signature(9,1). Hence

\[
2=r^2\ge L^2\mathbf1^TG_S^{-1}\mathbf1-\frac{2A^2}{C-2}>2,
\]

a contradiction in each case. Only two actual source vertices are needed; the proof does not infer a global diameter from samples.

This closes the no-low-center branch left open by the [uniform radial reduction](fresh_h9_h8_global_cap_gap_20261014_radial.md). The [separate polar-anchor gate](fresh_h9_h8_corner_anchors_20261014_polar_gate.md) treats the remaining low-center branch. The present certificate independently establishes only the all-high statement.

The covering hypothesis is material: these are the twelve frozen old-wall balls used in the raw attachment with the fixed source's perpendicular sidewalls. Additional nonperpendicular covering walls, replacement of old walls, and other compact H8 sources remain outside this theorem. This is not a nonexistence theorem for compact H9 Coxeter polytopes.
