# Two-angle repairs on nonelliptic H20 section supports

This separate bounded search found **two genuine finite-volume noncentral sections repairing two actual angles simultaneously**. Both still have another exact actual obtuse ridge, so neither is a Coxeter polytope. No general exclusion of nonelliptic three-wall supports follows.

The [independent standard-library/Fraction replay](../scripts/fresh_h20_bivariate_nonelliptic_verify_20261013.py) passes. It checks all declared angle-equation roots, the strict timelike interiors, both repaired actual ridges and the rejecting actual ridges. It also rejects a deliberately forged repaired-ridge point.

## Two constructed sections

All indices refer to the original one-based Borcherds H21 wall list.

| Section normal | Exact positive parameter | Repaired actual angles | Remaining actual defect |
|---|---|---|---|
| `r₁+s r₂₀−2s r₆₂` | `s=(√133−4)/36 ≈0.20923784985` | π/3 at 1/22; π/2 at 190/1 | Obtuse at 1/19 |
| `r₁+s r₄₈−2s r₆₁` | `s=(√494−2)/28 ≈0.72236109896` | π/3 at 1/22; π/2 at 42/1 | Obtuse at 1/4 |

Their ambient triple Gram matrices, in the displayed order, are

\[
 \begin{pmatrix}2&0&-2\\0&4&-4\\-2&-4&4\end{pmatrix},
 \qquad
 \begin{pmatrix}2&-2&-2\\-2&4&-2\\-2&-2&4\end{pmatrix}.
\]

Both spans are Lorentzian. Consequently these sections lie outside the [new elliptic-span obstruction](fresh_h20_bivariate_sections_20261013.md). Both normals have nonzero pairing with the original center, so the earlier central theorem does not reject them either.

Each section has an explicitly verified strict timelike interior point and a spacelike normal. Its finite volume follows from intersecting the finite-volume ambient Klein polytope with this genuine interior hyperplane. Every repaired and rejecting angle has an exact timelike ridge-point certificate; the verifier checks all 210 original inequalities and any redundant zero inequalities. The claimed failures therefore concern actual dihedral angles, not just raw pairwise Gram entries.

## Construction and bounded scope

Start with independent coefficients

\[
 e=r_i+s r_j-t r_k,\qquad s,t>0.
\]

For an auxiliary original wall `l` with `(r_l,r_k)<0`, impose the exact orthogonality equation `(e,r_l)=0`. It gives

\[
 t=a+bs,\qquad
 a=\frac{(r_l,r_i)}{(r_l,r_k)},\quad
 b=\frac{(r_l,r_j)}{(r_l,r_k)}.
\]

The auxiliary wall remains unchanged under projection, and its intersections with ambient walls orthogonal to it can furnish actual right-angled ridges. This is the first angle equation. The second equation sets a separately selected actual pair's squared cosine to `cos²(π/m)`, for `m=3,4,5,6`.

The scout explicitly excludes the line that would force every normal to be orthogonal to the old center. It selects twelve lines on twelve nonelliptic triple-orbit representatives, using the projected original center to ensure a viable nearby section and to identify an actual defective pair worth repairing. The line selection is numerical proposal generation; all resulting algebraic geometry is exact. There is no claim that the scout exhausts the 34 nonelliptic orbits, their sign choices, auxiliary walls or coefficient ratios.

The tested triple representatives were:

```
{1,2,62}, {1,3,62}, {1,20,54}, {1,20,55},
{1,20,62}, {1,20,79}, {1,20,84}, {1,20,91},
{1,20,111}, {1,20,143}, {1,20,190}, {1,48,61}.
```

The complete positive-root census for the declared equations has twelve roots. Two give the actual π/3 repairs above. The other ten solve the conjugate polynomial for the π/5 squared cosine but have the wrong real value; they are discarded by exact signs. There are no unresolved center or ridge scouts in this saved batch. An unsuccessful center scout would not have constituted a geometric exclusion.

The independent replay reconstructs all twelve lines and angle equations directly from the authoritative wall data. A squarefree product and rational Sturm count certify their complete positive-root denominators; disjoint isolating intervals match all archived roots. It then verifies each accepted section and each wrong-conjugate rejection. It does not rerun the line-selection heuristic.

## Artifacts and replay

- [Exact proposals and geometry certificates](../data/enumerations/fresh_h20_bivariate_nonelliptic_20261013.json).
- [Producer and deterministic scout](../scripts/fresh_h20_bivariate_nonelliptic_20261013.py).
- [Independent verifier](../scripts/fresh_h20_bivariate_nonelliptic_verify_20261013.py) and [PASS record](../data/enumerations/fresh_h20_bivariate_nonelliptic_verify_20261013.json).

```sh
python scripts/fresh_h20_bivariate_nonelliptic_verify_20261013.py
```

Both producer/verifier pairs in this round compile. This bounded result supplies exact noncentral algebraic constructions and identifies their remaining defects; it does not close the nonelliptic search space or produce a new Coxeter example.
