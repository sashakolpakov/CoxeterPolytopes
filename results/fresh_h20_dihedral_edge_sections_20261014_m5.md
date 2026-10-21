# Label 5: both dihedral-extension orientations are excluded

**Exact corollary:** every real transverse one-contact `A2 -> I2(5)`
section in the saved Borcherds H21 chamber is non-Coxeter. This supplements
the [all-label `m>=7` theorem](fresh_h20_dihedral_edge_sections_20261014.md).
The remaining labels in this one-contact family are **4 and 6**.

The corollary depends on the frozen
[coupled-H3 obstruction](fresh_h20_coupled_h3_obstruction_20261014.md)
and the earlier
[one-contact A3 obstruction and q5 pair table](fresh_h20_h4_face_closure_20261014.md).
Its [small exact checker](../scripts/fresh_h20_dihedral_edge_sections_20261014_m5.py)
and [certificate](fresh_h20_dihedral_edge_sections_20261014_m5.json)
independently verify the reduction to those results. They do not purport
to rerun the H3 irreducibility census independently.

Write the original norm-two pair as `a,b`, with `(a,b)=-1`, and normalize
`(e,a)=1`, `(e,b)=0`, `q=e²=(1+sqrt(5))/4`. In a hypothetical finite-volume
Coxeter section their actual `I2(5)` face `F` has dimension 18 and is compact:
no affine Coxeter diagram contains the label-5 edge. At a finite boundary
facet, a connected extension of its two simple roots is necessarily `H3`.
There are **two orientations**, and both must be considered.

If the new wall `l` attaches by label 3 to the unchanged end `b`, it is
orthogonal to the modified end. Thus, writing `k=(l,a)` and `d=l²`, its
projected norm satisfies `d-q k²=2(l,b)²`. Irrationality of `q` and the
literal integer pairing pools force `d=2`, `k=0`, `(l,b)=-1`. The old three
walls are `A3` with the original one-contact endpoint modification. The
earlier A3 theorem excludes this `H3` face.

If `l` attaches to the **modified** end `a`, it is orthogonal to `b`.
Put `g=(l,a)` and `x=(e,l)`. The exact label-3 equation is

\[
 (3+2q)x^2-8gqx+4g^2q^2-dq(2q-1)=0,
 \qquad x>gq.
\]

The old three-wall Gram is positive definite because it is the projected
spherical Gram plus a positive rank-one matrix. The literal source pools
leave exactly four possibilities:

| `d` | `g` | Old diagram | Applicable theorem |
|---|---|---|---|
| 2 | 0 | `A1+A2` | joined two-contact H3 |
| 4 | 0 | `A1+A2`, different isolated norm | joined two-contact H3 |
| 2 | −1 | `A3`, center `a` | connected two-contact H3 |
| 4 | −2 | `B3`, center `a` | connected two-contact H3 |

For the connected cases, normalize `e` to square one. The center contact
has square `A=1/q=2-1/(2cos²(pi/5))`, and the contact ratio is exactly `x`.
The displayed quadratic agrees identically with the appropriate first two
connected cases in the coupled-H3 certificate. For the disconnected cases,
`Y=1/(2q)` and `X=(1-Y)/(1+3Y)` give precisely its joined-H3 normalization.
Consequently neither orientation admits a boundary extension.

Every facet of `F` is therefore orthogonal to the dihedral pair. The frozen
`q5` table now applies without alteration: any finite neighbor of a modified
facet is right and has an unchanged rational source normal. If all facets
are unchanged, `F` itself is a compact rational reflection polytope; otherwise
a modified facet is such a compact polytope of dimension 17. In both cases
the original integral reflections preserve the full lattice in the rational
span. Meyer isotropy contradicts compactness exactly as in the preceding
face theorem. This is not an inference from rational coordinates alone.

```sh
python scripts/fresh_h20_dihedral_edge_sections_20261014_m5.py --verify
```

Scope: the original 210 inequalities and the prescribed one-contact source
`A2` modification are retained. General section normals, facet replacements,
and Coxeter subdivisions are not excluded by this corollary.
