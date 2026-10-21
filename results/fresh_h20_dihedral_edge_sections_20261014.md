# One-contact dihedral sections of the Borcherds H21 chamber

**Exact result:** none of the one-contact sections turning a norm-two source
`A2` edge component into `I2(m)`, for any integer `m >= 7`, is a finite-volume
Coxeter polytope. The argument covers every real transverse position, not a
finite label or parameter sample. Labels **4, 5, and 6 remain outside this
theorem**.

The separate [label-5 corollary](fresh_h20_dihedral_edge_sections_20261014_m5.md)
now closes label 5 using the newly certified coupled-H3 theorem. The
[G2 report](fresh_h20_dihedral_edge_sections_20261014_g2.md) gives three
complete continuous label-6 family exclusions and a broader numerical
census; labels 4 and 6 remain open in general.

This is an obstruction to this specified section construction, not to new
Coxeter polytopes in dimension 20. In particular, it does not address deleting
source inequalities, subdividing a section, or general section normals with
several nonzero pairings on the source edge.

## Source and normalization

The literal rational source has 210 walls in dimension 21, with 42 normals of
square 2 and 168 of square 4. The checker independently reconstructs its Gram
matrix and signature `(21,1)`, and checks that every reflection preserves the
coordinate lattice `Z^22`. Distinct old wall pairings are precisely

| Squared norms | Pairings |
|---|---|
| 2, 2 | −1, 0 |
| 2, 4 | −2, 0 |
| 4, 4 | −8, −6, −4, −2 |

Let `a,b` be a norm-two `A2` component of a positive-definite 20-wall source
edge. Normalize the spacelike section normal by
`(e,a)=1`, `(e,b)=0`, with zero pairings on the other 18 incident walls.
Write `e=w_a+z`, where `w_a` is the old dual weight and `z` belongs to the
orthogonal Lorentzian two-plane. Then `w_a^2=2/3`. Requiring the projected pair
to have angle `pi/m` gives

\[
 q=e^2=\frac{2C}{4C-1},\qquad C=\cos^2(\pi/m),\qquad z^2=q-\frac23.
\]

For `m >= 7`, `2/3 < q < 3/4`. The complete saved source walk contains 6,492
literal edge supports, including 2,524 `A2` components and hence 5,048 anchored
continuous families. The proof below is uniform and does not rely on testing
those families one at a time.

## The compact dihedral face

Suppose the section were a finite-volume Coxeter polytope. Its initial
`I2(m)` face `F` is an actual dimension-18 face: the initial 20 projected
normals are independent and positive definite, and all other source
inequalities are strict at the transverse edge point.

For `m >= 7`, `I2(m)` has no proper connected finite or affine Coxeter
extension. Therefore every bounding facet of `F` is orthogonal to both of
its dihedral normals. There can be no ideal vertex of `F`: the bounded
Euclidean link at such a vertex would place the label-`m` pair in an affine
component, which is impossible. Thus `F` is compact, and its orthogonal
bounding walls make it a compact Coxeter polytope in its own hyperbolic
space.

For an old wall `l` bounding `F`, put `k_l=(l,a)`. Orthogonality gives

\[
 (l,b)=0,\qquad (l,e)=qk_l,\qquad
 s_l=l-k_l e,\qquad s_l^2=d_l-qk_l^2.
\]

Positivity and the old pairing divisibilities leave four types:
`(d_l,k_l)=(2,0),(4,0),(2,-1),(4,-2)`. The first two are unchanged rational
roots; call the others modified short and modified long. For two face walls,

\[
 (s_l,s_j)=(l,j)-qk_lk_j.
\]

The exact checker handles every source-compatible pair of types on the
entire open interval `2/3<q<3/4`. Apart from right angles to unchanged roots,
all pairs involving a modified wall are divergent or have squared cosine
strictly between consecutive Coxeter values, with one exception: a modified
long wall and an unchanged long wall with old pairing −2 have squared cosine

\[
 \frac1{4(1-q)}.
\]

The other possible finite squared cosines lie in `(1/4,1/2)` or
`(3/4,4/5)`. The latter interval precedes the label-7 value, since
`cos^2(pi/7)>4/5`. These strict gaps exclude every integer Coxeter label.

## The sole exceptional label is 8

If the exceptional pair has finite label `n`, then `n>=7`, and substitution
of `q(m)` gives

\[
 \cos(2\pi/m)\cos(2\pi/n)=\frac12.
\]

Its sole integer solution for `m,n>=7` is `(8,8)`. By symmetry assume `m<=n`.
If `m>=8`, monotonicity gives a product at least `cos^2(pi/4)=1/2`, with
equality only at `(8,8)`. For `m=7`, all `n<=9` give a product below `1/2`,
and all `n>=10` give one above `1/2`.

The checker proves the strict bounds using the cubic polynomials for
`cos(2pi/7)` and `cos(2pi/9)`:

\[
 \frac{\sqrt5-1}{2}<\cos(2\pi/7)<\frac58,
 \qquad \cos(2\pi/9)<\frac45.
\]

Here `(5/8)(4/5)=1/2` and
`((sqrt5-1)/2) cos(pi/5)=1/2`. Polynomial derivative and exact endpoint-sign
checks certify the bounds; no numerical angle cutoff is used.

## Rational compact faces give the contradiction

If `F` has no modified facet, all its facet normals are unchanged rational
source roots. If `m!=8` and it has a modified facet, that facet is a compact
dimension-17 Coxeter polytope whose bounding normals are all unchanged
rational source roots, meeting the modified wall orthogonally.

For `m=8`, a modified short facet has the same property. A modified long
facet either also has only orthogonal neighbors, or it meets an unchanged
long facet at angle `pi/8`. Their intersection is a compact dimension-16
`I2(8)` face, whose bounding walls are all orthogonal to both. The exact
table shows that every right neighbor of a modified long wall is an
unchanged **short** source root. In particular, there is no overlooked
long-long right pair. This again produces a compact face with only original
rational root normals.

Each of these rational compact faces is impossible. Compactness forces its
facet normals to span its full Lorentzian space `U`; thus `U` is rational.
The lattice `U intersect Z^22` is full rank and preserved by its wall
reflections, since the literal old reflections preserve `Z^22`. A compact
Coxeter chamber would make its reflection group a cocompact lattice of
finite index in this integral orthogonal group. But the rational indefinite
form on `U`, of rank at least 17 here, is isotropic by Meyer's theorem; its
integral orthogonal group is nonuniform. A finite-index subgroup cannot be
cocompact. This uses an actual invariant lattice, not rational coordinates
alone.

## Remaining constructive cases and replay

The exact theorem leaves `m=4,5,6`. The apparent shortcut for `m=5` only
handles an `H3` extension attached to the unchanged end of the dihedral
pair. An extension attached to its modified end needs a separate argument;
it is not silently included in the exclusion above.

For `m=6`, the `G2` face can have ideal vertices. Its facets at finite points
are nevertheless orthogonal to `G2`. Since its initial 18 facets cannot
enclose a finite-volume dimension-18 face, at least one additional old wall
must be orthogonal to both projected `G2` walls. This gives a finite
necessary closure equation for each continuous source-edge family. A
separate numerical scout explores this equation and does not extend the
exact theorem to label 6.

Exact source and certificate:
[checker](../scripts/fresh_h20_dihedral_edge_sections_20261014.py),
[JSON](fresh_h20_dihedral_edge_sections_20261014.json).

```sh
python scripts/fresh_h20_dihedral_edge_sections_20261014.py --verify
```

The finite and affine diagram classification, compact-face argument, and
Meyer/arithmetic-group implication are mathematical premises explained
above. The saved replay checks the independent literal source arithmetic,
all uniform inequalities, the exceptional-label proof, and the edge counts.
An independent [source/theorem review and fresh saved replay](fresh_h20_dihedral_edge_sections_20261014_review.md)
also passed.
