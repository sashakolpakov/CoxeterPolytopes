# Fresh H9 construction: facet rigidity and changed Gram diagrams

**Outcome (2026-10-13): no new compact example.** The first part proves an
exact restriction on a proposed construction: a literal Bugaenko H8 facet
cannot meet any adjacent facet nonorthogonally. The second part changes the
old Gram diagram itself and performs a separate bounded numerical search
for compact polar truncations in H9. The numerical search is not an exact
nonexistence certificate.

## Exact rigidity of a literal Bugaenko facet

Suppose a hyperbolic Coxeter polytope has a facet isometric to the calibrated
eleven-wall Bugaenko H8 polytope. Write its normalized facet roots as `s_i`,
with square two and Gram matrix `G`. If the ambient adjacent wall at `s_i`
meets the distinguished facet at angle `beta_i=pi/m_i`, its ambient normal
can be written

\[
 r_i=\sin(\beta_i)s_i-\cos(\beta_i)e,
 \qquad(e,e)=2,\quad e\perp s_i.
\]

Thus

\[
 (r_i,r_j)=\sin\beta_i\sin\beta_j\,G_{ij}
             +2\cos\beta_i\cos\beta_j.                 \tag{1}
\]

When `G_ij=0`, acuteness forces at least one of `beta_i,beta_j` to equal
`pi/2`. When `G_ij=-1` and one endpoint is orthogonal to the distinguished
facet, the other must be orthogonal as well: otherwise (1) lies strictly
between `-1` and `0`, which is not a Coxeter pairing for square-two roots.

In the archived numbering, the label-three graph has components

```
{0}, {1,2,3,4,5,6,7,9,10}, {8}.
```

Nodes 1 and 3 in the large component are orthogonal. At least one is
therefore untilted, and the label-three edges propagate orthogonality to
all nine nodes of that component.

Each remaining endpoint, 0 or 8, has one label-five edge to an already
untilted neighbor. A nonright ambient angle `pi/m`, `m>=3`, would change
that edge pairing to `-phi*sin(pi/m)`. The exact squared magnitudes are

\[
 m=3:\quad\frac{3+3\phi}{4}\in(1,2),\qquad
 m=4:\quad\frac{1+\phi}{2}\in(1,2),
\]

and for all `m>=5`, monotonicity gives

\[
 0<\phi^2\sin^2(\pi/m)
 \leq\phi^2\sin^2(\pi/5)=\frac{2+\phi}{4}<1.
\]

None is a Coxeter pairing: the first nonzero square is 1 for label 3,
the next is 2 for label 4. Therefore all eleven adjacent ambient facets
meet the Bugaenko facet at right angles.

The normalization here is square-two roots: a finite Gram pairing is
`-2*cos(pi/k)`. The displayed squared pairings, including the JSON field
`endpoint_remaining_pairing_squares`, are therefore `4*cos(pi/k)^2`,
not normalized `cos(pi/k)^2`.

This rules out obtaining a new construction merely by tilting adjacent
walls while preserving that literal facet. It does **not** rule out a
compact H9 polytope with a right Bugaenko facet, or any other compact H9
polytope. No worldwide novelty claim is made for this local rigidity fact.

The exact checker reads the calibrated source, verifies its label-three
components and orthogonal pair, and performs every remaining comparison
over `Q(phi)`:

```sh
PYTHONPATH=src python3 scripts/fresh_bugaenko_facet_tilt_rigidity_20261013.py
```

Result: `PASS_EXACT_LITERAL_FACET_TILT_RIGIDITY`.

## A different construction: deform the Gram matrix and truncate

The subsequent search does not keep Bugaenko's side-wall Gram matrix.
Start with a nine-node chain, assign its eight edges labels from
`{3,4,5,6}`, attach two extra nodes to any of the 45 unordered position
pairs on the chain, assign their two attachment edges labels from the same
alphabet, and join the extra nodes by a dotted edge. There are

\[
 45\cdot4^{10}=47,185,920
\]

labelled templates before rank and signature restrictions.

Let `B` be the chain's `9 x 9` Gram matrix. Write the attachment columns
as `u=-a*e_p` and `v=-b*e_q`, where `a=2*cos(pi/l)` and similarly for `b`.
For nonsingular `B`, set

\[
 s_0=2-u^tB^{-1}u,\quad
 s_1=2-v^tB^{-1}v,\quad
 k=u^tB^{-1}v.
\]

The full 11-node Gram has rank ten precisely when the two-node Schur
complement has rank one. The two candidate dotted weights are

\[
 t=-k\pm\sqrt{s_0s_1}.                                \tag{2}
\]

The search requires `t>2`. It treats both ways to obtain signature
`(9,1,1)`: a Lorentzian `B` with a positive rank-one Schur complement,
or positive-definite `B` with a negative rank-one Schur complement.
Singular bases and zero diagonal Schur cases are outside this selected
construction and are not claimed impossible.

For a nondegenerate kernel vector `z`, a cone vertex is represented by a
wall-pairing vector supported on two indices `i,j` with opposite signs
of `z_i,z_j`. Its two nonzero entries can be taken as
`-|z_j|,-|z_i|`. These give all rays of
`{y<=0 : z.y=0}`. If `G^+` is the inverse of the full Gram on its image,
the vertex square is `y^t G^+ y`.

For every spacelike vertex, the search proposes the polar truncation
plane, whose normal is that vertex. Its pairings with old walls are
already the entries of `y`; incident old walls are automatically
orthogonal. It then tests the remaining old-wall angles and all
polar-polar angles. The numerical polar-angle test uses the nearest
integer Coxeter label and has no imposed maximum label. This does not
provide a theoretical bound on the original diagram's label alphabet.

**Original walls remain facets.** Let `w_i>0` be a Perron eigenvector for
the negative eigenvalue of the connected old Gram matrix, and put
`h=sum w_i*r_i`. Then `h` is timelike and satisfies all old inequalities
strictly. Each polar normal `v` is a feasible old cone ray, so
`(v,r_i)<=0` and `(h,v)=sum w_i*(r_i,v)<0`. Hence `h` also satisfies
every new polar inequality strictly. For each old wall, the projection

\[
 h_i=h-\frac{(h,r_i)}2r_i
\]

has square `h^2-(h,r_i)^2/2<0`, lies on that wall, and satisfies every
other old or polar inequality strictly: its pairing with a normal `s`
is `(h,s)-(h,r_i)*(r_i,s)/2 <= (h,s)<0`. Thus every original wall has
a relative-interior hyperbolic point after all polar cuts and remains
a facet. This argument does not require polar-polar acuteness.

**Remaining gate restriction.** The scan also treats every proposed
polar wall as a retained facet. It does not test whether some polar
inequalities become redundant after all cuts are added. A bad
polar/old pairing rejects that retained wall list; it does not exclude
a Coxeter polytope obtained from the same intersection after removing
redundant polar inequalities. This is an additional possible source of
false negatives in the exploratory search.

The first broad run gives:

| stage | count |
|---|---:|
| positive-definite or Lorentzian nonsingular chain bases | 2,609 |
| rank-ten Lorentz skeletons satisfying (2) and `t>2` | 113,099 |
| parked because the kernel has a zero entry | 2,465 |
| parked because at least one ray is numerically ideal | 43,290 |
| remaining skeletons rejected by a polar/old-wall angle | 67,344 |
| candidates passing the polar-angle gate | 0 |

The three outcome counts sum to 113,099. In particular, no compact
candidate proceeds to an exact Gram or volume certification. The separate
exact facet-rigidity statement above does not convert these numerical
counts into a theorem.

All rank, signature, ray, and angle decisions in this second experiment
use floating-point arithmetic with tolerance `1e-7`. Parked degenerate
and cusp-bearing systems, nonpolar truncations, other graph topologies,
other finite label alphabets, and reduced intersections after polar-cut
redundancy removal remain open. No claim of nonreflectivity
of an arithmetic host, or of nonexistence in dimension nine, follows.

Replay:

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_gram_polar_h9_20261013.py --mode complete_3456
```

Artifacts:

- [Exact rigidity checker](../scripts/fresh_bugaenko_facet_tilt_rigidity_20261013.py)
- [Exact rigidity certificate](../data/enumerations/fresh_bugaenko_facet_tilt_rigidity_20261013.json)
- [Numerical Gram and polar search](../scripts/fresh_gram_polar_h9_20261013.py)
- [Numerical transcript](../data/enumerations/fresh_gram_polar_h9_complete_3456_20261013.json)
