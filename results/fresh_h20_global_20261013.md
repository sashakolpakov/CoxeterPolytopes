# Every interior two-wall section of Borcherds' H21 chamber fails

The new construction test gives an exact obstruction, **not a new Coxeter polytope**:

> Let `P` be the archived 210-facet Borcherds chamber in H21. No hyperplane meeting the strict interior of `P`, with normal in the real span of at most two of its facet normals, cuts out a Coxeter polytope in H20.

The statement covers all **21,945 unordered wall pairs and every positive real parameter**, including irrational and transcendental parameters. All these interior sections have finite volume. Their failure is a genuine forbidden dihedral angle at an actual ridge, not a non-Coxeter entry between possibly redundant inequalities. It says nothing about sections whose normals require three or more walls, new tile assemblies, other ambient chambers, or general H20 Coxeter polytopes.

An independent verifier, using only Python's standard library and `Fraction`, passed the complete proof and rejected two deliberately corrupted certificates.

## Exact input and the construction family

The input is [the frozen Borcherds wall packet](../data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json), under `coxeter_input`. We reverse the archived quadratic-form sign, so spacelike normals have positive square. The 210 normals `r_i` comprise 42 roots of square 2 and 168 of square 4. This is Borcherds' Example 5; his chamber symmetry group has order 241,920. [Borcherds, *Automorphism groups of Lorentzian lattices*](https://math.berkeley.edu/~reb/papers/al/al.pdf).

The primitive sum `h` of the 210 normals satisfies

\[
 h^2=-14,\qquad (h,r_i)=-1\text{ for short walls},\qquad
 (h,r_i)=-8\text{ for long walls}.
\]

For each distinct pair we tested the full pencil

\[
 e(t)=r_i-t r_j,\qquad t>0,\qquad P_t=P\cap e(t)^\perp.
\]

Since all distinct ambient wall pairings are nonpositive, `e(t)` is spacelike. The section meets the strict interior for every positive `t`: let

\[
 z_i=h-\frac{(h,r_i)}{r_i^2}r_i,\quad
 z_j=h-\frac{(h,r_j)}{r_j^2}r_j,\quad
 A=-(r_i,z_j)>0,\quad B=-(r_j,z_i)>0.
\]

Each `z` is a future timelike relative-interior point of its indicated facet. The vector `A z_i+t B z_j` is future timelike, strictly satisfies every ambient inequality, and is orthogonal to `e(t)`.

Conversely, a nonzero same-sign combination of two outward wall normals pairs strictly with every strict interior point, so its perpendicular hyperplane cannot meet the interior. Thus the positive-parameter pencils cover all interior sections with normals supported on two walls.

**Finite volume is inherited here.** In the Klein model, the closure of the known finite-volume finite-facet chamber is a polytope with finitely many finite or ideal vertices. Its intersection with a hyperplane through an interior point is a full-dimensional H20 polytope. Strict convexity of the bounding sphere implies that any ideal vertex of the intersection is an inherited ideal vertex; new edge intersections lie inside the ball. Triangulating into finitely many finite/ideal H20 simplices proves finite volume. Equivalently, an inherited cusp has a bounded Euclidean section of the original cusp link. There is no assumption that the section avoids cusps.

This differs from the previous [facet gluing search](new_dimension_h20_facet_gluing_20261013.md) and [mixed tile/rotated-plane dictionaries](fresh_h20_mixed_20261013.md): these hyperplanes run through the chamber interior, and the entire positive parameter interval is checked. No facet-deletion search inside the old 580-, 767-, 1142-, or 1262-plane dictionaries is repeated.

## Complete symmetry reduction and central calibration

We lift generators of the 42-short-wall incidence graph to exact integral ambient isometries. Every lifted matrix is checked against the quadratic form, all 210 roots, and `h`. Their induced pair orbits cover every unordered pair. Completeness of the 21,945-pair coverage does not depend on trusting a claimed group order: the independent checker reconstructs the orbits from certified matrices.

The central section uses `e=(h,r_j)r_i-(h,r_i)r_j`. Its nine representatives already have exact strict timelike forbidden-ridge witnesses:

| Ambient pair, one based | Pair Gram off-diagonal | Pair count | Central forbidden squared cosine | Full positive-parameter proof |
|---|---:|---:|---:|---|
| 1, 2 (short/short) | 0 | 336 | 1/7 | Persistent triples; gcd `t²` |
| 1, 3 (short/short) | 0 | 420 | 1/7 | Persistent triples; gcd `t²` |
| 1, 20 (short/long) | 0 | 5,040 | 32/65 | Persistent triples; gcd `t²` |
| 1, 22 (short/short) | −1 | 105 | 3/11 | Persistent triples; gcd `1` |
| 1, 48 (short/long) | −2 | 2,016 | 81/163 | Persistent triples; gcd `1` |
| 20, 21 (long/long) | −2 | 3,360 | 3/10 | Persistent triples; gcd `1` |
| 20, 54 (long/long) | −4 | 3,780 | 2/7 | Actual ridge for every `t>0` |
| 20, 84 (long/long) | −6 | 6,048 | 5/18 | Actual ridge for every `t>0` |
| 20, 143 (long/long) | −8 | 840 | 1/121 | Actual ridge for `0<t≤1`, plus reciprocal symmetry |

The central projected inequalities have 209 distinct positive normal rays in the first five cases and 208 in the last four. These counts describe projected inequalities; the proof does not assume all of them define facets.

## Six elliptic pair types: persistent actual faces and incompatible angles

For each of the first six rows, the ambient pair has positive-definite Gram matrix. We select ambient positive-definite triples `(r_i,r_j,r_k)` and save an exact future timelike point orthogonal to precisely these three normals, with every other ambient wall pairing strictly negative. This is an actual codimension-three ambient face. It lies in every section `e(t)⊥`. In the section the first two projected normals satisfy `r'_i=t r'_j`, so the triple becomes an actual dihedral ridge for every positive `t`.

Write the triple Gram as

\[
 \begin{pmatrix}a&c&u\\c&b&v\\u&v&d\end{pmatrix},\qquad
 \Delta=ab-c^2>0.
\]

The squared cosine at the ridge is `N(t)/D(t)`, where

\[
 N(t)=\bigl((av-cu)+t(bu-cv)\bigr)^2,
\]

\[
 D(t)=\Delta\bigl(d(a-2ct+bt^2)-(u-tv)^2\bigr).
\]

Exact coefficient signs prove `D>0` and `0≤N/D<3/4` throughout `t>0`, with nonpositive outward normal pairing. Any Coxeter value in this interval must come from labels 2, 3, 4, or 5. Consequently a necessary condition is

\[
 N(4N-D)(2N-D)(4N-3D)(16N^2-12ND+D^2)=0.
\]

The unnecessary label-6 factor and the other conjugate of the label-5 squared cosine only enlarge the necessary zero set; they cannot cause a false exclusion. Taking the polynomial gcd over the saved persistent triples gives `t²` in the first three rows and `1` in the next three. Neither has a positive zero. No restriction on the field or algebraic degree of `t` is used.

## Three nonelliptic pair types: direct actual-ridge formulas

For these long-wall pairs a convenient strict center is

\[
 x(t)=(1+t)h+2r_i+2t r_j\in e(t)^\perp.
\]

The archive gives a polynomial timelike ridge point and all 210 ambient wall pairings for each family. Exactly the two indicated ridge inequalities vanish, and every other pairing is strictly negative over the stated domain. The two projected normals have positive-definite Gram matrix and negative mutual pairing. Thus these are actual acute dihedral angles between true facets near the certified points.

For pair **20,54**, the squared cosine is

\[
 \frac{(t+1)^2}{2(t^2+4t+2)}\in(1/4,1/2),\qquad t>0.
\]

For pair **20,84**, it is

\[
 \frac{t^2+3t+1}{2(t^2+6t+2)}\in(1/4,1/2),\qquad t>0.
\]

These open intervals contain no Coxeter squared cosine: their endpoints correspond to labels 3 and 4.

For pair **20,143**, the squared cosine is

\[
 \frac{t^4}{(t^2+8t+2)^2}\in(0,1/4),\qquad 0<t\le1,
\]

which lies strictly between labels 2 and 3. Every required strict polynomial sign is certified by nonnegative coefficients after the substitution `t=u/(1+u)`, plus an exact check at `t=1`. A saved word of five certified ambient symmetry generators swaps walls 20 and 143. It takes the section parameter `t` to `1/t`, proving the same exclusion for `t≥1`.

During exploration we also enumerated 28 positive algebraic candidates obtained from minimal-polynomial angle equations for labels 2 through 12 at the central defective pairs. They all had exact forbidden-ridge certificates; some are only conjugate roots of the desired angle equation. This bounded sample is retained for reproducibility but is superseded by the all-positive-parameter proof. The independent universal verifier does not depend on it.

## Why this mechanism was selected from the literature

Ordinary rolling alone does not supply a finite closing cover here. Alekseevski–Michor–Neretin identify the developed mirror chamber with a universal cover of the boundary cut along even seams and describe the associated deck group. Combined with the already archived nontrivial infinite-order hexagon holonomy, any finite cover retaining the same transitions still contains a positive power of that holonomy. Adding different seam transitions remains outside that argument. [*Rolling of Coxeter Polyhedra Along Mirrors*, Propositions 2.7–2.9 and Theorem 2.10](https://www.mat.univie.ac.at/~michor/coxeter.pdf).

The recent primary paper by Shimada was also checked. Its involution invariant hyperbolic lattices have ranks 18, 14, 10, or 2, hence hyperbolic dimensions 17, 13, 9, or 1; these do not provide an H20 or above-H21 seed. [Shimada, *Involutions in the affine Conway group*, Theorem 1.2 and Table 1.1, August 2026](https://arxiv.org/pdf/2608.29065).

Interior sections therefore offered a concrete different operation with exact existing input and automatic finite volume. The result now identifies a precise requirement for a further section search: its normal must involve at least three ambient walls. This report does not claim that such a wider search will succeed.

## Artifacts and replay

- [Central sections and exact symmetry matrices](../data/enumerations/fresh_h20_global_interior_sections_20261013.json); [generator/checker](../scripts/fresh_h20_global_interior_sections_20261013.py).
- [All-parameter certificates and exploratory algebraic samples](../data/enumerations/fresh_h20_global_section_parameters_20261013.json); [generator/checker](../scripts/fresh_h20_global_section_parameters_20261013.py).
- [Independent standard-library verifier](../scripts/fresh_h20_global_fraction_verify_20261013.py); [PASS record, including corruption rejection](../data/enumerations/fresh_h20_global_fraction_verify_20261013.json).

```sh
python scripts/fresh_h20_global_interior_sections_20261013.py --verify
python scripts/fresh_h20_global_section_parameters_20261013.py --verify
python scripts/fresh_h20_global_fraction_verify_20261013.py
```

The independent replay verifies the source geometry, lifted symmetries, complete pair-orbit denominator, nine central actual ridges, six universal gcd obstructions, three universal polynomial-ridge obstructions, and reciprocal symmetry. Two mutation checks reject a changed central ridge point and a changed parametric point coordinate. No numerical optimization or sampled sign decision enters the universal proof.
