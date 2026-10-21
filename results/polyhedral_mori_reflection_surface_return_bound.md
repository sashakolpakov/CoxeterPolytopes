# A surface-side reflection bound ready for `P -> surface -> P` transfer

**Status (2026-08-09).**  This note isolates a classical theorem of Nikulin
as an exact return engine for the inverse-carrier program.  It is already a
theorem about every *strong direct full-nef surface carrier*.  It becomes a
theorem about a surface-free Coxeter class only after a `P`-first compiler is
proved universal for that class.

## 1. The surface theorem

Let `X` be a smooth projective surface over `C` with Picard rank
`rho>=3` and finite polyhedral Mori cone.  Suppose every irreducible negative
curve spanning an extremal Mori ray defines an integral reflection of
`NS(X)`,

\[
 s_E(x)=x-\frac{2(x,E)}{E^2}E\in NS(X).
\]

Then

\[
 \boxed{\rho(X)\le 22.}
 \tag{1.1}
\]

Nikulin records two sharper uniform branches:

\[
 \begin{array}{c|c|c}
 \text{extremal negative curves}&\rho(X)&
 \dim \mathbf P(\operatorname{Nef}(X))\\ \hline
 E^2=-1\ \text{for all }E&\rho\le9&n\le8,\\
 E^2=-2\ \text{for all }E&\rho\le19&n\le18,\\
 \text{arbitrary integral reflecting }E&\rho\le22&n\le21.
 \end{array}
 \tag{1.2}
\]

The first line uses the classification of hyperbolic lattices reflective in
norm `-1`; the second uses the norm-`-2` reflective-lattice bound; the last is
the Vinberg--Esselmann maximal-rank bound quoted by Nikulin.  These are
finite-covolume bounds and therefore apply, a fortiori, when the nef
projectivization is compact.

Primary source: V. V. Nikulin, *A remark on algebraic surfaces with
polyhedral Mori cone*, Nagoya Math. J. 157 (2000), Example 1.3 and Example
1.4, [arXiv:math/9806047](https://arxiv.org/abs/math/9806047),
[journal DOI](https://doi.org/10.1017/S0027763000007194).

## 2. Exact polytope return statement

Let `C` be a class of marked arithmetic Coxeter data defined without naming
a surface.  Assume a finite `P`-first compiler has proved the following
universal statement for every actual `P in C`:

1. some compiled package is a smooth projective surface `X`;
2. the marked cone extracted from the package is the **entire** nef cone of
   `X`, and its normalized projectivization is exactly `P`;
3. every facet is dual to an extremal irreducible negative curve `E`; and
4. the corresponding Coxeter reflection preserves `NS(X)`.

Then (1.1) pushes back through that compiler and gives the unconditional
polytope theorem

\[
 \boxed{P\in C\Longrightarrow\dim P\le21.}
 \tag{2.1}
\]

If the compiler proves that every extremal curve has square `-2`, the returned
bound is `dim P<=18`; if every one has square `-1`, it is `dim P<=8`.
Because exact full-nef equality is assumed, the convex--Hodge theorem also
returns finite hyperbolic volume automatically, and compactness is equivalent
to absence of isotropic extremal nef rays.

For a **compact** target these rank bounds are superseded by a much sharper
arithmetic return.  The integral facet reflections generate a finite-volume
subgroup of `O^+(NS(X))`, hence a finite-index subgroup.  Compactness is
therefore equivalent to rational anisotropy of `NS(X)`.  Meyer makes every
indefinite rational lattice of rank at least five isotropic, so

\[
 \boxed{P\text{ compact}\Longrightarrow\dim P\le3.}
 \tag{2.2}
\]

This applies to every direct full-nef integral-reflection surface carrier,
independently of the self-intersections or genera of its wall curves.  See
the standalone
[`direct full-nef compactness theorem`](direct_full_nef_surface_compactness_obstruction.md).

This is the desired logical shape:

```text
surface-free P
  -> exhaustive strong full-nef carrier compiler
  -> Nikulin's AG/reflection classification bound on rho(X)
  -> dimension bound for every P in the original class.
```

The conclusions are not conditional on choosing a convenient surface *after*
seeing `P`; it is conditional only on proving the stated compiler universal.
That last step remains open for broad compact arithmetic classes.

## 3. Why Hodge index and polyhedrality alone cannot work

Nikulin's Example 1.2 supplies the decisive countercalibration.  For every
fixed `delta_E>=3` and arithmetic genus bound `p_E`, iterated blowups of a
ruled surface produce finite-polyhedral Mori cones with

\[
 \rho(X_k)=3+k
\]

unbounded.  Thus Hodge signature, finite generation of the Mori cone, bounded
negative self-intersection, and bounded curve genus do **not** by themselves
give an absolute Picard-rank or finite-volume polytope-dimension bound.  The
integral-reflection condition is the ingredient that activates the bound
(1.1).

These ruled-surface examples retain the pullback of the ruling fibre as a
nonzero nef isotropic class after every blowup.  Their nef projectivizations
are therefore noncompact.  More strongly, every smooth projective rational
surface of Picard rank at least two has a pulled-back ruling class and hence
noncompact full nef projectivization.  They do **not** disprove the
possibility that compactness plus weaker non-rational geometric hypotheses
has its own bound; they prove only that Hodge index and finite polyhedrality
without such an extra hypothesis are insufficient.

This also explains why a weak divisor-lattice carrier is too soft for the
inverse goal.  Without full cone equality, there may be additional extremal
curves not represented by Coxeter facets, and Nikulin's hypothesis cannot be
verified from the displayed wall subsystem.

## 4. Scope exclusions

The theorem does not presently bound:

- invariant or eigen-nef sections of a larger nef cone;
- a divisor subsystem which is not the complete nef/Mori cone;
- automorphic Heegner arrangements not promoted to a surface nef cone; or
- a polytope class for which some actual members might have no direct
  full-nef surface carrier.

Consequently (2.1) is not being reported as a universal bound for compact
arithmetic Coxeter polytopes.  Its role is sharper and concrete: it supplies
the surface-side classification theorem needed by the requested return
architecture.  The remaining theorem task is now cleanly separated:

```text
prove universal direct full-nef promotion for a non-singleton surface-free C,
or prove a finite exhaustive disjunction in which every non-full-nef branch
has its own comparably strong return bound.
```

## 5. Relation to the sharper rational branch

For a smooth rational surface with complete polyhedral nef cone, primitive
smooth-rational facet curves, and integral facet reflections, unimodularity
forces each wall square to be `-1` or `-2`; adjunction makes `-K_X` nef, and
`K_X^2+rho(X)=10` gives `rho<=10`, hence `dim P<=9` for finite-volume
members.  The pulled-back ruling class rules out compactness in every
positive dimension.  These independent arguments are frozen in
[`rational_surface_primitive_wall_dimension_bound.md`](rational_surface_primitive_wall_dimension_bound.md).
The Nikulin theorem is broader--it does not assume rationality or primitive
rational wall classes--and returns the finite-volume ceiling `21`; the Meyer
argument returns the compact ceiling `3` for the entire direct carrier class.
