# Exact automaticity tests for the Coble incidence compiler

**Status (2026-08-09): a real criterion and three sharp no-go results, but
not the missing strict non-singleton L4 theorem.**

The broad Coble--Vinberg inverse theorem reduces carrier existence to the
nonemptiness of finitely many rational-surface incidence stacks.  This note
tests the most natural ways to make that nonemptiness automatic.  The main
positive statement is an exact theorem for the ruled-grid branch:

> In the complete-projectivity `P1 x P1` grid family, finite volume holds if
> and only if the four-point projectivity group is `S4`, equivalently the
> characteristic is three.  The resulting chamber is exactly the known
> Type-MII `H17/48` chamber.

The other three conclusions are negative but useful.

1. A toric surface or Looijenga pair cannot directly be a strong Coble
   carrier because it carries an effective anticanonical divisor.
2. The smooth-component plane-sextic ansatz has at most one finite-volume
   exact type: the characteristic-three Type-VI `H14` row.
3. Tumarkin's `n+3` work does not classify the relevant one-boundary
   `H9/12` root chambers.  Both known chambers have three non-simple ideal
   vertices, exactly outside Roberts' later one-non-simple-vertex census.

The frozen calculation is in
[`certify_coble_incidence_automaticity_obstruction.py`](../scripts/certify_coble_incidence_automaticity_obstruction.py),
with [JSON](../data/enumerations/coble_incidence_automaticity_obstruction.json)
and [tests](../tests/test_coble_incidence_automaticity_obstruction.py).

## 1. A uniform rational-sextic strong-Coble factory

Work over an algebraically closed field of characteristic different from two,
unless a characteristic-zero subcase is explicitly mentioned.  Let

\[
                         6=d_1+\cdots+d_b                 \tag{1.1}
\]

be any partition of six.  Choose irreducible rational nodal plane curves
`C_i` of degrees `d_i`, generally enough that distinct components meet
transversely, all pairwise intersection points are distinct, and no such
point is a self-node.  Put

\[
 \delta_i=\frac{(d_i-1)(d_i-2)}2.
\]

Blow up every self-node and every pairwise intersection.  Rational nodal
plane curves of the required degrees exist; arbitrary-characteristic Severi
geometry is treated by Christ--He--Tyomkin,
[*On the Severi problem in arbitrary characteristic*](https://arxiv.org/abs/2005.04134).
The fixed proximity-diagram spaces themselves are honest smooth schemes by
Kleiman--Piene--Tyomkin,
[*Enriques diagrams, arbitrarily near points, and Hilbert schemes*](https://arxiv.org/abs/0905.2169).

The number of blowups is

\[
\begin{aligned}
 N&=\sum_i\delta_i+\sum_{i<j}d_id_j\\
  &=\frac{\sum_i d_i^2-18+2b}{2}
    +\frac{36-\sum_i d_i^2}{2}=9+b.             \tag{1.2}
\end{aligned}
\]

If `B_i` is the strict transform of `C_i`, then

\[
\begin{aligned}
 B_i^2
 &=d_i^2-4\delta_i-d_i(6-d_i)\\
 &=-4.                                                   \tag{1.3}
\end{aligned}
\]

Every center is a double point of the reduced sextic, so

\[
       B=\sum_iB_i=6H-2\sum_{p=1}^{9+b}E_p=-2K_S.        \tag{1.4}
\]

This already proves the strong-Coble linear-system assertions, without a
dimension count.  If `D in |-2K_S|`, then

\[
                         D.B_i=(-2K_S).B_i=-4.
\]

Hence every `B_i` is a component of `D`.  Their sum already has the full
class, so `D=B`.  If `A in |-K_S|`, then `2A=B`, impossible because `B` is
reduced.  Thus

\[
                 |-K_S|=\varnothing,\qquad |-2K_S|=\{B\}.  \tag{1.5}
\]

This is a genuine automatic, non-singleton algebraic surface factory: all
eleven partitions of six occur.  It still constructs only the boundary
skeleton, not all walls of an input finite-volume Coxeter chamber.

### AG pushback before any wall search

An exact member of `D_CobVin` has finite root-reflection index and therefore
finite automorphism group.  Kondō's exhaustive classification then requires

\[
                         b\in\{1,2,5,8\}.
\]

Since a plane sextic has at most six components, the eleven partitions
collapse to five:

```text
[6], [5,1], [4,2], [3,3], [2,1,1,1,1].
```

Their possible hyperbolic dimensions are only `H10`, `H11`, and `H14`.
This is already a finite search denominator derived from the plane geometry,
not from a list of the nine Coble answers.

Every branch is noncompact.  After the first point blowup, the pencil of
lines through that point is base-point-free.  Its pullback through every
later blowup is a primitive nef isotropic class.

## 2. Why smooth plane components cannot give the strict class

If every plane component is smooth and rational, then every `d_i` is one or
two.  The only partitions are

```text
[2,2,2], [2,2,1,1], [2,1,1,1,1], [1,1,1,1,1,1],
```

with boundary counts `3,4,5,6`.  Kondō's finite-automorphism theorem leaves
only `b=5`.  His classification then forces characteristic three and Type
VI.  Therefore:

> **Smooth-sextic no-go.**  Any finite-volume exact Coble--Vinberg carrier
> whose plane branch has a reduced sextic with smooth rational components
> is the Type-VI `H14` class.  In particular this natural ansatz has at most
> one abstract finite-volume return and cannot be the required
> non-singleton denominator.

Allowing tangencies does not alter this conclusion.  It replaces an ordinary
intersection by a multiplicity-two proximity chain and can force resolution
`(-2)`-curves, but it changes neither the degree partition nor `b`.  The known
Type-VI conic plus four tangent lines is the surviving special point.

There is also a useful expected-dimension obstruction to obtaining the rest
of a finite-index root stock generically.  A degree-`d` plane curve has
projective parameter dimension

\[
                         M_d=\frac{d(d+3)}2.
\]

At ordinary simple base points, a strict transform of square `-2` must pass
through `d^2+2` points.  The excess number of conditions is

\[
        d^2+2-M_d=\frac{d^2-3d+4}{2}>0                 \tag{2.1}
\]

for every integer `d`.  Thus no additional `(-2)`-curve is forced on the
general ordinary-point locus.  Finite-volume wall stocks necessarily live on
special tangency, collinearity, finite-field, or other determinantal strata.
This explains why the soft strong-Coble factory does not by itself solve the
polytope problem.

## 3. The complete ruled-grid family has an exact iff

Suppose a branch contracts to `P1 x P1`, and the boundary images are `r`
fibers of the first ruling and `s` fibers of the second.  Blow up their `rs`
grid points.  The boundary equality is

\[
 rF_1+sF_2-2\sum E_{ij}
   =-2K_S=4F_1+4F_2-2\sum E_{ij}.                       \tag{3.1}
\]

It forces

\[
                        r=s=4.                          \tag{3.2}
\]

So the anti-bicanonical condition itself fixes `b=8`, Picard rank 18, and
hyperbolic dimension 17.  There is no higher-`q` finite-field ladder hidden
in this construction.

Let the ordered grid coordinates be `a_1,...,a_4` and `b_1,...,b_4`.  A
smooth `(1,1)` curve through the matching

\[
                 (a_i,b_{\sigma(i)}),\quad 1\le i\le4,
\]

is the graph of the unique projectivity inducing `sigma`, if that
projectivity exists.  Modulo `PGL2`, an ordered four-point set is

\[
                         (\infty,0,1,\lambda).
\]

The exact cross-ratio stabilizers in characteristic different from two are:

| stratum | equation | complete projectivity group | graph roots |
|---|---|---:|---:|
| generic | none | `V4` | 4 |
| harmonic | `lambda in {-1,2,1/2}` | `D4` | 8 |
| equianharmonic | `lambda^2-lambda+1=0` | `A4` | 12 |
| characteristic three | `lambda=-1=2` | `S4` | 24 |

For each row, take the eight boundary walls, the sixteen grid exceptional
walls, and exactly the graph walls in the complete projectivity group.  The
same exact `cddexec_gmp` path used for the MII theorem gives:

| group | walls | finite rays | ideal rays | hyperideal rays | volume |
|---|---:|---:|---:|---:|---|
| `V4` | 28 | 2,000 | 1,178 | 2,500 | infinite |
| `D4` | 32 | 10,024 | 2,380 | 3,344 | infinite |
| `A4` | 36 | 5,360 | 3,002 | 5,340 | infinite |
| `S4` | 48 | 190,808 | 5,244 | 0 | finite, noncompact |

All retained walls remain irredundant: the full MII archive gives a strict
relative-interior witness for every wall, and deleting other inequalities
preserves that witness.  Every wall system spans the full rank 18.  The
hyperideal rays are therefore decisive, not an artifact of redundant
inequalities.

This proves the promised criterion:

\[
\boxed{
 \text{complete grid cone has finite volume}
 \iff \operatorname{Aut}_{\mathbf P^1}(\{a_i\})\cong S_4
 \iff \operatorname{char}k=3.
}                                                        \tag{3.3}
\]

The sole survivor is the known Type-MII `H17/48` chamber.  This is useful
theoretical leverage and an independent explanation of the characteristic,
but it is one survivor, so it does not satisfy strict non-singleton `Car`.

## 4. Tumarkin does not close the one-boundary denominator

The folded one-boundary Coble--Mukai chamber lies in `H9` and has twelve
root facets, so it is tempting to invoke an `n+3` classification.  The
published scope is narrower:

- Tumarkin's
  [noncompact `n+3` paper](https://arxiv.org/abs/math/0311272) proves
  `n<=16`, gives the unique `H16` example, and treats the pyramid lane; it
  does not give a complete noncompact list in all dimensions.
- His [compact paper](https://arxiv.org/abs/math/0406226) is exhaustive only
  for compact polytopes.
- Roberts'
  [later classification](https://arxiv.org/abs/1511.08451) treats exactly
  one non-simple vertex and explicitly leaves both the simple case and the
  nonpyramidal case with more than one non-simple vertex.

An exact cdd computation on the two known twelve-root chambers gives:

| chamber | finite vertices | simple ideal vertices | non-simple ideal vertices | active facets at a non-simple vertex |
|---|---:|---:|---:|---:|
| Type I | 14 | 6 | 3 | 10 |
| Type II | 38 | 14 | 3 | 10 |

The ambient dimension is nine, so ten active facets mean a non-simple
vertex.  Both examples have three such vertices.  They lie precisely in the
unclassified Roberts stratum.  Consequently there is no published
exhaustive list that can currently prove

\[
                  D_{1,\mathrm{split}}=R_{1,\mathrm{split}}.
\]

Calling Tumarkin's dimension bound a classification here would be a serious
error.

There is, however, a sharply bounded next computation: enumerate only the
`H9/12` Gale diagrams whose Coxeter Gram is simply laced at intersecting
facets, integral on `E10`, of rank ten and corank two, and has the required
three non-simple cusps.  Solve the determinant equations for the integral
dotted-edge entries, then apply the characteristic-wall and complete-cusp
filters of `D_1_split`.  This is much smaller than completing the general
noncompact `n+3` classification, but it has not yet been run.

## 5. Why toric and Looijenga pairs do not directly solve it

A toric surface has its toric boundary in `|-K|`.  A Looijenga pair is, by
definition, a rational surface together with an effective anticanonical
cycle.  Boundary blowups preserve such an effective anticanonical divisor.
Gross--Hacking--Keel's toric-model machinery is powerful precisely in this
log Calabi--Yau setting; see
[*Mirror symmetry for log Calabi--Yau surfaces I*](https://arxiv.org/abs/1106.4977).

A strong Coble surface instead satisfies `|-K|=empty` and has a reduced,
disconnected anti-bicanonical boundary.  Hence neither a toric surface nor a
Looijenga pair can be a direct output of the Coble incidence compiler.  A
quotient, double cover, or degeneration could still connect the theories,
but that extra operation must prove the strong-Coble and exact-wall gates; it
cannot be suppressed as formal toric realization.

## 6. A genuine global-realizability obstruction

The all-contraction compiler contains ordinary point-line incidence schemes
as special cases.  Rank-three matroid realization spaces satisfy
scheme-theoretic Mnëv--Sturmfels universality; see Lee--Vakil,
[*Mnëv--Sturmfels universality for schemes*](https://arxiv.org/abs/1202.3934),
and Vakil's
[*Murphy's Law in algebraic geometry*](https://arxiv.org/abs/math/0411469).
Thus incidence spaces can have arbitrary finite-type singularity behavior.

This does **not** make the compiler undecidable: each fixed branch is still
a finite elimination problem.  It does prove that local intersection
numbers, an expected-dimension count, or the mere representability of a
fixed Enriques diagram cannot force nonemptiness across arbitrary branches.
A strict theorem needs one of:

1. a globally realizable structural matroid/contact class;
2. a K3/period or fibration theorem that forces the required curves; or
3. an exhaustive restricted Coxeter census whose every survivor has an
   explicit incidence realization.

## 7. Exact next finite search

The grid lane is closed: it has one finite-volume survivor.  The two live
computations are now precise.

1. **`H9/12` one-boundary Gale census.**  Enumerate the restricted integral
   `E10` stratum described in Section 4.  This is the only current route to
   proving or refuting class-wide essential surjectivity for `D_1_split`.
2. **Plane rational-sextic contact census.**  Start with the five
   Kondo-surviving degree partitions

   ```text
   [6], [5,1], [4,2], [3,3], [2,1,1,1,1].
   ```

   Enumerate multiplicity-two contact/proximity diagrams, generate all
   resolution-forced negative curves, and then enumerate low-degree special
   curves using (2.1) as the first bound.  Solve each saturated incidence
   ideal before applying the exact Coble--Vinberg return theorem.

The second search is the first genuinely finite version of the proposed
“design the surface from the combinatorics” fishing expedition.  It covers
the possible `H10`, `H11`, and `H14` plane branches and excludes compact
output from the outset.

### What the `H9/12` census would have to prove

The strict completion must use the **whole** independently defined
`D_1_split` denominator.  Adding “three non-simple cusps” merely because both
known answers have three would be answer-selected unless that count is first
derived from the existing raw axioms.

An exhaustive full-denominator census would complete all five arrows only if
it proves that every integral Gram passing the raw volume, `E10`,
characteristic-wall, and complete-cusp conditions is isometric to one of the
two archived Grams.  Then:

```text
C:        D_1_split is independent and non-singleton;
Car:      each of its now-exhausted two rows has the explicit K3/Coble branch;
Return:   the existing one-boundary theorem returns literally the same P;
AG:       Kondō applies to every member, not merely a positive locus;
Pushback: every raw input is one of exactly two noncompact H10 polytopes.
```

That would be a strict L4 criterion with two positive survivors.  It would
not improve the global dimension bound because the denominator already fixes
dimension ten.  The exact missing arrow is a proof that the integral
Gale/Gram enumeration terminates: dotted-edge inner products and subgroup
index are not presently bounded by the twelve-facet hypothesis.  If the
census instead produces an additional raw row, that row must receive a
carrier or be excluded by a natural condition declared independently of
compiler nonemptiness; deleting it after a failed carrier test is forbidden.

## 8. Strict score

| gate | result |
|---|---|
| independent structural class | pass for sextic and grid skeletons |
| carrier for every boundary skeleton | pass |
| all listed walls effective for every finite-volume input | fail in the plane lane; pass only for the single `S4` grid survivor |
| exact same-`P` return | inherited whenever a full incidence branch is realized |
| AG pushback | `b in {1,2,5}` and dimensions `H10,H11,H14` for plane branches; grid finite volume iff MII/char 3 |
| strict non-singleton L4 | **not achieved** |

Replay the cheap certificate and tests:

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_incidence_automaticity_obstruction.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_incidence_automaticity_obstruction.py
```

Replay both exact cdd layers, including the 196,052-ray full grid cone:

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_incidence_automaticity_obstruction.py --recompute-cm-rays --recompute-cones
```
