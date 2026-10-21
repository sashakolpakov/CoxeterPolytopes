# A finite maximal-rank `P -> OG10 -> same P` compiler

## Outcome

Fix a real-quadratic integral realization of a **compact** Coxeter polytope
in `H10`.  If its restriction-of-scalars divisor lattice has rank `22`, then
the entire OG10 carrier question is a finite exact computation:

1. finitely many positive binary complements and discriminant gluings;
2. finitely many oriented rank-two periods;
3. the four published OG10 wall types; and
4. finitely many actual wall vectors whose mirrors can meet the compact input
   chamber.

The fourth assertion is the new closure point.  A fixed OG10 wall square has
infinitely many lattice representatives in general, so a finite list of wall
**types** was not previously a terminating compiler.  Compactness supplies a
positive-definite ellipsoid containing every wall vector relevant to the
input polytope.  Complete MBM/hidden-wall promotion is therefore finite.

This gives a literal necessary-and-sufficient

```text
fixed compact quadratic H10 datum P
  -> finite OG10 embedding/period/MBM package list
  -> exact test whether the closure of a marked OG10 Kähler slice
     projectivizes to the same closed polytope P.
```

It does not enumerate all quadratic `H10` Coxeter data, and it does not prove
that every compact arithmetic `H10` polytope must have an OG10 carrier.  It
therefore does not yet lower the universal carrier ceiling from `10` to `9`.
The gate-by-gate independent check is the
[hostile audit](og10_h10_maximal_return_compiler_hostile_audit.md).

There is, however, a clean uniform branch elimination.  If the quadratic
action is required to extend from `NS_Q` to a unital self-adjoint Hodge action
on all of `H2_Q`, the period plane needs two transcendental dimensions over
the field.  Maximal OG10 `H10` leaves only one.  Thus the **coherent full-`H2`**
branch is empty and has ceiling `H9`; see the
[coherent real-multiplication bound](coherent_full_h2_real_multiplication_bound.md).

## 1. The surface-free input class

Let `C_OG10^max,comp` consist of marked arithmetic data

\[
\mathcal D=(K,\mathcal O,M,h,\alpha;P,r_1,\ldots,r_f)
\]

with the following properties, stated before a variety is chosen.

- `K` is a real quadratic field, `O` an order, and `M` a projective
  `O`-module of rank `11`.
- The two real signatures of `h` are `(10,1)` and `(11,0)`.
- The trace form

  \[
  B(x,y)=\operatorname{Tr}_{K/\mathbf Q}(\alpha h(x,y))
  \]

  is even integral of signature `(1,21)`; write `N=(M,B)`.
- The distinguished eigenspace `E` has signature `(1,10)`, and the oriented
  halfspaces of the displayed field normals cut out an exactly certified
  compact Coxeter polytope `P` in `H(E)`.
- The normalized Coxeter Gram, irredundancy, compactness, and preservation of
  the marked `O`-module by every declared facet reflection are part of the
  input certificate.

This is a class of lattice-marked polytopes, not a class of surfaces.  It
deliberately fixes the integral realization.  Starting from an abstract
diagram requires a preceding exhaustive normal-scaling, order-module, and
saturation compiler; failure for one such realization is not failure of the
abstract polytope.

## 2. Finite primitive-embedding compiler

Put

\[
\Lambda=\Lambda_{OG10}
 =3U\oplus2E_8(-1)\oplus A_2(-1),
\qquad |\det\Lambda|=3.
\]

For a primitive embedding `j:N->Lambda`, its orthogonal complement `T` is a
positive even lattice of rank two.  Let

\[
H=\Lambda/(N\oplus T)
\]

be the discriminant gluing graph and put `k=|H|`.  Primitivity makes both
projections

\[
H\longrightarrow A_N,\qquad H\longrightarrow A_T
\]

injective.  The determinant formula is

\[
                   \det T={3k^2\over|\det N|}.          \tag{2.1}
\]

Consequently `k` divides both `|det N|` and `det T`.  There are only finitely
many possible `k`.  For each surviving determinant, Gauss reduction gives a
finite list of positive even binary forms

\[
T=\begin{pmatrix}2a&b\\b&2c\end{pmatrix},
\qquad 0\leq b\leq a\leq c,
\qquad4ac-b^2=\det T.                                  \tag{2.2}
\]

For each `T`, enumerate the finite isotropic subgroups of `A_N+A_T` with
injective projections.  Retain precisely those for which the residual
discriminant form is the form of `A2(-1)` and the glued even lattice is
isometric, with the required orientation, to `Lambda`.  Quotienting by
`O(N) x O(T)` gives the **unmarked-domain** primitive embedding orbits.  For
the present input, `N`, its `O`-action, and its facet data are marked, so one
must not quotient by an arbitrary element of `O(N)`: either quotient only by
`O(T)` (and by the stabilizer of all marked data in `O(N)`) or, more simply,
keep every graph.  The latter list is still finite and is the convention of
the exact marked compiler.

### The corank-two local filter

There is a cheap obstruction before discriminant-form enumeration:

\[
 \boxed{\ell_p(A_N)\leq2\quad(p\ne3),\qquad
        \ell_3(A_N)\leq3.}                              \tag{2.3}
\]

To prove it, reduce the primitive inclusion `N->Lambda` modulo `p`.  Its
image has codimension two.  An element of the radical of the restricted form
on `N/pN` is determined by its two pairings with a complementary vector
space, up to an element of `Rad(Lambda/pLambda)`.  Hence

\[
\dim\operatorname{Rad}(N/pN)
 \leq2+\dim\operatorname{Rad}(\Lambda/p\Lambda).
\]

The last dimension is zero away from `3` and one at `3`; the left side is
the `p`-length of `A_N`.  This proves (2.3).

If `3` does not divide `det N`, the conclusion is sharper.  At every prime
`p` dividing `D=|det N|`, put `a=v_p(D)` and `b=v_p(k)`.  Equation (2.1)
gives `v_p(det T)=2b-a`.  Injectivity into `A_N` gives `b<=a`, while
injectivity into `A_T` gives `b<=2b-a`, hence `b>=a`.  Thus `b=a` at every
such prime.  Since `3` does not divide `D`, no other prime can divide `k`, so
the gluing consumes the full discriminant group and

\[
               k=|\det N|,\qquad \det T=3|\det N|.      \tag{2.4}
\]

Moreover no primitive class in `N` can have ambient divisibility three.  If
one did, its reduction would lie in the one-dimensional radical of
`Lambda/3Lambda`; pairing with `N/3N` and the nondegeneracy of `N mod 3`
would make the class divisible by three already in `N`, a contradiction.
Thus only the `(-2,1)` and `(-4,1)` OG10 wall types survive in the
three-coprime branch.

## 3. Periods are finite and exact at maximal Picard rank

For every accepted embedding, `T` is positive definite of rank two.  Its two
orientations give the two marked period lines

\[
[\sigma]\subset T_\mathbf C,
\qquad(\sigma,\sigma)=0,
\qquad(\sigma,\bar\sigma)>0.
\]

The real and imaginary parts span all of `T_R`.  Hence the integral `(1,1)`
lattice is exactly `N`; there is no Noether--Lefschetz avoidance left to do.
Since `N` contains a positive integral class, period surjectivity and the
projectivity criterion give a projective marked OG10 variety.  Orientation
and monodromy-component choices are finite and are retained in the package.

Thus, at `H10`, the continuous moduli problem collapses to a zero-dimensional
marked period problem.  All remaining information is lattice and chamber
information.

## 4. Compactness makes the complete wall search finite

Mongardi--Onorati's OG10 wall-divisor classification gives exactly

\[
(-2,1),\quad(-4,1),\quad(-6,3),\quad(-24,3),             \tag{4.1}
\]

where the second entry is divisibility in `Lambda`, not merely in `N`.
The first and third types are the prime-exceptional/movable walls; all four
types cut Kähler chambers.

For an embedding package `j`, first take any exact strict interior point
`p_0 in P`.  Let `R` be a certified radius of a hyperbolic ball about `p_0`
containing the compact polytope.  It can be computed exactly from the finite
vertex list:

\[
 \cosh^2R=max_{v\text{ a vertex of }P}
 { (p_0,v)^2\over(p_0,p_0)(v,v)}.                        \tag{4.2}
\]

Let `delta in N` have one of the four squares `s` in (4.1), the required
ambient divisibility, and suppose the restriction of `delta^perp` meets
`P`.  Decompose `delta=delta_E+delta_-` according to the hyperbolic
eigenspace and its negative-definite conjugate.  A meeting mirror has
`delta_E^2<0`, and

\[
0<-\delta_E^2\leq-s.
\]

The distance formula from `p_0` to a hyperplane gives

\[
 { (p_0,\delta)^2\over(p_0,p_0)(-\delta_E^2)}
     \leq\sinh^2R.
\]

Now define the positive-definite form on `N_R`

\[
\widehat Q_{p_0}(\delta)=
 {2(p_0,\delta)^2\over(p_0,p_0)}-(\delta,\delta).
\]

Indeed, if `x=a p_0+z` with `z perpendicular p_0`, then
`Qhat_(p_0)(x)=a^2(p_0,p_0)-(z,z)>0` for `x!=0`.  Using
`delta^2=s` in the distance bound gives (4.3).

Every relevant wall vector satisfies the explicit ellipsoid bound

\[
 \boxed{\widehat Q_{p_0}(\delta)
       \leq(-s)(2\cosh^2R-1).}                          \tag{4.3}
\]

A lattice has only finitely many points in a bounded positive-definite
ellipsoid.  Exact algebraic-number comparison supplies a rational lower
eigenvalue bound and hence a finite coordinate box; Fincke--Pohst or an exact
branch-and-bound enumeration then lists every vector.  Equation (4.3) is
therefore a terminating **complete-wall** algorithm, not merely a finite
wall-type alphabet.

Regularity of the point later used to orient the walls is effective, rather
than an extra genericity assumption.  The projective `K`-module structure on
`N_Q` is faithful, so no nonzero rational vector lies wholly in just one
irrational real eigenspace.  In particular, no nonzero integral wall
functional vanishes identically on `E`.  After (4.3) has listed the finitely
many walls meeting `P`, choose an exact point `p` in the strict interior of
`P` off their finitely many mirrors.  If a broader input category does allow
a wall functional to vanish identically on `E`, that package is rejected:
its Kähler slice in `E` is empty.

This argument uses compactness essentially.  For a finite-volume noncompact
polytope, `R` is infinite and (4.3) gives no enumeration.  Compact and
finite-volume-noncompact inputs must remain separate compiler classes.

## 5. Necessary-and-sufficient return criterion

For each primitive embedding/period package, enumerate the finite set from
(4.3), retain the vectors with one of the exact ambient types (4.1), choose
the regular point `p` just described, and orient every restricted functional
to be positive at `p`.  The package returns the same closed chamber as `P`
if and only if:

1. every irredundant input facet is the restriction of at least one retained
   OG10 wall divisor; and
2. every retained OG10 wall functional is nonnegative on all of `P`.

Both conditions are finite exact linear/conic tests.  The first gives one
inclusion between the proposed cone and the complete wall chamber; the
second gives the reverse inclusion.  Hence they are necessary and
sufficient for equality.  Any wall omitted by (4.3) has no zero on the
connected set `P`; after orientation at `p`, it is positive throughout `P`
and imposes no further cut.

Here equality is deliberately a closure statement:

\[
 \mathbf P\bigl(\overline{\mathcal K_X}\cap E^+\bigr)=P,
 \qquad
 \mathbf P\bigl(\mathcal K_X\cap E^+\bigr)=\operatorname{int}P.
\]

It is not literally correct to identify the open Kähler cone with the closed
polytope.

The marked-model step uses more than the one-way assertion that a Kähler
cone is a wall-complement component.  Mongardi--Onorati's Theorem 5.5 says
that every OG10 Kähler cone is a component for the four numerical wall
types.  Conversely, any component of that complement contains a very
general class; Markman's Lemma 5.11 puts such a class in a Kähler-type
chamber, and Theorem 5.16 realizes every Kähler-type chamber by a marked
model in the same period fibre.  Since a Kähler-type chamber is itself a
wall-complement component, the two components coincide.  Onorati's
`Mon^2=O^+` theorem accounts for the retained orientation/marked-component
choices.  Thus a same-period projective marked model has the closure above.
The input Gram, labels, orientation, and integral field reflections were
retained from the polytope datum, so the recovered object is the same
labelled Coxeter polytope, not merely an isomorphic nerve or a larger divisor
cone.

For a prime-exceptional or movable-**closure** return, run the identical
compiler with only `(-2,1),(-6,3)`.  This computes the closure of the
fundamental exceptional/birational-Kähler chamber (and hence the closed
movable cone, under that convention); its interior is the corresponding
open chamber.  A claim at this level must not be silently promoted to a
Kähler/nef statement, where the `-4` and `-24` walls return.

This proves:

> **Maximal compact OG10 return theorem.**  Strong exact-Picard OG10 Kähler
> (respectively prime-exceptional/movable-closure) carrier existence is
> decidable by a finite exact algorithm for every datum in
> `C_OG10^max,comp`.  Every allowed primitive
> embedding, period orientation, OG10 wall, hidden cut, and same-chamber
> alternative occurs in its output.

## 6. Frozen rank-22 calibration

The executable front end replays the exact rank-22 lattices already in the
repository.

| lattice | `A_N` | determinant/glue lanes before discriminant anti-isometry | selected complement | current chamber result |
|---|---|---|---|---|
| golden orthogonal `H10` | `Z/2 + Z/22` | `k=44`, `det T=132`, four reduced binary forms | `<2>+<66>` | full norm-two movable-chamber closure has **infinite volume** |
| A2-coupled golden `H10` | `Z/2 + Z/10` | `k=20`, `det T=60`, four forms | `<2>+<30>` | every norm-two chamber has **infinite volume** via the complete fixed `H9` slice |
| golden norm-six `H10` | `(Z/6)^2` | `(k,det T)=(12,12),(36,108)`, respectively two and five forms | `<2>+<6>` in the partial `3`-glue lane | exact projective OG10 carrier with a genuine `(-6,3)` wall; complete chamber and volume open |
| Borcherds `D4 H21/210` calibration | `(Z/2)^2` | `k=4`, `det T=12`, two forms | `A2(2)` | weak OG10 carrier exists, but a hidden `(-4,1)` wall crosses the strict interior |

The last row is not a quadratic `H10` eigencarrier; it is included because it
is the existing exact hostile calibration of the same rank-22
complement/glue/hidden-wall compiler.

The binary-form counts use full `GL(2,Z)` lattice equivalence, not oriented
proper equivalence: `b` and `-b` are identified.  The executable audit also
tests every pair in each displayed reduced list by an exact finite integral
isometry search, so reduction-boundary duplicates are absent.  These are
only complement candidates **before** discriminant anti-isometry and gluing;
the table records one already constructed complement in each calibration and
does not claim that every determinant lane has a unique surviving embedding.

No compact `H10` survivor is present in this table.  The first two are
proved infinite-volume, the third has only a positive carrier-and-wall
certificate, and the fourth fails complete-cone promotion.

## 7. What would complete the dimension drop

There are now two sharply separated routes.

### Coherent route -- complete

If the allowed carrier category requires the field action to extend to a
unital self-adjoint action on the full K3-type `H2` Hodge structure, the
period-slot theorem eliminates OG10 `H10` uniformly.  Any jointly exhaustive
compiler into coherent K3/OG10 carriers returns `n<=9` (and the K3 coherent
ceiling is `H8`).

### NS-only route -- still open

The known rank-22 examples prove that an action on `NS` alone survives the
period and primitive-embedding gates.  To lower the ceiling for this broader
category one must do one of the following.

1. Prove a **finite P-first denominator theorem** for the intended compact
   arithmetic `H10` class: finitely many diagrams, normal scalings,
   saturations, and trace lattices.  Run the compiler above on every output;
   if every package fails the finite hidden-wall/compactness test, the branch
   is exhausted.
2. Prove a uniform OG10 theorem that every complete quadratic rank-22 MBM
   chamber has an infinite stabilizer or another noncompactness witness.
   The golden infinite-dihedral argument is one instance, not yet a general
   theorem.
3. Find a package passing the finite return test.  That would be an actual
   compact `H10` Coxeter polytope with an exact OG10 realization and would
   move the positive boundary instead of lowering the negative one.

The compiler makes route 1 finite once the P-first candidate denominator is
finite.  It does not itself prove that denominator theorem.

## 8. Reproduction and sources

Run

```sh
PYTHONPATH=src:. python3 \
  scripts/certify_og10_h10_maximal_return_compiler.py

PYTHONPATH=src:. pytest -q \
  tests/test_og10_h10_maximal_return_compiler.py
```

The current executable audit is the arithmetic **front end**: it checks the
Smith data, local length gate, determinant lanes, and complete
reduced-binary-form counts for all four calibration lattices.  It imports the
already frozen exact Gram certificates rather than retyping them.  The
discriminant anti-isometry and ellipsoid enumeration are finite algorithms
proved above, but are not implemented end-to-end by this particular script.

The algebraic-geometric inputs are Rapagnetta's OG10 Beauville lattice
([arXiv:math/0606409](https://arxiv.org/abs/math/0606409)), Onorati's full
OG10 monodromy group
([arXiv:2002.04129](https://arxiv.org/abs/2002.04129)), Markman's
Kähler-type and exceptional chamber theory
([arXiv:1101.4606](https://arxiv.org/abs/1101.4606)), and
Mongardi--Onorati's complete OG10 wall-divisor and Kähler-cone theorem
([arXiv:2010.12511](https://arxiv.org/abs/2010.12511)).
