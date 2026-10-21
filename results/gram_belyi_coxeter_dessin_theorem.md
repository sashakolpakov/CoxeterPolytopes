# A universal Gram--Belyi dessin and an exact field census at bounded degree

**Status (2026-08-10): proved general real-marked
`P -> dessin surface -> same P` association, with a canonical numerical
complexity and an exact high-dimensional lower-bound calibration.**  This
gives a literal Belyi
object for every arithmetic Coxeter chamber and proves that bounding its
degree gives a finite, effective list of chambers and Vinberg fields without
ABSW, Borel--Prasad, Prasad's volume formula, or Mostow rigidity.

It does **not** yet prove that reflection maximality bounds this degree.  The
frozen `H10` K3 redoubling family proves that the degree is necessarily
unbounded on finite-index descendants even when the dimension, rational
field, divisor lattice, root-norm alphabet and finite-angle alphabet stay
fixed.  Thus a surface theorem must genuinely control the Belyi degree; it
cannot divide by the cover degree or discard the arithmetic marking.

## 1. The invariant packet

Let `P` be an arithmetic finite-volume Coxeter chamber with ordered facets
`F_1,...,F_f` and outward primitive roots `r_i`.  Put

\[
 u_{ij}=\frac{(r_i,r_j)^2}{r_i^2r_j^2}\in\overline{\mathbf Q}_{\geq0}
 \qquad(i<j).
\tag{1.1}
\]

This is independent of every root rescaling.  It is zero for an orthogonal
pair, takes the familiar cyclotomic values for finite Coxeter labels, equals
one for a parallel pair, and is greater than one for an ultraparallel pair.
For `u_{ij}>0`, define

\[
                 \lambda_{ij}=\frac{u_{ij}}{1+u_{ij}}\in(0,1).
\tag{1.2}
\]

A normalized Belyi polynomial is a polynomial `B(x)` with

```text
B(0), B(1) in {0,1},       every finite critical value in {0,1}.
```

For an algebraic number `lambda`, let `h_Bel(lambda)` be the least degree of
a normalized Belyi polynomial sending `lambda` to zero or one.  Existence and
effectivity are the genus-zero marked-point case of Belyi's theorem; see
[Khadjavi](https://doi.org/10.1006/jnth.2001.2748) and
[Rodriguez](https://arxiv.org/abs/1104.2027).  There are only finitely many
normalized Belyi polynomials of bounded degree.

For each nonorthogonal pair retain the finite set of all minimal polynomials,
together with the four ordered source points

\[
                         (0,1,\infty,\lambda_{ij}).
\tag{1.3}
\]

Replace `B` by `1-B` when needed so that `B(lambda_{ij})=1`.  For an
orthogonal pair retain instead the discrete orthogonality bit and the fixed
degree-two packet

\[
 \lambda_{ij}=\tfrac12,\qquad B_{ij}(x)=4x(1-x).
\tag{1.4}
\]

The bit distinguishes this packet from a nonorthogonal pair having
`u_{ij}=1`, which also gives `lambda=1/2`.

The marked packet also retains the distinguished real embedding of its
coefficient algebra, the acute convention for the outward unit normals, and
the chosen Lorentzian time sheet.  These real data are essential: the bare
object over `Qbar` remembers the Galois orbit of the squared Gram entries,
not which conjugate is the original hyperbolic chamber.

Define the **Gram--Belyi degree**

\[
 \beta_{\rm GB}(P)=\sum_{i<j}d_{ij},\qquad
 d_{ij}=\begin{cases}
 2,&u_{ij}=0,\\
 h_{\rm Bel}(\lambda_{ij}),&u_{ij}>0.
 \end{cases}
\tag{1.5}
\]

This is invariant under ordered arithmetic isometry; minimizing over facet
orders, or simply quotienting the finite packet by `S_f`, gives the unmarked
invariant.  It is not a covolume or a renamed facet count: its individual
coordinates measure the arithmetic complexity of the normalized dotted Gram
entries.  Nevertheless it automatically sees the number of pairs.  A
degree-one normalized Belyi polynomial is `x` or `1-x`, while every
`lambda_{ij}` in (1.2) lies strictly between zero and one.  Hence every
`d_{ij}` is at least two and

\[
                     \boxed{\beta_{\rm GB}(P)\geq f(f-1).}
\tag{1.6}
\]

The standard finite labels have now been calibrated exactly except for one
honest interval:

    m=2 assigned bit-packet: 2;  m=3: 5;  m=4: 3;
    m=5: at least 29 (exact minimum open);  m=6: 7;
    parallel: 2;  dotted u=2: 3.

For \(m=5\),
\(\lambda=(7+2\sqrt5)/29\) satisfies
\(29\lambda^2-14\lambda+1=0\); its \(29\)-adic Newton polygon gives a
conjugate with valuation \(-1\), so Rodriguez forces degree at least \(29\).
A completely explicit composed Belyi polynomial supplies a finite upper
witness.  Thus, if \(N_m\) counts facet pairs with finite label \(m\),

\[
\beta_{\rm GB}(P)\ge
2\binom f2+3N_3+N_4+27N_5+5N_6,
\tag{1.7}
\]

before extra dotted-coordinate contributions.  See the
[standard-packet calibration](gram_belyi_standard_packet_calibration.md).

## 2. A genuine Belyi surface, not merely a name

The packets assemble into a finite algebraic carrier.  Take one normalized
copy `C_{ij}=P1` for every facet pair.  Add a backbone `C_0=P1` with map
`x -> x^M`, where `M=binom(f,2)`.  Its fibre over one consists of `M`
distinct roots of unity.  Glue the marked point `lambda_{ij}` on each
`C_{ij}` to a different point of that fibre.  Since both branches map to
one, the componentwise maps give a finite morphism

\[
                     \beta_P:C_P\longrightarrow\mathbf P^1
\tag{2.1}
\]

from a connected nodal curve, branched only above `0,1,infinity`.  The
normalization, conductor, facet-pair labels, orthogonality bits and ordered
points `(0,1,infinity,lambda)`, together with the distinguished real-place,
acute-sign, and time-sheet marking, are retained.  The total cover degree is

\[
             \deg\beta_P=M+\beta_{\rm GB}(P).
\tag{2.2}
\]

Finally put

\[
 S_P=C_P\times\mathbf P^1,
 \qquad (\beta_P,\mathrm{id}):S_P\longrightarrow
       \mathbf P^1\times\mathbf P^1.
\tag{2.3}
\]

This is a projective normal-crossing Belyi surface with a finite map whose
branch lies over the three vertical rulings, equipped with the labelled
normalization/conductor and marked vertical curves.  Minimal polynomials need
not be unique, so the canonical output is a **finite exhaustive stack** of
such marked surfaces, not a falsely canonical single surface.
The stack is naturally defined over a finite extension containing the Belyi,
marking, and backbone gluing data; that field is not asserted to equal the
Vinberg field.

This complexity is an ordinary polarized-surface invariant.  Let

\[
 \Phi_P=(\beta_P,\mathrm{id}),\qquad
 A_P=\Phi_P^*\mathcal O_{\mathbf P^1\times\mathbf P^1}(1,1).
\tag{2.4}
\]

The map has degree `M+beta_GB(P)`, `A_P` is ample Cartier, and the projection
formula gives

\[
 \boxed{A_P^2=2\bigl(M+\beta_{\rm GB}(P)\bigr),\qquad
 \beta_{\rm GB}(P)=\tfrac12A_P^2-M.}
\tag{2.5}
\]

Thus a bound on the polarized Hilbert polynomial of these **marked exact
carriers** bounds the dessin complexity directly.  This is the precise
surface attribute that can replace an arithmetic-volume bound.  Bounding an
unmarked underlying surface, or normalizing by an unbounded cover degree,
does not control (2.5).

### Exact return

On each normalized leaf, the marked four-point cross-ratio recovers
`lambda_{ij}`.  The discrete bit gives `u_{ij}=0` in the orthogonal case;
otherwise

\[
                         u_{ij}=\frac{\lambda_{ij}}{1-\lambda_{ij}}.
\tag{2.6}
\]

At the retained distinguished real embedding, the acute outward-normal
convention chooses the off-diagonal entry `-sqrt(u_ij)`, with zero in the
orthogonal case.  Thus the packet recovers the complete outward unit-normal
Gram matrix; the retained time sheet fixes the intended labelled halfspaces.
Its rank and inertia recover the hyperbolic dimension and the halfspaces
recover the same chamber.  This does **not** recover primitive integral root
lengths or an integral lattice unless those are added as a separate marking.
Literal
Vinberg edge and full parabolic-neighbourhood tests remain the exact terminal
finite-volume filter.  The surface is therefore faithful finite encoding,
not evidence inferred from incidence alone.

The fields are outputs.  The numbers `u_{ij}` recover their finite
coefficient algebra; adjoining their distinguished square roots recovers the
unit-normal Gram field, and cyclic Gram products recover the Vinberg ground
field.  No field is prelisted in the construction.

## 3. Dessin-native finiteness and field effectivity

### Theorem 3.1

For every integer `B`, there are only finitely many arithmetic finite-volume
Coxeter chambers, in all dimensions, with
`beta_GB(P)<=B`, up to hyperbolic isometry.  Their normalized Grams and
Vinberg fields are effectively enumerable from finite permutation data and
algebraic elimination.

### Proof

Equation (1.6) gives

\[
 f\leq\frac{1+\sqrt{1+4B}}2,
 \qquad n\leq f-1.
\tag{3.1}
\]

There are finitely many normalized Belyi polynomials of degree at most `B`,
and each has finitely many points above zero and one.  Hence there are
finitely many labelled pair packets, finitely many matrices `(u_ij)`, and
finitely many square-root/cyclic-product field packages.  For each candidate,
exact rank, Lorentz signature, Coxeter-angle, irredundancy and Vinberg
finite-volume predicates decide whether it returns a chamber.  A full
spanning ordered Gram determines its hyperplane configuration up to one
ambient isometry, so no separate rigidity theorem is needed.  Enumeration is
effective by listing finite monodromy permutation pairs/passports and solving
their finite coefficient equations.  This is precisely the usual finite
dessin mechanism.  QED.

This is the requested field-list mechanism: **a bound on a genuinely
algebro-geometric Belyi degree directly gives the finite list of fields**.
It does not import that list from arithmetic covolume finiteness.

By (2.5), the entirely geometric hypothesis `A_P^2<=V` lets one take
`B=floor(V/2)`.  In particular

\[
 f(f-1)\leq\beta_{\rm GB}(P)\leq\tfrac12V.
\tag{3.2}
\]

So a log-bounded or Hilbert-bounded family of these marked surfaces gives an
ABSW-style finiteness theorem directly.  The genuinely open step is not the
surface-to-field implication; it is deriving one uniform `V` from reflection
maximality or another intrinsic Coxeter condition.

There is also a deliberately crude closed field-degree bound.  A degree-`d`
Belyi cover has a transitive monodromy pair in `S_d x S_d`, so there are at
most `(d!)^2` ordered permutation pairs before simultaneous-conjugacy
quotienting.  For a polynomial, `infinity` is the unique totally ramified
pole, while the marked source points `0,1,lambda` lie in
`B^{-1}({0,1})`, a set of size at most `2d`.  Thus their ordered choices
number at most `(2d)^3`, which we harmlessly enlarge to `(2d)^4`.  The ordered frame
kills every source `PGL2` automorphism.  Consequently the number of geometric
points in the degree-`d` marked-packet scheme, and hence every residue-field
degree in that finite scheme, is at most `(2d)^4(d!)^2`.  If
`M=binom(f,2)` and the nonorthogonal leaf
degrees are `d_ij`, the Galois orbit of the complete lambda packet is
therefore at most

\[
              \prod_{u_{ij}>0}(2d_{ij})^4(d_{ij}!)^2.
\tag{3.3}
\]

The field generated by the `u_ij` has degree no larger than the packet count
in (3.3).  Adjoining
the distinguished square roots needed for the unit-normal Gram costs at most
`2^M`; the Vinberg field is a cyclic-product subfield.  Consequently

\[
 [k_{\rm Vin}(P):\mathbf Q]
 \leq
 2^M\prod_{u_{ij}>0}(2d_{ij})^4(d_{ij}!)^2
 \leq
 2^{\lfloor B/2\rfloor}
 \bigl((2B)^4(B!)^2\bigr)^{\lfloor B/2\rfloor}.
\tag{3.4}
\]

The final expression is intentionally enormous, but it is explicit and
dessin-native.  It is a field-degree consequence of the geometric degree
bound, not an input field list.

## 4. The high-dimensional lower bound is real

The frozen K3 divisor-lattice redoubling family gives finite-volume
noncompact polytopes `P_N` in `H10` with `N+12` facets and reflection index
`N`.  In its fixed lattice

\[
 F=U+E_8(-2)+\langle-2\rangle,
\]

let `a,b` be the two parallel doubling roots and `c=q_0` the first moving
root.  The even moving roots are

\[
 q_{2m}=c+2m^2a+2m(m-1)b.
\tag{4.1}
\]

Exact lattice arithmetic gives

\[
 c^2=q_{2m}^2=-4,
 \qquad (c,q_{2m})=4(m^2-1),
 \qquad u(c,q_{2m})=(m^2-1)^2.
\tag{4.2}
\]

Let `p` be any odd prime, take `m=p+1`, and choose the phase-`A` chamber
`P_{2(p+1),A}`.  Both `q_0` and `q_{2m}` are boundary facets.  Then

\[
 u=[p(p+2)]^2,
 \qquad \lambda=\frac{u}{u+1},
 \qquad v_p(\lambda)=2.
\tag{4.3}
\]

Rodriguez's main theorem says that a nonzero algebraic number with nonzero
`p`-adic valuation has normalized Belyi-polynomial height at least `p`.
Therefore

\[
 \boxed{\beta_{\rm GB}(P_{2(p+1),A})\geq p.}
\tag{4.4}
\]

So the Belyi degree is unbounded already in fixed `H10`, over `Q`, in one
fixed Picard lattice and with one fixed root-norm/divisibility alphabet.  The
machine certificate checks (4.2) symbolically from the archived roots and
checks (4.3) for every odd prime through 97.  The symbolic argument proves it
for every odd prime.

## 5. What this changes, and the remaining completion step

The general construction closes three pieces that were previously vague:

```text
every arithmetic P
  -> finite exhaustive marked Belyi-surface stack
  -> exact same unit-normal Gram and Vinberg field;

bounded Gram-Belyi degree
  -> finite explicit dessins
  -> finite explicit Coxeter chambers and field list.
```

The remaining hard theorem is now sharp:

\[
 \boxed{\text{reflection-maximal P (or bounded intrinsic surface data)}
        \Longrightarrow \beta_{\rm GB}(P)\leq B(n,\text{geometry}).}}
\tag{5.1}
\]

Maximality by itself has not been shown to imply (5.1).  The maximal
`Gamma16` chamber has no Wang/Borcherds correction, so a proof cannot simply
identify maximality with purity of a reflective modular form.  Conversely,
the `H10` theorem above proves that a class containing arbitrary finite-index
descendants cannot have a uniform bound.  A viable surface proof must bound
the **absolute marked** Belyi/Hilbert degree on a canonical carrier, or bound
the Coxeter cover index separately; normalized degree on an unbounded cover
cannot work.

## 6. Exact replay

Sections 1--3 are a mathematical construction and finiteness proof; the
current executable does not enumerate arbitrary Belyi passports or formally
verify that general argument.  Its exact machine scope is Section 4: the
symbolic `H10` Gram formula and its `p`-adic lower-bound calibration.

- producer:
  [`certify_gram_belyi_h10_lower_bound.py`](../scripts/certify_gram_belyi_h10_lower_bound.py);
- archive:
  [`gram_belyi_h10_lower_bound.json`](../data/enumerations/gram_belyi_h10_lower_bound.json);
- test:
  [`test_gram_belyi_h10_lower_bound.py`](../tests/test_gram_belyi_h10_lower_bound.py).

Run:

```bash
python3 scripts/certify_gram_belyi_h10_lower_bound.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_gram_belyi_h10_lower_bound.py
```

The replay uses exact integer/rational arithmetic and no cone enumeration.
