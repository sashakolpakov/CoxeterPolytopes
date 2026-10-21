# Golden normal scalings: a finite theorem and the exact remaining gap

**Status (2026-08-09): theorem, with hostile boundary examples.**  The
all-normal-scaling part of the golden K3 compiler is no longer open for the
class defined below.  Infinite powers of the fundamental unit do occur in
Cartan matrices, but they are vertex coboundaries: they do not produce new
wall lattices.  In particular the equal-norm compiler is exhaustive for
Bugaenko's compact `H8` chamber and for the exhaustive twelve compact
equal-root golden Lanner diagrams.

This closes a finite-enumeration problem.  It does **not** prove Missing
Lemma U: primitive K3 embeddability, survival of a saturation, and equality
with the complete visible root arrangement remain independent gates.

## 1. The surface-free class

Put

\[
F=\mathbf Q(\sqrt5),\qquad O=\mathbf Z[\phi],\qquad \phi^2=\phi+1.
\]

Let `C_35^conn` be the class of oriented crystallographic Coxeter
**realizations** with the following properties.

1. The facet normals span a nondegenerate `F`-space.  This is automatic for
   the full-dimensional compact polytopes of interest.
2. For chosen normals `r_i`, every Cartan entry

   \[
   a_{ij}=\frac{2(r_i,r_j)}{(r_i,r_i)}
   \]

   belongs to `O`, and the Cartan matrix is symmetrizable.  Equivalently, if
   `d_i=(r_i,r_i)/2`, then

   \[
   d_i a_{ij}=d_j a_{ji}.
   \tag{1.1}
   \]

3. The graph on the facets whose edges have Coxeter label three or five is
   connected.

No surface, K3 embedding, root completeness, or carrier-success condition
occurs in this definition.  A realization means the actual wall reflection
representation.  An angle-only Coxeter matrix with an infinity entry need
not specify such a realization; Section 6 gives an exact infinite-family
warning.

## 2. Unit-coboundary rigidity

### Theorem 2.1 -- unique golden normal scaling

For each realization in `C_35^conn`, all integral symmetrizable choices of
oriented facet normals form one class under independent positive unit
rescaling of the normals and one common homothety.  The class has a unique
normalized symmetric representative

\[
(r_i,r_i)=2,qquad
(r_i,r_j)=
\begin{cases}
-1,&m_{ij}=3,\\
-\phi,&m_{ij}=5,\\
0,&m_{ij}=2.
\end{cases}
\tag{2.1}
\]

Every other pair, including a dotted pair, is also symmetric and integral in
this normalization, although its value is part of the geometric realization
unless rank or rigidity forces it.

#### Proof

For an edge of label `m`, crystallographicity and Coxeter geometry give

\[
a_{ij}a_{ji}=4\cos^2(\pi/m)=
\begin{cases}
1,&m=3,\\
\phi^2,&m=5.
\end{cases}
\tag{2.2}
\]

Both products are square units of `O`.  Since both factors in (2.2) are
algebraic integers, each factor is a unit.  Equation (1.1) gives

\[
\frac{d_j}{d_i}=\frac{a_{ij}}{a_{ji}}
=\frac{a_{ij}^2}{a_{ij}a_{ji}},
\tag{2.3}
\]

which is the square of a unit.  Connectivity propagates this statement from
one base vertex to every vertex.  Hence there are units `u_i`, chosen
positive at the distinguished embedding, such that `u_i^2 d_i` is one
common value.  Replacing `r_i` by `u_i r_i` and dividing the form by that
common value gives square two for every root.

Now (1.1) makes the normalized Cartan matrix symmetric.  The oriented
Coxeter signs and (2.2) force `-1` and `-phi`, proving (2.1).

There is no hidden cycle parameter.  More explicitly, compare two integral
symmetrizable systems with the same oriented Cartan-entry ideals and write
`c_ij=a'_ij/a_ij`, a unit.  Around any cycle, symmetrizability gives

\[
\left(\prod c_{ij}\right)^2=1.
\]

All `c_ij` are positive at the distinguished place, so the product is `+1`,
not `-1`.  Thus `c` has trivial holonomy and is a vertex coboundary
`c_ij=u_j/u_i`.  This proves uniqueness, including arbitrarily large powers
of `phi`.  `square`

### Corollary 2.2 -- the wall span does not change

Let `R=sum_i O r_i`.  Vertex unit rescaling leaves every rank-one module
`O r_i` unchanged, hence leaves `R` unchanged.  A common rescaling changes
it only by the declared homothety.  Moreover `O=Z[phi]` is a PID, so a
finitely generated torsion-free wall span is free; there is no nonprincipal
Steinitz class hiding behind the choice of generators.

If there are more facets than the ambient rank, the conclusion is unchanged.
The relation module is

\[
O^s\cap\ker_F(G),
\]

which is saturated because the actual wall span is torsion-free.  Once the
normalized Gram `G` is fixed, its integral wall lattice is therefore fixed.

## 3. Finite saturation theorem

### Theorem 3.1 -- all normalized integral ambient lattices are finite

Let `R` be the full-rank wall lattice from Theorem 2.1 with its normalized
integral form.  Every integral `O`-lattice `M` in the same `F`-space that
contains the displayed wall roots satisfies

\[
R\subseteq M\subseteq R^\vee.
\tag{3.1}
\]

Consequently all field-integral, even, `O`-stable, and reflection-stable
overclosures are obtained by testing submodules of the finite module
`R^vee/R`.  For the slightly broader K3 convention in which only the
`Z`-trace form must be even integral, apply the same argument to the trace
lattice `N_R`: every candidate lies between `N_R` and `N_R^vee`, and
`O`-stability is a finite test on that discriminant group.  Thus the
equal-norm plus discriminant-module compiler is exhaustive for every
realization in `C_35^conn`.

The hypotheses in this sentence matter:

- the normals must span the ambient `F`-space;
- the ambient lattice must contain the wall roots and be integral for the
  normalized form;
- the common scale must be locked at root square two; and
- distinct geometric realizations of an infinity label must either be
  supplied or separately forced.

Without full span, arbitrary orthogonal lattice summands survive.  Without
integrality, (3.1) fails.  Without scale locking, Section 6 gives infinitely
many trace discriminants.

There is also a slightly broader finiteness statement.  If all nonzero
products `p_ij=a_ij a_ji` are fixed algebraic integers, the ideals
`(a_ij)` are divisors of `(p_ij)`, so there are finitely many ideal choices.
For each ideal choice the same cycle argument shows that all unit twists are
one vertex-unit class.  Hence fixed products give finitely many Cartan
classes.  The square-unit connectivity in Theorem 2.1 is stronger: it makes
that finite list a singleton and forces equal norm.

## 4. Bugaenko: the labelled matrix forces everything

Bugaenko's eleven-facet diagram has ten `m=3,5` edges, and they form a
spanning tree.  Theorem 2.1 therefore kills every unequal-root and infinite
unit twist.  There is one apparent ambiguity left in an angle-only matrix:
facets 9 and 10 are ultraparallel.  Write `x` for their symmetric normalized
Gram entry.

The wall span has rank nine.  Hence every ten-by-ten principal minor must
vanish.  Direct exact determinant expansion gives

\[
\begin{aligned}
D_0(x)&=(-9+7\phi)
\left(x-\frac{36+2\phi}{31}\right)(x+2\phi),\\
D_1(x)&=(-16+12\phi)
\left(x-\frac{4+2\phi}{5}\right)(x+2\phi).
\end{aligned}
\tag{4.1}
\]

The two secondary roots differ, so the only common zero is

\[
\boxed{x=-2\phi.}
\tag{4.2}
\]

The resulting Gram has rank nine and is exactly the archived Bugaenko Gram.
Thus the labelled Coxeter matrix **plus dimension eight** determines not only
the normal-scaling class but also its sole dotted distance.  This is stronger
than the fixed-Gram theorem previously used by the K3 compiler.

The first nine archived walls have coordinate determinant `-3+2phi`, of
field norm `-1`, so they are an `O`-basis of the ambient Bugaenko lattice.
Its field discriminant module has order four and is `O/(2)`.  Since `2` is
inert, this is a simple `O`-module.  Its only nonzero submodule would give the
full dual lattice, whose line generator has square `-1/(2phi)` and is not
integral.  Therefore Bugaenko has exactly one field-integral
`O`-overclosure: the wall lattice itself.

The broader trace-even check gives the same answer.  The trace lattice is
`<2>+<-2>+E8(-1)^2`, with discriminant group `(Z/2)^2`.  In its standard
basis multiplication by `phi` acts by

\[
\begin{pmatrix}0&1\\1&1\end{pmatrix}\pmod2,
\]

which cycles all three nonzero vectors.  Hence it preserves no order-two
line; in particular the unique isotropic line defining the ordinary
index-two even overclosure is not `O`-stable.  Thus the wall trace lattice is
also the only `O`-stable even trace overclosure.

The conclusion for Bugaenko is exact:

```text
labelled matrix + H8 rank
  -> one normalized Gram
  -> one wall lattice
  -> one field-integral O-saturation and one O-stable even trace saturation.
```

This completely closes its all-scaling and stable-saturation bookkeeping.
Its K3 primitive embedding is already known to pass, but the present theorem
does not turn such passage into a universal statement for other diagrams.

## 5. The twelve golden Lanner diagrams

Every frozen seed uses only labels `2,3,5`, has no dotted edge, and has a
connected `m=3,5` graph.  Its labelled matrix therefore directly gives the
entire normalized Gram (2.1), with one normal-scaling class.  The exhaustive
twelve-row census is

```text
H2: (2,5,5), (3,3,5), (3,5,5), (5,5,5)
H3: 3-5-3, 5-3-5, D4(5,3,3), C4(one 5), C4(opposite 5s)
H4: 5-3-3-3, 5-3-3-5, D5(terminal 5).
```

Their finite discriminant overclosure searches may have more than one row,
but no additional normal-scaling lane exists.

## 6. Exact hostile boundaries

### 6.1 Connected Coxeter diagram is not enough

Take two walls with one infinity edge and put

\[
p=3+\phi,qquad
A=\begin{pmatrix}2&-1\\-p&2\end{pmatrix},\qquad
D=\operatorname{diag}(p,1).
\]

Then

\[
G=DA=\begin{pmatrix}2p&-p\\-p&2\end{pmatrix}.
\tag{6.1}

At the identity embedding, `p>4` and `G` has signature `(1,1)`; at the
conjugate embedding, `0<p'<4` and `G'` is positive definite.  The two
reflections bound a compact arithmetic interval in `H1`.  But `N_F/Q(p)=11`,
so `p` is not a square in `F`.  Equal root length would require a square root
of `p`.  Thus an equal-norm compiler is false when the square-unit graph is
not connected, even though the Coxeter diagram is connected and the Cartan
matrix is integral and symmetrizable.

### 6.2 An infinity label alone has infinitely many arithmetic realizations

For every `k>=1`, put

\[
p_k=3+\phi^{2k},\qquad
A_k=\begin{pmatrix}2&-1\\-p_k&2\end{pmatrix}.
\tag{6.2}

At the identity place `p_k>4`, while at the conjugate place
`3<p_k'=3+phi^{-2k}<4`.  Hence (6.1) with `p_k` gives a compact admissible
arithmetic `H1` reflection interval.  The products `p_k` are distinct and
are invariant under all normal rescalings.  All members nevertheless have
the same angle-only Coxeter matrix, consisting of one infinity edge.

Therefore no theorem may claim a finite realization list from an
angle-only matrix with dotted edges.  It must use the actual arithmetic
Gram, a rigidity theorem, or rank equations such as (4.1).

### 6.3 Common form homothety is harmless geometrically but not for K3 data

For any normalized wall form `h_0`, the forms `h_m=m h_0`, with `m` a
positive integer, have the same Cartan matrix, wall reflections, and
polytope.  If the trace scalar `alpha` is held fixed, however, a rank-`r`
`O`-lattice has `Z`-rank `2r` and its trace determinant is multiplied by
`m^(2r)`.  These are infinitely many different discriminant forms.

Thus common homothety can be quotiented in the geometric scaling theorem only
if the AG compiler locks root square two, or if it simultaneously replaces
`alpha` by `alpha/m` so that the trace form is unchanged.  Quotienting it
silently while keeping a fixed K3 intersection form would be invalid.

## 7. What has and has not been completed

For `C_35^conn`, the route is now

```text
P and its wall reflection representation
  -> unique equal-norm wall lattice
  -> finite reflection-stable saturation list
  -> finite K3 embedding and hidden-root tests.
```

The first two arrows are the theorem here.  They solve the normal-scaling
clause of Missing Lemma U for this independently defined class.  They do not
supply the existential statement that at least one row passes the last
arrow.  Therefore:

- a failed finite list now excludes a carrier of the normalized full-trace
  species for that `P`;
- it does not exclude `P` as a Coxeter polytope unless a separate theorem
  proves that every member of the class must possess such a carrier; and
- the K3 Picard ceiling remains a carrier ceiling, not a general polytope
  ceiling.

## Replay

The exact certificate checks the unit-exponent cancellation, all twelve
Lanner graphs, Bugaenko's spanning tree, both determinant factorizations,
the forced dotted entry, the unit-basis wall span, its unique integral
`O`-overclosure, the nonsquare bridge, and the first twelve members of the
infinite dotted family.

```sh
PYTHONPATH=.:src python3 \
  scripts/certify_golden_normal_scaling_universality.py --check-archive
PYTHONPATH=.:src pytest -q \
  tests/test_golden_normal_scaling_universality.py
```

Artifacts:

- `scripts/certify_golden_normal_scaling_universality.py`;
- `data/enumerations/golden_normal_scaling_universality.json`;
- `tests/test_golden_normal_scaling_universality.py`.
