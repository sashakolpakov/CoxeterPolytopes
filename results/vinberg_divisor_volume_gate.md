# Vinberg's volume gate for divisor chambers

## Outcome

The algebro-geometric track has not been abandoned.  It now has a hard,
exact decision gate:

1. algebraic geometry supplies a Lorentzian divisor lattice and a **complete**
   finite wall configuration;
2. Vinberg's diagram criterion decides compactness or finite covolume from
   that complete configuration; and
3. on a K3 surface, every ideal-vertex test has a direct elliptic-fibration
   interpretation, with rank defect measured by the Mordell--Weil group.

The gate has been run directly on complete K3 divisor diagrams through
hyperbolic dimension 18.  These polytopes are known, not newly discovered,
but the exact all-subdiagram certificates and the local Mordell--Weil
dictionary make them useful calibration objects for inverse design.

| divisor chamber | dimension | facets | elliptic edges | finite vertices | cusps | decision |
|---|---:|---:|---:|---:|---:|---|
| Bugaenko/K3 eigen-nef | 8 | 11 | 164 | 41 | 0 | compact |
| quotient-nef, $\rho=18$ cover | 9 | 10 | 45 | 8 | 2 | finite volume, noncompact |
| quotient-nef, $\rho=19$ cover | 10 | 13 | 208 | 34 | 4 | finite volume, noncompact |
| full K3 nef, $\rho=15$ | 14 | 21 | 2676 | 328 | 15 | finite volume, noncompact |
| full K3 nef, $\rho=16$ | 15 | 19 | 1563 | 186 | 7 | finite volume, noncompact |
| full K3 nef, $\rho=17$ | 16 | 19 | 1087 | 125 | 3 | finite volume, noncompact |
| full K3 nef, $\rho=18$ | 17 | 19 | 729 | 80 | 2 | finite volume, noncompact |
| full K3 nef, $\rho=19$ | 18 | 24 | 5160 | 525 | 10 | finite volume, noncompact |

Every listed edge has exactly two certified endpoints.  Every connected
affine subdiagram has a full rank-$(n-1)$ parabolic zero-neighborhood; there
are respectively 34, 13, 6, 3, and 22 such connected-affine audits in the
five full K3 cases, with no failure.

### Sharp theorem for full complex K3 nef chambers

For a complex projective K3 surface $X$, let $P_X$ be the full nef chamber
inside the projectivized positive cone of $NS(X)$.  The K3 Weyl-chamber
theorem, arithmetic finite covolume of $O^+(NS(X))$, and global Torelli give

$$
\operatorname{Vol}_{\mathrm{hyp}}(P_X)<\infty
\Longleftrightarrow
[O^+(NS(X)):W_{(-2)}(X)]<\infty
\Longleftrightarrow
|\operatorname{Aut}(X)|<\infty.
$$

Equivalently, the nef/effective cone is rational polyhedral and the Cox ring
is finitely generated.  Artebani--Hausen--Laface's
[Theorems 2.7, 2.11, and 2.12](https://arxiv.org/abs/0901.0369) classify the
eligible Néron--Severi lattices; there is no rank-20 case.  Therefore

$$
\boxed{\dim P_X=\rho(X)-1\leq18.}
$$

The bound is attained: the audited rank-19 surface with
$NS(X)=U\oplus E_8(-1)^2\oplus A_1(-1)$ has a complete 24-wall nef chamber
of finite volume in $\mathbb H^{18}$.  This is a sharp boundary only for
**full nef chambers of complex projective K3 surfaces**.  It says nothing
about invariant or irrational eigen-nef slices, Picard-lattice reflection
groups enlarged by non-divisorial roots, supersingular K3 surfaces, or
arbitrary hyperbolic Coxeter polytopes.

The reusable checker is
[`vinberg_volume.py`](../src/coxeter_polytopes/verify/vinberg_volume.py), the
five expensive high-rank replays are in
[`certify_k3_high_rank_vinberg_volume.py`](../scripts/certify_k3_high_rank_vinberg_volume.py),
and their frozen output is
[`k3_high_rank_nef_vinberg_volume.json`](../data/known_examples/k3_high_rank_nef_vinberg_volume.json).
Each high-rank row carries a canonical SHA-256 digest of the complete finite
vertex list, ideal neighborhoods and component decompositions, every edge
endpoint extension, and every full parabolic-neighborhood audit; the archive
therefore freezes the proof data beyond the displayed counts without storing
multi-megabyte lists.

## The exact criterion, in its two useful forms

Let $G=(g_{ij})_{i,j\in I}$ be the indecomposable Gram matrix of a
finite-sided acute Coxeter chamber in $\mathbb H^n$: $g_{ii}>0$,
$g_{ij}\leq0$, and
$\operatorname{inertia}(G)=(n,1,|I|-n-1)$.  The following is Vinberg's
criterion, stated in [Vinberg, Section 4](https://www.mathnet.ru/links/df899e7a8aeed4a3951f34fd2ab8427c/rm2140_eng.pdf)
and in a modern diagrammatic form in
[Heckman, Theorem 5.12](https://www.math.ru.nl/~heckman/CoxeterGroups.pdf#page=106).

### Full-neighborhood/Lannér form

For a connected affine subdiagram $J$, put

$$
Z(J)=\{i\in I\setminus J:g_{ij}=0\text{ for every }j\in J\},
\qquad N(J)=J\sqcup Z(J).
$$

The chamber has finite volume if and only if both conditions hold:

1. for every connected affine $J$, the **full** diagram on $N(J)$ is a
   disjoint union of $m\geq1$ connected affine diagrams and

   $$|N(J)|=n+m-1,$$

   equivalently its Gram rank is $n-1$; and
2. for every proper Lannér subdiagram $J$, the corresponding face cone obeys
   $C_{N(J)}\cap V'=\varnothing$ in the restricted dual Lorentz space.

Compactness is the same statement with no affine subdiagram at all and with
the Lannér condition retained.  The easy sufficient screen
$|N(J)|=n+1$ implies the Lannér condition, but is not its definition and is
not necessary in general.

The word **full** is essential.  It is wrong to select an affine subset and
ignore a wall orthogonal to all of it.  For example, if
$N(J)=\widetilde A_1\sqcup A_1$, the affine part alone looks like a cusp but
the finite $A_1$ contaminant makes the horospherical link noncompact.  The
new adversarial regression test rejects exactly this case.

### Edge-endpoint form

The equivalent form used for the complete finite diagrams is
[Guglielmetti, p. 761, Propositions 1.9--1.10](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/A117FD150E108CACCACF9EF00C61A23C/S1461157015000273a.pdf/div-class-title-coxiter-computing-invariants-of-hyperbolic-coxeter-groups-div.pdf#page=8):

- an elliptic subdiagram on $n-1$ nodes is an edge;
- an elliptic subdiagram on $n$ nodes is a finite endpoint; and
- a componentwise-affine diagram of Coxeter rank $n-1$ is an ideal
  endpoint.

The literal first bullet for finite covolume is that the diagram contain at
least one elliptic rank-$n$ subdiagram or one parabolic rank-$(n-1)$
subdiagram--that is, at least one finite or ideal vertex.  Every elliptic
rank-$(n-1)$ edge must then extend in exactly two ways to an elliptic rank-$n$
or parabolic rank-$(n-1)$ subdiagram.  Here **all** parabolic subdiagrams of
the required rank are counted, not merely a preselected list of maximal
parabolics.  For compactness the diagram must have an elliptic rank-$n$
vertex and both extensions of every edge must be elliptic.

At a non-simple cusp with $m$ affine components, the full parabolic diagram
has $n-1+m$ nodes.  An incident elliptic edge is obtained by deleting one
node from **each** component.  The whole parabolic neighborhood counts as
one endpoint; counting its added nodes or components separately gives the
wrong answer.  Three of the four quotient-nef $H^{10}$ cusps exercise this
reducible case.

Because Propositions 1.9--1.10 are themselves if-and-only-if criteria, this
complete edge-extension audit does not need a second, independent Lannér
calculation.  In Heckman's alternative Theorem 5.12 formulation the Lannér
cone condition remains an independent requirement; the easy equation
$|N(J)|=n+1$ is only sufficient for that requirement.  The two formulations
must not be mixed by checking only affine neighborhoods and silently
dropping Lannér diagrams.

For the quotient-nef $\mathbb H^9$ and $\mathbb H^{10}$ diagrams, a separate
literal all-subset replay now enumerates every elliptic edge, every finite
vertex, and every parabolic rank-$(n-1)$ subdiagram.  Its endpoint-type
histograms are respectively
$\{(2,0):28,(1,1):16,(0,2):1\}$ and
$\{(2,0):137,(1,1):66,(0,2):5\}$, and its endpoint sets agree exactly with
the optimized full-neighborhood routine.  Those counts become a proof only
after the complete, irredundant, acute Coxeter Gram hypotheses above have
been established; counts from an arbitrary wall sample would not suffice.

## The wall-completeness gate

There are two legitimate modes and one illegitimate one.

1. **A complete algebraic wall list.**  If a theorem gives all irreducible
   negative curves (or all projected facet orbits), the criterion directly
   decides its chamber.  This is the mode used for the five full K3 nef cones
   and the two quotient-nef cones.
2. **An exactly replayed Vinberg accepted-root prefix.**  If the current
   outer chamber passes the finite-volume criterion, Vinberg's stopping
   theorem proves that no later accepted simple wall exists.  Thus finite
   volume closes the list rather than presupposing it.  The generic join
   [`verify_vinberg_prefix_stopping`](../src/coxeter_polytopes/verify/vinberg.py)
   checks equality of the replayed accepted roots and certificate walls,
   exact extreme-ray geometry, and the diagram criterion.  It supports both
   compact and cofinite noncompact stopping; Bugaenko is the compact
   calibration.
3. **An arbitrary bounded wall sample or height prefix.**  A failed
   finite-volume test is inconclusive: later walls may cut off the observed
   hyperideal escape.  Conversely, if the sampled inequalities are proved to
   be genuine, coherently oriented walls containing the final chamber and
   their outer polyhedron is exactly certified finite-volume (or compact),
   monotonicity gives the same one-sided conclusion for the final chamber.
   What incompleteness does *not* justify is treating the sample's vertices,
   cusps, or facet list as final.  The low-level public API lacks the extra
   geometric data needed to distinguish these cases, so it conservatively
   returns `undecidable_incomplete_simple_roots` and emits no prefix vertex
   counts; the dedicated Vinberg-prefix verifier supplies the stopping step
   when its stronger hypotheses hold.

This distinction removes the earlier ambiguity: a finite-volume **verified
Vinberg prefix** is a stopping certificate, while a finite numerical sample
of walls is not.

## K3 interpretation of the cusp test

Let $X$ be a projective K3 surface of Picard rank $\rho$, and suppose here
that we are testing a finite-sided **full nef chamber with its complete curve
wall list**.  Use the intersection form on $NS(X)$, of signature
$(1,\rho-1)$ by the Hodge index theorem.  The Hodge theorem supplies
hyperbolicity for free, but it says nothing about the edge-endpoint design
restrictions.

Let $F$ be a primitive nef isotropic class on the boundary of the nef
chamber.  It defines a genus-one pencil.  Every active $(-2)$-wall
$C$ satisfies $F\cdot C=0$, so $C$ is vertical; the active curves are the
components of reducible fibers.  Their diagrams are affine ADE, and after
deleting one component from each fiber the finite root lattices $R_v$ have
total rank equal to the Coxeter rank of the full parabolic cusp diagram.
Vinberg's local finite-volume condition becomes

$$
\sum_v\operatorname{rank}R_v=\rho-2.
$$

For an elliptic fibration with a section, Shioda--Tate says

$$
\rho=2+\sum_v\operatorname{rank}R_v
      +\operatorname{rank}MW(X/F).
$$

Hence

$$
\boxed{\text{full-rank Vinberg cusp}\iff
       \operatorname{rank}MW(X/F)=0.}
$$

For a genus-one fibration without a section, the statement is made for its
Jacobian (equivalently the frame lattice); the section hypothesis must not
be silently assumed.  A rank defect is therefore not mysterious
hyperbolic leakage: it is the translation direction supplied by a
positive-rank Mordell--Weil group.

In this full, complete-nef setting this gives the desired classical algebraic
reading:

- compactness means that the full nef-chamber boundary has no primitive nef
  isotropic ray, hence no genus-one pencil;
- local finite covolume at every cusp means that every boundary elliptic
  fibration has full affine ADE rank, equivalently Mordell--Weil rank zero;
- the remaining global condition says that every negative-definite
  rank-$(\rho-2)$ curve configuration defining an edge closes at exactly two
  contractible ADE or genus-one-fiber endpoints; and
- for the **full** K3 nef cone, finite hyperbolic volume is also equivalent
  to finite automorphism group, rational polyhedrality, and finite generation
  of the Cox ring.  Hyperbolic volume itself is not a Chern number: all K3
  surfaces have the same elementary Chern numbers, while their nef-chamber
  behavior varies.

## Rational-surface interpretation and the positive/negative carrier criterion

The same endpoint gate has a useful form on a smooth rational surface, even
though its isotropic pencils need not have genus one.  Suppose the complete
Mori cone is generated by irreducible negative curves `C_i`, their primitive
reflections are integral, and their Gram is acute Lorentzian.  For a
primitive nef isotropic ray `F`, every connected affine component `A_j` of
the complete zero-neighbourhood has a positive null cycle

$$
D_j=\sum_{i\in A_j}a_iC_i=m_jF.
$$

Thus Vinberg's full parabolic-rank condition is visible in the divisor
algebra: the orthogonal negative curves must be a union of affine
configurations of total rank `rho-2`, and their null cycles must span the
cusp ray.  If two components occur, suitable multiples of `D_1,D_2` are
distinct effective linearly equivalent divisors.  They are supported on
orthogonal, disjoint curve configurations and have intersection zero, so
they generate a base-point-free pencil.  When `m_1=m_2=1`, the primitive
class `F` itself is semiample, with fibre genus

$$
g(F)=1+\frac{K_X\cdot F}{2}.
$$

This is exact in the new rational `P12`--`P16` carrier families.  Independent
dual-cone enumeration reconstructs all 28 cusps; every affine null cycle is
already primitive and equals `F`.  Sixteen multi-component cusps, five line
rulings, and one quartic pencil give the first 22 certified semiample classes.
For five of the six one-component residues, exact complete-ideal jet matrices
have nullity two and provide a second plane section with the assigned
infinitely-near multiplicities, no exceptional component, and no common
plane factor with the affine wall cycle.  Their square-zero pencils are
therefore base-point-free.  The last genus-three primitive system has nullity
one, as do its second and third multiples, while the fourth multiple has
nullity two.  Its stored coprime degree-16 section has exact multiplicities
`(4^16)`, so `|4F|` is base-point-free and the semiample index is exactly
four.  Thus all 28 cusps are semiample.  See the
[`cusp-fibration certificate`](rational_i1n_cusp_fibrations.md).  The later
[`symbolic cusp-family certificate`](rational_i1n_symbolic_cusp_family.md)
proves this `27+1` split throughout the nonempty carrier-open family.  Its
paired exact rank/section gates have gcd one, so no further parameter is
excluded by cusp mobility.

The complementary **negative carrier criterion** is surface-MMP based.  If
the proposed complete facet cone contains `K`-negative square-minus-one
curves, contract every such current facet image, updating all classes by

$$
f^*f_*D=D+(D\cdot e)e.
$$

Exhaust the finite contraction tree through every rank-one and minimal ruled
endpoint, reconstruct `K_X` at each endpoint, and quotient by the integral
weighted-diagram group.  A proposed carrier survives only if at least one
canonical orbit admits all required irreducible curves.  Applied to the
classical `P17`, `P18`, and `P19` facet cones, this gives respectively 8, 24,
and 120 canonical markings and reduces every one to the same Abel
contradiction.  In inverse design this theorem is now a pruning criterion:
it rules out a whole facet cone on smooth rational surfaces, rather than
encouraging further Cremona relabellings of the same impossible design.

Together the two statements give a genuine primal--dual surface test:

1. contraction/MMP data decide whether the prescribed negative curves can
   exist on any smooth rational carrier;
2. the complete Mori cone and integral angle gate produce the Coxeter
   chamber; and
3. negative-definite and affine curve configurations implement Vinberg's
   finite/ideal endpoint test, with affine null cycles giving the actual
   fibrations when they move.

### Fixed and eigen-slices

For $M\subset NS(X)$ of signature $(1,n)$, such as an invariant or eigen
lattice, one must use the **visible projected walls** in
$\operatorname{Nef}(X)\cap\mathbb P(M_\mathbb R)$.
Their normals may be orbit sums or projections rather than individual
$(-2)$-classes.  At an isotropic $F\in M$, Vinberg requires the active
projected diagram to have full rank $n-1$ and to equal its full
zero-neighborhood inside the complete projected wall list.  It is not valid
to insert the rank of a root system in the whole $NS(X)$ into this test.
Only when the visible walls are identified with actual fiber-component
orbits does the projected rank defect inherit the Mordell--Weil
interpretation.

## Status of the incomplete carriers

- The golden OG10 $H^{10}$ search (1,144 recorded walls) and its direct
  golden $H^9$ fixed search (71 walls) do not have terminated finite simple
  root lists.  Their finite-diagram Vinberg status is therefore undecidable.
- The intrinsically complete golden $H^9/H^{10}$ mirror chambers have an
  independent infinite-volume proof from an infinite chamber stabilizer.
  That global theorem remains valid, but it is not a volume inference from
  either finite prefix.
- The H11 and discriminant-19 Vinberg runs are bounded accepted-root
  prefixes whose current outer cones fail the stopping gate; later roots are
  still possible.
- Published higher-dimensional Alexeev/Borcherds examples cannot be run
  through this local checker until their complete wall Gram matrices are
  present in the repository.  A source-level facet count is not a substitute
  for the matrix.

The practical inverse-design target is now precise: design a K3 or quotient
surface whose complete curve-orbit diagram satisfies the two-endpoint rule,
and whose every boundary genus-one fibration has Mordell--Weil rank zero.
Hodge index supplies the Lorentzian ambient space; the fiber and curve
configurations supply the hard combinatorics.
