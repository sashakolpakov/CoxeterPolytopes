# Positive breakthroughs frozen on 2026-08-07

Last strengthened: 2026-08-10.

## Result in one line

TOOT: a complete Coble nef chamber now supplies two first-generation and
seventeen second-generation additional abstract finite-volume Coxeter
polytopes, together with infinite construction sequences in `H14`.  The 35
physical `(-1)`-divisor doubling walls of the
characteristic-three Type-VI `H14/60` seed form two weighted-Gram orbits of
sizes 15 and 20.  Their doubles are the nonisomorphic chambers
`COBLE-VI3-D2A-H14-100` and `COBLE-VI3-D2B-H14-101`, with exact
finite/ideal/hyperideal counts `54,318/1,993/0` and `55,130/2,043/0`.
All facets are irredundant; all physical wall vectors are primitive integral
divisor classes spanning `NS=I_(1,14)`; all reflections are integral; and
their normalized Grams agree entrywise with the arithmetic certificates.
Each double contains an explicit parallel pair among its 68 doubling walls,
so Allcock's redoubling theorem gives a finite-volume noncompact chamber at
every positive subgroup index and hence infinitely many pairwise
nonisometric descendants.  The chambers and families are new to this
project ledger, not claimed worldwide-new.  Doubling crosses the seed nef
wall, so no nef/effectivity/surface-automorphism statement is made for the
descendants.

The 68 doubling walls of `H14/100` and `H14/101` split into respectively
7 and 10 orbits under their intrinsic undecorated diagram groups, of orders
16 and 12.  Doubling one representative of every orbit exhausts the marked
second generation and gives 17 chambers at reflection index four in the
seed, with facet histogram

```text
3*H14/180 + 4*H14/181 + 4*H14/182 + 6*H14/183.
```

Distinct exact normalized-Gram color-refinement invariants prove that the 17
are pairwise nonisometric within this census.  Complete union-ray transport
gives, for facet counts 180 through 183, the four exact
finite/ideal/hyperideal profiles
`107,264/3,925/0`, `108,076/3,975/0`, `108,888/4,025/0`, and
`109,700/4,075/0`.  Literal Vinberg endpoint and affine-neighborhood audits
cover 1,847,848 finite vertices, 68,225 ideal vertices, and 16,221,212
rank-13 elliptic edges across all 17 rows, with zero bad endpoints, invalid
cusps, incomplete zero-neighborhoods, or hyperideal rays.  Every child has
134 doubling walls and an archived parallel pair, so redoubling continues.
These are 17 new-to-project abstract chambers; `H14` already had Allcock's
infinite families, and worldwide priority remains open.  See the
[`Type-VI single-doubling theorem`](coble_type_vi_h14_single_doublings.md)
and the
[`second-generation census`](coble_type_vi_h14_second_generation_doublings.md).

TOOT: two marked infinite `H14` rays now have exact all-length combinatorics.
For the common W12/W20 parallel doubling pair in `H14/100` and `H14/101`,
the dihedral translation satisfies `(T-I)^3=0`.  Solving every possible
wall-collision equation for all `z in Z`, rather than scanning a bounded
number of descendants, proves

```text
f_100(N)=81N+19-floor(N/2),
f_101(N)=82N+19-floor(N/2),   N>=1.
```

Both chambers at length `N` have index and volume multiple `2N` over the
common `H14/60` seed.  Unequal lengths are distinguished by volume, and at
equal length the two facet counts differ by `N`, so all chambers across both
rays are pairwise nonisometric.  The `N=2` coordinate stocks agree exactly
with second-generation rows A02/180 and B04/182.  Every path retains the full
arithmetic lattice and full physical `NS=I_(1,14)` divisor lattice.  This is
a new-to-project exact classification of two marked rays whose existence was
already supplied by Allcock; worldwide priority is not claimed, and the
post-nef-crossing effectivity warning remains.  The audit is deliberately
OOM-safe: 15-by-15 exact reflection transport, all-integer quadratic
collision solving, and small sparse wall unions, with no `cdd`, ray
enumeration, or CoxIter.  See the
[`parallel-path theorem`](coble_type_vi_h14_parallel_path_families.md).

The inverse lane now has a class-wide scoped theorem, but deliberately no
TOOT because it fails the nonthin pushback gate.  For the surface-free class
of twelve-root `E10` nodot presentations with connected subcubic support,
pairings in `{0,1,2}`, at most three pairing-two pairs, finite root volume,
and complete characteristic-boundary/cusp lift, the exhaustive count is
`10,456,309 -> 1,038 -> 837 -> 2`.  Exactly Type I and II survive and every
input obtains an exact complex strong-Coble same-`P` return, `H10/22` or
`H10/30`.  This closes essential surjectivity on a natural non-singleton raw
class.  It remains `L3`, not `L4`, because the P-side cone enumeration has
already proved the two-type classification before Kondō is invoked.
Classical `P10` supplies the sharp boundary: it passes the broad raw
characteristic/`E10`/cusp package but has ten root facets and therefore no
exact complex one-boundary finite-Aut strong-Coble full-nef carrier, while
Kondō's exhaustive Types I and II have twelve.  Hence
`R_CobVin(C) proper_subset D_CobVin`; no positive-characteristic or
other-carrier no-go is claimed.  See the
[`complete subcubic theorem`](coble_nodot_maxdegree3_closure.md).

Positive calibration, deliberately without a TOOT: D'Urso's published
Cremona fundamental cone for `Bl_n(P2)` at very general points has been
replayed exactly and is Coxeter precisely for `n=10,11,13`.  The three
finite/ideal/hyperideal splits are `17+2+0`, `26+2+0`, and `44+2+0`.  The
first two are Vinberg `P10,P11`; the `H13/15` output is Tumarkin
`13-tum04_13_01`, not the nonisomorphic Vinberg-`P13`/Tumarkin-`02` row with
`43+3+0`.  This adds a new project-level algebraic interpretation of a known
abstract chamber.  It is only a Cremona sector in the `K`-negative
nef-effective cone, not a full-nef claim, and its square-`-4`, divisibility-one
canonical mirror is nonintegral on raw `NS` in dimension 13.  The AG-leverage
delta is zero, since the divisor angle formula itself selects the Coxeter
dimensions.  See the
[`D'Urso transfer and lineage theorem`](durso_cremona_nef_coxeter_transfer.md).

TOOT: the algebraic carrier program now has a general finite-volume bridge.
For every smooth projective surface over an algebraically closed field of
arbitrary characteristic and Picard rank at least three, a complete rational-
polyhedral nef cone has finite-volume hyperbolic projectivization;
compactness is equivalent to the absence of isotropic extremal nef rays.  The
proof is elementary but decisive: Hodge index puts every nef ray in the
closed Lorentz cone, so the Klein polytope has only finite or ideal vertices,
and finite/ideal simplices have finite volume in dimension at least two.  The
rank-two exception is sharp (`P1 x P1` gives an infinite-length ideal
interval).  This does not provide Coxeter angles, reflection integrality, or
full-nef equality for a proposed divisor subsystem.  Once those gates are
proved, however, rational polyhedrality is a classical algebro-geometric
certificate of finite volume; after the reflection gate identifies a
fundamental chamber, it is also a finite-covolume certificate.  Isotropic
extremal nef rays are precisely its ideal vertices, hence its cusps in the
Coxeter setting.
See the
[`projective-surface nef-volume theorem`](projective_surface_nef_hyperbolic_volume_theorem.md).

TOOT (structural, classical input): the inverse program now has an explicit
surface-side dimension-bound engine.  Nikulin's theorem says that a smooth
complex projective surface with finite-polyhedral Mori cone, all of whose
extremal negative curves define integral Neron--Severi reflections, has
`rho<=22`; the all-`(-2)` and all-`(-1)` branches give `rho<=19` and
`rho<=9`.  Hence an exhaustive direct full-nef carrier theorem would return
polytope bounds `n<=21`, `n<=18`, or `n<=8` to its original surface-free
class.  His unbounded-Picard-rank polyhedral examples after dropping the
reflection condition also prove that Hodge index plus polyhedrality alone is
too soft for the finite-volume problem; their nef isotropic ruling class makes
them noncompact, so they do not close the compact case.  This is not yet a
universal Coxeter bound--the missing result is precisely carrier
universality--but it supplies the previously missing AG classification step
of the desired `P -> surface -> P` architecture.  See
the [`reflection-surface return theorem`](polyhedral_mori_reflection_surface_return_bound.md).

TOOT: the strict inverse lane now has a completed theorem-shaped prototype,
not just a carrier construction.  Start with the surface-free class of
marked complete finite-covolume norm-two chambers of even integral lattices
having a rational lattice Weyl vector.  Exact integral, full-root, Vinberg,
and Weyl-equation gates force a K3 Picard/nef carrier; the converse K3
criterion is finite automorphism group plus a rational ample ray of equal
degree on every smooth rational curve.  Gritsenko--Nikulin's classification
therefore transfers back to the polytopes: exactly 59 integral lattice
carriers, hyperbolic dimension at most 17, and exactly five compact carrier
classes with compact dimension at most 3.  The positive endpoint is the
known noncompact `U+2E8` chamber in `H17`, with 19 facets and a two-dimensional
K3 moduli locus.  Bugaenko's compact `H8` chamber is rigorously outside this
class: its generalized Weyl equations have ranks `9/10` and residual
`4-8phi`, of field norm `-80`.  Thus the result is not a universal compact
bound, but it is the first exact `P -> surface -> classification -> theorem
on P` benchmark for the broader program.  See the
[`polytope-to-K3 Weyl-vector transfer theorem`](polytope_to_k3_weyl_vector_transfer.md).

TOOT: the whole **direct full-nef surface** lane now has a sharp compact
return theorem.  If a `P`-first compiler identifies the entire nef cone of a
smooth projective surface with a finite-volume Coxeter chamber and every
facet reflection is integral on `NS`, the reflection group is finite-index
in `O^+(NS)`.  Compactness is therefore rational anisotropy; Meyer forces
`dim(P)<=3`.  This is sharp: the K3 Weyl-vector classification returns the
anisotropic rank-four carrier `S_(4,28)` and its compact `H3` chamber.  The
rational branch is stronger still--every rational surface of Picard rank at
least two has a pulled-back nef isotropic ruling class, so its full nef cone
is never compact.  Hence every compact example in `H4+` must use a proper
eigen/folded/automorphic section or a non-direct carrier.  See the
[`direct full-nef compactness theorem`](direct_full_nef_surface_compactness_obstruction.md).

TOOT: the Weyl-vector condition has now been removed on a genuine
non-singleton, fixed-Gram totally-real class.  From a compact field-even
norm-two `Z[phi]` Coxeter realization alone, a finite compiler enumerates all
`O`-stable trace-even overclosures, K3 complement/glue orbits, and all visible
`(-2)` shells.  Primitive embedding plus equality with the complete restricted
root arrangement is an exact iff; a surviving package returns the **same
polytope** as a K3 eigen-nef section.  K3 embedding is automatic in
`H2--H4`, sufficient through `H8` under `ell(A_N)<=18-2n`, and Bugaenko meets
the boundary equality.  The exhaustive twelve-member compact golden Lanner
class certifies a non-singleton low-dimensional positive class; Bugaenko is the explicit
compact `H8` endpoint.  For each fixed `H9` trace lattice the carrier gate is
a finite positive-binary-complement problem, while no `O`-stable integral
full-trace quadratic K3 carrier exists in `H10+`.  That last ceiling is not a
polytope no-go.  Normal scaling is now separately solved on the connected
golden `m=3,5` class, but automatic K3 success or exhaustive alternative
carriers and the rank-twenty classification remain the exact Missing Lemma U.
An independent hostile audit confirms the result with this strict scope.  See the
[`fixed-Gram golden K3 theorem`](polytope_first_golden_k3_leverage.md) and
[`hostile audit`](polytope_first_golden_k3_leverage_hostile_audit.md).

TOOT: the apparent infinite golden normal-scaling problem collapses exactly.
For every full-rank crystallographic symmetrizable `Z[phi]` realization whose
`m=3,5` graph is connected, edge products `1,phi^2` force all root-length
ratios to be unit squares; every unit twist is a vertex coboundary, so the
equal-norm wall span and finite discriminant-overclosure compiler are
exhaustive.  Bugaenko is stronger: two rank-nine minors force its only dotted
entry to be `-2phi`, its first nine walls are an `O`-basis, and its integral
and even-trace saturation lists are singletons.  Hostile examples show the
scope is sharp: an infinity edge can have a nonsquare product, and one bare
infinity diagram has infinitely many arithmetic dotted realizations.  This
closes scaling/full-Gram bookkeeping, not carrier success or a polytope
bound.  See the [`golden scaling theorem`](golden_normal_scaling_universality.md)
and [`hostile audit`](golden_normal_scaling_universality_hostile_audit.md).

TOOT: the equal-root golden Lanner class is now exhaustive, and all twelve
members are literal complete chambers rather than existential subchambers.
A permutation-complete diagram search gives `4/5/3` orbits in `H2/H3/H4`;
Lanner's theorem closes `H5+`.  Exhaustive two-place support and sign-pattern
enumeration finds every norm-two mirror meeting each simplex.  The contact
counts are `9,8,10,12 / 25,27,26,28,44 / 64,114,108`, and every row has
exactly zero mixed-sign/interior-cutting root.  Deleted-wall Gram tests make
all vertices finite.  Their trace lattices pass the stable primitive K3
embedding gate, and the eigen-nef construction returns exactly the same
twelve polytopes in period dimensions `14/12/10`.  Four exact K3 returns were
missing from the earlier handpicked list.  These are new algebro-geometric
realizations of classical abstract Lanner simplices, not worldwide-new
abstract polytopes.  See the
[`exhaustive twelve-chamber theorem and replay`](golden_lanner_exhaustive_census.md).

TOOT: the inverse lane now has a genuinely uniform compact-arithmetic
`P ->` algebraic-geometry theorem, although not yet a bound.  From the even
primitive totally-real wall lattice of any marked compact arithmetic `P`, the
coordinate-free determinant line produces finitely many reflection-stable
`(2,n)` gluings and proper orthogonal Shimura stacks.  The canonical real
slice pulls facet Heegner divisors back to the exact mirrors; oriented
halfspace tests recover actual faces; marked dihedral inertia of order
`2m_ij` recovers Coxeter labels; and compactness makes hidden sheets in a
fixed Heegner shell finitely enumerable.  Bruinier's vector-valued cusp
pairing is an exact finite iff for the Whittaker/Borcherds principal-divisor
carrier.  This does **not** yet transfer to a polytope bound: automatic
correction is false in general, a uniform high-weight obstruction estimate is
missing, and facet counts require a further sheet/degree inequality.  See the
[`determinant-line transfer theorem`](polytope_to_automorphic_divisor_transfer_theorem.md).

TOOT: the first compact instance of that polytope-first theorem is now exact.
From Bugaenko's marked `H8/11` wall lattice alone, the compiler constructs
`<2phi>+E8(-1)+<4phi>` of signatures `(2,8)/(0,10)`.  Its discriminant
module is `O/(2)+O/(4)` of order 64 and has no nonzero isotropic element, so
there is exactly one even determinant-line glue package.  All eleven facet
roots are primitive `(-2)` vectors of divisibility one.  Exact closure of the
finite parabolic root systems at all 41 vertices gives the complete
real-slice shell

```text
320 = 11 facets + 309 lower-face-only contacts + 0 interior cuts.
```

The zero proves that the open chamber is not refined.  The 309 prove that
the full numerical Heegner type is not a facet-only divisor on the closed
chamber.  The separate Whittaker gate is now exact too: all eleven facets
collapse to the one standard full-group column `(phi/sqrt(5),0)`, and the
certified diagonal estimate `|E|<0.742149649308345<1` proves that its pure
norm-two principal part is obstructed.  The complete first compensator shell
is now exact too: across all 64 components the global minimum is one
three-component unit orbit of index `phi/(4sqrt(5))` and norm `1/80`.  Its
special divisor is nonempty but wholly disjoint from the real slice, and the
target-plus-shell `2x2` Poincare matrix is strictly diagonally dominant.
Hence this complete support has no target-preserving kernel.  Larger and
multi-shell supports remain open because no intrinsic support bound is known.
This is a new project-level AG arrangement realization and exact carrier
obstruction for the known compact polytope, not a new abstract example or a
transferred bound.  It also proves that a universal inverse theorem must
allow controlled compensators or an exhaustive non-Borcherds branch.  See
the [`Bugaenko determinant-line calibration`](bugaenko_h8_determinant_line_ag.md),
the [`Whittaker gate`](bugaenko_h8_whittaker_gate.md), and the
[`first compensator theorem`](bugaenko_h8_whittaker_first_compensator.md).

TOOT: the requested `P -> surface` direction now has its first exact
classification-sensitive middle theorem and continuous positive example.
For a fixed lattice-combinatorial Looijenga decoration extracted from a
rational arithmetic Coxeter datum, the intended nodal walls define a subtorus
of `Hom(D^perp,C*)`; deleting all hidden-root hypertori is an exact
necessary-and-sufficient character criterion, including the finite torsion
check.  A lifted Weyl branch makes the intended simple roots effective, and
independent equality with `C_D++` exhausts every exceptional and boundary
wall.  The first complete calibration starts with eight walls in `I_(1,4)`
and reconstructs a one-dimensional marked Looijenga-pair family `C*\{1}`.
The eight curves are exactly the complete nef walls; all reflections are
integral; all nine extreme rays are certified (`5` timelike, `4` null).
Hence this is a full-nef finite-volume noncompact arithmetic `H4` carrier
family.  It is a project-new inverse-moduli theorem and certificate, not a
worldwide novelty claim for the abstract polytope or weak del Pezzo surface.
See the
[`Looijenga inverse-carrier theorem`](looijenga_h4_inverse_carrier.md),
[`exact ledger`](../data/enumerations/looijenga_h4_inverse_carrier.json), and
[`certifier`](../scripts/certify_looijenga_h4_inverse_carrier.py).

TOOT: the fixed-decoration calibration has now become a genuinely variable
surface-free compiler.  For every primitive full-span acute crystallographic
wall system in `I_(1,4)` having at most eight square-minus-one/two walls, the
`C_triangle(4,8)` algorithm reconstructs the canonical class, exhausts all
165 anticanonical triangles and all period/effectivity branches, and decides
triangular full-nef carrier existence.  Its nonempty locus has exactly four
weighted Gram types `T6a,T6b,T7,T8`, with `6,6,7,8` facets and exact
finite/ideal profiles `5+2,4+3,5+3,5+4`; hence all four are finite-volume
noncompact `H4` Coxeter polytopes.  `T6a,T6b,T7` are project-new complete-nef
carriers.  Their abstract diagrams are not new: exact matches are Tumarkin
families `07_02`, `07_01` and Im Hof item 13/Ma--Zheng `P_(2,1)`.
All four chambers are nevertheless linked by a new exact arithmetic theorem:
they are the mod-two preimages of `W(A2+A1),W(A2),W(A1+A1),W(A1)` in
`W(A4)=S5`, of indices `10,20,30,60`.  `T8` worldwide abstract priority
remains unresolved.  The theorem is exhaustive for this carrier species,
not universal for all carriers, so rejection is not an abstract-polytope
nonexistence result.  See the
[`triangular compiler theorem`](triangular_looijenga_h4_compiler.md),
[`lineage/congruence theorem`](triangular_looijenga_h4_worldwide_lineage.md),
[`certificate`](../data/enumerations/triangular_looijenga_h4_compiler.json),
and [`certifier`](../scripts/certify_triangular_looijenga_h4_compiler.py).

**Directional Coble update.**  On `surface -> polytope`, Kondō's complex
one-boundary Types I and II have complete full-Neron--Severi nef chambers in
`H10`.  Type I has 22 facets and exact rays `315+38+0`; Type II has 30 facets
and `1333+207+0`.  Both pass literal Vinberg audits with zero failures.  Type
II is the known `D16-D5-H10-30` chamber with a new complete-nef Coble carrier.
Type I is a project-new exact chamber row distinct from the known Roulleau
`H10/22` chamber; worldwide literal priority remains open.

TOOT: the one-boundary Coble factory now has a third and substantially larger
`H10` output.  Kondō's unique characteristic-five Type-VII surface has the
complete full-nef wall split

```text
1 boundary (-4) + 20 roots (-2) + 51 cusp lifts (-1) = 72 walls.
```

The one-boundary Halphen correspondence proves that the 51 classes `E=f-b`
are all the `(-1)` curves, rather than merely candidate walls.  Exact cdd
enumerates `18,761` finite, `1,731` ideal, and zero hyperideal rays, all on
one time sheet.  Exhaustive rank-nine/rank-ten Coxeter scans recover all
`143,565` edges and all finite endpoints, and all 1,731 ideal supports have
full affine rank nine; there are zero Vinberg failures.  Hence
`COBLE-VII-p5-b1-H10-72` is an arithmetic finite-volume noncompact Coxeter
polytope and the complete nef chamber, not a partial divisor cone.  The
surface, Shimada's twenty-wall `20C=20D` `H9` chamber, Mukai's abstract
51-fibration census on a different Type-VII Coble surface, the ambient
`I_(1,10)` class, and the repeated index `777,093,120` are known.  The
project-new result is the characteristic-five one-boundary cusp-to-curve
lift and resulting complete `H10/72` face/volume certificate.  No checked
source or pinned catalogue prints the literal chamber, but worldwide
abstract priority remains open.  This is strictly `SURFACE -> POLYTOPE` and
scores `N/A/L0` for inverse leverage.  See the
[full-nef theorem](coble_type_vii_char5_one_boundary_full_ns_h10.md) and
[hostile lineage audit](coble_type_vii_char5_one_boundary_full_ns_h10_lineage_audit.md).

TOOT: the one-boundary factory has a fourth exact `H10` output and a new
largest facet count in that lane.  Kondō's characteristic-five Type-VI
surface has the complete negative-curve split

```text
1 boundary (-4) + 20 roots (-2) + 57 cusp lifts (-1) = 78 walls.
```

The labelled `Gamma_VI` roots span `E10` unimodularly, giving an integral
bridge to Shimada's known `20E` chamber without an overlattice choice.  The
one-boundary Halphen theorem makes the 57 classes `E=f-b` all the `(-1)`
curves.  Exact cdd gives `21,661` finite, `1,797` ideal, and zero hyperideal
rays.  Exhaustive scans classify 26,769,085 rank-nine and 40,619,993
rank-ten compatible masks, recover all 166,320 edges, and certify every one
of the 1,797 disconnected cusp supports with zero failures.  Thus
`COBLE-VI-p5-b1-H10-78` is a complete arithmetic finite-volume noncompact
nef chamber.  Its exact invariants are `chi=-527/6480`, index
`906,608,640` in `W(P10)`, and volume `1054*pi^5/382725`.  The surface,
`20E` seed, ambient lattice, and repeated index are known.  The project-new
positive content is the complete characteristic-five cusp-to-curve lift and
78-wall full-nef face certificate; worldwide literal priority is open.  This
is strictly `SURFACE -> POLYTOPE`, inverse score `N/A/L0`.  See the
[full-nef theorem](coble_type_vi_char5_one_boundary_full_ns_h10.md) and
[hostile lineage audit](coble_type_vi_char5_one_boundary_full_ns_h10_lineage_audit.md).

TOOT: the same surface-first lane now reaches `H11`.  Kondō Type-I with two
boundary components has a complete 26-wall full-nef chamber whose curve split
is `2*(-4)+10*(-2)+14*(-1)`.  Exact double description gives `508` finite,
`42` ideal, and zero hyperideal rays; all 3,641 elliptic edges and 106 affine
neighbourhood audits pass.  The integral arithmetic reflection lattice is
the canonical doubled extended Picard lattice
`U+E8(-2)+<-2>^2`.  This is a project-new exact complete realization with
worldwide abstract priority open, and it is explicitly only
`SURFACE -> POLYTOPE`.  See the
[two-boundary Coble theorem](coble_type_i_two_boundary_full_ns_h11.md).

TOOT, repaired and genuinely inverse: the same known `H11/26` chamber can now
be fed in only as its odd-unimodular Gram and 26 oriented walls.  Exhaustive
contraction gives 52 `P2` endpoints and exactly four intrinsic line--quintic
branches.  Their compiler-generated `23x21` jet matrix has rank-drop locus
`a^2-b^2=0`; the reduced horizontal component constructs over `Q` a line and
an absolutely irreducible rational quintic with two `A5` singularities.
Blowing up the forced `(3,2,3,3)` proximity chains gives a strong
two-boundary Coble surface.  Explicit interpolation constructs all ten
cusp-lift curves omitted by the original root seed, while Riemann--Roch gives
an independent effectivity proof.  Fresh exact replay returns all 26 walls
and `508+42+0` rays.  Therefore

```text
raw H11/26 P -> plane determinantal locus -> strong Coble surface -> same P
```

is proved for this input.  The original 16-of-26 handoff failed as written;
the repaired primary and two helper-free hostile audits pass.  The project
novelty is this line--quintic inverse/birational realization, not the known
abstract polytope or Coble class, and it is one-`P` `L3`, not variable-class
`L4`.  See the [repaired inverse theorem](coble_type_i_two_boundary_plane_inverse.md),
[AG hostile audit](coble_type_i_two_boundary_plane_inverse_hostile_audit.md),
and [raw hostile replay](coble_type_i_two_boundary_plane_inverse_raw_hostile_audit.md).

TOOT: the two-boundary lane now closes Type MI and reconnects it exactly to
an existing positive Fermat-quartic chamber.  Kondō's characteristic-three
Type-MI surface has forty complete Coble--Mukai roots and a CM chamber with
`6,210+247+0` rays.  Exactly 82 cusps lift over both boundaries, producing
164 additional exceptional curves.  Its complete negative-curve split is

```text
2*(-4) + 30*(-2) + 174*(-1) = 206 walls in H11.
```

Independent GMP double description gives `327,566` finite, `18,354` ideal,
and zero hyperideal rays.  Literal Vinberg audits 2,688,957 elliptic edges
and 60,458 connected affine components with no bad endpoint or parabolic
neighborhood.  More strongly, an explicit 206-facet permutation and a
determinant-one integral isometry identify all 42,436 Gram entries with the
already certified project chamber `F3-C2-H11-206`.  Kondō's double cover is
the supersingular Fermat quartic, so this is a project-new complete Type-MI
full-nef realization and quotient/K3 same-polytope theorem for a known
abstract polytope.  The lattice and maximal group are the known Alexeev
`(12,10,1)` class; worldwide priority for the surface interpretation is
open.  This remains `SURFACE -> POLYTOPE`, inverse score `N/A/L0`.  See the
[`Type-MI theorem`](coble_type_mi_full_ns_h11.md).

TOOT: the final finite-automorphism Coble template is now closed.  For
Kondō's characteristic-three Type-V surface, the folded twenty-root Gram
intrinsically singles out its four `(-1)`-roots as the unique four-clique of
pairing two.  The visible divisor lattice has determinant four and exactly
one integral index-two overclosure, namely `NS=I_(1,11)`.  Seven of the 29
primitive CM cusps lift over both boundaries and the other 22 over neither,
so the complete negative-curve split is

```text
2*(-4) + 16*(-2) + 18*(-1) = 36 walls in H11.
```

Exact cdd gives `6,356` finite vertices, `432` ideal vertices, and zero
hyperideal rays; every facet is irredundant.  Literal Vinberg audits all
47,929 elliptic edges (`25,372 FF + 19,172 FI + 3,385 II`) and all 1,154
connected affine-component neighbourhoods with zero failures.  All 36
reflections act integrally on the known doubled extended Picard lattice
`U+E8(-2)+<-2>^2`.  This is a project-new exact full-nef realization;
the surface, folded `H9/20` chamber, arithmetic lattice, and maximal group
are known, while worldwide priority for the literal `H11/36` chamber is
open.  It is strictly `SURFACE -> POLYTOPE`, inverse score `N/A/L0`.  With
this row, all nine classified finite-automorphism Coble surfaces have exact
full-`NS` polytope certificates.  See the
[`Type-V theorem`](coble_type_v_two_boundary_full_ns_h11.md).

Separately, on `polytope -> surface`, a no-import compiler starts from the
pre-existing D16-D5 wall lattice and uniquely recovers the characteristic
boundary, `E10` saturation, twelve roots, and every one of the seventeen cusp
lifts before invoking Kondō existence.  Across six local odd-unimodular
`H10` chamber types, 24 presentations, and 154 characteristic candidates,
Type II is the unique positive one-boundary finite-Aut Coble carrier.  This is
a scoped inverse recognition/classification theorem for known polytopes.  It
is neither a new abstract example nor a polytope nonexistence theorem.  See
the [full-NS Coble theorem](coble_one_boundary_full_ns_h10_theorem.md) and the
[independent P-first audit](polytope_first_coble_h10_hostile_audit.md).

TOOT (scoped P-first criterion, not a new abstract example): the complete
fifteen-edge, maximum-degree-three census is now frozen through three
pairing-two edges.  The surface-free denominator consists of connected
twelve-node supports with exactly fifteen edges, at most three pairing-two
edges, and pairing one on every other support edge.  Across all four exact
shards there are 3,306,816 unordered marked trials, 238 target Lorentzian
presentations, 211 weighted root classes, and 180 even-unimodular `E10`
branches.  Exact cdd/Vinberg proves

```text
finite volume iff the weighted root class is Type II.
```

The positive class has two integral markings; complete cusp lifting sends
both to the same known Type-II `H10/30` full chamber and the frozen
split-K3/Coble compiler returns the same polytope.  Every other one of the
210 weighted root classes has explicit hyperideal rays.  Timelike strict-
interior witnesses are certified even for cones with no finite extreme ray,
so the exclusion does not rely on the equal-weight ray-sum heuristic.  This
is a positive/negative Coxeter classification inside the stated finite lane,
not AG pushback, all nodot, strict L4, a new polytope, or abstract
nonexistence after allowing degree above three, four-plus pairing-two edges,
other edge counts, or other root alphabets.  See the
[`three-pairing-two theorem`](coble_nodot_15edge_three_pair2_lane.md).

TOOT: two previously surface-first Coble rows now have structural
`POLYTOPE -> finite geometry -> SURFACE -> same P` returns.

For the raw `H11/206` chamber, the two boundaries and ten intrinsic
`(1,1)`-profile `(-1)` walls contract unimodularly to `P1 x P1`.  The thirty
roots recover the Steiner/Möbius design `S(3,4,10)`.  The prime-power
Frobenius factory has

```text
B^2=-(q-1)^2, C^2=1-q, K.C=q-3,
B+B'=-2K+(q-3)(F+G),
```

so it is a strong-Coble `(-4)/(-2)` package exactly at `q=3`; the general
root/cusp compiler returns the exact `327,566+18,354+0` Type-MI chamber.

TOOT: the same Frobenius/Möbius factory has an exact `q=2` companion.  The
five graph intersections over `F4` blow up to the degree-three del Pezzo
whose anticanonical model is the smooth Hermitian cubic
`z0^2*z3+z0*z3^2+z1^2*z2+z1*z2^2=0`.  Its 27 explicitly realized curves
exhaust all exceptional classes, and a unimodular marking sends them exactly
to the 27 standard normals of the classical right-angled `P^6`, dual to
Gosset `2_21`.  Exact cdd is `72+27+0`; all 648 Vinberg edges and 27 cusp
neighborhoods pass.  In the full prime-power factory a circle has square
`-(q-1)` and divisibility one, so its reflection is integral iff `q-1 | 2`:
the only endpoints are `q=2` and `q=3`.  This is a project-new explicit
presentation and construction criterion, not a new cubic surface, abstract
polytope, `q>3` nonexistence theorem, or class-wide inverse payoff.  See the
[`q=2` Hermitian-cubic theorem](mobius_q2_cubic_p6.md).
The independent [hostile audit](mobius_q2_cubic_p6_hostile_audit.md) imports
neither the production constructor nor the project Vinberg verifier and
rebuilds the exact `72+27+0` rays, 648 edge endpoints, and 135 full affine
zero-neighborhoods with no defect.

TOOT: the known Type-II `H10/30` chamber now has a project-new strict
polytope-first plane realization.  From only the raw rank-eleven lattice and
thirty walls, all 573 sequential contraction endpoints are enumerated and
exactly two unimodular `7+3` plane branches survive.  Their incidences force
the `B3` arrangement and the Choi--Lam sextic; over `Q/C` the latter is
absolutely irreducible with `3 A3 + 4 A1`, and its prescribed resolution is a
strong Coble surface with `h0(-K)=0`, `h0(-2K)=1`.  All thirty walls transport
exactly.  An independent hostile program gives
`1540=1333 finite+207 ideal+0 hyperideal` rays and the same
`1 boundary+12 roots+17 cusp lifts` marked chamber.  This is a new complete
realization of a known finite-volume noncompact polytope--not a new equation,
abstract polytope, dimension, compact example, or class-wide inverse theorem.
It supplies a serious finite contraction/incidence/Hilbert compiler for the
future `P -> all carriers -> AG -> P` theorem.  See the
[plane inverse theorem](coble_type_ii_b3_plane_inverse.md) and
[no-import hostile audit](coble_type_ii_b3_plane_inverse_hostile_audit.md).

For the raw `H10/78` chamber, a unique ten-wall contraction frame gives the
ten edges of `K5`; its roots are exactly ten triangle lines and ten six-edge
cut conics.  Over `F5`, the modular quotient
`{sum x_i=0}/<all-ones>` realizes the points `[e_i-e_j]`.  Cubic evaluation
has zero kernel, while the sextic double-point system has a one-dimensional
kernel.  Exact factorization and saturated-Jacobian calculation prove that
its generator is geometrically integral with precisely ten ordinary nodes;
all 90 root directions avoid the boundary branches.  The blowup is a strong
Coble surface.  The raw exceptional basis has determinant one, all twenty
root vectors equal the constructed plane classes, and 120 colored incidence
isomorphisms certify the marked handoff; the same-P return is the known
`21,661+1,797+0` chamber.

The `H11/206` and `H10/78` returns are project-new P-first `L3` calibrations and positive complete
realizations of known abstract polytopes, not new abstract examples or
class-wide AG bounds.  A normalized arbitrary-input finite-geometry
recognizer now closes the extraction gate for the polar/chord, grid, Möbius,
Desargues, and Type-VII tetrahedral cores without archive or stored-Gram
lookup.  Its strict score
is still L3 because the rigid core grammar already fixes the dimensions and
Kondō gives no nonthin new pushback.  See the [seven-row audit](coble_remaining_seven_contraction_design_audit.md),
[Möbius hostile audit](coble_four_point_constructible_incidence_hostile_audit.md),
the [generic recognizer](coble_finite_geometry_raw_recognizer.md), and
[Desargues inverse theorem](coble_type_vi_char5_desargues_inverse.md).
The independent
[Desargues hostile audit](coble_type_vi_char5_desargues_inverse_hostile_audit.md)
closes the formerly profile-only handoff by reconstructing the unique
unimodular frame, all 120 colored markings, exact root coordinates, and the
full 78-wall transport.  Its verdict remains deliberately `L3`, not a
classwide theorem.

TOOT: the raw characteristic-five Type-VII `H10/72` chamber now has a third
strictly directional status in addition to its surface-first full-nef
certificate.  Starting only with the 72 walls, it has exactly five intrinsic
unimodular plane frames, each forcing a `4+6` tetrahedral design.  Over `F5`,
six conic points and four matching poles support twelve lines, four smooth
conics, and four geometrically integral nodal quartics.  The unique boundary
sextic has four `A1` nodes and six `A2` cusps (Jacobian length 16, radical the
ten centers), and its strict transform is smooth.  Every one of 120 colored
markings returns all twenty root vectors exactly; every wall transports
integrally; and the root cone has `450` finite, `51` ideal, and zero
hyperideal rays, with the 51 cusps exactly `{E+b}`.  The full same-P return is
the known `18,761+1,731+0` chamber.  This is a project-new explicit
realization and P-first L3 calibration of a known polytope, not a new abstract
polytope, new dimension, or class-wide AG theorem.  See the
[`Type-VII tetrahedral inverse`](coble_type_vii_char5_tetrahedral_inverse.md).

Boundary update (deliberately **no TOOT**): the finite-geometry extension
theorem is now rigorous after replacing the false base-point-free Bertini
step by resolved base clusters, relative closed immersion, length-one and
length-two evaluation, and explicit Hilbert bad-locus certificates.  For
every witnessed branch that independently passes the uniform Coble return,
Kondō forces boundary counts `{1,2,5,8}` and dimensions
`{10,11,14,17}`; characteristic five leaves `{1,2}` and `{10,11}`.  But
the class is indexed by existence of the coordinate/Hilbert witness rather
than forced by raw P-data, and the nine-row census found zero nonzero
extensions.  This is a corrected sufficient theorem and a precise search
boundary, not a new example, an exhaustive carrier theorem, or L4.  See the
[`corrected extension theorem`](coble_corrected_extension_pushback.md).

TOOT: the inverse-carrier track now has its first **strong surface
classification-to-polytope theorem**.  The full nef chambers of complex
Enriques surfaces with finite automorphism group have exactly seven
congruence types in `H9`: two with 12 facets and five with 20 facets.  Exact
double description gives 2,214 finite vertices, 294 cusps, and no hyperideal
ray; an independent literal Vinberg audit checks 13,182 elliptic edges and
584 full affine neighborhoods with zero failures.  Conversely, each of the
seven surface-free primitive normalized Gram types has such a full-nef
Enriques carrier.  These are strong realizations of known Kondō--Shimada
configurations—not claimed as new abstract polytopes—but the exact
seven-image transfer and complete face/cusp census are a positive project
theorem.  See the
[`Enriques seven-image transfer theorem`](enriques_finite_automorphism_nef_transfer_theorem.md)
and its independent
[`raw-data hostile audit`](enriques_finite_aut_nef_chambers_hostile_audit.md).

TOOT: the direct K3 full-root lane now gives a higher-dimensional unmarked
classification/finiteness theorem.  The ambient surface-free class
`K3Root[11,20]` consists of primitively K3-embeddable even lattices with their
complete square-minus-two root chambers; finite reflection index is not an
ambient hypothesis.  The pre-existing Nikulin--Vinberg reflective-lattice
classification, as imported by the surface sources, identifies its
finite-volume sublocus `K3RootFV[11,20]` as exactly nineteen types in ranks
11--19 and none in rank 20.  The facet lists are
`H10:12,14,22,90`, `H11:13,14,19,59`, `H12:14,17,39`,
`H13:15,20,27`, and `H14` through `H18`: `21,19,19,19,24`;
there is no `H19` type.  All nineteen are pairwise distinct by dimension and
facet count, finite-volume, and noncompact.  Periods, exact Weyl recovery, and
full-root cone completeness give every ambient K3 carrier.  On the classified
surface strata, finite automorphism group, polyhedral effective cone, finitely
generated Cox ring, finite reflection index, and finite hyperbolic covolume
are equivalent, so volume is an output rather than an input.  The
`N`-polarized coarse moduli spaces are irreducible and
the nonempty exact-Picard loci have dimensions `20-rho`; each lattice has one
primitive embedding orbit under both full `O(Lambda_K3)` and oriented
`O^+(Lambda_K3)`.  All nineteen complete wall systems now pass exact literal
Vinberg endpoint/full-neighbourhood audit: 117,845 elliptic edges, 14,098
finite vertices, 1,156 cusps, and 4,255 affine audits, with zero failures.
The fixed-component chamber orbit and oriented arithmetic period/change-of-
marking images are also classified: pointwise `N`-polarization has `T`-image
`tilde O^+(T)`, while the unmarked ambient chamber normalizer has `T`/period
image `O^+(T)`.  This is an image statement, not a group isomorphism or an
independent geometric-monodromy theorem.  These are known published
polytopes, not nineteen new examples.  The positive project result is the
completed exact carrier/return package, all-row literal volume certification,
and oriented arithmetic calibration.  It is **not strict `L4`**: the
nineteen-type exhaustion and no-`H19` boundary are inherited from the
reflective-lattice classification, not obtained by independent AG knowledge
and pushed back to a new polytope theorem.  Strict M3c remains open because the ambient class is not the
bounded raw `IC0` specification, explicit presentations of the nineteen
infinite arithmetic groups are absent, and geometric monodromy remains
conditional.  See the
[`high-rank K3 transfer theorem`](k3_high_rank_full_root_classification_transfer.md)
with its [`machine manifest`](../data/enumerations/k3_high_rank_full_root_classification.json),
the [`rank-11--14 literal audit`](k3_rank11_14_literal_vinberg_volume.md), and
the [`oriented arithmetic theorem`](k3_high_rank_oriented_monodromy.md).  The
[`strict L4 hostile audit`](k3_high_rank_l4_promotion_hostile_audit.md)
records the source-dependency obstruction and the infinite-volume carrier
control.

TOOT: the first Enriques downstairs compiler is now exact.  Type I has
27,648 possible two-sheet anti-Grams from its curve-intersection alphabet,
but only 18 sheet-switch orbits.  Exact inertia leaves one rank-nine and one
rank-ten negative-semidefinite orbit.  The rank-nine orbit is precisely the
explicit geometric `E8(-2)+<-4>` package; half-sum glue forces an integral
anti-invariant `(-2)` vector in the saturation of the rank-ten residue for
all 512 sheet assignments.  Therefore saturated rootlessness forces the
geometric orbit without a Picard-rank hypothesis.  The surviving known
Type-I surface independently passes the complete `L1`--`L4`, full-nef, and
Vinberg gates.  This is not a new abstract polytope, but it is the first
finite downstairs-to-upstairs uniqueness criterion in the surface-softening
program.  See the
[`Type-I anti-Gram compiler`](enriques_type_i_downstairs_anti_gram_census.md),
its [hostile audit](enriques_type_i_downstairs_anti_gram_census_hostile_audit.md),
and the [explicit positive package](enriques_type_i_explicit_lift_package.md).

TOOT: the second Enriques downstairs compiler is now exact and hostile-audited.
Type II has 32,768 labelled two-sheet anti-Grams but only sixteen switching
classes.  Exact signature leaves two rank-nine candidates; their complete
half-sum congruence lattices are `D9(-2)` and `<-2>^9`.  The cubic class has
eighteen anti-invariant `(-2)` roots, whereas `D9(-2)` is rootless and is the
known very-general Type-II carrier.  An independent implementation reproduces
the orbit census and both isometry classes and checks all 73,728 transports of
the 36 signed root witnesses.  It also shows that binary witnesses depend on
the sheet gauge, so full-lattice rootlessness is the exact criterion.  This
recovers a known surface and known `H9/12` polytope, but it is a new finite
downstairs-to-upstairs uniqueness theorem.  See the
[`Type-II anti-Gram compiler`](enriques_type_ii_downstairs_anti_gram_census.md)
and its
[`hostile audit`](enriques_type_ii_downstairs_anti_gram_census_hostile_audit.md).

TOOT: the K3 inverse-carrier track now has its first exact Gram-only
datum compiler.  Starting solely from the pinned `T2/C01 H9/14` 14-facet
Gram, independent implementations enumerate all 18 minimal even
crystallographic scalings and all `16+5` reflection-stable even overclosures.
The long-class, rootless-sheet, and K3-complement gates force the unique fixed
lattice `U+E8(-2)`.  Explicit graph glue constructs a rank-18 `(18,4,0)` K3
lattice with involution and 24 square-minus-two source roots folding exactly
to the input.  Exhausting all `2^10` sheet orders gives four integral marked
isometries to the archived `C01` action; the hostile witness checks the form,
involution, all ordered roots and folded sums, primitive eigensublattices, and
the action on all 48 ambient walls.  Thus datum-level `IC0/R` are closed and
the existing strong promotion transports complete-cone and Vinberg `C/V`
(`90` finite and `20` ideal vertices).  This is a new inverse-carrier
criterion and realization of a known project polytope, not a new abstract
polytope or class-wide universality theorem; `U` remains open.  See the
[`surface-free unfolding pilot`](h9_14_t2_surface_unfolding_pilot.md),
[`Gram-only hostile audit`](a18_c01_t2_lattice_precompiler_hostile_audit.md),
and [`marked-bridge hostile audit`](a18_c01_t2_marked_bridge_hostile_audit.md).

TOOT: the same arithmetic compiler now reaches a second nonisomorphic
`H9/14` diagram.  In the pinned primary six-marking/four-type corpus, the
surface-free `D34` conditions retain exactly T1 and T2.  A corrected HNF root-
lattice quotient exhausts 36 scalings and 28 stable even overclosures, and
both diagrams uniquely force `U+E8(-2)`.  T1 yields a new weak rank-18 K3/
Nikulin package of type `(18,2,0)`, determinant `-4`, complement `U+U(2)`,
and 24 acute square-minus-two roots folding exactly to its known 14-facet
Gram.  The roots span index two, recorded explicitly.  This is a positive new
algebraic carrier lineage for the known T1 polytope, not a new polytope or a
class-wide theorem: the denominator is corpus-relative, formal `U` is open,
and the strong T1 promotion is in fact impossible for this compiled carrier.
Five explicit exchanged-root hidden cuts cross the proposed cone; the first
has exact ray pairings `-6,+6`, so K3 Riemann--Roch forces a missing nef
inequality.  This leaves the weak carrier intact and excludes no different
algebraic T1 realization.  See the
[`D34` two-diagram census](h9_14_t1_t2_d34_unfolding_census.md), its
[`hostile replay`](h9_14_t1_t2_d34_unfolding_census_hostile_audit.md), the
[`T1` obstruction](h9_14_t1_strong_k3_nef_obstruction.md), and the general
[`exchanged-root criterion`](invariant_nef_exchanged_root_obstruction.md).

TOOT: the D34 compiler is now genuinely variable.  From two unlabelled
fourteen-vertex supports it exhausts 37,180 raw eight-`m4` labellings, 3,610
support orbits, and 748 scaled-dual integral Gram lanes.  Exact rank/inertia
leaves precisely one diagram on each support, deriving T2 and T1; all 28
stable even overclosures and both weak K3 packages are reconstructed without
access to the old Gram corpus.  Fresh literal Vinberg gives `90+20` and
`112+15` finite/ideal endpoints.  This closes arithmetic and weak-surface
universality for a nontrivial raw class.  Better, the oriented `D16+`
subclass has one marked Gram orbit and one viable strong package, so the exact
T2 bridge and invariant-nef theorem close marked `U/R/C/V`.  Its adjacent
`D17+` class has a nonempty T1 polytope but an empty minimal strong stratum,
making the exchanged-root cut a genuine criterion boundary.  These are new
transfer/classification theorems and reconstructions of known polytopes, not
new abstract diagrams.  See the
[`variable cut denominator`](h9_14_d34_cut_denominator.md), its independent
[`no-import hostile replay`](h9_14_d34_cut_s16_s17_hostile_audit.md), and the
[`S16/S17 transfer theorem`](h9_14_s16_strong_transfer_and_s17_countercalibration.md).
The hostile replay independently certifies the raw support/orbit/748-lane
arithmetic and full `V`; the overclosure/K3/Torelli/nef and T1 rank-18
obstruction layers remain transported from hash-pinned certificates.

The discriminant-eleven K3 fishing lane has a new exact positive checkpoint.
Five escape-directed `V4` orbits add twenty effective square-minus-two
divisor walls to the previous 309-wall subsystem, giving **329** pairwise
Coxeter walls.  A full independent replay reconstructs all 35 galleries and
2,687 Shimada-atlas transitions and classifies every one of the 53,956
pairings.  One exact hyperideal ray still survives.  Thus this is genuine
algebro-geometric construction progress and a sealed resumable search state,
but it is not counted as a new finite-volume polytope.  See the
[`329-wall directed K3 search`](shimada_disc11_directed_329.md).

The project now has **twenty-nine distinct K3/divisor chamber-realization
slots**.  The previous twenty-one consist of sixteen earlier special-period,
complete-nef-fixed-section, or mixed-divisor realizations plus the already
known finite-volume noncompact `H10/16`, `H10/45`, and `H11/29` polytopes,
the `A15_11` `H10/21`, and the `SS2-2A` `H13/50` chamber.  The new rank-ten
factory adds three ambient section/cone orbits and three `H9/14` diagrams.
The first two are apparently new after a bounded lineage audit; the third is
exactly Bredon's `P_star`, with a new project-level K3 divisor-lattice
realization.  All three pass exact cone, time-sheet, integral-reflection, and
literal Vinberg gates and are finite-volume noncompact.  They are weak
full-Picard divisor carriers at the factory stage; subsequently `T2/C01` and
`T3/C07` were promoted to genuine automorphism-fixed, full invariant-nef
cones.  The first is carried by a symplectic Nikulin involution with
coinvariant `E8(-2)`; the second by a purely non-symplectic order-three
automorphism with coinvariant `-(A2 tensor D4)`.  Both are therefore nef
cones of normal quotient surfaces under finite pullback.  `T1` remains weak.
See the [`H9/14` factory](alexeev_k3_h9_14_factory.md) and the
[`C01/C07 strong promotion`](a18_c01_c07_k3_nef_promotions.md).

TOOT: the known `T3=P_star` row now has a genuinely polytope-first surface
lineage.  Bredon's four-arm Coxeter diagram forces two integral long/short
orientations.  One has Gram lattice `U+A2(-1)^4`, determinant `-81`, Smith
`1^6 3^4`, and an index-one 24-root `D4`-triality unfolding equal to the
archived rank-18 `(18,4,0)` K3 model; period/Torelli and hidden-root checks
return exactly the full invariant-nef and normal-quotient-nef cone.  The
reversed unfold has rank 28 and inertia `(25,3,4)`, so K3 geometry really
selects the surviving arithmetic branch.  Exact cone/Vinberg replay gives
`128+6+0` vertices and zero failures.  This is a new project-level
`P -> K3 -> same P` realization of Bredon's published-new abstract polytope,
not a new polytope or arithmetic group.  Its strict score is only `L3`:
Bredon's one-`m6` and three-`m6` triangles obstruct a uniform integral
triality unfolding across the 24-member ADEG class, and no new AG-derived
class-wide polytope theorem follows.  See the
[`P_star` polytope-first audit](bredon_pstar_h9_14_polytope_first_k3_triality_audit.md).

TOOT: the ADEG stock search now has eight additional exact polytope-first
complete-nef divisor-lattice returns.  A primary-source compiler adds 17
Roulleau rank-5--8 complete K3 nef stocks to the former 38-stock local
denominator.  Five new matches join the earlier three:

```text
N2-01/o1 -> [4]+D4(-1)
N3-01/o0 -> U+D4(-1)
N3-04/o0 -> U(2)+D4(-1)
N3-05/o1 -> U+A2(-1)^3
N5-02/o0 -> U+A2(-1)+D4(-1)
N5-01/o0 -> U+D4(-1)^2
N6-01/o1 -> U+D4(-1)+A3(-1)
N7-02/o0 -> U+D4(-1)^2.
```

Every weighted-Gram match extends unimodularly, preserves a strict invariant
ample class, uses a published complete nef-wall stock, and folds by orbit
sums to exactly the input ADEG Gram.  The residual 30-package ledger is
`P=8,R=15,O=7`.  These are positive inverse calibrations of known polytopes,
not new abstract chambers or the requested class-wide AG pushback theorem;
all eight transported trialities act nontrivially on discriminant groups, so
no K3 automorphism is claimed.  See the
[rank-5--8 stock compiler](roulleau_rank5_8_matching_full_nef_stocks.md),
[eight-promotion theorem](bredon_adeg24_remaining_k3_stock_promotions.md),
and [30-package ledger](bredon_adeg24_remaining_k3_stock_promotion_ledger.md).

Five further slots come from the complete parked fixed-rank-ten denominator.
Together with `T1`, that denominator has exactly six rows, and all six give
finite-volume noncompact `H9` sections.  The five new rows have 12, 16, 13,
15, and 15 facets.  They are pairwise nonisomorphic and distinct from the
factory and doubling-cube diagrams.  The 12-facet row is exactly Im Hof's
classical `9-imhof85_03`, now with an Alexeev/K3 divisor-lattice extraction;
the other four are bounded-corpus unmatched, with worldwide priority open.
All five are weak full-Picard divisor carriers.  Four also retain their
original archived lattice involution on genuine symplectic Nikulin K3 action
carriers; this does not make their walls effective or their cones full-nef.
An independent hostile replay confirms every ray and
Vinberg gate and shows that the five unique maximal stable even overclosures
belong to known 2-elementary arithmetic classes; hence there is no new
maximal arithmetic group claim.  See the
[`parked rank-ten H9 theorem`](alexeev_parked_rank10_h9_examples.md).
The maximal envelopes are explicitly
`U+D4(-1)^2`, `U(2)+D4(-1)^2`, and `U+E8(-2)`, with their classical
ten-dimensional nonsymplectic-K3-involution moduli and fixed loci recorded
in the
[`overclosure/K3 lineage theorem`](alexeev_parked_rank10_h9_lineage_overclosures.md).
That maximal-envelope interpretation is conditional.  A separate exact
rank-18 completion now promotes the first four original involutions to
genuine Nikulin actions with rootless coinvariant `E8(-2)`; the fifth is
obstructed only for that completion route.  See the
[`parked K3-action theorem`](alexeev_parked_h9_k3_action_carriers.md).

Four parked rows are positive amplifiers.  One complete `C2^3` cube adds
three finite diagram types with 14--16 facets.  Two parallel-pair seeds give
infinite path-times-cube families with facet/index formulas
`4N+13, 2N` and `3N+14, 4N`.  A fourth seed has two commuting parallel pairs
and gives, for every `m,n>=1`, an `H9` polytope with `mn+m+n+12` facets and
volume multiplier `mn`.  Six representatives pass direct GMP and literal
Vinberg audits.  A no-import hostile companion proves the group normal forms,
wall stabilizers, collision exclusions, and facet formulas for all integer
parameters, rather than extrapolating from the finite representatives.  These
descendants stay in the same five carrier slots, so the twenty-nine-slot
count does not increase.  See the
[`H9` product-redoubling theorem](alexeev_parked_h9_doubling_products.md) and
its [all-parameter hostile audit](alexeev_parked_h9_doubling_products_hostile_audit.md).

The full-Picard Jacobian--Kummer `H12/172` carrier is now a second major
amplifier.  Its complete census has 124 doubling walls.  Three single
doubles have 296, 317, and 326 facets.  W31/W36 generate `D_infinity`, and
for every `N>=1` their finite path has exactly `124N+48` facets and volume
multiplier `N`; a separate no-import symbolic audit proves the formula for
all integral translation exponents.  Independently, the orthogonal walls
W31/W58/W59/W68 form `C2^4`.  Its 15 nonempty marked subsets contain at least
seven nonisomorphic types, distinguished by facet counts
`296,317,523,559,941,1004,1827`, with volume multipliers `2,4,8,16`.
All are finite-volume noncompact and remain weak full-Picard divisor carriers
on the same seven-dimensional K3 family; no nef/effective/action claim is
made.  Worldwide priority is open.  See the
[`H12` path theorem](jacobian_kummer_h12_single_doublings.md),
[`all-parameter hostile audit`](jacobian_kummer_h12_parallel_redoubling_hostile_audit.md),
and [`C2^4` cube theorem](jacobian_kummer_h12_doubling_cube.md).

A still larger `H12` product is now production-certified.  Eight exact
parallel-pair-times-three-cube configurations give faithful
`D_infinity x C2^3` actions.  For the canonical W68/W63 times
W31/W58/W59 configuration, every `N>=1` gives a finite-volume noncompact
polytope with exactly `886N+55` facets and volume multiplier `8N` on the same
weak K3 carrier.  A no-import audit solves 3,526,884 signed wall equations
over all integer translations, proves full Picard span, and shows the eight
listed configurations plus both phases form one isometry class for fixed
`N`.  Direct exact reconstruction through `N=8` reaches 7,143 facets.  See
the [`H12` path-times-cube theorem](jacobian_kummer_h12_path_cube_family.md)
and its
[`all-parameter hostile audit`](jacobian_kummer_h12_path_cube_hostile_audit.md).

TOOT: the same stock contains exactly two direct `D_infinity^3`
configurations, conjugate under the seed chamber symmetry.  For every
`N1,N2,N3>=1` their finite box has volume multiplier `N1*N2*N3` and

```text
facets = [177*N1*N2*N3 + 29*(N1*N2+N1*N3+N2*N3)
          + 13*(N1+N2+N3) + 40 + epsilon]/2.
```

Here `epsilon` is phase-refined: if at least one `Ni` is odd it equals one in
every phase; if all three are even it equals zero or two according to the
endpoint-orientation parity.  Thus `(2,2,2)` has 941 or 942 facets.  The
certificate solves all 440,896 signed wall equations over `Z^3`, replays 240
finite boundaries, proves full Picard span for all eight phase cores, and
shows that no `D_infinity^4` or extra orthogonal `C2` factor exists inside
this exact 124-wall stock.  This maximum-product statement is scoped to
direct products of parallel pairs and orthogonal factors; arbitrary
redoubling trees remain open.  The family is finite-volume noncompact and a
weak divisor-reflection carrier only.  Worldwide priority is open.  See the
[`H12` triple-product theorem](jacobian_kummer_h12_triple_parallel_family.md).

TOOT: the unrestricted finite-subtree amplifier is substantially larger.
Independent reconstruction from Dolgachev--Shimada's raw 316-wall data finds
a unique maximum W61--W172 of **112 pairwise disjoint doubling walls** and a
unique maximum W77--W172 of **96 pairwise ultraparallel doubling walls**.
The full maximal-clique size censuses are

```text
disjoint:       36^9, 55^36, 74^42, 93^12, 112^1,
ultraparallel:  32^32, 57^12, 61^16, 96^1.
```

Separate conflict matchings prove the upper bounds `124-12=112` and
`124-28=96`.  Allcock 4.1 turns every finite subtree of the resulting
112-regular Cayley tree into a finite-volume noncompact `H12` Coxeter
polytope with exact index/volume multiplier `|V(T)|`; Allcock 4.3 supplies
exponentially many nonisometric choices.  The fixed 60-wall exterior core
W1--W60 has Smith diagonal `1^13`, so every descendant boundary spans the
same full Picard lattice.  Allcock's abundance mechanism and prior `H12`
existence are classical; the exact surface-derived maxima and carrier
inheritance are the project gain, with worldwide priority open.  See the
[`H12` maximal-disjoint theorem](jacobian_kummer_h12_maximal_disjoint_doubling_family.md)
and its [independent raw-Gram hostile audit](jacobian_kummer_h12_maximal_disjoint_doubling_family_hostile_audit.md).

TOOT: the first complete non-product slice of that 112-regular factory is
now explicit.  Reconstructing the 384-element folded chamber action gives
two orbits on the 112 generators and exactly 35 orbits on all 6,216
unordered generator pairs.  Every one-edge double was rebuilt: they are
exactly the already recorded 317- and 326-facet types.  Every three-vertex
colored path was then classified through its orbit representative, producing
**35 pairwise nonisomorphic finite-volume noncompact `H12` polytopes**:

```text
4 types with 462 facets,
7 types with 471 facets,
24 types with 480 facets.
```

The nonisomorphism proof uses only the abstract Coxeter diagrams: the sorted
multisets of exact normalized edge-label degree profiles at the facets are
all distinct.  They remain distinct even after parallel and ultraparallel
pairs are collapsed to the ordinary Coxeter-matrix entry `infinity`, so fine
dotted-edge distances play no role.  Direct reflected-root, Coxeter-pair,
integrality, strict-interior, and Smith-`1^13` checks retain the same weak
full-Picard K3 carrier.
The facet counts `462/471/480` are new to this project and absent from the
pinned CoxIter lanes; Allcock's `H12` abundance is classical and worldwide
priority for these exact diagrams remains open.  See the
[`H12` complete small-subtree census](jacobian_kummer_h12_small_subtree_census.md).

TOOT: the entire four-vertex layer is now exact.  Among 924,112 unrooted
properly colored trees before the seed quotient, the 384-element action gives
1,964 path orbits and 680 star orbits.  Direct integral reconstruction proves
that all **2,644** resulting ordinary Coxeter matrices are pairwise
nonisomorphic finite-volume noncompact `H12` polytopes, with index and volume
multiplier four and facet histogram

```text
607: 23,   616: 170,   625: 844,   634: 1,607.
```

An independent raw-source audit reconstructs the seed/action, orbit counts,
all boundary and relation hashes, and the Smith-`1^13` carrier core.  Degree
profiles distinguish 2,535 types; a permutation-stable weighted two-walk
invariant separates every one of the remaining 102 collision groups.  The
focused hostile suite passes in full.  These are weak divisor-reflection
carriers on the common K3 family, not nef/effective/action chambers.  Allcock's
exponential abundance is classical and worldwide priority for these exact
types remains open.  See the
[`four-vertex census`](jacobian_kummer_h12_four_vertex_subtree_census.md) and
its [`hostile audit`](jacobian_kummer_h12_four_vertex_subtree_hostile_audit.md).

TOOT: a first five-chamber pass adds five more explicit positive `H12`
examples.  Properly colored paths with `B=0,1,2,3,4` have respectively
`752,761,770,779,788` facets.  Direct raw-seed reconstruction verifies the
five chamber transforms, four internal mirrors, primitive integral boundary
roots, strict timelike interiors, complete ordinary Coxeter relations, and
the fixed Smith-`1^13` carrier core.  Allcock gives index and volume
multiplier five; all five polytopes are finite-volume noncompact and pairwise
nonisomorphic by facet count.  This is a pilot, not a full five-vertex census
or a proof of the displayed facet formula for all colorings.  The common K3
carrier is weak, and worldwide priority remains open.  See the
[`five-chamber pilot`](jacobian_kummer_h12_five_chamber_pilot.md).

TOOT: the complete five-vertex `K1,4` layer contains **16,872 pairwise
nonisomorphic** finite-volume noncompact `H12` Coxeter polytopes.  Burnside and
explicit orbit cover agree on 16,872 seed-action orbits among 6,210,820 raw
proper colorings.  Every exact boundary passes ordinary-Coxeter, primitive-
integral, strict-interior, and Smith-`1^13` gates; Allcock gives index and
volume multiplier five.  Exact degree profiles isolate 16,355 types, and
canonical color refinement separates every member of the 498 residual
collision groups.  The formula `f=752+9B` holds on the whole layer, with

```text
752^19, 761^169, 770^1582, 779^6117, 788^8985.
```

The full five-vertex color denominator is 417,981 orbits—201,747 paths,
16,872 stars, and 199,362 forks.  The later theorems below close the entire
path census, star--`P5` comparison, and all five fork layers.  Thus all 417,981
rows are classified and no five-vertex denominator remains.
The five headers were already visible in the pilot; the complete 16,872-type
star theorem is the new positive content.  The common K3 carrier remains weak
and worldwide priority is open.  See the
[`five-vertex star census`](jacobian_kummer_h12_five_vertex_star_census.md)
and [`all-shape denominator`](jacobian_kummer_h12_five_vertex_tree_denominators.md).

TOOT: the first complete five-vertex `P5` tranche contains **20,693 pairwise
nonisomorphic** finite-volume noncompact `H12` Coxeter polytopes.  The seven
path-reversal shards with at most two `B` edges cover 7,562,616 raw proper
colorings.  Exact degree profiles leave 629 collision groups and 1,293
members; canonical color refinement separates every residue.  Primitive
integral reflections, strict interiors, the fixed Smith-`1^13` core, and
Allcock's index/volume multiplier five pass throughout.  The facet histogram
is

```text
752^142, 761^2091, 770^18460.
```

The exact increment beyond the five earlier path pilots is 20,690 types.  The
next paragraphs record the completed `ABBB`, `BABB`, and `BBBB` shards,
closing the complete 201,747-row `P5` denominator.  At this historical
checkpoint star-vs-`B<=2` SHA screens were negative but were not a pure exact
cross-shape theorem; the later literal replay and complete-fork theorem close
the full 417,981-row denominator.  The common K3 carrier
remains weak and worldwide priority is open.  See the
[`P5`, `B<=2` census](jacobian_kummer_h12_five_vertex_p5_b_le_2_census.md).

TOOT: the complete `ABBB` path shard adds **36,247 pairwise nonisomorphic**
finite-volume noncompact `H12` Coxeter polytopes, every one with 779 facets
and reflection-index/volume multiplier five.  The exact reduction is
`13,862,400` raw proper colorings to 36,247 seed/path-reversal orbits; all
2,068 members of the 1,023 initial invariant-collision groups are separated
by exact discrete refinement.  Facet count proves exact separation from the
completed `B<=2` types, giving 56,940 completed `P5` types in total and an
`ABBB` increment of 36,246 beyond its prior pilot.  The following `BABB` and
`BBBB` theorems raise the total again and close the complete `P5` census.
The carrier is weak full-Picard only and worldwide priority is open.  See the
[complete `ABBB` census](jacobian_kummer_h12_five_vertex_p5_abbb_census.md).

TOOT: the complete streaming `BABB` shard adds **36,648 pairwise
nonisomorphic** finite-volume noncompact `H12` Coxeter polytopes.  Its
14,008,320 raw proper colorings reduce to 36,648 seed/path-reversal orbits;
all 4,015 members of the 1,950 exact degree-collision groups are replayed and
separated by discrete canonical refinement.  Every row has 779 facets,
primitive integral reflections, a strict timelike interior, Smith `1^13`, and
reflection-index/volume multiplier five.  Digest screening has logical force
only for inequality here: there are zero common degree digests against all
36,247 `ABBB` rows and all 6,117 same-header star rows, so no literal equality
candidate exists.  At this stage the completed `P5` stock had **93,588** types
and the exact star+`ABBB`+`BABB` union had **89,767** types; the following
`BBBB` theorem closes the path census.  This is new-to-project positive output on a weak
full-Picard K3 carrier, not a nef/action theorem, and worldwide priority is
open.  See the
[complete `BABB` census](jacobian_kummer_h12_five_vertex_p5_babb_streaming_census.md).

TOOT: the final streaming `BBBB` shard closes the complete five-vertex `P5`
census.  Its 41,154,000 raw proper colorings reduce by explicit cover and
Burnside to **108,159** seed/path-reversal orbits.  Every reconstructed
788-facet boundary passes the primitive/integral reflection, strict-interior,
Smith-`1^13`, ordinary-Coxeter, and Allcock index/volume-multiplier-five
gates.  Exact degree profiles leave 6,366 literal collision groups containing
13,590 rows; exact canonical refinement separates them all.  Against the
8,985 same-header star rows there is one common literal degree-profile pair,
and exact refinement separates that sole pair.  Hence all **201,747** `P5`
rows are pairwise nonisomorphic, with facet histogram
`752^142,761^2091,770^18460,779^72895,788^108159`.  The exact
star+`ABBB`+`BABB`+`BBBB` high-`B` union has 197,926 types; star versus
`B<=2` is closed by the later literal-byte theorem below.  At that checkpoint
only the 106,385-row fork `B=4` matrix remainder remained; the complete
five-vertex theorem below has since closed it.  These are weak
full-Picard K3 divisor carriers, not nef/action chambers; worldwide priority
remains open.  See the
[complete `P5` census](jacobian_kummer_h12_five_vertex_p5_complete_streaming_census.md).

TOOT: the last all-`B` fork layer contributes **106,385 pairwise
nonisomorphic** finite-volume noncompact `H12` Coxeter polytopes, all with 788
facets and reflection-index/volume multiplier five.  Its 40,720,800 raw
proper colourings form exactly 106,385 seed-action orbits.  Exact refinement
separates all 9,410 members of 4,559 internal collision groups and all 5,536
literal equal-degree pairs against 117,144 same-facet predecessors (`P5 BBBB`
108,159 and star `B=4` 8,985).  A hostile terminal replay aligns every
archived predecessor colour/digest pair, rehashes all 106,385 compact rows,
corrects the discrete-refinement diagnostic to 13,945, and passes archive
equality.

This closes the fork at **199,362 pairwise nonisomorphic** types and gives the
complete three-shape theorem

```text
201,747 complete P5 + 16,872 complete star + 199,362 complete fork
= 417,981 pairwise nonisomorphic H12 polytopes.
```

No five-vertex denominator remains.  This is new-to-project positive output
on the same weak seven-dimensional full-Picard K3 divisor carrier.  Worldwide
novelty/priority and full-nef/effectivity/action remain open.  See the
[complete five-vertex census](jacobian_kummer_h12_five_vertex_three_shape_complete_census.md).

Historical superseded milestone: the complete degree-three-fork `B=3` layer
added **72,758 pairwise
nonisomorphic** finite-volume noncompact `H12` Coxeter polytopes, all with
779 facets and reflection-index/volume multiplier five.  Exact Burnside and
explicit orbit-cover counts are

```text
BB|AB 36,247 + BA|BB 18,055 + AB|BB 18,456 = 72,758.
```

All primitive/integral, strict-interior, Smith-`1^13`, ordinary-Coxeter, and
Allcock gates pass.  Exact refinement separates all 4,043 members of 1,968
internal degree-collision groups.  The complete same-header cross-shape audit
uses 79,012 predecessor rows (`P5 ABBB` 36,247, `P5 BABB` 36,648, star
`B=3` 6,117); literal replay finds 2,690 equal-degree pairs and refinement
separates all 2,690.  A hostile terminal replay additionally pins the exact
production/import hashes, verifies every archived predecessor colour/digest
alignment, and rehashes all 72,758 compact rows.  Consequently

```text
201,747 complete P5 + 15,102 star B>=3 + 92,977 fork B<=3
```

is an exact **309,826-type** pairwise-nonisomorphic union.  At that checkpoint
the complete fork `B<=3` census had **92,977** types and only `B=4` with
**106,385** rows remained.
The subsequent exact comparator rebuilds the 1,770 formerly omitted star
`B<=2` rows and all 20,693 same-header `P5` rows.  Literal degree-profile byte
membership finds zero matches.  Hence complete star plus complete `P5`
contains **218,619 pairwise nonisomorphic** types, and

```text
201,747 complete P5 + 16,872 complete star + 92,977 fork B<=3
```

is an exact **311,596-type** three-shape union.  This historical checkpoint is
superseded by the complete 417,981-type theorem above.  It was new-to-project
positive output on a weak K3 divisor carrier;
worldwide priority, full-nef/effectivity/action, and worldwide novelty remain
open.  See the
[fork `B=3` census](jacobian_kummer_h12_five_vertex_fork_b3_complete_census.md)
and the [complete low-`B` star--`P5` separation](jacobian_kummer_h12_star_p5_b_le_2_exact_separation.md).

Historical milestone retained: the first complete degree-three-fork matrix
layer added **20,219
pairwise nonisomorphic** finite-volume noncompact `H12` Coxeter polytopes:
116 with 752 facets, 2,042 with 761 facets, and 18,061 with 770 facets.
Every exact matrix passes the primitive/integral, strict-interior,
Smith-`1^13`, ordinary-Coxeter, and Allcock index-five gates.  Exact
refinement separates all 637 members of the 317 internal degree-collision
groups and all 299 literal equal-degree pairs against completed same-header
star/`P5` rows.  Hence

```text
201,747 complete P5 + 15,102 star B>=3 + 20,219 fork B<=2
```

is an exact **237,068-type** pairwise-nonisomorphic union.  The independent
twelve-pattern denominator partitions all 199,362 fork orbits as
`116+2,042+18,061+72,758+106,385`, leaving at that checkpoint the exact
179,143-row `B=3,4` matrix census; the stronger theorem above has since closed
`B=3`.  This is new-to-project positive output on a weak K3 divisor
carrier; worldwide priority, full-nef/effectivity/action, and the omitted
small-header star--`P5` comparison remain open.  See the
[fork `B<=2` census](jacobian_kummer_h12_five_vertex_fork_b_le_2_complete_census.md).

TOOT (theorem-facing, not a new polytope): every marked arithmetic compact or
finite-volume Coxeter chamber now has a literal real-marked Gram--Belyi
surface.  The scaling-invariant squared Gram entries are encoded as marked
points on minimal normalized-Belyi leaves, glued to one Belyi backbone; the
product with `P1` returns the complete geometric unit-normal Gram and same
chamber.  Its natural ample divisor satisfies

```text
A_P^2 = 2 (binom(f,2) + beta_GB(P)),   beta_GB(P) >= f(f-1).
```

Hence bounded absolute Hilbert degree gives a finite effective chamber and
Vinberg-field list without ABSW, Prasad, or Mostow.  The exact `H10`
prime-index redoubling family proves `beta_GB>=p` for every odd prime `p`, so
fixed dimension, field, Picard lattice and root alphabet do not supply the
missing bound.  The result is a hostile-audited universal `P -> surface ->
same geometric P` theorem with genuine conditional AG pushback; the open
completion is maximality or another intrinsic condition forcing bounded
`A_P^2`.  It does not recover primitive integral root lengths without an
extra lattice marking.  See the
[Gram--Belyi theorem](gram_belyi_coxeter_dessin_theorem.md) and
[hostile audit](gram_belyi_coxeter_dessin_theorem_hostile_audit.md).

The first exact local calibration makes the invariant concrete:
\(m=3,4,6,\parallel\) contribute \(5,3,7,2\), while \(m=5\) contributes
at least \(29\).  The \(m=5\) lower bound is the exact \(29\)-adic obstruction
for \(\lambda=(7+2\sqrt5)/29\); a factored composed Belyi polynomial gives a
finite, explicitly certified upper endpoint, while the exact minimum remains
open.  The resulting global weighted-edge bound is
\[
\beta_{\rm GB}\ge2\binom f2+3N_3+N_4+27N_5+5N_6
\]
before dotted excess.  See the
[standard-packet calibration](gram_belyi_standard_packet_calibration.md).

TOOT: the complete star layer and `ABBB` shard are now separated by a literal
invariant, not a hash assumption.  Their only common facet header is 779:
the certifier rebuilds all 6,117 star rows and 36,247 `ABBB` rows as sorted
7,790-byte facet degree-profile strings and finds zero exact intersection.
All other 10,755 star rows have different facet counts.  Since each source
layer already has an exact internal nonisomorphism proof, their union contains
exactly **53,119 pairwise nonisomorphic** finite-volume noncompact `H12`
Coxeter polytopes, all of reflection-index/volume multiplier five.  This is
project-scoped positive accounting on a weak K3 carrier; it does not establish
worldwide priority and it does not separate the star layer from `B<=2`.  The
full independent archive replay passed all 42,364 common-header rows.  See the
[exact star--`ABBB` theorem](jacobian_kummer_h12_star_abbb_exact_separation.md).

The factory also produces nine new-to-project abstract `H9` entries.  Each
of `T1` and `T2` has four mutually orthogonal doubling walls; the 30 nonempty
marked subsets collapse to nine pairwise nonisomorphic labelled diagrams
with 15--18 facets and exact volume multipliers `2,4,8,16`.  Every type is
finite-volume and noncompact.  A hostile audit independently replayed the
iterated convex gluing and all exact cone/Vinberg gates.  The diagrams are
absent from the pinned ordinary and doubling CoxIter corpora, while worldwide
priority remains open.  They inherit weak divisor carriers on the same
saturated fixed lattices; the doubled `T2` chambers are not nef cones because
the construction crosses nef walls.  See the
[`H9 doubling-cube theorem`](alexeev_h9_doubling_cubes.md).

Of an earlier five-addition batch, the first two are effective finite-index
divisor subsystems on order-four K3 surfaces.  The third is an explicit full-rank
mixed divisor-reflection subsystem in `U+D8(-1)+A1(-1)^2`.  The fourth has
full Picard carrier `U+E8(-2)+<-2>` on a very general K3 and a named
invariant-divisor occurrence on Shimada's `X_2`.  The fifth has a scaled-dual
full Picard carrier `U+D8(-1)+D4(-1)` on a complex K3.  None of these five is
claimed to be the full nef chamber of its ambient Picard lattice.

The examples now support a structural positive theorem rather than only a
carrier census.  For every finite-sided spherical face of Conway's chamber
whose saturated root complement passes the exact K3 discriminant gate, the
sign-reversed complement is the full Néron--Severi lattice of a projective K3
surface, every induced Coxeter wall is an integral divisor hyperplane after
primitive rescaling, and the imported cone contains an ample class.  All 21
types in Allcock's published ladder from `H2` through `H19` pass, with

\[
 \dim_{\mathbf C}{\cal M}_{K3}+\dim_{\mathbf R}P=19.
\]

Congruent rooted finite redoublings stay on the same marked carrier family.
At the zero-dimensional H19 endpoints, the carrier lattices specialize
exactly to Vinberg's classical determinant-three and determinant-four K3
Picard lattices `U+E8(-1)^2+A2(-1)` and
`U+E8(-1)^2+A1(-1)^2`; Borcherds's E6 and D6/D7 complement
calculations give the historical identification.
This is the first frozen forward `polytope => surface` universality theorem
for a class defined without surfaces.  The root-complement/K3 mechanism has
classical Borcherds--Kondō lineage; the uniform ladder, exact gate, 21-row
certificate, dimension identity, and redoubling synthesis are the project
contribution.  Outside the separately integral `D6+D6` case, it is a
divisor-hyperplane rather than an integral-reflection, effective-curve, or
full-nef theorem.  See the
[Conway-face K3 inverse-carrier theorem](conway_face_k3_inverse_carrier_theorem.md).

The Enriques calibration has also become a complete variable-surface theorem.
For an arbitrary complex Enriques surface `Y`, Kondō's reflection fundamental
domain and Torelli index theorem give

\[
 \operatorname{vol}_{\rm hyp}P(\operatorname{Nef}Y)<\infty
 \Longleftrightarrow |\operatorname{Aut}Y|<\infty.
\]

Artebani--Hausen--Laface identify the same condition with finite generation
of `Cox(Y)`, while Kondō--Martin force exactly the seven complete nodal
types.  The independent exact audits therefore classify **all** full-nef
finite-volume Enriques outputs: two 12-facet and five 20-facet noncompact
`H9` chambers, with no compact, eighth, or 14-facet type.  This is an
algebra-to-covolume theorem rather than a new abstract-polytope claim; see the
[full classification](enriques_full_nef_finite_volume_classification.md).

Several rows have already gained further positive realizations.  The
`A15_11 H10/21` carrier now comes from an explicit
nine-parameter quadratic base change with `I2+22I1` and exact very-general
`NS=U+E8(-2)+<-2>`.  The previously special-period `A18-C04 H10/28`
chamber has a second full-Picard carrier: after its unique stable odd closure,
doubled duality identifies all 28 oriented mirrors with divisor reflections
on `U+D8(-1)+A1(-1)`.  This second carrier is now a nine-parameter
equation-level family with a square-1408 ample class inside all 28 walls.
Likewise, the `H13/50` equation-level carrier now has forty explicit effective
reducible mirror divisors and a primitive square-12970 ample class inside all
50 walls.  The `H11/29` chamber has an arithmetically different carrier on
the characteristic-three Fermat fixed divisor lattice, and `A18-C03 H11/28`
has a second full-Picard complex-K3 carrier obtained by explicit primitive
K3-lattice glue.  The two previously different `A15-C2` and `SS2-V4`
realizations of the known `H10/34` chamber now converge on one full-Picard
K3 carrier `U+D4(-1)+A1(-1)^5`; all 34 mirrors generate Picard, and a
square-2178 ample divisor lies in their common interior.  These strengthen
existing slots rather than changing the count.  A second common-carrier
result places the 123-wall Fermat-`V4` and 84-wall Jacobian--Kummer `H10`
chambers on one nine-dimensional complex-K3 family with
`NS=U+E8(-2)+<-4>`; every one of the 207 reflections preserves the full
Picard lattice.  The `H11/28` scaled-dual carrier is now equation-level too:
an explicit eight-moduli family has `I0*+6I2+6I1`, exact very-general
`NS=U+D4(-1)+A1(-1)^6`, trivial Mordell--Weil group, and a determinant-
minus-one marking of all 28 walls.

The 37-wall rational-involution chamber on Shimada's `X_2` has now joined
the full-Picard equation-level list as an enrichment of its existing slot.
A pinned determinant-one isometry transports all 37 primitive mirrors to
the geometric basis of the same nine-parameter quadratic-base-change family
used for `A15_11 H10/21`.  On a very-general member,
`NS=U+E8(-2)+<-2>` is the full Picard lattice; the 37 mirrors span it with
index one, every reflection is integral, and a square-12 class orients every
halfspace with minimum pairing one.  Literal Vinberg remains
`3236+261+0`.  Thus one explicit family carries distinct `H10/21` and
`H10/37` chambers.  This does not change the twenty-one-slot count and is
not a full-nef, automorphism-reflection, or abstract-polytope novelty claim;
see the [`H10/37` full-Picard carrier](x2_h10_37_ue8m2a1_full_picard_carrier.md).

The known full `K3-S3-C2-H10/13` chamber now closes the same equation-family
triangle.  A newly frozen determinant-one `GL(11,Z)` isometry identifies its
archived saturated fixed Gram with the geometric
`F=U+E8(-2)+<-2>` basis.  All 13 primitive mirrors transport integrally,
span full Picard with index one, and split as
`2*(-2,div1)+2*(-2,div2)+9*(-4,div2)`.  A square-358 class is strict on all
walls.  Fresh literal, optimized, and exact-cone replays give 34 finite
vertices, four cusps, no hyperideal rays, and 208 good edges.  Thus one
explicit family carries the distinct 13-, 21-, and 37-wall `H10` chambers.
This strengthens an existing slot without changing the count and is not a
full-nef, generic-`S3`, reflection-automorphism, or abstract-novelty claim;
see the
[`H10/13` full-Picard carrier](k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.md).

That 13-wall chamber is no longer only a carrier calibration: it generates an
infinite positive Coxeter family.  Its one-based walls 10 and 13 are the
only doubling walls, and they are parallel, with raw product square `4` but
normalized squared ratio exactly `1`.  Allcock's finite-path construction
therefore gives, for every `N>=1`, a finite-volume noncompact `H10` polytope
with `N+12` facets, reflection index `N`, and volume `N` times the seed
volume.  Odd `N` has one path-reversal type; every even `N` has two
nonisomorphic phases.  Closed polynomial root formulas keep all walls in the
same `F=U+E8(-2)+<-2>` Picard lattice.  The two first 14-wall descendants
also pass independent exact Vinberg/extreme-ray audits with `41+5+0` and
`52+6+0`.  This is new to the project abstract-polytope ledger, with
worldwide priority open; the 49-record pinned H10 CoxIter nonmatch is only
bounded evidence because it misses the known seed as well.  See the
[K3-S3 redoubling-family theorem](k3_s3_h10_redoubling_family.md).

The common family now carries a structural compiler, not just three isolated
examples.  Every genuine finite-sided oriented `H10` Coxeter polyhedron with
primitive negative integral mirrors in
`F=U+E8(-2)+<-2>` is a weak full-Picard carrier on a marked very-general
member.  The discriminant exponent makes
`(-2,div1)`, `(-2,div2)`, and `(-4,div2)` the exhaustive negative mirror
alphabet.  The full 13-wall group is nonuniform, while the 21- and 37-wall
groups have exact indices `6720` and `241920`; consequently every
finite-volume chamber in this carrier class is noncompact.  This adds a
project-level common-moduli criterion, not an abstract polytope or full-nef
slot.  See the
[common-moduli carrier theorem](ue8m2a1_common_k3_moduli_carrier.md).

The generic-Hessian transposition `H10/27` chamber has now crossed the same
full-Picard threshold.  Its saturated rank-11 lattice `S_H` has determinant
256, Smith diagonal `1^4 2^6 4`, and exact discriminant form
`<3/4>+u(2)^2+v(2)`.  A full order-48 discriminant anti-isometry glues the
rank-16 Hessian Picard Gram to `U+U(2)+A2(-2)`, producing a displayed even
unimodular signature-`(3,19)` lattice in which `S_H` embeds primitively.
Hence a nine-dimensional very-general K3 family has full `NS=S_H`.  All 27
primitive mirrors generate Picard with index one and act integrally.  Fresh
Vinberg and GMP replays agree on 895 finite vertices, 67 cusps, 6,221
elliptic edges, zero failures, and zero hyperideal rays.  This strengthens an
existing Hessian/divisor slot without changing the twenty-one-slot count;
it is period-level, not equation-level, full-nef, surface-automorphism, or
abstract-polytope novelty.  See the
[`GHESS H10/27` carrier](ghess_h10_27_full_picard_k3_carrier.md).

The two high `A20` sections have now crossed the full-Picard threshold
together.  Their saturated rank-13 and rank-12 wall lattices embed
primitively in one marked K3 lattice, giving seven- and eight-dimensional
very-general algebraic K3 families with full `NS` equal to the respective
carrier.  The 19 and 20 primitive mirrors generate Picard with index one.
For each reflection group, an exhaustive discriminant-module audit finds
stable integral subgroup orders `[1,2]` but stable even orders only `[1]`:
the unique index-two overclosure is odd.  Thus the wall group forces the
unique even K3 carrier denominator.  Independent Vinberg and GMP gates give
`527+25+0` for `H12/19` and `311+29+0` for `H11/20`, with no edge or cusp
failures.  This strengthens two existing special-automorphism slots without
changing the count; the generic families do not inherit those automorphisms,
and no equation-level, full-nef, or abstract-polytope claim is made.  See the
[`A20` paired full-Picard theorem](a20_pair_full_picard_k3_carriers.md).

The known 17-wall `K3-D4 H12` triality chamber has likewise gained a
very-general full-Picard carrier without changing the slot count.  Its
saturated fixed lattice is explicitly
`U+E8(-1)+A2(-1)+A1(-1)`, and a visible unit minor embeds it primitively in
the K3 lattice.  On the resulting seven-dimensional family the 17 mirrors
span full Picard with index one and split as
`13*(-2,div1)+2*(-2,div2)+2*(-6,div3)`.  Literal and optimized volume lanes
agree on 93 finite vertices, three cusps, zero hyperideal rays, and 623 good
edges.  The generic carrier does not realize triality or the reflections as
surface automorphisms and is not a full-nef theorem; see the
[`K3-D4 H12/17` carrier](k3_d4_h12_full_picard_k3_carrier.md).

There is also a carrier theorem outside the twenty-one K3/divisor slot count.
Every admissible totally-real arithmetic Coxeter datum has a weak
semialgebraic divisor-eigenslice carrier on a rational blow-up: an exact
four-square construction embeds its trace lattice, after common square
rescaling, in `I_(1,N)`, and an integral divisor-product correspondence cuts
out the distinguished Lorentzian eigenspace.  More intrinsically, the trace
space is the image of a rational divisor-product Chow idempotent and the whole
coefficient field acts on that Tate-motive summand under exact correspondence
composition; a chosen generator becomes integral after scaling.  The compact
calibration is
fully explicit: Bugaenko's known `H8/11` chamber occurs on
`Bl_{p_0,...,p_18}(P^2_Q)`, with 11 integral divisor cuts and a 40-term golden
correspondence.  This is a new project-level algebraic realization and a
universality theorem, not a new abstract polytope or a nef/automorphism
claim; see the
[rational-surface eigen-carrier theorem](rational_surface_eigencarrier_universality.md).

Three further known `H11` chambers now have characteristic-zero
eight-dimensional full-Picard K3 families: Fermat `C2-H11/206`, Fermat
`V4-H11/328`, and characteristic-two `3B-H11/22`.  Their 556 walls all act
integrally and generate the respective Picard lattices.  The 328-wall case
uses the mandatory unique stable index-four closure, repairing 36
nonintegral source reflections.

The 172-wall Jacobian--Kummer `H12` chamber now has a seven-dimensional
full-Picard K3 family with `NS=U(2)+D5(-1)^2+<-4>`.  Its walls generate
Picard, and the transported square-eight ample class is strict on all 172.

The last three audited machine-Gram rows form a uniform nonstable `H10`
series: the known 30-, 47-, and 196-wall chambers have very-general full
Picard lattices `S_a=<2a>+<-2>^10`, `a=1,2,3`.  Explicit full-discriminant
glue embeds them in the K3 lattice, all 273 mirrors generate Picard, and both
volume lanes are green.

The known `SS2-2C H14/98` chamber has also acquired a stronger positive
enrichment without increasing this count.  After its unique
reflection-stable odd index-two closure `M`, the doubled dual `M^vee(2)` is
explicitly the actual Picard lattice `(15,9,1)` of a Hilbert square.  All 98
oriented halfspaces and reflections transport, and the primitive mirrors
generate the full Picard lattice.  Alexeev already published the target
98-wall chamber, so this is a new coordinate-level algebraic carrier bridge,
not a new abstract chamber or geometric-chamber slot.

The earlier positive enrichments remain frozen: the rank-12 `H10/12` slot
has an explicit four-parameter order-four Weierstrass carrier, the `H10/15` slot has an
independent explicit `D8+A1^2` carrier, and the `H10/38` slot has an explicit
`D4^2+A1^2` carrier, completing the rank-12 trilogy.  The known full `H17/19`
chamber has an explicit very-general two-parameter Inose carrier with exact
Néron--Severi lattice and a complete Vinberg audit.  A codimension-one
collision in the same pencil adds an `I2` fiber and gives an explicit
very-general carrier for the adjacent full `H18/24` chamber.  At rank 20,
two explicit presentations of Vinberg's `X3` surface carry the complete
36-mirror mixed-divisor chamber in `H19`.  Independently of the new
single-surface divisor-subsystem carriers, transverse ADE enhancements in the
explicit `D8+A1^2` and `D4^2+A1^2` families realize the `H10/16` and `H10/45`
actions as divisor monodromy.  This parallel transport is not asserted to be
an automorphism of one returned fiber.

There is a stronger common-carrier enrichment.  One actual member of the
`U+E8(-1)+A1(-1)^2` family, with its one order-four row-39 automorphism,
simultaneously carries all five rank-12 sections `H10/12`, `H10/15`,
`H10/16`, `H10/38`, and `H10/45`.  The complete chamber-group censuses have
rank-11 involution classes `[2]`, `[2,3]`, and `[4,12]`; their exact
Smith/discriminant fingerprints force the five archived actions.  `H10/12`
is the full-NS nef section on this surface, while the other four are
effective finite-index divisor-subsystem sections.

This count is geometric, not a claim of twenty-one new abstract polytopes.
All five new slots were already certified abstract Coxeter polytopes in the
project ledger.  Their new gain is project-level algebro-geometric
realization.  Several literal subgroup chambers remain worldwide-priority
candidates, but the present denominator does not prove such priority.

## The nine newly closed chamber realizations

| chamber | polytope | effective K3 action | exact positive gain |
|---|---:|---|---|
| `R12-E8A1x2-C2-H10-12` | `H10/12` | BH `0.4.4.1`, purely non-symplectic `C4` | fixed-section realization from a complete nef chamber, now with an explicit four-parameter Weierstrass family |
| `R12-D8A1x2-C2a-H10-15` | `H10/15` | BH `0.4.4.2`, purely non-symplectic `C4` | genuine realization of the archived 15-wall fold, now with an explicit four-parameter `I4*+2I2+10I1` family |
| `R12-D4x2A1x2-p67-H10-38` | `H10/38` | BH `0.4.4.6`, with BH `1.4.4.2` enrichment | genuine realization, now with an explicit four-parameter `2I0*+2I2+8I1` family and also admitting an effective `C4 x C2` action |
| `R14-E8A1x4-C2D-H11-19` | `H11/19` | BH `0.4.3.4`, purely non-symplectic `C4` | genuine double-transposition realization |
| `R14-D8D4-T2-H11-16` | `H11/16` | BH `0.8.1.3`, purely non-symplectic `C8` | genuine Type-2 realization |
| `R14-D8D4-T4-H10-14` | `H10/14` | BH `0.16.0.3`, purely non-symplectic `C16` | genuine Type-4 realization |
| `R14-E8D4-C3-H11-13` | `H11/13` | BH `0.6.3.1`, purely non-symplectic `C6` | genuine triality realization, with an explicit Weierstrass family |
| `R14-E8D4-C2-H12-14` | `H12/14` | BH `0.16.0.1`, purely non-symplectic `C16` | genuine transposition realization of Vinberg `P12` |
| `R16-UE8D6-C2-H14-17` | `H14/17` | BH `0.4.2.2`, extended by `1.4.2.1` | genuine order-four realization of Kondo's known Vinberg `P14` chamber, now with an explicit two-parameter Weierstrass family |

The dimension/facet distribution is therefore

`4 x H10 + 3 x H11 + 1 x H12 + 1 x H14`.

Every row was checked through the complete wall system and the literal
Vinberg endpoint and maximal-parabolic gates.  Signature alone was never
used as a volume certificate.

## Explicit algebro-geometric inverse design

For very general pairwise-distinct nonzero parameters `1,a,b,c`, the
Weierstrass surface

\[
  X_{a,b,c}:\quad y^2=x^3+t^3(t-1)(t-a)(t-b)(t-c)
\]

is a K3 surface with

\[
  \operatorname{NS}(X)=U+E_8(-1)+D_4(-1).
\]

If `omega` is a primitive cube root of unity, then

\[
  \zeta(x,y,t)=(\omega x,-y,t)
\]

has order six.  The fibers are `II* + I0* + 4 II`; the `II*` fiber supplies
`E8`, and the action cyclically permutes the three outer components of the
`I0*` fiber, giving literal `D4` triality.  The very-general Picard-rank and
saturation arguments are exact, so this is an equation-level realization of
BH `0.6.3.1` and of the complete `H11/13` fixed nef chamber.

This is the requested algebraic-object-first mechanism in a usable form:
encode the divisor roots as reducible Kodaira fibers, encode the desired
diagram symmetry by a fiber/base automorphism, and then prove generic
Neron--Severi exactness before applying the chamber and Vinberg gates.

The same explicit family now carries a second chamber by a different
mechanism.  The `A18_4_0_C06` `H13/31` wall lattice embeds with index four in
`U+E8(-1)+D4(-1)`, all 31 facet reflections extend integrally, and primitive
reduction gives `26*(-2,div1)+5*(-4,div2)`.  Thus the family is also an
equation-level mixed divisor-reflection carrier of that already counted
chamber.  The complete Vinberg replay has 3,433 finite vertices, 238 cusps,
28,766 elliptic edges, 529 affine-neighbourhood audits, and zero failures.
The five primitive square-minus-four facets and possible extra Picard roots
prevent a full-nef claim, and no claim is made that the order-six family
automorphism preserves this 31-wall chamber.  See the
[`H13/31` equation-level carrier theorem](a18_c06_ue8d4_divisor_overlattice_carrier.md).

That carrier is now nef-framed and ample-interior.  A signed 235-reflection
Weyl transport makes all 26 square-minus-two mirrors effective: five are
published irreducible curves and 21 are exact reducible divisors.  The five
square-minus-four vectors remain separate.  With `H` the primitive
square-1540 big-nef carrier class, `D506+7H` is primitive ample of square
90666 and lies strictly inside all 31 walls; seven is minimal in this
integral family.  See the
[`H13/31` nef/effectivity theorem](a18_c06_ue8d4_nef_effectivity.md).

The same inverse design now reaches `H14`.  On the explicit nonempty
parameter open recorded in the theorem, the very general surface

\[
 y^2=x^3+t^2(t^2-3)x+t^3(2+bt^2+ct^4)
\]

has fibers `I2*+II*+6I1`, hence
`NS=U+E8(-1)+D6(-1)`.  The order-four action
`(x,y,t)->(-x,i*y,-t)` exchanges the two `D6` spinor components.  Its fixed
complete-nef section is the known finite-volume noncompact `H14/17`
Vinberg chamber.  The component action has two independent proofs: a local
resolution exchanges the two tangent branches, and the identity alternative
violates the Artebani--Sarti Lefschetz congruence.

It also now reaches the rank-12 `H10/12` row directly.  Put `s=t^2` and

```text
A(s)=-3+alpha(s-1)+beta(s-1)^2,
C(s)=2-(alpha+1)(s-1)+gamma(s-1)^2+delta(s-1)^3.
```

For a very general parameter quadruple in the certified open, the minimal
resolution of

```text
y^2=x^3+A(t^2)x+tC(t^2)
```

has fibers `II*+2I2+10I1`, exact
`NS=U+E8(-1)+A1(-1)^2`, and order-four automorphism
`(x,y,t)->(-x,i*y,-t)` exchanging the two `A1` roots.  The fixed lattice,
rootless coinvariant, discriminant action, and index-two glue are exactly BH
`0.4.4.1`; hence its fixed complete-nef chamber is the finite-volume
noncompact `H10/12` polytope.

The adjacent `H10/15` row has an independent explicit inverse design.  With
the four normalized polynomials recorded in the theorem, the family

```text
y^2=x^3+P(t^2)x+tQ(t^2)
```

has fibers `I4*+2I2+10I1`, very-general
`NS=U+D8(-1)+A1(-1)^2`, and the same coordinate order-four action.  A local
equivariant branch calculation is essential here: it proves that all three
`D8` branch germs are fixed separately, so the `D8` spinor nodes are not
swapped, while the two `I2/A1` components are exchanged.  Mordell--Weil
torsion is excluded by the sharp correction bound
`2+1/2+1/2<4`.  The resulting fixed/coinvariant/glue data are exactly BH
`0.4.4.2`, giving the complete 15-facet chamber.

The rank-12 trilogy is completed by

```text
y^2=x^3+t^2P(t^2)x+t^3Q(t^2),
P(s)=-3+alpha(s-1)+beta(s-1)^2,
Q(s)=2-alpha(s-1)+gamma(s-1)^2+delta(s-1)^3.
```

Very generally the fibers are `2I0*+2I2+8I1` and
`NS=U+D4(-1)^2+A1(-1)^2`.  Invariant Tate-root coordinates at zero and
infinity prove both `D4` graphs are fixed vertexwise; the `I2/A1` pair is
exchanged.  The height correction bound `1+1+1/2+1/2<4` excludes torsion.
The exact BH `0.4.4.6` handoff yields the 38-facet chamber, with 526 finite
vertices, 141 cusps, 4,518 elliptic edges, and zero bad audits.

The common generic-Picard step is now certified rather than inferred from a
parameter count.  In every family the normalized coefficients recover all
four parameters with identity Jacobian; the residual base and short-
Weierstrass stabilizers are finite; and forgetting the marked fiber and
section has zero-dimensional fibers because their integral divisor classes
are locally constant and a K3 surface has no infinitesimal automorphisms.
The `D8` and both `D4` component actions are likewise derived in the actual
equivariant blowup charts, rather than imported from a normal-form label.

Two further rank-12 chambers are now algebraic, with a deliberately different
scope from the preceding complete-nef statements.  For `H10/16`, the same
explicit `U+E8(-1)+A1(-1)^2` order-four family contains a stable index-two
divisor sublattice

```text
M = U+D8(-1)+A1(-1)^2
```

whose full lattice-action pair is integrally conjugate to the archived
`C2b` pair.  Its 19-root subsystem has the certified 16-facet fixed chamber.
The inverse-overextension census is exhaustive: among all 67 subspaces of
the discriminant group, there is exactly one nonzero action-stable even
isotropic line, and adjoining it produces the effective BH `0.4.4.1`
Picard action.  The literal Vinberg audit gives 70 finite vertices, six
cusps, 426 elliptic rank-nine edges, 14 affine-neighbourhood audits, and no
failures.

For `H10/45`, two actual order-four K3 carriers are certified.  The explicit
`U+D8(-1)+A1(-1)^2` family contains the archived
`U+D4(-1)^2+A1(-1)^2` lattice with index two, while the explicit
`U+E8(-1)+A1(-1)^2` family contains it with index four.  The full integral
action, fixed/coinvariant glue, and all discriminant classes are transported,
not merely their ranks or determinant histograms.  Exhausting the relevant
64-element discriminant group gives exactly seven admissible index-two and
three admissible index-four extensions, landing respectively in the two
displayed Picard lattices.  The inherited complete Vinberg audit has 1,048
finite vertices, 93 cusps, 7,332 elliptic edges, 302 affine-neighbourhood
audits, and no failures.

In both cases the ambient algebraic coinvariant is the rootless lattice
`<-4>`.  A generic invariant ample class therefore avoids every extra-root
hyperplane, and K3 Riemann--Roch chooses the effective signs of all subsystem
roots simultaneously.  Extra roots of the saturated Picard lattice may
refine these chambers; neither subsystem is called the full nef chamber.

There is also an independent family-monodromy realization of each action.
For `H10/16`, the transverse slice
`(r,alpha,beta,gamma)=(1,2,1,-2/3+z)` has terminal `D8` branches satisfying
`3q^2+2q+z+1/3=0`; a meridian around `z=0` exchanges the two spinor
components at the `D8 -> D10` enhancement.  For `H10/45`, the slice
`(alpha,beta,gamma,delta)=(1,1,1,2+z)` has Tate cubic
`rho^3-3rho-z+2`, whose discriminant is `27z(4-z)`; the meridian exchanges
two outer `D4` components at the `D4 -> E6` enhancement.  Explicit unimodular
changes of basis conjugate both component braids to the archived integral
actions.  These are algebraic parallel-transport statements, not
single-fiber automorphisms.  See the
[`H10/16` divisor-overlattice theorem](r12_h10_16_order4_divisor_overlattice_carrier.md),
[`H10/45` divisor-sublattice theorem](r12_h10_45_bh_divisor_sublattice_carrier.md),
and [independent monodromy theorem](r12_h10_16_h10_45_monodromy_carriers.md).

More rigidly, the same explicit BH `0.4.4.1` surface and the same order-four
automorphism carry all five rank-12 chambers at once.  In the pinned Picard
marking their divisor-lattice indices are `1,2,2,4,4`.  A single averaged
invariant ample class selects a stable root chamber in each lattice.  The
exhaustive complete chamber-group censuses leave respectively only classes
`[2]`, `[2,3]`, and `[4,12]` at fixed rank 11, and the exact fixed-lattice
fingerprints uniquely select `H10/12`, `H10/15`, `H10/16`, `H10/38`, and
`H10/45`.  This avoids incompatible target-specific Weyl conjugators.  See
the [five-chamber simultaneous-carrier theorem](rank12_five_chamber_simultaneous_carrier.md).

The same explicit `U+D8(-1)+A1(-1)^2` family now carries a structurally
different positive in `H11`.  Let `M` be the determinant-256 rank-12 wall
lattice of the Alexeev-dual 29-facet chamber.  Its complete
reflection-stable overgroup census has a unique nontrivial even survivor,
of index four.  An explicit integral isometry identifies that overgroup with
the very-general Picard lattice of the BH `0.4.4.2` family, and a matrix
`I` with `det(I)=-4` verifies `I^T Gram(NS) I=Gram(M)`.  All 29 facet
reflections extend integrally.  After primitive reduction in `NS`, the
mirrors are exactly `16*(-2,div1)+1*(-2,div2)+12*(-4,div2)`; the six
archived square-minus-eight presentation vectors become twice primitive
square-minus-two vectors.  This is a mixed divisor-reflection subsystem,
not a full nef cone or an automorphism-fixed section.  See the
[`H11/29` divisor-overlattice theorem](alexeev_h11_29_ud8_divisor_overlattice_carrier.md).

The complex carrier is now nef-framed as well.  After 87 exact curve
reflections, all seventeen square-minus-two mirrors are effective: one is
the irreducible published curve `A16` and sixteen are explicit reducible
divisors.  The twelve square-minus-four mirrors remain separate.  The
primitive class `H+D136` is ample of square 44258 and is strict on all 29
imported walls with minimum pairing 30.  See the
[`H11/29` nef/effectivity theorem](a16_h11_29_ud8a1x2_nef_effectivity.md).

The dual `A16_6` model now gives a second integral carrier of this same
29-wall chamber.  Its unique stable odd index-eight closure `M` has

```text
M^vee(2) = U + E8(-2) + <-2>^2,
```

with an explicit determinant-one isometry to the actual saturated
involution-fixed divisor lattice on the supersingular Fermat quartic in
characteristic three.  All 29 oriented walls transport, with primitive split
`5*(-2,div1)+7*(-2,div2)+17*(-4,div2)`, and generate the target lattice with
index one.  Vinberg remains `3876+253+0`.  This determinant-1024 carrier is
arithmetically different from the determinant-16 complex-K3 carrier, but it
is not a new abstract polytope, complete nef chamber, or automorphism-fixed
section; see the
[Fermat scaled-dual theorem](a16_h11_29_fermat_char3_scaled_dual_carrier.md).

The `A15_11_1_C02` `H10/21` row closes the formerly arithmetic-only slot.
The unique even reflection-stable overclosure of its wall lattice has index
four and is explicitly `F=U+E8(-2)+<-2>`.  On a very general `F`-polarized
K3, `NS=F` and all 21 facet reflections are integral on the full Picard
lattice.  On Shimada's named `X_2`, the same lattice is the saturated
invariant divisor lattice of a rational involution; all 21 vectors remain
algebraic and fixed, but exactly eleven reflections extend to rank-20
`NS(X_2)` and the ten primitive square-minus-four reflections do not.  The
literal Vinberg test and an independent GMP double-description replay agree
on 607 finite, 50 ideal, and zero hyperideal rays.  This is a new K3/divisor
realization of the known chamber, not a full-nef, surface-automorphism, or
abstract-novelty claim.  See the
[`A15_11 H10/21` two-tier carrier](a15_11_h10_21_x2_k3_carrier.md).

The same full Picard lattice now has an explicit geometric factory.  For a
normalized nine-parameter rational elliptic surface with twelve nodal fibers,
quadratic base change ramified over one `I1` and one smooth fiber gives a K3
with `I2+22I1`.  Its invariant divisors are exactly
`U+E8(-2)+<-2>`: `E8(-2)` is the doubled pullback of the rational surface's
Mordell--Weil lattice, and `<-2>` is the new `I2` component.  The family has
the full nine moduli of `(11,9,1)` pairs, hence this is the very-general full
`NS`; all 21 reflections therefore attach to an explicit equation.  See the
[`H10/21` base-change family](a15_11_ue8m2a1_basechange_family.md).

That same family now carries the `X2-RAT-C2-H10/37` chamber as well.  The
archived determinant-one standard-to-fixed isometry transports all 37
primitive normals from `NS(X_2)^iota` to the displayed geometric basis.
Their coordinate matrix has Smith diagonal `1^11`, so they generate the
very-general full Picard lattice, and all 37 reflections remain integral.
The transported square-12 ample class from `X_2` is a strict positive-cone
interior with minimum wall pairing one.  The finite-volume statement remains
the independent literal Vinberg audit `3236` finite, `261` ideal, zero bad
edges and zero bad affine neighbourhoods.  This is an equation-level
full-Picard deformation realization of the known chamber, not a full-nef or
surface-automorphism claim; see the
[`H10/37` full-Picard carrier](x2_h10_37_ue8m2a1_full_picard_carrier.md).

It now carries the known `K3-S3-C2-H10/13` full reflection chamber too.  The
explicit determinant-one bridge is constructed by splitting a displayed
integral isotropic `U` and solving the remaining rank-nine definite
isometry.  The 13 transported roots generate `F`, every reflection is
integral, and a square-358 interior has minimum wall pairing two.  A fresh
all-subset Vinberg/CoxIter replay gives `34+4+0`, with all 208 rank-nine
edges two-ended and all seven affine-neighborhood audits good.  Kondo's
special `S3` quotient geometry is not transferred to the very-general
family, and no full-nef claim is made; see the
[`H10/13` full-Picard carrier](k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.md).

The `A18-C04 H10/28` chamber also admits a carrier unrelated to its original
special Picard-rank-19 fixed lattice.  If `M` is its unique stable odd
index-two closure, then `M^vee(2)` is explicitly
`U+D8(-1)+A1(-1)`.  Exact positive halfspace transport gives
`16*(-2,div1)+8*(-2,div2)+4*(-4,div2)`, all 28 reflections integral, and
index-one generation of the target Picard lattice.  This is a second K3
divisor realization of the known chamber, not a full-nef or new-polytope
claim; see the [`H10/28` scaled-dual carrier](a18_c04_ud8a1_scaled_dual_carrier.md).

The carrier is now attached to an explicit nine-parameter elliptic K3 family
with `I4*+I2+12I1`, trivial Mordell--Weil group, and exact very-general
`NS=U+D8(-1)+A1(-1)`.  Deterministic reduction through 226 reflections in
Roulleau's 14 complete-nef walls gives all 28 primitive mirrors in the
geometric fiber basis with maximum absolute coefficient six.  Their
square-508 sum is big and nef and remains strict for the imported chamber.
Adding the published square-208 ample class gives a square-1408 ample class
whose minimum pairings are one on the complete nef walls and seven on the
imported walls.  Hence the imported chamber meets the ample cone explicitly.
This is an equation-level ample-intersection realization, not a full-nef or
automorphism claim; see the
[`H10/28` equation-level theorem](a18_c04_ud8a1_weierstrass_carrier.md).

The neighbouring `A18-C03 H11/28` chamber has a second full-Picard K3
carrier as well.  The doubled dual of its stable odd index-four closure is

```text
U + D4(-1) + A1(-1)^6,       (r,a,delta)=(12,8,1).
```

An explicit anti-isometry of all 256 discriminant classes glues this lattice
primitively to `U+U(2)+A1(-1)^6` inside the K3 lattice.  Generic periods give
an eight-dimensional family with this exact full Néron--Severi lattice.  All
28 oriented mirrors transport, generate it with index one, and retain the
literal `1252+64+0` Vinberg audit.  This is a new full-Picard realization of
the known chamber, arithmetically different from its special-period
Picard-rank-20 carrier; it is not a complete-nef or automorphism claim.  See
the [`H11/28` scaled-dual theorem](a18_c03_ud4a1x6_scaled_dual_carrier.md).
This carrier now has an explicit eight-moduli Jacobian elliptic model.  The
identity `4A^3+27B^2=S^2 t^6 H` gives `I0*+6I2+6I1` on a certified nonempty
open; tangent transversality proves that the moduli image has dimension
eight, and a witness Groebner basis `[1]` excludes global two-torsion.  Hence
the fiber trivial lattice is saturated and equals very-general Picard.  An
exact determinant-minus-one fiber marking carries all 28 roots and their
strict square-608 interior; see the
[`H11/28` equation-level theorem](a18_c03_ud4a1x6_weierstrass_carrier.md).

The known `A15-C2 = SS2-V4 H10/34` chamber has an unusually strong common
carrier.  Exhaustive reflection-stable closure on the Alexeev side and
scaled duality on the characteristic-two side both give
`U+D4(-1)+A1(-1)^5`, with an exact oriented identification of all 34 mirrors
and index-one generation of the full Picard lattice.  In Roulleau's actual
90-wall K3 stock, five Weyl reflections put the mirror sum in the closed nef
cone; adding the published square-18 ample class gives a square-2178 ample
divisor with minimum imported pairing four.  Literal Vinberg replay has
`1774` finite vertices, `146` cusps, `12792` elliptic edges, and zero bad
audits.  The 34 walls are not asserted to be the complete 90-wall nef cone;
see the [common-carrier theorem](h10_34_ud4a1x5_k3_carrier.md).

Two much larger known `H10` chambers share a second full-Picard family.  Both
`F3-V4-H10/123` and `JK-H10-0007/84` have integral carrier lattice
`S=U+E8(-2)+<-4>`.  The explicit primitive embedding
`S -> Lambda_K3` has complement `U+E8(-2)+<4>`, so periods give a
nine-dimensional family with exact very-general `NS=S`.  All 123 and 84
oriented walls reconstruct, every reflection is integral on `S`, and each
wall span has exact index two.  Their literal Vinberg replays respectively
give `37451+1732` and `15698+862` finite/ideal vertices, `246610` and
`107233` elliptic edges, and zero failures.  These are mixed-divisor
carriers, not full-nef or automorphism realizations; see the
[shared K3-carrier theorem](f3_v4_jk_h10_shared_k3_carrier.md).

That shared family now yields a second common-moduli compiler.  Self-pairing
divisibility first excludes `(-2,div4)` from the exponent-four upper bound;
the full discriminant form then excludes `(-2,div2)` and `(-8,div4)`.
Thus `(-2,div1)`, `(-4,div2)`, and `(-4,div4)` are the exact primitive
negative mirror alphabet.
The complete 13-wall `W(S)` chamber has three cusps; the 123- and 84-wall
groups have exact indices `65,802,240` and `16,450,560`.  Thus every genuine
finite-volume Coxeter fundamental polyhedron integrally realized in `S` is
noncompact.  This is a period-level carrier-class criterion, not an explicit
equation family or a new abstract-polytope slot; see the
[`S` common-moduli theorem](ue8m2a4_common_k3_moduli_carrier.md).

The carrier census also closes three `H11` rows at once.  The full Picard
lattices are respectively `U+E8(-2)+<-2>^2`, the unique minimal even stable
index-four closure of the Fermat-`V4` fixed lattice (determinant 48), and
`U+A2(-1)^5`.  Each has rank 12, a primitive K3 embedding, and hence an
eight-dimensional very-general family.  All 206, 328, and 22 wall sets have
index-one span and integral reflections.  Exact volume remains independent:
the first and third use literal Vinberg audits with zero failures; the
328-wall case uses a complete exact cone enumeration with `225408` finite,
`11554` ideal, and zero hyperideal rays.  See the
[three-carrier theorem](three_h11_full_picard_k3_carriers.md).

The same census reaches `H12`: the 172-wall Jacobian--Kummer chamber has
rank-13 full Picard carrier `U(2)+D5(-1)^2+<-4>`.  Its discriminant length is
five, so `9>=5+2` gives a primitive K3 embedding and a seven-dimensional
very-general family.  All walls generate Picard, every reflection is
integral, and the square-eight published ample class has minimum imported
pairing one.  The literal Vinberg audit (`910316` edges,
`113274` finite vertices, `6214` cusps) and complete GMP cone replay agree
with zero failures and zero hyperideal rays.  This promotes a known invariant
divisor sublattice to a full-Picard deformation carrier; see the
[`JK-H12/172` theorem](jacobian_kummer_h12_full_picard_k3_carrier.md).

The smaller full `K3-D4 H12/17` chamber now has a second carrier tier too.
The saturated triality-fixed lattice admits a determinant-one decomposition
as `U+E8(-1)+A2(-1)+A1(-1)` and an explicit primitive K3 embedding.  Generic
periods give a seven-dimensional family with exactly this full Picard
lattice.  All 17 mirrors generate it with index one; the two norm-minus-six,
divisibility-three mirrors remain mixed divisor cuts.  The independent
literal and optimized gates give `93+3+0`, 623 edges, and zero failures.
This generic-period construction neither realizes triality as a surface
automorphism nor identifies the chamber with the full nef cone; see the
[`H12/17` carrier theorem](k3_d4_h12_full_picard_k3_carrier.md).

The current machine-ingested residue then closes uniformly.  The chambers
`D16-D5-H10/30`, `SS2-2D-H10/47`, and `F3-C2-H10/196` have Picard lattices

```text
S_a=<2a>+<-2>^10,  a=1,2,3.
```

Every easy stable-range test is inconclusive because the discriminant length
and complement rank are both 11.  Instead, explicit complements
`T_a=<2>^2+<-2a>+<-2>^8` and full discriminant anti-isometries glue each
`S_a` primitively into the K3 lattice.  The resulting nine-dimensional
families carry all 30, 47, and 196 walls; every reflection is integral and
each wall span has index one.  Vinberg and exact GMP counts agree row by row:
`1333+207`, `3136+216`, and `212824+10560` finite/ideal vertices, with zero
failures and zero hyperideal rays.  See the
[uniform three-carrier theorem](three_h10_uniform_full_picard_k3_carriers.md).

The `SS2-2A H13/50` chamber now has a second geometric life in
characteristic zero.  For its determinant-1,024 source lattice `M`, the
scaled dual `N=M^vee(2)` is explicitly
`U+D8(-1)+D4(-1)`, the full Picard lattice of genuine complex K3 surfaces.
The contragredient map transports every oriented halfspace with positive
scale, exchanges the primitive types
`40*(-4,div2)+10*(-2,div1)` with
`40*(-2,div1)+10*(-4,div2)`, and the 50 target normals generate all of `N`.
Every reflection is integral.  Vinberg's 391,614-edge audit and an
independent exact-cone replay agree on 45,954 finite, 2,598 ideal, and zero
hyperideal rays.  This is a scaled-dual complex-K3 mixed-divisor carrier,
not an overclosure, full-nef chamber, automorphism-fixed chamber, or new
abstract polytope.  See the
[`H13/50` scaled-dual carrier](ss2_2a_h13_50_ud8d4_divisor_carrier.md).

That carrier is now equation-level.  For six parameters `a,b,c,d,e,f`, put

```text
P(t)=-3t^4+6t^3+a*t^2+b*t+c,
Q(t)=2t^6-6t^5+(3-a)t^4+(1+a-b)t^3+d*t^2+e*t+f,
y^2=x^3+t^2*P(t)*x+t^3*Q(t).
```

On a certified nonempty open, the minimal resolution has
`I4*+I0*+8I1`, trivial Mordell--Weil group, and exact very-general
`NS=U+D8(-1)+D4(-1)`.  An explicit unimodular 14-by-14 marking transports
all 50 mirrors into the literal fiber-divisor basis, where they still
generate the full Picard lattice.  Clingher--Malmendier classify this
six-dimensional `P'_14` family and standard fibration; the new gain is the
equation-to-fifty-wall attachment, not priority for the family or polytope.
See the
[`H13/50` equation-level carrier](ss2_2a_h13_50_ud8d4_weierstrass_carrier.md).

This carrier now has a simultaneous exact nef/effectivity framing.  After
changing the marking from `X` to `-X`, a 60-letter Weyl word moves its
square-2012 strict interior into Roulleau's complete closed nef chamber.
All forty transported square-minus-two mirrors have nonnegative connected
decompositions in the complete 20-curve stock.  None is one of those
irreducible curves, so all forty are effective reducible divisors.  The ten
square-minus-four mirrors remain algebraic reflection vectors without an
effectivity claim.  The primitive class

```text
B2 = 2*v_nef + D154
```

has square 12970, degree at least one on every complete-nef curve, and
pairing at least four on every imported wall.  It is therefore an explicit
ample class strictly inside the 50-wall chamber.  See the
[`H13/50` nef/effectivity theorem](ss2_2a_h13_50_ud8d4_nef_effectivity.md).

Finally, a full-chamber rather than fixed-section enrichment reaches `H17`.
For very general `(A,B)` on

```text
A(A^3-(B+1)^2)(A^3-(B-1)^2) != 0,
```

the classical Inose pencil

```text
y^2=x^3-3At^4x+t^5(t^2-2Bt+1)
```

has fibers `II*+II*+4I1` and exact unimodular
`NS=U+E8(-1)^2`.  Its complete 19-wall nef chamber is the known
finite-volume noncompact `H17/19` chamber.  The literal Vinberg audit has
729 elliptic edges, 80 finite vertices, two cusps, three affine-neighborhood
audits, and no failures.  This is a new project-level equation-to-chamber
bridge, not a novelty claim for the Inose equation or abstract chamber.

The collision `A=s^2`, `B=1-s^3` in that pencil gives a second full-chamber
rung:

```text
y^2=x^3-3s^2t^4x+t^5(t^2-2(1-s^3)t+1).
```

For very general `s^3 != 0,1`, the fibers are
`II*+II*+I2+2I1` and the exact Picard lattice is
`U+E8(-1)^2+A1(-1)`.  Its complete 24-wall nef chamber lies in `H18`; the
Vinberg audit has 5,160 elliptic edges, 525 finite vertices, ten cusps, 22
affine-neighborhood audits, and no failures.  This is a standalone explicit
carrier theorem for the known rank-19 Vinberg--Kondo chamber.

The rank-20 endpoint is a mixed-divisor rather than a full-nef-wall result.
The extremal models

```text
y^2=x^3+t^5(t-1)^2,
y^2=x^3+(t^3-1)^4
```

both have exact `NS=U+2E8(-1)+A2(-1)`.  The first realizes the lattice by
`II*+II*+IV`; the second arises from an order-three quotient of
`E_omega x E_omega`, with a visible `3`-torsion section producing the
index-three primitive glue from `U+E6(-1)^3`.  All 36 Vinberg mirrors are
integral divisor classes: 24 square-minus-two rational-curve walls and 12
square-minus-six divisibility-three cuts.  The complete `H19` chamber has
87,046 elliptic edges, 8,346 finite vertices, 107 cusps, and no bad edge.

## A named rational member and cyclic `C12` enhancement

The symmetric order-six subfamily supplies a further positive refinement.
The rational K3 surface

\[
 y^2=x^3+t^3(t^2-1)(t^2-4)
\]

has geometric Neron--Severi lattice exactly
`U+E8(-1)+D4(-1)`.  This is proved for the named member, not only very
generally: direct counts over `F_(7^n)` for `n=1,2,3,4`, the order-six
eigenspace pairing, an elementary exclusion of all cyclotomic Frobenius
roots, Tate, and good-reduction specialization give Picard rank 14.

The surface carries a cyclic `C12` action.  Its order-four subgroup is the
identity on `NS`, while the order-six and order-twelve generators induce
the same literal `D4` triality.  Exhaustive replay uniquely identifies BH
rows `0.4.3.1`, `0.6.3.1`, and `0.12.1.1`.  Thus one named rational surface
simultaneously gives equation-level order-4, order-6, and order-12
realizations of the known full `H13/15` nef chamber and the finite-volume
noncompact `H11/13` triality chamber.  The construction extends to a genuine
one-parameter `C12` subfamily.

## Full-group and simultaneous-carrier enrichments

The same exact search also found noncyclic enhancements that do not increase
the number of abstract chambers but strengthen the geometry:

- BH `1.4.3.1` gives `C4 x C2` on the `H11/19` carrier;
- BH `1.8.1.1` gives `C8 x C2` on the `H11/16` carrier;
- BH `1.16.0.1` gives one `C16 x C2` K3 surface that simultaneously realizes
  `H11/16` and `H10/14` through two subactions;
- BH `1.4.4.2` enriches the rank-12 `H10/38` carrier to `C4 x C2`; and
- BH `1.4.2.1` enriches the rank-16 `H14/17` carrier to `C4 x C2`.

The rank-16 case is also a useful theorem about search design: an earlier
order-two lift obstruction remains correct, but the same Picard involution is
realized by an order-four K3 automorphism with square acting trivially on
`NS` and as `-1` on the transcendental lattice.  Thus a failed order-two glue
test is a boundary, not a reason to discard all higher-order lifts.

## Reusable positive criteria

The successful search has compiled eleven construction and promotion criteria.

1. **Effective-action intersection.**  Start with a complete divisor chamber,
   enumerate its finite wall actions, reconstruct `NS` from every effective
   K3-action row, and require exact integral agreement of the source action,
   fixed and coinvariant lattices, discriminant forms, primitive glue,
   invariant Kahler chamber, and literal Vinberg volume tests.
2. **Elliptic fiber-root inverse design.**  Realize `U` by a Jacobian elliptic
   fibration and the negative root summands by reducible Kodaira fibers.  Pick
   a CM or base automorphism with the desired component action; then prove
   Euler/minimality, exact generic `NS`, Mordell--Weil and saturation gates,
   chamber transport, and Vinberg completeness.  A Dynkin-edge/component
   match alone is necessary but not sufficient.
3. **Unimodular fiber-root shortcut.**  When the forced trivial lattice is
   already unimodular, as for `U+E8(-1)^2`, primitivity, torsion, glue, and
   same-rank saturation collapse automatically.  One then needs only a
   full-dimensional moduli argument for generic Picard exactness before the
   complete-nef and Vinberg gates.  This makes high-rank explicit carriers
   unusually tractable.
4. **Visible torsion-glue and mixed divisors.**  A quotient or elliptic
   torsion section can construct the primitive closure of a reducible-fiber
   stock explicitly, as the diagonal order-three glue in `E6^3` does for
   `X3`.  Once the divisor lattice is exact, allow reflective divisor classes
   of other norms only with their divisibility recorded and with a complete
   mixed-root Vinberg audit; do not mislabel them as irreducible curve walls.
5. **Finite-index divisor carrier.**  If an effective finite-order K3 action
   on `N=NS(X)` preserves a full-rank finite-index root sublattice `M`, first
   identify the complete integral lattice-action pair, including its
   fixed/coinvariant discriminant glue.  If the ambient coinvariant is
   rootless, choose a generic invariant ample class off all extra-root
   hyperplanes.  The `M`-chamber containing it is invariant, and
   Riemann--Roch makes its chosen square-minus-two walls effective.  Extra
   roots of `N` may refine the chamber, so the output is a divisor-subsystem
   chamber unless equality with the full nef cone is proved separately.
   Stable even isotropic subgroups of `A_M` make the possible ambient Picard
   lattices a finite exact search.  For several stable sublattices in one
   Picard lattice, average one ample class once and use the complete
   chamber-group census to identify every induced action; this is the
   simultaneous finite-index carrier criterion.
6. **Resolved-ADE monodromy carrier.**  In an algebraic family with exact
   generic `NS`, take a transverse ADE enhancement and compute the local root
   cover in resolution coordinates.  A meridian can permute exceptional
   divisor components and hence induce the desired integral wall action.
   Exact pair conjugacy and the full Vinberg endpoint/parabolic audit then
   produce a geometric Coxeter carrier.  The output is parallel transport in
   a family and must not be relabelled as a regular automorphism of one fiber.
7. **Reflection-stable Picard overclosure.**  For a Coxeter wall lattice
   `M`, exhaust the subgroups of `A_M` invariant under every facet
   reflection and retain the even overgroups of K3 signature.  An explicit
   isometry from a surviving overgroup `N` to the Néron--Severi lattice of an
   actual surface turns the same hyperplanes into algebraic divisor mirrors,
   provided every primitive ambient reflection is integral.  Primitive-
   reduce the embedded normals before recording their squares and
   divisibilities.  This produces a mixed divisor-reflection subsystem
   without requiring a surface automorphism or a full-nef interpretation.
   More generally, if `N` embeds primitively in the K3 lattice, period
   surjectivity gives a very general `N`-polarized K3 with `NS=N`; this alone
   upgrades the arithmetic chamber to a full-rank divisor carrier, with
   square-minus-two signs oriented effectively by Weyl framing and
   Riemann--Roch.  The `F3-V4-H10/123` and `JK-H10/84` pair illustrates that
   the mirrors need not generate Picard: their common wall spans have index
   two, yet every reflection acts on the full primitively embedded Picard
   lattice.
8. **Ambient-root nef-promotion gate.**  For a full-dimensional pointed fixed
   cone whose displayed facets are genuine ambient square-minus-two walls and
   whose fixed space is not trapped in a root hyperplane, full-nef promotion
   is equivalent to the absence of any other ambient root hyperplane through
   the relative interior.  On a complete extreme-ray list, every ambient
   root must have one weak sign and must not vanish identically.  This turns
   extra Picard roots and non-root primitive facets into finite exact design
   filters before an equation-level surface is constructed.
9. **Scaled-dual contragredient carrier.**  If a wall lattice `M` is awkward
   geometrically, test an integral scaling of `M^vee`.  For
   `N=M^vee(c)`, transport a source normal `r` to its primitive covector
   `p=Gram(M)r/div_M(r)` and prove the oriented halfspace identity with a
   positive scale, not merely conjugacy of reflections.  If `N` is an actual
   Picard lattice and the primitive target reflections are integral, the
   same polytope becomes a divisor-reflection chamber on a different
   algebraic carrier.  Primitive norms can exchange long and short mirrors,
   as `H13/50` changes `40*(-4)+10*(-2)` into `40*(-2)+10*(-4)`.
   Before dualizing, enumerate reflection-stable closures as well: for
   `H14/98`, the source lattice itself has the wrong discriminant, while its
   unique stable odd index-two closure `M` has
   `M^vee(2)=Pic(Hilb^2(X))`.  Exact positive transport then gives the target
   split `14*(-2,div1)+70*(-2,div2)+14*(-4,div2)` and index-one generation by
   the 98 mirrors.
10. **Ample-intersection promotion.**  Do not confuse failure of full-nef
    equality with failure of algebraic relevance.  For an oriented imported
    cone `C`, its interior meets the ample cone exactly when it meets the nef
    cone in a positive class.  Indeed, if a positive nef `H` is strict on all
    imported walls, then `mH+A` is ample and remains strict for all sufficiently
    large `m`, for any ample `A`.  When the complete nef chamber has finitely
    many known walls, finding `H` is a finite exact system of weak nef-wall
    inequalities and strict imported-wall inequalities.  The `H10/28`
    equation-level carrier passes this gate with square-508 `H`, and `H+A`
    already gives the square-1408 ample witness at `m=1`.
11. **Effective-mirror Weyl descent.**  Once a complete finite stock of
    irreducible `(-2)`-curves and one ample height are known, correct the
    positive-cone sign, fire any curve having negative pairing, and record
    the exact reflection word.  Ample height strictly decreases, so the
    integral descent terminates.  Apply the same word to every imported
    mirror.  For each square-minus-two orientation, a second negative-wall
    descent gives an explicit nonnegative curve decomposition; comparison
    with the complete irreducible stock decides irreducibility rather than
    guessing it.  Keep other mirror norms separate.  For `H13/50` this yields
    40 effective connected decompositions, proves all 40 reducible, and then
    the ample-intersection criterion gives the primitive square-12970 class
    `2*v_nef+D154` inside all 50 walls.  The same mechanism is exceptionally
    short for `H10/34`: five reflections make the imported sum big and nef,
    after which one published ample class produces the square-2178 strict
    witness.

Negative computations are now retained only as finite prefilters for these
positive constructions or as independently stated frontier theorems.

The direct primitive generic-action denominator is exhausted: all 587 pure
and 3,561 full BH v2.1 rows reduce to 32 exact ambient-lattice hits, and every
Picard-image element and possible `V4/S3` subgroup in those hits misses the
target discriminant forms under equality of the Picard lattice.  This is a
useful boundary, not a contradiction to the positive completions: those place
the target root lattices nonprimitively inside larger effective Picard rows.

The first proposed special-period replacement is now also decided and used
constructively.  For `H10/16`, the complement `U+U(2)+D6(-1)` has maximal
totally singular subspace dimension six modulo two, whereas a primitive
`E7(-2)` would inject a seven-dimensional totally singular subspace.  Hence
the primitive symplectic Nikulin route is impossible, and an order-two
nonsymplectic lift with the same fixed lattice is excluded by its `Z/4`
discriminant factor.  The positive completion deliberately changes both
features: it uses an order-four automorphism and places the archived lattice
nonprimitively, with index two, in `U+E8(-1)+A1(-1)^2`.  The old no-go is
therefore retained as the boundary that predicts the unique stable even
overextension; see the
[sharp obstruction](r12_h10_16_nikulin_obstruction.md) and
[positive completion](r12_h10_16_order4_divisor_overlattice_carrier.md).

For `H10/45`, primitive `E7(-2)` embeddings in `U(2)^2+D6(-1)` do exist,
but the primitive symplectic lane nevertheless fails universally.
Primitivity and integrality force glue index eight; the mod-two glue then
makes every fixed discriminant class have even quadratic value, with
histogram `(16,0,16,0)`, while the required Picard transposition has
`(8,8,8,8)`.  The positive completion again leaves that lane: order-four
actions on the two explicit Picard overgroups contain the target divisor
lattice with indices two and four.  The exhaustive extension census is now a
positive carrier criterion rather than an unresolved search instruction;
see the [`H10/45` boundary theorem](r12_h10_45_e7_embedding_search_boundary.md)
and [positive completion](r12_h10_45_bh_divisor_sublattice_carrier.md).

## TOOT: the toric \(3A_2\) endpoint and its 216-fold cubic gluing

There is now a complete positive theorem for smooth projective toric
surfaces, not merely an incidence carrier.  The complete nef cone is an
integral Coxeter chamber exactly when every invariant boundary curve has
square at least \(-2\), equivalently when \(-K\) is nef.  Reflexive-polygon
duality and an exact fan replay then give twelve hyperbolic carrier types,
all finite-volume noncompact, and the sharp ceiling \(H^6\).

The unique endpoint is a nine-wall chamber with cyclic labels
\((4,4,3)^3\) and exact ray split \(18+3+0\).  Its abstract polytope is
Im Hof's known dimension-six Napier-cycle example, so this is not a new
worldwide Coxeter polytope.  The positive gain is a rigid full-nef
realization: the surface is the minimal crepant resolution of
\(XYZ=W^3\), with three \(A_2\) exceptional chains.  Their reflections form
\(W(A_2)^3\); the three remaining wall orbits are exactly all 27 cubic-line
classes.  Hence
\[
  W_{\mathrm{tor}}=\Gamma_2^6\rtimes W(A_2)^3,
  \qquad [W_{\mathrm{tor}}:\Gamma_2^6]=216,
\]
and the classical smooth-cubic \(P^6\) chamber is tiled by exactly 216 toric
endpoint chambers.  A separate exact audit proves chamber containment by
expressing every cubic-line root as a \(0/1\) sum of toric boundary classes;
the smoothing interpretation is made on the product local
simultaneous-resolution marking cover, not asserted as the full monodromy of
one pencil.

The inverse direction is a terminating twelve-row marked
polytope-to-fan/divisor recognizer, rated L3 rather than strict L4.  Its
accepted denominator is the already toric-completable locus, so the
independent broad P-only pre/post AG survivor delta is zero.  See the
[smooth-toric theorem](toric_full_nef_coxeter_factory.md), [direct hostile
audit](toric_full_nef_coxeter_factory_hostile_audit.md), and [congruence and
smoothing audit](toric_h6_congruence_smoothing_audit.md).

## TOOT: a Prasad-free maximal-reflection/Borcherds counterexample

The complete `U+E8(-1)+E7(-1)` K3 nef chamber is exactly Tumarkin's
finite-volume noncompact `Gamma16`, with 19 facets, `125+3` vertices and
1,087 edges.  A direct Felikson--Tumarkin facet-rank argument proves that it
is maximal among **all discrete reflection groups**: an overgroup chamber
would have 17, 18 or 19 facets, while the complete simplex and `n+2`
classifications remove the first two and Tumarkin uniqueness makes the last
the same chamber.  This uses no Prasad/ABSW volume machinery.

At the same time determinant/discriminant-form arithmetic proves that the
lattice is maximal integral and `W_all=W_2`, while Wang's Theorem 7.2 proves
that it has no automorphic correction.  Thus maximal arithmetic reflection
does not imply Borcherds correction.  It is not a maximal arithmetic group:
the exact diagram involution has Lorentz eigendimensions `(9,8)`, not
`(16,1)`, and gives the nonreflection extension
`Theta16=Gamma16 semidirect C2`.

The other seventeen Wang-negative rows all have explicit same-lattice
`(-4,div 2)` or `(-6,div 3)` reflection enlargements, so the experiment is
closed `17+1+0`.  Dolgachev's decorated genus-two model also makes the AG
failure sharp: the same maximal chamber occurs over a three-dimensional K3
family, so the Gram image is a point but the carrier locus is not.  This is
a new project theorem/cross-identification, not a new abstract polytope.  See
the [counterexample theorem](wang_h16_maximal_reflection_no_correction.md)
and [complete ledger](wang18_same_lattice_mixed_root_ledger.md).

## Frozen proof packages

- [smooth-toric full-nef theorem](toric_full_nef_coxeter_factory.md)
- [smooth-toric hostile audit](toric_full_nef_coxeter_factory_hostile_audit.md)
- [toric \(3A_2\) congruence/smoothing audit](toric_h6_congruence_smoothing_audit.md)
- [rank-12 three-lift theorem](roulleau_rank12_bh_k3_lifts.md)
- [rank-14 `E8+D4` theorem](roulleau_ue8d4_generic_ns_k3_automorphisms.md)
- [explicit order-six Weierstrass family](ue8d4_order6_weierstrass_family.md)
- [named rational `C12` member](ue8d4_named_c12_member.md)
- [rank-14 `D8+D4` two-lift theorem](roulleau_ud8d4_bh_k3_lifts.md)
- [rank-14 double-transposition theorem](roulleau_r14_c2d_special_k3_order4_lift.md)
- [rank-14 full-group enrichments](roulleau_bh_full_group_enrichments.md)
- [rank-16 `H14/17` theorem](kondo_rho16_h14_bh_order4_lift.md)
- [explicit rank-16 order-four family](ue8d6_order4_weierstrass_family.md)
- [explicit rank-12 order-four family](ue8a1x2_order4_weierstrass_family.md)
- [explicit `D8+A1^2` rank-12 order-four family](ud8a1x2_order4_weierstrass_family.md)
- [explicit `D4^2+A1^2` rank-12 order-four family](ud4x2a1x2_order4_weierstrass_family.md)
- [`H10/16` Nikulin obstruction and higher-order completion criterion](r12_h10_16_nikulin_obstruction.md)
- [`H10/45` primitive-glue obstruction and higher-order completion criterion](r12_h10_45_e7_embedding_search_boundary.md)
- [`H10/16` order-four divisor-overlattice carrier](r12_h10_16_order4_divisor_overlattice_carrier.md)
  ([certificate](../data/enumerations/r12_h10_16_order4_divisor_overlattice_carrier.json))
- [`H10/45` order-four divisor-sublattice carriers](r12_h10_45_bh_divisor_sublattice_carrier.md)
  ([certificate](../data/enumerations/r12_h10_45_bh_divisor_sublattice_carrier.json))
- [`H10/16` and `H10/45` ADE-monodromy carriers](r12_h10_16_h10_45_monodromy_carriers.md)
  ([certificate](../data/enumerations/r12_h10_16_h10_45_monodromy_carriers.json))
- [one row-39 K3 action carrying all five rank-12 `H10` chambers](rank12_five_chamber_simultaneous_carrier.md)
  ([certificate](../data/enumerations/rank12_five_chamber_simultaneous_carrier.json))
- [`H11/29` mixed divisor-overlattice K3 carrier](alexeev_h11_29_ud8_divisor_overlattice_carrier.md)
  ([certificate](../data/enumerations/alexeev_h11_29_ud8_divisor_overlattice_carrier.json))
- [`H11/29` exact nef/effectivity and ample-intersection framing](a16_h11_29_ud8a1x2_nef_effectivity.md)
  ([certificate](../data/enumerations/a16_h11_29_ud8a1x2_nef_effectivity.json))
- [`H11/29` characteristic-three Fermat scaled-dual carrier](a16_h11_29_fermat_char3_scaled_dual_carrier.md)
  ([certificate](../data/enumerations/a16_h11_29_fermat_char3_scaled_dual_carrier.json))
- [`A18-C03 H11/28` scaled-dual full-Picard K3 carrier](a18_c03_ud4a1x6_scaled_dual_carrier.md)
  ([certificate](../data/enumerations/a18_c03_ud4a1x6_scaled_dual_carrier.json))
- [`A18-C03 H11/28` explicit eight-moduli Weierstrass carrier](a18_c03_ud4a1x6_weierstrass_carrier.md)
  ([certificate](../data/enumerations/a18_c03_ud4a1x6_weierstrass_carrier.json))
- [`A15-C2 = SS2-V4 H10/34` common full-Picard K3 carrier](h10_34_ud4a1x5_k3_carrier.md)
  ([certificate](../data/enumerations/h10_34_ud4a1x5_k3_carrier.json))
- [`F3-V4-H10/123` and `JK-H10/84` shared full-Picard K3 carrier](f3_v4_jk_h10_shared_k3_carrier.md)
  ([certificate](../data/enumerations/f3_v4_jk_h10_shared_k3_carrier.json))
- [`U+E8(-2)+<-4>` common period-moduli carrier compiler](ue8m2a4_common_k3_moduli_carrier.md)
  ([certificate](../data/enumerations/ue8m2a4_common_k3_moduli_carrier.json))
- [explicit Type-I Enriques/K3 good lift package](enriques_type_i_explicit_lift_package.md)
  ([certificate](../data/enumerations/enriques_type_i_explicit_lift_package.json))
- [Type-I downstairs anti-Gram compiler](enriques_type_i_downstairs_anti_gram_census.md)
  ([certificate](../data/enumerations/enriques_type_i_downstairs_anti_gram_census.json))
- [independent Type-I explicit-package and anti-Gram hostile audits](enriques_type_i_explicit_lift_package_hostile_audit.md)
  ([census audit](enriques_type_i_downstairs_anti_gram_census_hostile_audit.md))
- [three `H11` full-Picard K3 carriers with 556 walls](three_h11_full_picard_k3_carriers.md)
  ([certificate](../data/enumerations/three_h11_full_picard_k3_carriers.json))
- [`JK-H12/172` full-Picard K3 deformation carrier](jacobian_kummer_h12_full_picard_k3_carrier.md)
  ([certificate](../data/enumerations/jacobian_kummer_h12_full_picard_k3_carrier.json))
- [three `JK-H12` single doubles and the `124N+48` path family](jacobian_kummer_h12_single_doublings.md)
  ([certificate](../data/enumerations/jacobian_kummer_h12_single_doublings.json))
- [independent all-parameter hostile audit of the `JK-H12` path family](jacobian_kummer_h12_parallel_redoubling_hostile_audit.md)
  ([certificate](../data/enumerations/jacobian_kummer_h12_parallel_redoubling_hostile_audit.json))
- [four-wall `JK-H12` doubling cube through 1,827 facets](jacobian_kummer_h12_doubling_cube.md)
  ([certificate](../data/enumerations/jacobian_kummer_h12_doubling_cube.json))
- [independent raw-data hostile audit of the `JK-H12` cube](jacobian_kummer_h12_doubling_cube_hostile_audit.md)
  ([certificate](../data/enumerations/jacobian_kummer_h12_doubling_cube_hostile_audit.json))
- [`JK-H12` `D_infinity x C2^3` path-times-cube family](jacobian_kummer_h12_path_cube_family.md)
  ([certificate](../data/enumerations/jacobian_kummer_h12_path_cube_family.json))
- [independent all-parameter hostile audit of the `JK-H12` path-times-cube family](jacobian_kummer_h12_path_cube_hostile_audit.md)
  ([certificate](../data/enumerations/jacobian_kummer_h12_path_cube_hostile_audit.json))
- [`JK-H12` phase-refined `D_infinity^3` family and scoped maximum-product theorem](jacobian_kummer_h12_triple_parallel_family.md)
  ([certificate](../data/enumerations/jacobian_kummer_h12_triple_parallel_family.json))
- [`JK-H12` complete one-edge and three-vertex finite-subtree census](jacobian_kummer_h12_small_subtree_census.md)
  ([certificate](../data/enumerations/jacobian_kummer_h12_small_subtree_census.json))
- [`X2-RAT-C2-H10/37` equation-level full-Picard K3 carrier](x2_h10_37_ue8m2a1_full_picard_carrier.md)
  ([certificate](../data/enumerations/x2_h10_37_ue8m2a1_full_picard_carrier.json))
- [`K3-S3-C2-H10/13` equation-level full-Picard K3 carrier](k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.md)
  ([certificate](../data/enumerations/k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.json))
- [five finite-volume `H9` chambers from the complete parked rank-ten stock](alexeev_parked_rank10_h9_examples.md)
  ([certificate](../data/enumerations/alexeev_parked_rank10_h9_examples.json))
- [independent geometry, arithmetic-overclosure, and lineage audit of the five parked chambers](alexeev_parked_rank10_h9_examples_hostile_lineage.md)
  ([certificate](../data/enumerations/alexeev_parked_rank10_h9_examples_hostile_lineage.json))
- [primary-source arithmetic-overclosure and K3-moduli lineage of the five parked chambers](alexeev_parked_rank10_h9_lineage_overclosures.md)
  ([certificate](../data/enumerations/alexeev_parked_rank10_h9_lineage_overclosures.json))
- [genuine Nikulin action carriers for four parked `H9` chambers](alexeev_parked_h9_k3_action_carriers.md)
  ([certificate](../data/enumerations/alexeev_parked_h9_k3_action_carriers.json))
- [three finite cube types and three infinite `H9` product-redoubling families](alexeev_parked_h9_doubling_products.md)
  ([certificate](../data/enumerations/alexeev_parked_h9_doubling_products.json))
- [independent all-parameter hostile audit of the `H9` product families](alexeev_parked_h9_doubling_products_hostile_audit.md)
  ([certificate](../data/enumerations/alexeev_parked_h9_doubling_products_hostile_audit.json))
- [infinite `K3-S3-H10-P_N` finite-path redoubling family](k3_s3_h10_redoubling_family.md)
  ([certificate](../data/enumerations/k3_s3_h10_redoubling_family.json))
- [`U+E8(-2)+<-2>` common-moduli carrier compiler](ue8m2a1_common_k3_moduli_carrier.md)
  ([certificate](../data/enumerations/ue8m2a1_common_k3_moduli_carrier.json))
- [`GHESS H10/27` very-general full-Picard K3 carrier](ghess_h10_27_full_picard_k3_carrier.md)
  ([certificate](../data/enumerations/ghess_h10_27_full_picard_k3_carrier.json))
- [`K3-D4-H12/17` very-general full-Picard K3 carrier](k3_d4_h12_full_picard_k3_carrier.md)
  ([certificate](../data/enumerations/k3_d4_h12_full_picard_k3_carrier.json))
- [`A20 H12/19` and `H11/20` paired full-Picard K3 carriers](a20_pair_full_picard_k3_carriers.md)
  ([certificate](../data/enumerations/a20_pair_full_picard_k3_carriers.json))
- [remaining full-Picard carrier ingestion queue](full_picard_carrier_ingestion_queue_2026_08_08.md)
- [uniform 30-/47-/196-wall `H10` full-Picard K3 carriers](three_h10_uniform_full_picard_k3_carriers.md)
  ([certificate](../data/enumerations/three_h10_uniform_full_picard_k3_carriers.json))
- [`H13/31` equation-level `U+E8+D4` mixed divisor carrier](a18_c06_ue8d4_divisor_overlattice_carrier.md)
  ([certificate](../data/enumerations/a18_c06_ue8d4_divisor_overlattice_carrier.json))
- [`H13/31` exact nef/effectivity and ample-intersection framing](a18_c06_ue8d4_nef_effectivity.md)
  ([certificate](../data/enumerations/a18_c06_ue8d4_nef_effectivity.json))
- [`A15_11 H10/21` generic/named K3 divisor carrier](a15_11_h10_21_x2_k3_carrier.md)
  ([certificate](../data/enumerations/a15_11_h10_21_x2_k3_carrier.json))
- [`A15_11 H10/21` explicit quadratic-base-change family](a15_11_ue8m2a1_basechange_family.md)
  ([certificate](../data/enumerations/a15_11_ue8m2a1_basechange_family.json))
- [`A18-C04 H10/28` scaled-dual `U+D8+A1` K3 carrier](a18_c04_ud8a1_scaled_dual_carrier.md)
  ([certificate](../data/enumerations/a18_c04_ud8a1_scaled_dual_carrier.json))
- [explicit nine-parameter `U+D8+A1` Weierstrass family](ud8a1_weierstrass_family.md)
  ([certificate](../data/enumerations/ud8a1_weierstrass_family.json))
- [`A18-C04 H10/28` equation-level ample-intersection carrier](a18_c04_ud8a1_weierstrass_carrier.md)
  ([certificate](../data/enumerations/a18_c04_ud8a1_weierstrass_carrier.json))
- [`SS2-2A H13/50` scaled-dual complex-K3 divisor carrier](ss2_2a_h13_50_ud8d4_divisor_carrier.md)
  ([certificate](../data/enumerations/ss2_2a_h13_50_ud8d4_divisor_carrier.json))
- [`SS2-2A H13/50` explicit six-parameter Weierstrass carrier](ss2_2a_h13_50_ud8d4_weierstrass_carrier.md)
  ([certificate](../data/enumerations/ss2_2a_h13_50_ud8d4_weierstrass_carrier.json))
- [`SS2-2A H13/50` exact nef/effectivity and ample-intersection framing](ss2_2a_h13_50_ud8d4_nef_effectivity.md)
  ([certificate](../data/enumerations/ss2_2a_h13_50_ud8d4_nef_effectivity.json))
- [`SS2-2C H14/98` scaled-dual Hilbert-square carrier enrichment](ss2_2c_h14_98_hilb2_scaled_dual_carrier.md)
  ([certificate](../data/enumerations/ss2_2c_h14_98_hilb2_scaled_dual_carrier.json))
- [exact K3 full-nef promotion gate for `H10/16`, `H10/45`, and `H11/29`](k3_nef_promotion_gate_h10_h11.md)
  ([certificate](../data/enumerations/k3_nef_promotion_gate_h10_h11.json))
- [explicit rank-18 Inose `H17/19` family](ue8x2_inose_h17_family.md)
- [explicit rank-19 Inose-collision `H18/24` family](ue8x2a1_inose_h18_family.md)
- [two explicit `X3` mixed-divisor `H19/36` carriers](vinberg_x3_h19_explicit_carriers.md)
- [master lineage ledger](higher_dimensional_polytope_lineage_ledger.md)

The earlier focused packages and the two sharp rank-12 carrier boundaries
have pinned test and external-source replays.  The new positive carrier
packages have separate regression suites and exact archive replays.  The
corrected explicit-family payload
has SHA-256
`ef90365b2a3772e97e75da3d470d4d40f96ec61ff68bcdd0a401c2bf174cdfd8`.
The named-member and rank-16 explicit-family payloads are respectively
`39336a4a893b4d36f3db875f7f6af2f8d20b449568b83bec684e213734675d47`
and `445992677f5e8a3a143ee9c0693db4d23cb83ebecbf52a1fdacf47bc1ae13207`.
The rank-12 explicit-family and rank-18 Inose-family payloads are
`65a9a423a9d593294cde87e282b7b386433c47bc1d170d399c38c3ab26645904`
and `b02465173f6f73acbc743159aaabd3bc0a6f8784fb4ecbfdd8b95e731494bdaa`.
The rank-19 Inose-collision payload is
`6c4507a6eef69adb8ada547be7733986cf8df4ee5705662deec3d7413de0fc0f`.
The `D8+A1^2` rank-12 family and explicit `X3`-carrier payloads are
`284c419db6d040698fe4509535aa309694bc77d76481c668f8ca91b77655d710`
and `307e5d3b1ba8dc3118efb4f7a8a981b44ed5e118155fe6cb1d854ec3d8c28025`.
The `D4^2+A1^2` rank-12 family payload is
`f5c9378f61d62c42e6a23702ebc200eaadbf6556a2ddfe18229e19bc3be2df78`.
The `H10/16` divisor-overlattice, `H10/45` divisor-sublattice, and joint
monodromy payloads are respectively
`53a4d077584ea71528e496eace2e24ac0f699f520118dd2c14f6900e455db829`,
`412808ca6aafde96a7a13ddba04b7bb2ae1ea787f75b52617557c544063a9acb`,
and
`ce457d3a0b1718e73391eb68ec7f7ca0aadb8f2e1f32c4b3eb12a470c3aef9c0`.
The simultaneous five-chamber archive SHA-256 is
`325c306d262a28dbbb9bf52c8c5d6ddad79fc0eedf247c2e04aa661e2ec8cc80`.
The `H11/29` divisor-overlattice payload is
`24a73754ddfa9f1254f4a6f7a83575dedc19f50e54e1f96d67576b5c52d8634a`.
Its nef/effectivity and ample-intersection payload is
`401f828906e17445002e467707ca8e762d6b3c1ff1eb07a44d334b7c39c3f0d5`.
The `H11/29` Fermat scaled-dual carrier payload is
`eb5aff25a3e5573355ef6fe11f0ab93a2e62ec42035f6bd9d1783f625ec5cd0e`.
The `A18-C03 H11/28` scaled-dual K3 carrier payload is
`5225b9eeb5d48c585ae793abea8fe71f06ef6bd3a16d235a6d3f6859707900b2`.
Its equation-level Weierstrass carrier payload is
`ff33f2177c55c50eea24e7efb99dbaf4536e269ef64c3a20c3ee7a428b25c222`.
The common `A15-C2 = SS2-V4 H10/34` K3 carrier payload is
`55b990d3002e787f7199ca6577e2e39a05ccf713ca05365e08fafd8212395e88`.
The shared `F3-V4-H10/123` and `JK-H10/84` K3 carrier payload is
`39574ca56421155388ec734076cf728a49f81c94bdfcdbc7d472b104e5014f08`.
The `U+E8(-2)+<-4>` common period-moduli compiler payload is
`87dc623443c6386bb562f1fa969603768cd04f2f4bcf279f428cf07adfcf552c`.
The joint three-`H11` full-Picard carrier payload is
`58cf2623ba476ec5ee0d19a173864fa5955a9e7a49f8e3948b0a4d87831252a2`.
The `JK-H12/172` full-Picard carrier payload is
`07a31cf36fd1e3ca7d548dd89b976d1d6919455243087ec3250c12ea463a1420`.
The `JK-H12` single-doubling/path-family payload is
`01eb89baef2cedcf78d3a20dd1c55931db9a9a2ca83bc5f4b7c1a1329cbdca98`.
Its independent all-parameter hostile archive is
`ac7ae51a5df8615ad54d871152f71e0c0b10d53ae417872a31cd8cd5e514636d`.
The four-wall `JK-H12` cube archive is
`0267e3d56b5b0acaa7015e4cadac1d357d387532831e89485c1846eae28be3f1`.
Its independent raw-data hostile archive is
`b0e94b5486e47fccc9586e10c3c67515dc5726e770295ab286d1394f2df62a33`.
The `JK-H12` path-times-cube production archive is
`5bcd971e1c1e99ca7486f07387c6ab410f0fba99a1ad4f43c3bde511082a6024`.
Its independent all-parameter hostile archive is
`5c81e4c311e96a64a86ec377d10aa40b54dcd94f685432fee11d9bc9f8ea5663`.
The phase-refined `JK-H12` `D_infinity^3` production archive is
`bf6fa667f224df5117b97a4590a1f673b0504d3144dfd573219e177cec487c5a`.
The complete `JK-H12` one-edge/three-vertex finite-subtree census archive is
`796575ebbaaf8ec0ba3054911397d83834d8939e1eb4c20cd24aae1b56aee56d`.
The `X2-RAT-C2-H10/37` equation-level full-Picard carrier payload is
`50c48a6ef7dcaaa4719b21d91f303a6d1983d5ab15cd998e4832d498a7ae4e08`.
The `K3-D4-H12/17` very-general full-Picard carrier payload is
`d41e18b550d600d42ec75a7e068073e8bbd51366af975c25e6efdc798ef62550`.
The `K3-S3-C2-H10/13` equation-level full-Picard carrier payload is
`c033a684aa04b2c6d0bfe80d1b5393f0c52f5a9c9c1039140ca674ea9cc67a34`.
The five parked rank-ten `H9` chamber certificate has archive SHA-256
`e181b8c021fe0f37d23b53b33e4e3900559502d0f8c982a4a9e9cbcf0866e430`.
Its independent hostile-lineage archive has SHA-256
`79a0c1abb05d464b7a041cea68c17a9c3461122db2d94e47739bb0993db676b4`.
The primary-source overclosure/K3-lineage archive has SHA-256
`972f1eb50940046c90c3770ff356a5a923477dec2e87d1d7acbffd85bf785bcd`.
The four-row parked Nikulin-action carrier archive has SHA-256
`ba7ba78766d5e00bca012b6c5591ebb02572a07a9a2ef160888ad1bc57a41159`.
The parked-seed `H9` product-family archive has SHA-256
`10e6323895bacd2d18f897f33ea5f9d4bb1ee4eedad061d545aa5302b447c694`.
Its independent hostile all-parameter archive has SHA-256
`29366fe58e52aa2492693704134b78967e4340ca972e5f56d7d7f69972d646bd`.
The infinite `K3-S3-H10-P_N` redoubling-family payload is
`11da6289847d134ffe56081604fbf6df8779f48f3984fd2a4a0e2aee4b8309f9`.
The `U+E8(-2)+<-2>` common-moduli compiler payload is
`608f81a0fdd709134206e6cfe9e8a7de1338231b03f9ffc6e0e36c8c8e97c09a`.
The `GHESS H10/27` very-general full-Picard carrier payload is
`f8e2f2e128be4d212370b8b99dc298e0725da8d0b8e560decbb03db32019d036`.
The paired `A20 H12/19` and `H11/20` full-Picard carrier payload is
`2d6c68b8aed6fbe7c39b527d94a23a1b7e466504d3a6097b13536f9edc96cbfe`.
The uniform three-`H10` full-Picard carrier payload is
`aaac645cdc0d918dc7d5ea50a7944690bb4a1705b089820d9c37a16432c63d32`.
The `H13/31` equation-level divisor-overlattice payload is
`d34d7a9b36dd9d8476db53e253530635d2e51d0d4199e12e147c46b385eebf05`.
Its nef/effectivity and ample-intersection payload is
`f8b9abf7fe2589a5044159d3347ca17524396543c6e4e1b73be3d67b61406a5d`.
The exact K3 nef-promotion-gate payload is
`c14290b4bd93e8976482a5404ff3188f43ff1be3a61da7f786ac4eb8e622dd8f`.
The `A15_11 H10/21` K3 carrier payload is
`a2905e05f15d5538d7f5896314bac4f52228f79a0ba2609e2ea13e99ff2e963c`.
Its equation-level quadratic-base-change payload is
`83d4d64c5a4003c296efc99c55557bd7ec4170fe26d578b9e950c02edff45f6d`.
The `A18-C04 H10/28` scaled-dual K3 carrier payload is
`e173d05f70a9dbe47c2a3d4499a11cae38e5584cee8e49a1e4a0ffcd6174a61a`.
The standalone `U+D8+A1` equation-family payload is
`f10e1204fb4d15b3664c9a525c53363cab954cc95ce746c8b412ece6c83a30c0`.
The equation-to-`H10/28` ample-intersection payload is
`c97f519070d1915e7e6afce1544d266ada350b834a0abd9d289ba5847dca424e`.
The `SS2-2A H13/50` scaled-dual K3 carrier payload is
`6f27ddf2ee4cedb8b5d3aaf2d59bdf2369056179635812c29091bd7129a221a7`.
Its equation-level Weierstrass carrier payload is
`8b6a2d15e9a289b28c0be5b84f06786b808bdc45979207a29053cdf008b9ad54`.
Its nef/effectivity and ample-intersection payload is
`335f5e79e9efdab56ca4c3727e946bc0d3289d77d41a3862a57774151374df0c`.
The `SS2-2C H14/98` Hilbert-square carrier payload is
`4ab9ae96351b015ce694e2e072605d4db8fb7eb33497678aba4db5613313e8e8`.
