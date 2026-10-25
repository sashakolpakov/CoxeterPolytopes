# Native5 and native7+: an explicit compact Coxeter codimension-two face

23 October2026. **Complete MANUAL candidate, pending ENTIRE different-author
review and root acceptance.** No scientific worker, import, enumeration,
numerical/symbolic process, optimizer, test, compiler, Lean, TeX, build or
Git operation supplied the proof. Earlier mathematical sources, papers,
rules, guides and global57/18/16 bounds remain unchanged.

This is a direct special-good-face argument. It may reproduce a known
Borcherds-type face case; no historical novelty is claimed. The Coxeter
property is PROVED below for these particular faces, not assumed for an
arbitrary face. There is no arithmetic hypothesis or arithmeticity-of-
faces inference. Root requested the residual far-row8 finite5+ branch and
suggested using the published compact29 input to keep this restriction
independent of the stronger internal compact18 proof. I derived the
explicit facet classes/projection calculation. I contributed to older
local geometry/financing sources; that does not provide this new source's
outside review.

## 1. Exact theorem, inputs and nonclaims

Let P be a finite-sided finite-volume real hyperbolic Coxeter polytope
in H^n, n>=4. Suppose two of its ACTUAL facets meet at a finite vertex
with finite native label

    m=5 or m>=7.

Then their actual codimension-two face is a COMPACT Coxeter polytope
of dimension n-2. Consequently the PUBLISHED compact ceiling29 gives
n<=31. In particular a hypothetical finite-volume Coxeter57-polytope
has no such bond at a finite vertex. Its ideal links also have no such
bond, directly by affine classification.

More generally, conditional on ANY valid compact ceiling N, the same
geometric proof gives n<=N+2 in the presence of this bond. The accepted
internal compact18 theorem thus gives the stronger CONDITIONAL n<=20
for this bond-bearing class. The finite57 application uses ONLY the
published29 input:55>29. It does not need that stronger internal result.
Neither statement bounds all finite-volume dimensions by31 or20: the
other native labels2/3/4/6 remain. No dimension56/global364 cap, root
realization, new example, ledger debit or mathematical execution follows.

All literal bindings match at this source freeze.

|ID|Source|SHA256|Restricted interface|
|---|---|---|---|
|C|manuscripts/finite_volume57.tex|c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0|Actual all-wall acuteness/native products, finite simplicity, finite/affine Coxeter classification and finite-volume compactification/cusp links; no numerical57 conclusion used|
|V|manuscripts/compact_farthest_vertex_bounds.tex|530d964e1957981af6c3a8f7bfc68331908e30b7a52d0f98cb264393d64a94f1|Section1's explicit PUBLISHED compact29 starting theorem and its Vin84 bibliographic identification only; not the paper's stronger23 proof|
|A|docs/manuscript_assumptions.md|a5b578d3304e9adc7f185f22fe4281ed96fb08c1a6404e104bae344fd086f2cc|PublishedInputs.compactCeiling29 is the general compact Vinberg theorem, not an arithmetic-only premise|
|G|results/compact_dimension18_global_composition_and_whole_d19_exclusion_candidate_20261022.md|9a8558607afaac0889b1d359a361ec1ca42942e2dc4b725f320bf126bb382ed7|OPTIONAL separate internal compact18 interface for the final N+2 corollary only|
|GA|results/compact_dimension18_global_written_acceptance_20261022.json|f332b581b1d20eceef2af2e28905be616a711f388407ecfa5e8112acfde2de7a|Accepted exact G/global18 scope; unnecessary for the finite57 exclusion|
|F|results/finite57_far_row8_isolated_intervening_three_four_parallel_actual_cusp_financing_candidate_20261023.md|be66da2027adc88ea30f2e40e5a71c413962b8383d424997495449c0b31bad79|ONLY exact original far-row8/isolated-u roots, for the final residual application; not any new charge or contact5+ theorem|
|FA|results/review_finite57_far_row8_isolated_intervening_cusp_financing_root_acceptance_20261023.md|e36d02ae40efd815642062fdca44fef5587e4f325455ef5b1017820fac71a64d|Accepted exact F scope and different-author review|

The published external input is E.B. Vinberg, *The absence of
crystallographic groups of reflections in Lobachevsky spaces of large
dimension*, Trudy Moskov. Mat. Obshch.47 (1984),68--102; English
translation, Trans. Moscow Math. Soc. (1985),75--112. V records the
primary bibliographic link https://www.mathnet.ru/eng/mmo444 and states
the general compact n<=29 conclusion. This source uses the already
registered theorem as an EXTERNAL premise; it is not a new primary-text
scope audit or a reproof of Vinberg's bound.

## 2. The actual face exists and is compact

Normalize the two actual outward normals c,s to square1, with
c.s=-k, k=cos(pi/m). Their Gram B=[[1,-k],[-k,1]] is positive
definite. Because the pair occurs at an ACTUAL finite vertex, its face

    F=P intersect c-perpendicular intersect s-perpendicular

has dimension exactly n-2 and nonempty hyperbolic relative interior.
This follows also from finite simplicity at that vertex; an abstract
elliptic pair is not being asserted to create an ambient face.

The closure of a finite-sided finite-volume hyperbolic polytope has
only finite or ideal vertices, and its ideal boundary consists of its
finitely many ideal vertices. A face's compactification has vertices
among those of P. If an ideal vertex belonged to F, both actual c,s
would be incident there. Its whole diagram is a product of irreducible
affine diagrams. Within those diagrams the finite bonds are3/4/6;
the affine-A1 parallel pair does not create a finite m bond. Between
whole factors the product is0. Thus no native5 or native7+ pair can
be incident at any ideal vertex.

Hence F has NO ideal vertex or other ideal boundary point. Equivalently,
in the Klein compactification it is the convex hull of its finite
vertices. It is therefore compact. It is a finite-sided face, since P
has finitely many facets. There is no assumption that a finite-volume
polytope is globally compact, or that its cusps are simple.

## 3. Extraction of EVERY intrinsic facet and its finite rank-three list

Every facet of compact F has a finite relative interior and is cut by
an ACTUAL original facet f distinct from c,s. At a finite vertex of
that facet the ambient chamber is simple. Thus F has precisely its
two defining ambient facets, and each intrinsic facet comes from one
additional ambient wall. The three selected normals c,s,f are
independent and positive definite: their walls meet in hyperbolic space
at a codimension-three finite face. In particular their Schur residual
is strictly positive, not a null or indefinite projected normal.

Use unit f and put beta_f=(-f.c,-f.s), entrywise nonnegative. Its
intrinsic outward normal before normalization is the exact orthogonal
projection

    f_F=f+(c,s)*B^(-1)*beta_f,
    delta_f=f_F²=1-beta_f^T B^(-1) beta_f>0.       (1)

For two such facets f,g the exact product is

    f_F.g_F=f.g-beta_f^T B^(-1) beta_g.          (2)

These are metric orthogonal projections, not asserted actual ambient
roots or a transported lattice. Their restrictions give exactly the
original inequalities on the supporting H^(n-2), so these are ALL
intrinsic facets of F. Finite side count and nonempty interior were
established before examining their angles.

The COMPLETE finite rank3 list is A3,B3,H3 or disconnected factors.
No finite irreducible rank3 contains an m>=7 bond. For m5, its only
irreducible extension is H3: one endpoint gains a simple3 contact,
and the other has0. Two contacts give a cycle; a native4 or larger
attachment to a5 bond is not positive definite. Consequently:

* At m>=7, EVERY f is orthogonal to BOTH c,s.
* At m5, EVERY f is either orthogonal to both, or has exactly one
  simple3 contact to c or s and0 to the other.

There is no incomplete small-label cutoff, foreign block or unbounded
parameter in this classification. It is applied at finite hyperbolic
facet interiors, not at an ideal vertex or to an unproved Coxeter face.

## 4. Native7+: the intrinsic face is directly Coxeter

For m>=7 all beta_f=0. Formula(1) gives f_F=f and delta_f=1.
All intrinsic facet angles therefore have their ORIGINAL native
Coxeter values; disjoint hyperplanes remain disjoint/parallel. In
particular every actual adjacent pair of facets of F has a native
dihedral. The actual compact finite-sided F is a Coxeter polytope.

No finite-volume or compact dimension theorem was applied to F before
this explicit normal/angle verification.

## 5. Native5: all projected angles, including both H3 orientations

Let phi=(1+sqrt5)/2 and k=phi/2. For an H3-type facet with one
simple3 contact, beta is(1/2,0) or(0,1/2). Put

    Q=1/[4(1-k²)]=(5+sqrt5)/10,
    delta=1-Q=(5-sqrt5)/10>0.                    (3)

Orthogonal facets have projected squared norm1; every H3-type facet
has projected squared norm delta. Every ORIGINAL distinct wall pair
has nonpositive product, either0, -cos(pi/l) for finite l>=3, or a
parallel/dotted magnitude at least1. We now keep all those cases.

For two orthogonal-type facets formula(2) is unchanged. For one
orthogonal and one H3 facet, their normalized intrinsic product is
f.g/sqrt(delta). If f.g0 it gives a right angle. If f.g=-1/2 it is
exactly -cos(pi/10), since

    1/(4delta)=(5+sqrt5)/8=cos²(pi/10).           (4)

For EVERY other nonzero native contact, its original magnitude is at
least1/sqrt2; because delta<1/2 the new magnitude exceeds1. Original
parallel or dotted products likewise give a magnitude greater than1.
Thus the only new finite angle in this mixed pair is pi/10; all other
nonzero cases are disjoint hyperplanes, not unspecified acute angles.

For two H3-type facets at the SAME endpoint, the subtracted term in
(2) is Q. Their normalized product is at most

    -Q/delta=-phi²<-1.                          (5)

For OPPOSITE endpoint orientations, the subtracted term is kQ. Their
normalized product is at most

    -kQ/delta=-phi³/2<-1.                       (6)

The inequalities use the ACTUAL original f.g<=0, so no original
finite/dotted mutual branch is lost. They cover both endpoint roles
and every original wall defining an intrinsic facet. Such pairs do
not meet in F. They do not create an additional finite face angle.

We have proved all intrinsic adjacent-facet angles native, with only
possible newly projected label10 in(4). Hence compact F is a genuine
Coxeter polytope also at m5. This is the full special-good-face proof;
the intrinsic claim was not imported as arbitrary face inheritance.

## 6. Dimension consequence and the exact finite57 residual

Apply the EXTERNAL published compact ceiling29 to the compact Coxeter
F just proved. Its dimension is n-2, so n-2<=29 and n<=31. Thus at
n57 the hypothesized native5 or native7+ incident pair is impossible.
The argument contains no finite57 charge, ideal reserve, curvature
payment or new source-debt removal. Its numerical57 application does
not depend on the optional internal compact18 input G/GA.

The same completed geometric proof conditionally gives n<=N+2 for
ANY compact ceiling N. With accepted internal18 it gives n<=20 for
this BOND-BEARING class alone. It does not give20 for all finite-volume
polytopes. All earlier global ceilings and frozen editions are unchanged.

For the far-row8 application keep F's actual isolated finite neighbour u,
with m11/17, and let alpha=-c.s in[1,2) be finite. In unit normalization
put k=alpha/2 and H=-u.c>0. Its positive projection

    q0=u+[H/(1-k²)]c+[kH/(1-k²)]s

is future, zero on c,s, and strictly negative at the old short a and
every other unselected original wall, exactly by positive actual-wall
coefficients. It retains A_(m+4)+the two inert Ds. Those rank55 roots
plus c,s form a full positive rank57 forest. Thus its normalization is
an ACTUAL finite vertex on the c,s face, for BOTH m11 and17; no new
financing theorem or firstness claim is assumed to supply it.

At n57 any such label5 or>=7 is therefore excluded by this lemma.
After the separately accepted3/4
and parallel rules, the remaining finite5+ contact is therefore ONLY
native6, in BOTH affected-D17 and affected-D23 cases. Existing generic
finite-contact financing theorems remain valid as written; no frozen
source has to be changed to use this further restriction after review.

Native6 is deliberately NOT excluded. Its positive finite rank3 facet
extensions are orthogonal, so its corresponding intrinsic face can be
Coxeter, but its two roots CAN occur in an affine G2 cusp. The crucial
compactness step then fails. A finite-volume native6 face is not supplied
a compact dimension bound. Native2/3/4, strictly dotted contacts,
nonisolated endpoints and full global financing/coverage likewise require
their own arguments. No universal364 cap or56 conclusion is claimed.
