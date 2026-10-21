# Independent global quadrilateral review over Q(√2)

The complete quadrilateral reduction passes independent arithmetic and a hostile proof review. Together with the separately reviewed weight accounting, it supports the research theorem candidate **n ≤ 20 for compact arithmetic hyperbolic Coxeter polytopes with ground field Q(√2)**. This is not a statement about arbitrary compact polytopes, other arithmetic fields, or finite-volume noncompact chambers. No priority or current-best claim is made.

The [standalone checker](../scripts/fresh_sqrt2_quad_exhaustion_20261015_review.py) and [certificate](fresh_sqrt2_quad_exhaustion_20261015_review.json) bind the final four primary source/data pairs, the [complete Lannér catalogue](fresh_sqrt2_lanner_catalogue_20261015.md) through their recursive manifests, the [independently reviewed dotted-tail theorem](fresh_sqrt2_dotted_tail_20261015_review.md), and the [separate independent weight review](fresh_sqrt2_weight_bound_20261015_review.md). The [primary manuscript](fresh_sqrt2_quad_exhaustion_20261015.md) supplies the complete proposed argument. Its source/data are bound; its review-status prose is not an arithmetic dependency.

## Independently reconstructed finite calculations

The checker imports no producer module, project arithmetic, SymPy, NetworkX, numerical eigensolver, or interval sign engine. It uses a standalone Fraction field with basis

\[
1,s,t,st,u,su,tu,stu,\qquad s^2=2,\ t^2=3,\ u^2=2+s.
\]

Successive quadratic norms give inverses, and recursive exact comparisons give signs. Matrix inertia is computed by direct congruence, including a change of basis when every remaining diagonal vanishes. This differs from the primary algebraic-field corner-Schur implementation.

The independently ordered traversal reconstructs the 51 marked minimal-piece representatives and all 6,010 white-node identifications. It tries free contacts in a different order and tests finite diagrams using the classification by path labels, branching degree and arm lengths. It reproduces the **exact set of 115 labelled union supports**, rather than only its cardinality. A separate exhaustive permutation canonicalization, preserving the four black nodes and opposite relation, obtains **52 orbits**; all 115 saved full-pair bijections are checked individually.

Direct cycle products in the radical field and full-matrix physical/conjugate inertia reproduce the filters: **40 two-negative-direction exclusions, one cyclic-field exclusion, nine finite-label survivors and two dotted families**. Every one-extra-white contact domain is regenerated independently. All **44** finite-core patterns and all four surviving patterns agree. For rows rejected by a conjugate kernel, the independent calculation instead checks that the full conjugate extension actually has a negative direction.

The **18** dotted structural profiles are also regenerated. Four are rejected by literal positive affine (4,3,3) kernels, two by a (3,3,4) Lannér triangle orthogonal to the dotted pair, and twelve have contacts only at the dotted endpoints. The dotted distance is never discretized. The two-negative-direction argument for the orthogonal pieces works for every strict dotted distance.

The checker further verifies all four complete marked closure maps, the 16 pendant-diagonal cases, the D5 short-arm continuation through E6/E7/E8 and its first affine obstruction, and the core-37 determinant identity. For the latter, four exact evaluations establish the entire polynomial because the determinant is separately affine in the two variable diagonal entries. These are algebraic identity checks, not sampled sign evidence.

## Why the initial union enumeration is exhaustive

The full white support W of a compact quadrilateral has spherical Gram. Each W together with two adjacent black nodes is an actual finite corner. Each opposite extension is indefinite: a positive-definite extension would force the opposite facets to meet, and an affine extension is forbidden in a compact Coxeter chamber. Minimal indefinite subdiagrams therefore contain both opposite black nodes and belong to the complete arithmetic Lannér catalogue.

All contacts on the union of these two pieces are specified or enumerated. Apart from the already specified opposite pairs, every pair lies in an actual finite corner, so the finite alphabet applies. Shared white nodes, extra cross contacts and all abstract quadrilateral relabellings are included. Two orthogonal Lannér pieces are impossible already as a principal subdiagram, regardless of later white attachments. Requiring no marked path of length at most two is inherited by principal restrictions and is a valid necessary condition for the exceptional case.

No inference that the initial union equals the full star is used. The next steps supply the missing arbitrary-attachment argument.

## Arbitrarily many external white nodes

**Seven sealed finite cores.** Each individual external white node must pass the complete single-node domain. For cores 5, 9, 17, 32, 41, 47 and 50, every nonzero contact is excluded. Thus every external node is orthogonal to the entire core. An indirect path cannot evade this conclusion, since it would need a first edge from a core node. In particular the exceptional finite corner components used by the weights are entire components, not proper subdiagrams.

**Core 14.** Adding only orthogonal positive white blocks cannot turn its nonsingular Gram into the corank-one full quadrilateral-star Gram. Hence some external contact exists. Its unique allowed pattern gives sealed core 17, by the checked full-pair map. Every remaining external node is then orthogonal to this larger core.

Here the full-star rank assertion is geometric: W has rank n−2 and its Lorentzian complement is three-dimensional. The four projected sides of a compact hyperbolic quadrilateral span this complement. Otherwise their common orthogonal direction would make the intersection unbounded or degenerate. The n+2 original facet normals consequently have rank n+1 and corank one. Orthogonal external positive blocks add equally to size and rank.

**Core 37.** The only external contacts are 3-bonds to endpoints 0 and/or 2. An actual corner containing either endpoint already contains D5. Any connected finite diagram containing this D5 is a simply-laced D/E tree with a single branching vertex. An external component cannot touch the same endpoint twice, which would create a cycle. Two components at that endpoint would create a second branching vertex. Nor can the external component itself branch, for that branch remains visible in a corner retaining the attached endpoint. Hence a nontrivial attachment extends the specified short arm by a simple path. The arm classification permits only lengths one, two or three; length four already contains affine E8, so this is an all-size argument.

A component joining the endpoints is therefore a path of one, two or three white nodes. The checked maps identify the resulting cores with sealed 41, 47 or 50. If no component joins the endpoints, the remaining paths reduce the unit diagonals to p,q in {1,3/4,2/3,5/8}. The exact rank condition is

\[
(p+q-4pq)+(1-p-q)\sqrt2=0.
\]

Its irrational coefficient forces p+q=1, contradicting p+q≥5/4. Thus no further case is hidden behind a chain of nodes initially orthogonal to the core.

**Dotted cores.** The complete 18-profile calculation first proves that all external contacts occur only at the two dotted endpoints a,y. If an external component contacts both, eliminate all external whites simultaneously. Their Gram is positive definite and acute, so its inverse is entrywise nonnegative, and is strictly positive within each connected component. In unit normalization every inverse diagonal is at least one. Both reduced endpoint diagonals therefore satisfy 0<p,q≤3/4, and the effective norm-two opposite magnitude increases strictly: r_eff>r>2. Positivity of p,q follows from the corresponding actual finite corners. Other pendant components can only strengthen these bounds.

Elimination creates no other core contacts. The independently reviewed rank identities thus still apply. The single-contact identity gives

\[
r_{\rm eff}^2=q(4p-1+\sqrt2)
\le\tfrac34(2+\sqrt2)<4.
\]

For the double-contact identity, with t=(1−√2)/2<0, its only positive root is

\[
r_{\rm eff}=-t+\sqrt{(2p-t)(2q-t)}
\le\tfrac12+\sqrt2<2.
\]

Both contradict a dotted connection. These inequalities allow arbitrary finite connecting components and all allowed endpoint labels; there is no pendant assumption in this exclusion. With connecting components excluded, the remaining components attach to just one endpoint, and the previously reviewed all-size pendant theorem applies. Orthogonal extra white components do not affect the determinant conditions or the affected corner components.

## Dimension scope and verification

The closure above proves that every compact Q(√2)-arithmetic quadrilateral either has a marked distance-two corner or belongs to one of the sealed finite cores or three reviewed tail survivors. It supplies exactly the entire-component premise required by the independent weight review. That review reconstructs the literal weights, proves the all-rank vertex upper bound 2n−2, and checks the strict h-vector contradiction for every n≥21. Its arithmetic is inherited here and is not claimed as a second full weight replay. The triangle argument was also source-reviewed: the minimal Lannér piece contains all three black nodes, and either its connected triangle or its C4/C5 gaps supply two short marked paths.

The published compactness, spherical-diagram, arithmetic-conjugacy and h-vector inputs remain explicit mathematical inputs. This review does not extend them to noncompact chambers. Compactness is essential to the absence of affine subdiagrams and to all four finite corners.

Independent generation and saved replay pass. Recursive source manifests, Python syntax and local report links pass. No remaining mathematical discrepancy was found in the stated compact Q(√2) classification or its combination with the separately reviewed weights. Frozen predecessor artifacts were not changed.
