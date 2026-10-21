# Independent proof and source review: simple-cusp bound61

**PASS: no mathematical or source-level gap found in the stated scope.** The conclusion is that a finite-sided, finite-volume hyperbolic Coxeter polytope whose ideal vertices are all simple has dimension at most61. It permits finite vertices and imposes no arithmetic hypothesis. It does not change the unrestricted163 bound or the compact-facet127 branch.

This review read the [successor source](../scripts/fresh_simple_cusp_affine_exclusion_20261017.py), its [saved data](fresh_simple_cusp_affine_exclusion_20261017.json), the draft proof, and the inherited distance15 proof. It did not independently repeat the producer's 7,521 marked-graph calculations, 245 inherited graph checks or quasi-Lanner census. The producer reports successful generation and fresh saved replay. The review independently checked the geometric implications, classification scope, graph estimates and all-rank threshold algebra.

Bound frozen inputs:

* Source SHA-256: `a710c0f90911ec582ba6133838468466eaabebbd1b2f2982335d4e344bcc7c5c`.
* Data SHA-256: `01687a2e7e9896d0f077efad28f1633d61f7f14b12c20813358e823f03bbf34a`.

The producer report was still receiving verification-status wording at review time, so its hash is deliberately not bound here.

## Affine subdiagrams and actual ideal vertices

For a connected affine principal diagram, its positive kernel gives a normal combination \(z\) with \(z^2=0\). It is not the zero vector: its product with a common strict future timelike interior point is strictly negative. This also fixes its future orientation. Its products with all other outward wall normals are nonpositive by acuteness, so it is a feasible null ray.

Finite volume matters at the next step. The projective closure is a finite convex polytope contained in the closed Klein ball. Strict convexity of that ball makes every feasible boundary-sphere point extreme; therefore this ray is an actual ideal vertex. It cannot be merely an unattained direction or an interior point of a larger projective face. All roots in the affine diagram belong to its incident-wall diagram.

At a simple ideal vertex the Euclidean link is one simplex. Its full incident diagram is irreducible affine of order \(n\); every proper principal subdiagram is positive definite. The original singular affine diagram must consequently be the whole incident diagram. This argument includes affine subdiagrams consisting of white roots or containing all marked black roots.

The primary [Potyagailo–Vinberg paper, printed pp.64 and69](https://ems.press/content/serial-article-files/43032) supplies the local link structure and finite-volume convex-hull convention. There is also direct support in [Felikson–Tumarkin, Section1, PDF pp.3–4](https://arxiv.org/pdf/math/0502413): parabolic subdiagrams extend into ideal-vertex diagrams, and a simple ideal vertex has connected parabolic diagram. Their condition(iii) follows from simplicity; their different condition(ii), excluding elliptic diagrams of order \(n\), uses the all-ideal hypothesis. The present mixed-vertex argument does not import condition(ii).

## Minimal cores and payments

The inherited classification of inclusion-minimal indefinite diagrams gives Lanner or quasi-Lanner order at most10. A non-Lanner minimal core has a connected affine proper subdiagram of order at most9. This contradicts the affine order \(n\) result when \(n\ge10\). Thus every minimal indefinite core used in the argument is Lanner of order at most5. Smaller ambient dimensions already satisfy the final bound.

For an actual simplex face, deleting any black side root leaves roots at a finite or ideal vertex, hence a positive-semidefinite block. Every indefinite subdiagram therefore contains all black roots. The usual minimal-core argument then gives simplex dimension at most4, without assuming the induced simplex has Coxeter angles.

For a triangle, a spanning tree of its minimal core has two paths joining pairs of black nodes without the third in their interiors; each has length at most \(5-2=3\). For a quadrilateral, opposite-side blocks are indefinite: a positive-definite block would produce a finite common point, and an affine block a feasible ideal common point, contradicting disjointness of the **closed** opposite polygon edges. Each opposite minimal core contains both marked roots and has radius at most3 from that pair. The two cores must overlap or have a joining edge, since two disjoint orthogonal indefinite blocks would give two negative directions. A shortest path between their black sets has only white interior nodes and length at most7. It belongs to an actual corner diagram.

Consequently distance7 corner weights pay at least2 on triangles and at least1 on quadrilaterals, including ideal corners. This does not assume all polygon vertices are finite.

## All ranks and strict threshold

The inherited path, cycle, single-fork and double-fork formulas give at most \(N\) unordered pairs at each fixed positive graph distance. In particular the affine double-fork family is included in every rank. The small exceptional finite/affine diagrams have order at most9, so all their pairs total at most \(4N\le7N\). The affine D4 four-leaf star is treated as an exception. Thus every finite rank-\(n\) forest and every simple affine order-\(n\) vertex costs at most \(7n\). Finite sampled ranks merely confirm these explicit formulas; they do not stand in for a tail proof.

The strict simple-polytope lower bound is

\[
b_n=\begin{cases}n(n-6)/8&n\text{ even},\\(n-1)(n-5)/8&n\text{ odd}.\end{cases}
\]

Independent substitution gives

\[
8(b_{62+2t}-7(62+2t))=4t^2+124t,
\qquad
8(b_{63+2t}-7(63+2t))=4t^2+128t+68.
\]

Both are nonnegative for every integer \(t\ge0\). Equality at62 is excluded because the lower payment inequality is strict. Hence \(n\le61\).

The checker uses exact arithmetic, preserves source manifests and implements the stated checks. Its finite threshold loop confirms polynomial identities whose displayed algebra supplies the unbounded tail. This is an internal proof/source review, not external refereeing or formal proof-assistant verification. No novelty or best-published-bound assertion is supported or needed. The published simple **all-ideal** bound8 has a different hypothesis.
