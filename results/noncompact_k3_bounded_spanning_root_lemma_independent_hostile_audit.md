# Independent hostile audit: finite-volume K3 chambers have a bounded spanning root packet

**Verdict (2026-08-11): PASS in the exact finite-sided, finite-volume,
full-Picard/full-nef scope.** I reconstructed the proof without using the
owner certifier. For every \(3\leq\rho\leq20\), the closed star of a facet of
minimal visual aperture spans the rank-\(\rho\) Lorentz space and has all
normalized off-diagonal pairings strictly below \(14\). Literal K3
square-\((-2)\) facet curves therefore contain a rational basis

\[
 E_1,\ldots,E_\rho,
 \qquad 0\leq E_i.E_j\leq13\quad(i\neq j).
\]

The finite-host and one-complete-Weyl-chamber-orbit conclusions follow. No
compactness, cusp simplicity, affine-diagram classification, reflective-
lattice classification, or reflection maximality is hidden in the proof.

I found one implicit geometric sublemma, not a defect: the statement that
ultraparallel facet shadows are disjoint uses the irredundancy already
contained in the word *facet*. In the Klein ball, two noncrossing shadows are
either disjoint or nested; nesting would make one of their supporting
half-space inequalities redundant. Section 2 below supplies this argument
explicitly.

## 1. Frozen owner artifacts

This independent audit pins all five owner artifacts:

- `results/noncompact_k3_bounded_spanning_root_lemma.md`, SHA-256
  `de3b78a9070dc7035a8ecf51c5aa697b0736ad1559f0db8867b6629eb3cc4446`;
- `results/noncompact_k3_bounded_spanning_root_lemma_hostile_scope_audit.md`,
  SHA-256
  `cfaf5b5fd10c97e14f28105b99c379c8f9c391c6c89ffd9875ceefd1426ff851`;
- `scripts/certify_noncompact_k3_bounded_spanning_root_lemma.py`, SHA-256
  `895a0c1fe0e10d38f5516ab917154bde0d3a8d2eef89d1ad567260abaf2e5ffb`;
- `data/enumerations/noncompact_k3_bounded_spanning_root_lemma.json`,
  SHA-256
  `5c5ffa217c15aebb189102c3a92937c5e20fe5e465a00caaac44780fa9761874`;
- `tests/test_noncompact_k3_bounded_spanning_root_lemma.py`, SHA-256
  `bf4fdfd82d1c1f7529cb99ed7a62a9c11822d534755f304273eccc373707fd93`.

## 2. The visual-gap estimate is valid in exactly its asserted branch

Fix an interior observer \(O\). The visual shadow of a supporting hyperplane
is a round spherical cap of full aperture \(0<\theta<\pi\). In the Klein
model it is the set of directions whose rays meet the hyperplane before
leaving the ball.

There is a small irredundancy point to check. If two hyperplanes do not meet
in the open ball, their boundary spheres do not cross. The two observer-side
shadow caps are therefore disjoint or one is nested in the other. In the
nested case, the corresponding cut-off half-ball is contained in the other
cut-off half-ball. Equivalently, one inward half-space contains the other
inward half-space throughout the Klein ball. The weaker inequality cannot
support a facet of the irredundant polyhedron. Thus distinct ultraparallel
**facet** supports have disjoint visual shadows. Asymptotic shadows may touch
in one ideal direction and belong to the separate \(b=2\) branch.

Now let \(a=\theta_F/2\) and \(b=\theta_G/2\) be the half-apertures of two
ultraparallel supports, and let \(c>0\) be the angular gap between their
shadows. Their center directions are separated by \(a+b+c\). In the
observer-frame Lorentz normalization, the absolute unit-normal pairing is

\[
 p(F,G)=
 \frac{\cos a\cos b-\cos(a+b+c)}{\sin a\sin b}.
\]

The normalized Coxeter pairing is \(2p(F,G)\). Product-to-sum gives the exact
identity

\[
 2p(F,G)
 =4\frac{\sin(a+c/2)\sin(b+c/2)}{\sin a\sin b}-2,
\]

which is precisely the owner formula after substituting
\(a=\theta_F/2\) and \(b=\theta_G/2\). This formula is used only for a
positive gap. If the hyperplanes intersect in hyperbolic space, the direct
normal formula gives \(b=2|\cos\alpha|\leq2\); if they are asymptotic, it
gives \(b=2\).

Let \(E\) have minimal aperture and suppose both closed facets \(F,G\) meet
\(E\). Choose directions through points of
\(\overline F\cap\overline E\) and
\(\overline G\cap\overline E\). They lie in the closed shadow of \(E\), so

\[
 0<c\leq\theta_E\leq\min\{\theta_F,\theta_G\}.
\]

This remains true if either chosen incidence occurs only at infinity. For
\(0\leq c\leq\theta<\pi\),

\[
 \frac{\sin((\theta+c)/2)}{\sin(\theta/2)}<2.
\]

Indeed, for \(\theta\leq\pi/3\) it is at most
\(2\cos(\theta/2)<2\); for \(\theta>\pi/3\) it is at most
\(1/\sin(\theta/2)<2\). Both factors in the identity are therefore below
two, and

\[
                         b(F,G)<14.
\]

This proves the claimed bound for every pair in the projective closed star.
Compactness and simplicity play no role in it.

## 3. A complete nonsimple cusp packet spans \(u^\perp\)

Let \([u]\) be any ideal boundary point, not provisionally assumed to be a
vertex. Move it to infinity in upper half-space. Because there are finitely
many facets, a single sufficiently high horosphere lies above every bounded
hemisphere belonging to a nonincident support. Above that height the chamber
is exactly

\[
                         C\times(T,\infty),
\]

where \(C\subset\mathbf R^{\rho-2}\) is the full-dimensional Euclidean
polyhedron cut out by all incident vertical walls. Its volume is

\[
 \operatorname{vol}_{\mathrm E}(C)
 \int_T^\infty y^{-(\rho-1)}\,dy.
\]

The integral is finite and positive for \(\rho\geq3\). An unbounded
full-dimensional finite-sided Euclidean polyhedron has infinite Euclidean
volume, so finite hyperbolic volume forces \(C\) to be bounded.

Choose a null partner \(v\) with \((u,v)=1\) and negative-definite
\(A=(\mathbf Ru+\mathbf Rv)^\perp\). An incident normal has the unique form

\[
                         r_i=a_i u+w_i,
 \qquad w_i\in A.
\]

After orienting inequalities, the link is

\[
                  C=\{z:a_i+(z,w_i)>0\text{ in its interior}\}.
\]

Boundedness says that the recession cone is zero. Hence the \(w_i\) span
\(A\) and positively span in every direction; in particular there are
\(\lambda_i>0\) with

\[
                         \sum_i\lambda_iw_i=0.
\]

At an interior point \(z_0\), summing the strict inequalities gives

\[
 \sum_i\lambda_i a_i
 =\sum_i\lambda_i\bigl(a_i+(z_0,w_i)\bigr)>0.
\]

Thus the same dependence of the lifted normals is a nonzero multiple of
\(u\). Their span contains \(u\), and subtracting the \(a_i u\) terms then
shows that it contains all of \(A\). Consequently

\[
                         \operatorname{span}\{r_i\}
                         =\mathbf Ru+A=u^\perp.
\]

No simplex assumption occurred. The independent executable additionally
replays bounded box links with \(2(\rho-2)\) incident square-\((-2)\) lifts
in every rank \(3\) through \(20\); these links are nonsimple from rank four
onward and still have packet rank \(\rho-1\).

## 4. Two extreme points of one facet span the Lorentz space

The Klein closure \(Q=\overline E\) is a compact convex set of dimension
\(\rho-2\geq1\), so it has two distinct extreme points \(x,y\).

If \(x\) is inside the ball, it is a finite vertex. Were its active facet
normals not to span \(x^\perp\), their common kernel would have a tangent
direction independent of \(x\). A sufficiently small displacement along
that direction stays timelike, preserves every active equality, and keeps
all inactive inequalities strict. Projectivizing would give a segment of
\(Q\) through \(x\), contradicting extremality.

If \(x\) lies on the ideal sphere, Section 3 applies before any vertex
assumption: all facets incident at \(x\) span \(x^\perp\). Their common
projective intersection is consequently the single point \([x]\), so \(x\)
is a genuine finite-volume ideal vertex, possibly nonsimple. The same
argument applies to \(y\).

Every facet through either point lies in the closed star of \(E\). Since
\(x,y\) are distinct projective points,

\[
 (x^\perp+y^\perp)^\perp
 =\mathbf Rx\cap\mathbf Ry=0.
\]

Nondegeneracy of the ambient Lorentz form gives

\[
                         x^\perp+y^\perp=L_{\mathbf R}.
\]

Thus the closed star spans. Crucially, the visual estimate of Section 2
applies to **all** pairs in that star, including cross-pairs from the two
vertex packets. A duplicate of one packet would remain rank \(\rho-1\);
the second projectively distinct extreme point is exactly the transverse
repair.

## 5. Physical K3 roots give finite hosts

In the asserted full-Picard/full-nef lane, a facet normal is the literal
class of an irreducible smooth rational K3 curve. Hence its square is
\(-2\), and two distinct facet curves have nonnegative integral
intersection. For such physical normals,

\[
 \frac{2|(E_i,E_j)|}{\sqrt{|E_i^2E_j^2|}}=E_i.E_j.
\]

Choose \(\rho\) linearly independent star roots. Integral vectors that span
over \(\mathbf R\) also span over \(\mathbf Q\), and the strict normalized
bound becomes

\[
                         0\leq E_i.E_j\leq13.
\]

Their lattice \(R=\sum_i\mathbf ZE_i\) has a nonsingular Gram matrix with
diagonal \(-2\) and every off-diagonal entry in
\(\{0,1,\ldots,13\}\). For fixed rank there are at most
\(14^{\binom\rho2}\) such matrices before signature and determinant
filters. Since \(R\) has full rank in the integral Picard lattice \(L\),

\[
                          R\subseteq L\subseteq R^*.
\]

The quotient \(R^*/R\) has order \(|\det R|\), so it has only finitely many
subgroups and therefore only finitely many intermediate integral
overlattices. The finite rank range \(3\leq\rho\leq20\) proves host
finiteness. Projected roots would not justify the equality between the
normalized and integral pairings; the literal same-vector hypothesis is
essential.

## 6. Complete root chambers form one Weyl orbit per oriented host

The square-\((-2)\) mirror arrangement of a fixed integral hyperbolic
lattice is locally finite. Indeed, if a root mirror meets a fixed compact
set \(K\) in hyperbolic space, at a meeting point the root has fixed length
in the negative-definite tangent metric. Those tangent norms are uniformly
equivalent to a fixed Euclidean norm as the point ranges over \(K\). The root
therefore lies in a bounded subset of the discrete lattice, containing only
finitely many lattice vectors.

Choose interior points in two complete mirror chambers. A generic path
between them has compact image, avoids codimension-two intersections, and
crosses only finitely many mirrors. Reflection in each crossed
square-\((-2)\) root is the integral isometry

\[
                         s_r(x)=x+(x,r)r.
\]

It preserves the chosen positive-cone component, fixes its wall, and swaps
the two adjacent chambers. Successive crossed-wall reflections carry the
first chamber to the second. Thus the complete root chambers are one Weyl
group orbit for each oriented lattice.

This last statement needs completeness of the full Picard root arrangement.
It does not need reflection maximality. A selected subarrangement could have
several unrelated chambers and is outside the theorem.

## 7. Exact boundary and reproducibility

The theorem begins only after a finite-sided, finite-volume, exact
full-Picard/full-nef K3 carrier has been returned. It does not construct that
carrier, establish canonical/log/MMP fidelity, or prove exhaustive Property
4 for every raw maximal arithmetic polytope. It also says nothing about an
infinitely-faced chamber whose aperture infimum need not be attained.

The independent replay is
[audit_noncompact_k3_bounded_spanning_root_lemma_independent_hostile.py](../scripts/audit_noncompact_k3_bounded_spanning_root_lemma_independent_hostile.py),
with regressions in
[test_noncompact_k3_bounded_spanning_root_lemma_independent_hostile.py](../tests/test_noncompact_k3_bounded_spanning_root_lemma_independent_hostile.py).
It imports no owner code. It pins the five owner artifacts, proves the sine
identity symbolically, samples the sharp analytic interface, replays complete
nonsimple cusp packets and all four finite/ideal two-vertex combinations in
every K3 rank, checks the integral Gram box and finite discriminant-group
interface, and verifies sample root reflections as integral involutive
isometries.

The exact audited conclusion is

\[
 \boxed{\text{finite volume plus literal full-nef K3 return forces a
 bounded full-spanning physical root packet with }B=13.}
\]
