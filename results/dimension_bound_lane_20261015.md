# Dimension bounds: hypotheses, exact inequalities, and next proof targets

This lane targets dimension bounds for classes of hyperbolic Coxeter polytopes. It is separate from construction and from obstructions tied to one lattice or retained chamber. The first completed result is a derived necessary inequality for **all compact Coxeter polytopes**, and stronger consequences for a precisely defined compact combinatorial class. These compact vertex-rank inequalities themselves do **not** improve the unrestricted bounds 29 or 995; the separate Prokhorov refinement is discussed below. No publication-priority claim is made for the derived inequalities.

The [standalone exact checker](../scripts/dimension_bound_vertex_components_20261015.py) and [certificate](dimension_bound_vertex_components_20261015.json) verify the scalar identities, exact partition maxima, threshold tables, and quantitative consequences below. The geometric proof is given here; a finite table is not substituted for a universal proof.

## Hypothesis matrix and primary references

All dimensions below are hyperbolic dimensions. “Finite volume” includes compact polytopes unless explicitly qualified. These are upper bounds, not assertions that the endpoint dimension is attained.

| Class | Established upper bound | Primary reference / qualification |
|---|---:|---|
| Compact Coxeter polytopes, no arithmetic assumption | 29 | [Vinberg, 1984](https://www.mathnet.ru/eng/mmo444) |
| Finite-volume Coxeter polytopes, no arithmetic assumption | 995 | [Prokhorov, 1986/1987](https://www.mathnet.ru/eng/im1484); the compact case is already smaller |
| Arithmetic finite-volume reflection groups | 29 | Combines the compact bound 29 with the noncompact arithmetic bound 21 |
| Noncompact arithmetic reflection groups | 21, with dimension 20 excluded | [Esselmann, 1996](https://doi.org/10.1006/jnth.1996.0141); the theorem concerns Lorentz lattices of rank at most 22 |
| Compact right-angled polytopes | 4 | [Potyagailo–Vinberg, §2](https://ems.press/content/serial-article-files/43032), giving the Nikulin proof |
| Compact 3-free Coxeter polytopes | 12 | [Alexandrov, Theorem C and §4](https://arxiv.org/pdf/2203.07248) |
| Simple ideal Coxeter polytopes | 8 | [Felikson–Tumarkin, 2008](https://doi.org/10.4213/im852); both “simple” and “ideal” are required |
| Compact polytopes with $d+k$ facets, $k=5,6,7,8,9,10$ | $9,12,15,18,22,26$, respectively | [Burcroff](https://arxiv.org/pdf/2201.03437); these facet-excess restrictions are not unrestricted dimension bounds |

The literature check includes current primary research through the session date, including [McLeod's October 2026 revision on $d+4$ facets](https://arxiv.org/abs/2608.02894). A significant historical qualification emerged during the proof audit: [Nikulin's ICM 1986 report, printed p. 659](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf), announces a general bound below 300, with an outline but without the complete constants and affine-cusp case analysis. Thus the table records the standard fully developed bounds, not a claim that no stronger announcement exists. The known general construction records remain compact dimension 8 and noncompact dimension 21, as also recorded in [recent construction work](https://arxiv.org/pdf/2607.14715).

The [earlier literature audit](literature_frontier_recheck_20261013.md) and [noncrystallographic face gate](fresh_noncrystallographic_face_gate_20261014.md) remain unchanged. In particular, an arbitrary face need not be Coxeter. [Allcock's Theorems 2.1–2.2](https://arxiv.org/pdf/0903.0138) give a useful sufficient spherical-support condition, not a license to replace every face by a lower-dimensional Coxeter polytope.

## A compact vertex-component inequality

Let $P\subset\mathbb H^d$ be a compact finite-sided Coxeter polytope, $d\ge4$. It is simple. At a vertex $v$, let the irreducible components of its spherical Coxeter diagram have ranks

\[
r_1(v),\ldots,r_{c(v)}(v),\qquad \sum_i r_i(v)=d,
\]

and put $B(v)=\sum_i\binom{r_i(v)}2$. Write $f_i$ for the number of $i$-faces.

Here $B(v)$ counts **nonright intrinsic polygon angles at $v$**, not diagram edges. Indeed, the $d$ edge directions at a simple vertex are the negatives of the dual incident normals. Their scalar products are the entries of the inverse positive-definite incident Gram matrix. Each irreducible positive-definite matrix with nonpositive off-diagonal entries has a strictly positive inverse within its block; different blocks have zero inverse entries. Thus two edge directions form a nonright angle exactly when their released nodes belong to the same component.

A compact hyperbolic triangle with angles at most $\pi/2$ has at least two nonright angles; a quadrilateral has at least one. Every polygon with at least five sides satisfies the following inequality without an angle count. Consequently, for every two-face $F$,

\[
\#\{\text{sides of }F\}+\#\{\text{nonright angles of }F\}\ge5.
\]

Double counting vertex–two-face incidences gives

\[
5f_2\le\binom d2f_0+\sum_v B(v). \tag{1}
\]

The strict Nikulin average inequality, in the explicit form reproduced in [Potyagailo–Vinberg §2](https://ems.press/content/serial-article-files/43032), is

\[
\frac{\binom d2 f_0}{f_2}<U_d,
\qquad
U_d=\begin{cases}4(d-1)/(d-2)&d\text{ even},\\4d/(d-1)&d\text{ odd}.\end{cases}
\]

Combining it with (1) proves the necessary condition

\[
\boxed{\frac1{f_0}\sum_v B(v)>T_d},\qquad
T_d=\begin{cases}d(d-6)/8&d\text{ even},\(d-1)(d-5)/8&d\text{ odd}.\end{cases} \tag{2}
\]

This proof uses neither arithmeticity nor a bound on Coxeter labels or the number of facets. The earlier [global-bound audit, §5.4](ag_global_bound_and_transition_constant_hostile_audit.md) used the rank-at-most-two instance of the same counting mechanism. The general rank formulation here is a derived use of that classical mechanism.

### Consequence for bounded irreducible ranks

Suppose every spherical component at every vertex has rank at most $r\ge1$. If $d=qr+s$, $0\le s<r$, then

\[
B(v)\le C_{d,r}:=q\binom r2+\binom s2. \tag{3}
\]

This is the exact partition maximum: moving one unit from a smaller positive part into a larger part below $r$ increases the sum of binomial coefficients. Iterating leaves only full $r$-parts and one remainder.

Equations (2)–(3) imply

\[
\boxed{d\le4r}. \tag{4}
\]

The first excluded dimension needs the exact maximum, not a rounded estimate:

\[
C_{4r+1,r}=2r(r-1)=T_{4r+1}.
\]

For even $d\ge4r+2$, the weaker bound $B(v)\le d(r-1)/2$ suffices because the difference from $T_d$ is $d(d-4r-2)/8\ge0$. For odd $d=4r+3+t\ge4r+3$, that difference equals $((4r+3+t)(1+t)+5)/8>0$. These formulas prove the entire tail, independently of any finite table.

A distribution version is also useful. If $p$ is the fraction of vertices having a component of rank greater than $r$, and $T_d\ge C_{d,r}$, then

\[
\boxed{p>\frac{T_d-C_{d,r}}{\binom d2-C_{d,r}}}. \tag{5}
\]

For example, a compact $\mathbb H^{20}$ Coxeter polytope would require more than $5/36$ of its vertices to have a component of rank at least 3. In dimension 29, more than $4/27$ must have a component of rank at least 4. These are necessities, not existence assertions.

## A stronger charge for the compact 3-free class

**Definition used here:** $P$ is *3-free* if every collection of facets with empty intersection contains a disjoint pair. Equivalently, every minimal missing face of the facet nerve has size two; for a simple polytope this is the flag condition on its dual simplicial boundary. This is Alexandrov's hypothesis. It is stronger than merely having no triangular two-faces, and it is not merely absence of a missing face of size three.

Here is the quadrilateral charge behind [Alexandrov, Lemma 4.3](https://arxiv.org/pdf/2203.07248), retaining the actual diagram-edge count. A quadrilateral has four side facets in addition to its common support. Opposite side facets must be globally disjoint: otherwise all pairs in the corresponding empty family would meet, contrary to 3-freeness. If all four adjacent side-facet pairs were right, their Gram matrix would split into two indefinite opposite-pair blocks, impossible in a Lorentz space. Therefore some adjacent pair gives a direct nonright edge in the spherical diagram at its corner. Charge the quadrilateral to that corner and that edge. Simplicity makes the charge injective, since the corner and released facet pair determine the two-face. Spherical diagrams are forests, so

\[
q_4\le\sum_v |E(S_v)|=\sum_v(d-c(v)). \tag{6}
\]

Here $q_4$ denotes the number of quadrilateral two-faces, distinct from the number $f_4$ of four-dimensional faces. The same flag hypothesis excludes triangular two-faces. Thus $5f_2\le\binom d2f_0+q_4$, and Nikulin's strict inequality gives $q_4/f_0>T_d$. With (6),

\[
\boxed{\frac1{f_0}\sum_v c(v)<d-T_d}. \tag{7}
\]

The crucial distinction is that (6) charges **direct diagram edges**; the count $B(v)$ in (1) instead concerns all pairs in a common component. They must not be interchanged.

In dimension 12, (7) says the average number of components is strictly below three. Some vertex therefore has at most two components and a component of rank at least six. Combining (7) with $c(v)\ge\lceil d/r\rceil$, and with the general bound above, gives:

| Maximum component rank at every vertex | General compact bound from (4), capped at 29 | Compact 3-free bound from (7) |
|---:|---:|---:|
| 1 | 4 | 4 |
| 2 | 8 | 8 |
| 3 | 12 | 10 |
| 4 | 16 | 11 |
| 5 | 20 | 11 |
| 6 | 24 | 12 |
| 7 | 28 | 12 |

In particular, the bounds 11 and 10 are conditional class bounds. They do not replace 29 or 995. Without the rank restriction, (7) recovers Alexandrov's upper bound 12.

## A concrete next optimization, with the missing premise exposed

One finite target is the compact 3-free dimension-12 case. Nikulin forces $q_4>9f_0$. A proof that universally $q_4\le9f_0$ in that class would exclude dimension 12. The existing charge only gives $q_4\le\sum_v(12-c(v))\le11f_0$, so it does not do this.

The next proposed local catalogue concerns spherical vertex diagrams of rank 12 with one or two components, marked by their direct edges and possible quadrilateral charges. A necessary linear program would use vertex-type frequencies, marked-edge capacities, four-corner balance for quadrilaterals, and independently proved forbidden local attachments. Any unrestricted $I_2(m)$ families must remain symbolic or be covered by a proved all-label reduction. An arbitrary label cutoff would invalidate a dimension conclusion.

The [dimension-independent Schur-inertia ledger](fresh_h9_lorentz_inertia_20261015.md) supplies some reusable forbidden configurations. Its positive determinant plus negative residual diagonal implies a second negative direction, in every ambient dimension. Its negative-determinant exceptions can have signature $(10,1)$, so they are only rank obstructions below dimension ten. Neither category proves that a generic high-dimensional polytope must contain the saved local pattern.

Before an infeasible linear program could become a theorem, one must prove that **every** compact 3-free polytope induces a feasible point of that program. In particular, its charged quadrilaterals must satisfy the literal spherical-face and attachment premises of each imported obstruction. A fixed Gale nerve or a source-specific configuration supplies no such global occurrence theorem. No such complete linear program or occurrence theorem is claimed here.

## The unrestricted finite-volume lane

The parallel investigation recovered Prokhorov's primary weighted-face inequalities, including the treatment of nonsimple ideal vertices. Keeping an exact additive deficit in its spherical-diagram pair count gives an independently checked [refinement to 990](dimension_bound_prokhorov_20261015.md). This is separate from the compact inequalities above and is not called a best-known bound in view of Nikulin's below-300 announcement. The next broad target is to reconstruct that stronger narrow-face argument with complete affine-cusp coverage and explicit constants. Compact simplicity, (1), and (6) cannot simply be transferred to arbitrary cusped polytopes.

Arithmetic, compact, right-angled, flag, and facet-count hypotheses are legitimate separate targets. Any reduction between them must be proved. A nonreflective lattice, a failed cap family, or a fixed-source obstruction does not exclude a hyperbolic dimension.

## Verification and scope

The checker performs 2,813 exact partition dynamic-programming comparisons, 100 first-excluded odd-boundary checks, rational identity checks, and the saved class/distribution tables. Its universal conclusions rely on the displayed geometric and algebraic proofs, not on testing only dimensions up to 100. It makes no census of Coxeter diagrams or polytopes.

```sh
python scripts/dimension_bound_vertex_components_20261015.py --verify
```

The source and JSON are frozen; generation and saved replay pass. The [independent partition arithmetic](dimension_bound_vertex_components_20261015_review.md) and the [separate primary-source/proof audit](dimension_bound_prokhorov_20261015_proof_review.md) agree with the displayed arguments. The latter explicitly checks the full flag hypothesis and direct-diagram-edge charge. Neither review upgrades these conditional compact results to an unrestricted dimension theorem.
