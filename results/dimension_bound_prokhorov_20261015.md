# An exact refinement of the Prokhorov weighting estimate

The calculation below derives the upper bound **990** for finite-sided finite-volume hyperbolic Coxeter polytopes by keeping an additive term discarded in the classical weighted-diagram estimate. It changes neither the assigned weights nor the geometric hypotheses. The complete exact graph/scalar certificate and separate arithmetic/source audits pass. This is an explicitly checked refinement of Prokhorov's proof, **not a claim to the best known bound or publication priority**.

Historical qualification: [Nikulin's ICM 1986 report, printed pp. 654 and 659](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf), defines its class to include all finite-volume reflection groups, then announces $\dim\mathbb H<300$. Its outline invokes a narrow face, combinatorial inequalities, and a cusp estimate. The displayed example treats the affine-$A$ cusp and gives $n<140$; the report does not supply the full numerical case analysis establishing 300. We neither dismiss that announcement nor silently promote our weaker 990 result above it. Recovering a complete auditable proof of the announced stronger bound is the next broad-dimension target.

This is an unrestricted dimension argument: it imposes no arithmeticity, facet count, fixed lattice, prescribed source, simplicity at ideal vertices, or bound on finite Coxeter labels. Compact polytopes already satisfy the stronger bound 29. The [broader dimension-lane report](dimension_bound_lane_20261015.md) records the established literature and separately scoped compact-class inequalities.

## Published geometric input, unchanged

Use [Prokhorov, Math. USSR-Izv. 28 (1987), 401–411](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&option_lang=eng&paperid=1484&what=fullteng), printed pp. 405–410, Lemmas 1–2 and Proposition 5. For each edge $e$, its spherical diagram $T_e$ has rank $n-1$. A three-face containing $e$ selects an unordered node pair. The assigned weight is

\[
w(i,j)=\begin{cases}1&1\le d(i,j)\le7,\\1/3&8\le d(i,j)\le15,\\0&d(i,j)>15\text{ or disconnected}.\end{cases}
\]

Proposition 5 supplies nonnegative weights whose sum on every good three-face with $k$ facets is at least $7-k$. “Bad” means triangular bipyramid. Lemma 2 states that, if the sum at each edge is at most $C(n-1)$, with $C>0$, then

\[
n<96C+68. \tag{1}
\]

These statements include cusps: finite vertices are simple, ideal links are products of simplices, and Khovanskii's edge-simple inequalities apply. We retain those published results without changing the good/bad-face classification or its coverage. The original estimate takes $C=29/3$, obtaining $n<996$.

The landing-page abstract loses an inequality symbol; the paper's Theorem C explicitly excludes dimensions greater than 995. Metadata is not evidence for a different bound.

## Exact weighted pair maximum

For a spherical diagram $T$ of rank $N$, let $f_k(T)$ count unordered pairs of distinct nodes at graph distance at most $k$. The sum of the unchanged weights is

\[
W(T)=\frac{2f_7(T)+f_{15}(T)}3.
\]

We prove

\[
\boxed{W(T)\le\frac{29N-173}{3}\quad(N\ge17).} \tag{2}
\]

The finite Coxeter classification says that an irreducible spherical diagram, after forgetting edge labels, is a path, $D_r$, or one of $E_6,E_7,E_8$. In particular, the unbounded family $I_2(m)$ is represented by the same two-node path for every finite $m\ge3$. This is not an angle cutoff.

For a path with $r\ge17$ nodes,

\[
f_k=kr-\frac{k(k+1)}2\quad(k=7,15),
\qquad 29r-3W=176.
\]

For $D_r$, view its tree as a center with arms of lengths $1,1,r-3$. Counting center–arm, same-arm, and different-arm pairs gives

\[
f_k=kr-\frac{k(k+1)}2+1\quad(k=7,15),
\qquad 29r-3W=173.
\]

The condition $r\ge17$ includes both distances and is essential at the short end of the formula. There are only finitely many smaller unlabelled graphs. The certificate counts each of them and verifies

\[
17\bigl(29r-3W(T)\bigr)\ge173r\quad(1\le r\le16). \tag{3}
\]

For a disconnected diagram, weights add across components and no cross-component pair contributes. If there is a component of rank at least 17, it already contributes deficit at least 173, and every other deficit is positive. If all components have smaller rank, summing (3) gives deficit at least $173N/17\ge173$. This proves (2). Equality is attained by the graph $D_N$; its geometric occurrence as an edge diagram is not needed for the upper bound.

There is a second proof adequate for the dimension conclusion that avoids even the finite small-graph count. For any component of rank $r\le16$, simply use $W\le\binom r2$, so

\[
29r-3W\ge\frac{r(61-3r)}2\ge\frac{13r}2.
\]

If all components are small and $N\ge27$, their total deficit exceeds 173. If a component is large, its explicit path or $D$-formula supplies 173 by itself. Thus (2) for $N\ge27$, which is more than sufficient below, requires no small-rank enumeration.

## Substitution into the published geometric lemma

Suppose a polytope existed in dimension $n\ge991$. Apply (2) to each edge diagram with $N=n-1$, and set

\[
C_n=\frac{29}{3}-\frac{173}{3(n-1)}>0.
\]

For this fixed dimension, $C_n$ is a single constant valid at every edge. It is therefore a valid input to (1); the lemma does not require one constant to be used simultaneously in all dimensions. We obtain

\[
n<96C_n+68=996-\frac{5536}{n-1},
\]

equivalently

\[
n^2-997n+6532<0. \tag{4}
\]

But for $n=991+t$, $t\ge0$, the left side is

\[
t^2+985t+586>0,
\]

a contradiction. Thus the derived upper bound is $n\le990$. At $n=990$, the polynomial is $-398$; this argument does not exclude dimension 990.

## Independent incidence-algebra check

For additional error detection, the checker reconstructs the calculation before the simplification (1). Write $A_{31},A_{32}$ for the strict Nikulin–Khovanskii upper bounds for mean edges and mean two-faces of a three-face, and $p$ for the strict bad-face proportion bound. Their parity forms are

\[
\begin{array}{c|cc}
&n\text{ even}&n\text{ odd}\\\hline
A_{31}&12(n-1)/(n-4)&12n/(n-3)\\
A_{32}&6(n-2)/(n-4)&6(n-1)/(n-3)\\
p&(3n+6)/(4(n-3))&(3n+9)/(4(n-2)).
\end{array}
\]

If $S$ is the total weight and $f_3$ the number of three-faces, the incidence identities give

\[
7-A_{32}-p<\frac{S}{f_3}
<\frac{2W_{n-1}A_{31}}{(n-1)(n-2)},
\qquad W_{n-1}=\frac{29n-202}{3}.
\]

The left lower bound minus the right upper bound simplifies to

\[
\begin{cases}
\dfrac{n^3-1000n^2+9604n-19824}{4(n-4)(n-3)(n-2)},&n\text{ even},\\[5pt]
\dfrac{n^3-997n^2+6679n-147}{4(n-3)(n-2)(n-1)},&n\text{ odd}.
\end{cases}
\]

Both are positive in the excluded range. Substitution $n=992+t$ into the even numerator gives coefficients $(1634832,977796,1976,1)$ in increasing powers of $t$; substitution $n=991+t$ into the odd numerator gives $(726256,976868,1976,1)$. These positive coefficients prove the entire tails. This agrees with the simpler argument (4) and checks the distinction between $n$ and edge rank $n-1$.

## Coverage, verification, and remaining work

The [primary checker](../scripts/dimension_bound_prokhorov_20261015.py) uses only integer arithmetic and `Fraction`. The [JSON certificate](dimension_bound_prokhorov_20261015.json) records all small connected graph counts, all-rank path/$D$ formulas, the alternate rank-27 proof, exact rank-partition checks through 100, and both scalar calculations. Its finite loops check formulas; the all-rank arguments above establish their unrestricted scope.

The imported geometric facts are the published weighting construction and Lemma 2, rather than a fresh computer verification of all quasi-Lannér configurations. Independent review must check that import, including nonsimple ideal vertices, the unordered-pair correspondence, and weights at distances 7 and 15. Neither a numerical sample nor a fixed-source obstruction is involved. The graph estimate itself is label independent and covers disconnected diagrams.

Generation, saved replay, and compilation pass. The [separate graph/scalar replay](dimension_bound_prokhorov_20261015_review.json) uses graph BFS and independent rational calculations without primary imports. The [independent source/proof audit](dimension_bound_prokhorov_20261015_proof_review.md) confirms the unchanged edge and cusp hypotheses and records the historical qualification above. These checks do not re-prove the published quasi-Lannér classification.

The immediate next broad target is to reconstruct Nikulin's narrow-face argument, including all affine cusp-component partitions and the compact narrow-face alternative. Its functions and constants must be made explicit before any new numerical conclusion. A second option is genuinely geometric: change the distance weights or improve the bad-face estimate while proving every good-face lower bound remains valid. Merely optimizing an unjustified finite catalogue would not yield another dimension theorem.

```sh
python scripts/dimension_bound_prokhorov_20261015.py --verify
```
