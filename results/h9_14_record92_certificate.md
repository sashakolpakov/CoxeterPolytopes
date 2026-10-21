# Solver-independent exclusion of H9/14 record 92

> **Current status.** Later record-93/54/23/74/46/37/34 certificates, together with
> the multiwedge realizations of records 43, 44, 47, and 35, give the current
> partition \(121+6+15\); see `results/h9_14_facet_completion.md`.

## Result

**Record 92 is independently impossible.**  Its graph6 string is

```text
M?????????W?o?w??
```

and its missing edges are

\[
\{0,11\},\{1,11\},\{0,12\},\{1,12\},
\{0,13\},\{1,13\},\{2,13\}.
\]

Equivalently, the graph is \(K_{2,3}\), on parts \(\{0,1\}\) and
\(\{11,12,13\}\), together with the edge \(\{2,13\}\).

The proof is solver-independent.  It uses a vertex-link \(g_4\) double count,
Dehn--Sommerville, and blocker duality for the minimal-nonface clutter.  The
archived Z3 result is only corroboration.

This promotes record 92 from theorem-undecided to certified impossible.
**Historical stage snapshot:** the partition of the 142 graph shadows at the
time became

\[
114\text{ independently impossible}
\; +\;2\text{ explicit ordinary-polytopal selections}
\; +\;26\text{ undecided}.
\]

The sibling record-93 certificate subsequently advanced this to
\(115+2+25\), and the then-later record-54 certificate to \(116+2+24\).
These are proof-history snapshots; the current status is stated above.

It is one graph-shadow exclusion, not the complete H9/14 non-existence
theorem.

## 1. The positive-\(a\) schedules are impossible

Let \(K\) be the boundary of the hypothetical simplicial 9-polytope and let
\(a\) be its number of missing four-faces.  The three archived schedules have
\(a=0,1,2\).

For every vertex \(v\), its link is the boundary of a nonsimplex simplicial
8-polytope.  Every missing face of a link is controlled by a global missing
face, so the Coxeter--Lanner order bound five is inherited by the link.  If
\(g_4(\operatorname{lk}v)=0\), the generalized lower-bound equality theorem
would give a nontrivial 3-stacked triangulation and hence a missing link face
of order at least six.  Therefore

\[
g_4(\operatorname{lk}v)\geq1
\qquad\text{for all fourteen vertices.}
\]

For an 8-polytope boundary,

\[
g_4=f_3-6f_2+21f_1-56f_0+126.
\]

Incidence double-counting over all vertex links gives

\[
\begin{aligned}
\sum_v g_4(\operatorname{lk}v)
 &=5f_4(K)-24f_3(K)+63f_2(K)-112f_1(K)+14\cdot126\\
 &=18-6a.
\end{aligned}
\]

The left side is at least 14.  Thus \(18-6a\geq14\), and the integer
\(a\in\{0,1,2\}\) is forced to be

\[
\boxed{a=0}.
\]

In particular, the \(a=1\) and \(a=2\) branches have sums 12 and 6 and are
independently contradictory.

## 2. The remaining minimal-nonface clutter

The graph has 1,008 independent five-sets.  The \(a=0\) schedule requires
\(f_4=1005\), so exactly three independent five-sets are missing faces.  Call
them \(B_1,B_2,B_3\).  Because there are no missing four-faces, they are
minimal.  The Lanner order bound excludes any further minimal nonface of
order greater than five.  Hence the entire minimal-nonface clutter is

\[
\mathcal M=E(H)\cup\{B_1,B_2,B_3\}.
\]

The complements of maximal faces of a simplicial complex are precisely the
inclusion-minimal transversals of its minimal-nonface clutter.  A pure
eight-dimensional complex on 14 vertices must therefore have all minimal
transversals of order five.  The \(a=0\) schedule additionally requires

\[
f_8=78,
\]

so it needs exactly 78 such five-transversals.

## 3. Hit-pattern classification

Put \(I=\{3,\ldots,10\}\).  The three inclusion-minimal vertex covers of
the missing-edge graph are

\[
C_1=\{0,1,2\},\qquad
C_2=\{0,1,13\},\qquad
C_3=\{11,12,13\}.
\]

For an independent five-set \(B\), define

\[
p(B)=\{i:B\cap C_i\neq\varnothing\}.
\]

Direct independence gives exactly the following six types.

| \(p(B)\) | structure | count |
|---|---|---:|
| \(\varnothing\) | \(B\in\binom I5\) | 56 |
| \(\{1\}\) | \(B=\{2\}\cup Q,\ Q\in\binom I4\) | 70 |
| \(\{1,2\}\) | \(B\subset\{0,1,2\}\cup I\), meeting \(\{0,1\}\) | 336 |
| \(\{3\}\) | \(B\subset\{11,12\}\cup I\), meeting \(\{11,12\}\) | 196 |
| \(\{1,3\}\) | \(B=\{2\}\cup Q\), where \(Q\) meets \(\{11,12\}\) | 140 |
| \(\{2,3\}\) | \(B=\{13\}\cup Q,\ Q\in\binom{\{11,12\}\cup I}4\) | 210 |

Types \(\{2\}\) and \(\{1,2,3\}\) are forbidden by independence.

For each \(i\), let

\[
D_i=\{B_j:B_j\cap C_i=\varnothing\}.
\]

If \(D_i\) is empty, then \(C_i\) is already a transversal of order three.
If the members of \(D_i\) have a common vertex \(x\), then
\(C_i\cup\{x\}\) is a transversal of order at most four.  Thus absence of
a small transversal forces every \(D_i\) to have at least two members with
empty intersection.  In particular, the three hit patterns are pairwise
disjoint.

The compatible pattern multisets are

\[
(0,0,0),\quad (0,0,p),\quad
(0,\{1\},\{3\}),\quad
(0,\{1\},\{2,3\}),\quad
(0,\{1,2\},\{3\}),
\]

with \(p\) any of the five nonzero allowed patterns.  The middle cases fail:

- in \((0,0,p)\), the two zero-pattern five-subsets of the eight-set \(I\)
  cannot be disjoint;
- in the next two cases, a zero-pattern \(5I\) set cannot be disjoint from
  the four \(I\)-vertices in a type-\(\{1\}\) set.

Only two families remain.

## 4. The two surviving families have the wrong facet count

### Family I: pattern \((0,0,0)\)

All three \(B_j\) are five-subsets of \(I\), and the no-small-transversal
condition is

\[
B_1\cap B_2\cap B_3=\varnothing.
\]

Let \(R_j=I\setminus B_j\).  The three \(R_j\) are 3-sets whose union is
\(I\).  Their total incidence is nine over eight points, so exactly one point
occurs twice.  Consequently their three families of 2-subsets are disjoint.

There are

\[
8\binom73\cdot3=840
\]

unordered families: choose the repeated point, the 3-set not containing it,
and one of the three pairings of the remaining four points.

A graph-cover five-set meeting all three \(B_j\) must consist of one of the
three minimal active covers \(C_i\) and an isolate pair.  Of the
\(\binom82=28\) isolate pairs, exactly
\(3\binom32=9\) lie in some \(R_j\).  Thus every family has

\[
3(28-9)=\boxed{57}
\]

five-transversals.

### Family II: pattern \((0,\{1,2\},\{3\})\)

Choose \(B_0\in\binom I5\) and put \(R=I\setminus B_0\).  The two required
empty intersections force

\[
B_3=R\cup\{11,12\},
\]

while \(B_{12}\) is any five of the six elements
\(R\cup\{0,1,2\}\).  Hence there are

\[
\binom85\cdot6=336
\]

families.

The missing-edge graph has respectively 3, 8, and 6 active vertex covers of
orders 3, 4, and 5.  It therefore has

\[
3\binom82+8\binom81+6=154
\]

five-element vertex covers.  Among these, the numbers disjoint from
\(B_0,B_{12},B_3\) are respectively

\[
39=3\binom32+8\cdot3+6,
\qquad
15=\binom62,
\qquad
25=2\binom52+5.
\]

These three forbidden classes are pairwise disjoint.  Thus every family has

\[
154-39-15-25=\boxed{75}
\]

five-transversals.

Neither 57 nor 75 is the required 78.  If a wrong-sized minimal transversal
exists, purity already fails; if none exists, the facet count fails.  This
exhausts the \(a=0\) branch and proves record 92 impossible.

## 5. Independent verification and backend corroboration

The pure-Python verifier independently reconstructs:

- the missing-edge graph and its three minimal covers;
- the three Dehn--Sommerville schedules and the link-\(g_4\) sums
  \((18,12,6)\);
- all 1,008 independent five-sets and the six hit-pattern classes;
- all 6,773,256 pairwise-pattern-compatible unordered triples;
- the 1,176 triples with no transversal of order at most four; and
- the exact histogram \(\{57:840,75:336\}\).

The verifier uses no optional solver.  The earlier schedule-split Z3 archive
returned UNSAT for \(a=0,1,2\); the new manifest records that only as untrusted
corroboration.

Artifacts:

- `src/coxeter_polytopes/combinatorics/h9_14_record92.py`;
- `scripts/certify_h9_14_record92.py`;
- `data/enumerations/h9_14_record92_certificate.json`; and
- `tests/test_h9_14_record92_certificate.py`.

Reproduce with

```text
PYTHONPATH=src python3 scripts/certify_h9_14_record92.py --verify-only
PYTHONPATH=src pytest -q tests/test_h9_14_record92_certificate.py
```
