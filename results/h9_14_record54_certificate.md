# Solver-independent exclusion of H9/14 record 54

> **Current status.** Later record-23, record-74, record-46, record-37, and record-34
> certificates,
> together with the multiwedge realizations of records 43, 44, 47, and 35,
> give the current partition \(121+6+15\); see
> `results/h9_14_facet_completion.md`.

## Result

**Record 54 is independently impossible.** Its graph6 string is

```text
M?????????O?W?W??
```

and its missing-edge graph is

\[
H=K_2[\{0,11\}]\sqcup
K_{2,2}[\{1,2\},\{12,13\}]\sqcup8K_1.
\]

Thus its missing edges are

\[
\{0,11\},\quad
\{1,12\},\quad\{1,13\},\quad
\{2,12\},\quad\{2,13\}.
\]

The proof uses integer face arithmetic, strict vertex-link \(g_4\), blocker
duality, and a finite four-row incidence classification. It neither imports
nor trusts a SAT/SMT backend. The archived Z3 result is only corroboration.

**Historical stage snapshot.** Together with the independent certificates
for records 92 and 93, this certificate promoted the theorem-facing partition
of the original 142 graph shadows to

\[
116\text{ independently impossible}
\; +\;2\text{ explicit ordinary-polytopal selections}
\; +\;24\text{ theorem-undecided}.
\]

This excludes one graph shadow. It is not the complete compact H9/14
non-existence theorem.

## 1. Strict vertex-link \(g_4\) leaves only two schedules

Let \(K\) be the boundary of the hypothetical simplicial 9-polytope, and let
\(a\) be its number of missing faces of order four. The four archived
Dehn--Sommerville schedules are as follows; \(f_j\) counts faces of dimension
\(j\).

| \(a\) | \(f_3\) | \(f_4\) | \(f_8\) | \(\sum_v g_4(\operatorname{lk}_K v)\) |
|---:|---:|---:|---:|---:|
| 0 | 718 | 1,144 | 96 | 24 |
| 1 | 717 | 1,138 | 94 | 18 |
| 2 | 716 | 1,132 | 92 | 12 |
| 3 | 715 | 1,126 | 90 | 6 |

Every vertex link is the boundary of a nonsimplex simplicial 8-polytope. A
simplex link is impossible because its boundary has a missing face of order
nine, beyond the Coxeter--Lanner bound five. The generalized lower-bound
equality theorem similarly makes \(g_4=0\) impossible: a nontrivial
3-stacked link would have a missing face beyond that bound. Consequently

\[
g_4(\operatorname{lk}_K v)\geq1
\qquad\text{for every one of the fourteen vertices.}
\]

For an 8-polytope boundary,

\[
g_4=f_3-6f_2+21f_1-56f_0+126.
\]

Double-counting faces containing a vertex therefore gives

\[
\begin{aligned}
\sum_v g_4(\operatorname{lk}_K v)
 &=5f_4-24f_3+63f_2-112f_1+14\cdot126\\
 &=24-6a.
\end{aligned}
\]

The left side is at least 14, so

\[
\boxed{a\in\{0,1\}}.
\]

The rest of the proof excludes those two branches.

## 2. The common blocker test

Put

\[
U=\{0,1,2,11,12,13\},
\qquad I=\{3,4,\ldots,10\}.
\]

The missing-edge graph has exactly four inclusion-minimal vertex covers:

\[
\begin{array}{ll}
C_{0L}=\{0,1,2\},&C_{0R}=\{0,12,13\},\\
C_{1L}=\{11,1,2\},&C_{1R}=\{11,12,13\}.
\end{array}
\]

Let \(\mathcal M\) be the clutter of all minimal nonfaces. Complements of
maximal faces are precisely the inclusion-minimal transversals of
\(\mathcal M\). Since the target is pure of dimension eight, all such
transversals must have order five, and their number must equal \(f_8\).

Every transversal of the five missing edges contains one of the four covers
\(C\) above. Given the additional minimal nonfaces \(\mathcal F\), define

\[
\mathcal D_C=\{F\in\mathcal F:C\cap F=\varnothing\}.
\]

There is a transversal of order at most four extending \(C\) exactly when
either \(\mathcal D_C\) is empty or its members have a common vertex. Thus a
necessary and sufficient condition for avoiding this immediate purity
obstruction is

\[
\mathcal D_C\ne\varnothing,
\qquad
\bigcap_{F\in\mathcal D_C}F=\varnothing
\quad\text{for all four }C.
\tag{2.1}
\]

Once (2.1) holds, every five-set transversal is inclusion-minimal. Hence if
their number differs from the scheduled \(f_8\), then either a larger minimal
transversal exists, contradicting purity, or all minimal transversals have
order five and the facet count is wrong. In either case the branch is
impossible.

## 3. Exact local \(g_4\) capacities

First pretend that every independent set through order five is a face. The
numbers of such faces of orders two through five containing a fixed vertex,
and the resulting baseline vertex-link \(g_4\), are

| vertex type | orders 2, 3, 4, 5 | baseline \(\beta_v\) |
|---|---|---:|
| \(v=0,11\) | \((12,62,184,350)\) | 2 |
| \(v=1,2,12,13\) | \((11,54,156,294)\) | 2 |
| \(v\in I\) | \((13,73,235,483)\) | 4 |

Thus

\[
\beta_v=
\begin{cases}
2,&v\in U,\\
4,&v\in I.
\end{cases}
\tag{3.1}
\]

A missing face of order five containing \(v\) subtracts one from the local
\(g_4\). A missing face \(A\) of order four containing \(v\) adds six, while
each independent five-set containing \(A\) and \(v\) subtracts one. These
observations give all capacities used below directly from the face formula;
no local-link solver is involved.

## 4. The \(a=0\) branch: an exact four-row classification

There are 1,148 independent five-sets in \(H\), whereas the schedule requires
\(f_4=1,144\). With no missing four-face, exactly four independent
five-sets are minimal nonfaces. Call them

\[
B_0,B_1,B_2,B_3.
\]

Let \(d_v=|\{i:v\in B_i\}|\). Equation (3.1) gives

\[
g_4(\operatorname{lk}_K v)=\beta_v-d_v.
\]

Strict positivity therefore forces

\[
d_v\leq1\quad(v\in U),
\qquad
d_v\leq3\quad(v\in I).
\tag{4.1}
\]

In particular, the four active parts

\[
A_i=B_i\cap U
\]

are pairwise disjoint, and no isolate lies in all four \(B_i\).

### 4.1 Completeness of the compression

An independent subset of \(U\) has at most three elements. There are exactly
21 possible active parts, distributed by size as

\[
1,\ 6,\ 10,\ 4
\qquad\text{in sizes }0,1,2,3.
\]

For a cover \(C\), put

\[
M_C=\{i:A_i\cap C=\varnothing\}.
\]

Because distinct nonempty active parts are disjoint, condition (2.1) first
requires \(|M_C|\geq2\). For each isolate \(x\), encode its row membership by

\[
S(x)=\{i:x\in B_i\}\subseteq\{0,1,2,3\}.
\]

Then the empty-intersection part of (2.1) says precisely that no column
\(S(x)\) may contain any \(M_C\).

Consequently every admissible family, and only such a family, is encoded by
an unordered pairwise-disjoint quadruple \((A_i)\), followed by nonnegative
integers \(n_S\), \(S\subseteq[4]\), satisfying

\[
\begin{aligned}
&\sum_{S\subseteq[4]}n_S=8,\\
&\sum_{S\ni i}n_S=5-|A_i|\qquad(0\leq i<4),\\
&n_S=0\qquad\text{if }M_C\subseteq S
   \text{ for at least one minimal graph cover }C.
\end{aligned}
\tag{4.2}
\]

The last line includes \(n_{[4]}=0\), as required by (4.1). When two active
parts are both empty, equality of their isolate rows is explicitly excluded,
because the four minimal nonfaces must be distinct.

Independent enumeration of all \(21^4\) ordered active quadruples, followed
by sorting, reconfirms that exactly 235 active-part multisets obey the
pairwise-disjoint and \(|M_C|\geq2\) conditions. Enumerating the integer
solutions of (4.2) leaves 1,422 row-incidence states.

This compression loses no labelled family: expand a solution into the eight
isolate columns, assign the eight labels in all
\(8!/\prod_S n_S!\) ways, and quotient permutations of rows having the same
empty active part. The resulting labelled unordered-family count is
2,229,164.

### 4.2 The 159-orbit canonicalization

There is also a much smaller independently checkable quotient. Represent a
state by its ordered active masks \((A_0,A_1,A_2,A_3)\) and the multiset of
its eight four-bit isolate columns.

Act on this representation by

\[
\operatorname{Aut}(H|_U)\times S_4,
\qquad
|\operatorname{Aut}(H|_U)|
=|\operatorname{Aut}(K_2)|\,|\operatorname{Aut}(K_{2,2})|
=2\cdot8=16.
\]

The first factor acts on every active mask. If \(\pi\in S_4\) reorders the
rows, a column \(s\) is sent to the column \(s'\) whose bit \(j\) is bit
\(\pi(j)\) of \(s\). The \(S_8\)-action on the isolates has already been
quotiented by sorting the eight transformed column values. The canonical key
is the lexicographically least pair

\[
\bigl((A'_0,A'_1,A'_2,A'_3),
      \operatorname{sort}(s'_1,\ldots,s'_8)\bigr)
\]

over all \(16\cdot24\) transformations. There are exactly

\[
\boxed{159}
\]

distinct canonical keys.

### 4.3 None has the required 96 facets

The missing-edge graph has 198 five-element vertex covers. For each family
the exact number relevant to blocker duality is

\[
N(B)=\left|\left\{
T\in\binom{[14]}5:
T\cap e\ne\varnothing\ (e\in E(H)),\quad
T\cap B_i\ne\varnothing\ (0\leq i<4)
\right\}\right|.
\tag{4.3}
\]

The complete labelled histogram, together with its 159-orbit quotient, is

| \(N(B)\) | labelled unordered families | canonical orbits |
|---:|---:|---:|
| 60 | 25,200 | 2 |
| 62 | 26,880 | 2 |
| 64 | 75,600 | 8 |
| 65 | 13,440 | 1 |
| 66 | 174,720 | 6 |
| 67 | 16,800 | 2 |
| 68 | 414,120 | 29 |
| 69 | 107,520 | 2 |
| 70 | 449,400 | 21 |
| 72 | 287,420 | 21 |
| 74 | 248,640 | 9 |
| 76 | 270,340 | 27 |
| 82 | 3,360 | 1 |
| 84 | 59,360 | 12 |
| 86 | 25,760 | 3 |
| 88 | 20,160 | 2 |
| 90 | 4,200 | 6 |
| 92 | 1,120 | 1 |
| 95 | 224 | 1 |
| 98 | 4,900 | 3 |
| **total** | **2,229,164** | **159** |

The schedule requires \(f_8=96\), which does not occur. Families excluded
before the table have a transversal of order at most four; every family in
the table has the wrong number of order-five transversals. This proves the
\(a=0\) branch impossible.

## 5. The \(a=1\) branch collapses to 840 isolate triples

Let \(A\) be the unique missing four-face, and define

\[
c(A)=|\{w\notin A:A\cup\{w\}\text{ is independent in }H\}|.
\]

The schedule has ten independent nonfaces of order five. Exactly \(c(A)\)
of them contain \(A\), so the number \(b\) of minimal missing five-faces is

\[
b=10-c(A).
\tag{5.1}
\]

For a vertex \(v\), let

\[
s_A(v)=
\begin{cases}
c(A),&v\in A,\\
1,&v\notin A\text{ and }A\cup\{v\}\text{ is independent},\\
0,&\text{otherwise},
\end{cases}
\]

and let \(d_v\) count the minimal missing five-faces containing \(v\). The
exact local equation is

\[
g_4(\operatorname{lk}_K v)
=\beta_v+6\,[v\in A]-s_A(v)-d_v.
\tag{5.2}
\]

Write \(P\) for the \(K_2\) component of \(H|_U\), and write \(Q_1\) for
one vertex and \(Q_2\) for one complete side of the \(K_{2,2}\) component.
Up to \(\operatorname{Aut}(H)\), the 718 independent choices for \(A\) form
exactly the following six types.

| \(A\cap U\) | \(|A\cap I|\) | labelled choices | \(c(A)\) | \(b\) | strict-link result |
|---|---:|---:|---:|---:|---|
| \(\varnothing\) | 4 | 70 | 10 | 0 | an isolate of \(A\) has pre-\(d\) \(g_4=0\) |
| one endpoint of \(P\) | 3 | 112 | 9 | 1 | that endpoint has pre-\(d\) \(g_4=-1\) |
| \(Q_1\) | 3 | 224 | 8 | 2 | that vertex has pre-\(d\) \(g_4=0\) |
| \(Q_2\) | 2 | 56 | 8 | 2 | those vertices have pre-\(d\) \(g_4=0\) |
| \(P_1\cup Q_1\) | 2 | 224 | 7 | 3 | survives |
| \(P_1\cup Q_2\) | 1 | 32 | 7 | 3 | survives |
| **total** |  | **718** |  |  |  |

Thus only the last two rows remain. Normalize the chosen endpoint and chosen
side so that

\[
A=\{0,1,3,4\}
\quad\text{or}\quad
A=\{0,1,2,3\}.
\]

In either case (5.2) gives capacity zero on \(\{0,1,2\}\), capacity one on
the opposite minimal graph cover

\[
D=\{11,12,13\},
\]

and capacity two on every isolate. Hence the three minimal missing
five-faces may initially use vertices of \(D\), each at most once.

The blocker condition removes those apparent possibilities. Since
\(C_{0R}=\{0,12,13\}\) already meets \(A\), any use of 12 or 13 either leaves
only one of the three missing five-faces disjoint from \(C_{0R}\), or leaves
two. In the latter case the two missed faces avoid 12 and 13; at most one can
use 11, so their isolate parts have sizes at least four and five. They must
intersect inside the eight-set \(I\), violating (2.1). Therefore 12 and 13
occur in none of the three faces.

Now use \(C_{1L}=\{11,1,2\}\), which also meets \(A\). If one missing
five-face used 11, the other two would be five-subsets of \(I\), and hence
would intersect. Again (2.1) fails. Thus 11 cannot occur either.

Consequently the three minimal missing five-faces are, necessarily,

\[
B_1,B_2,B_3\in\binom I5,
\qquad B_1\cap B_2\cap B_3=\varnothing.
\tag{5.3}
\]

For example,

\[
\begin{aligned}
B_1&=\{3,4,5,6,7\},\\
B_2&=\{3,4,5,8,9\},\\
B_3&=\{6,7,8,9,10\}
\end{aligned}
\]

passes this reduced blocker test.

Put \(R_i=I\setminus B_i\). The three \(R_i\) are 3-sets with union \(I\).
Their total incidence is nine on eight points, so exactly one point belongs
to two of the \(R_i\), and every other point belongs to one. The three
families \(\binom{R_i}{2}\) are disjoint. There are exactly

\[
8\binom73\cdot3=840
\tag{5.4}
\]

unordered triples (5.3) for each fixed \(A\).

An isolate pair hits all three \(B_i\) exactly when it is contained in none
of the \(R_i\). There are therefore

\[
\binom82-3\binom32=19
\tag{5.5}
\]

allowed pairs. A five-element graph cover with four or five active vertices
has at most one isolate and cannot hit all three \(B_i\), because their
triple intersection is empty. Hence every five-transversal consists of one
of the four minimal active graph covers and one of these 19 pairs.

Three active covers already meet the active part of \(A\). The opposite
cover \(D\) meets \(A\) exactly when its isolate pair meets \(J=A\cap I\).
If \(L\) is the 19-edge graph of allowed isolate pairs, the facet count is
therefore

\[
N=3\cdot19+|\{e\in E(L):e\cap J\ne\varnothing\}|.
\tag{5.6}
\]

The point repeated among the \(R_i\) has degree three in \(L\); each other
point has degree five. Formula (5.6) gives the complete per-\(A\) histogram:

| \(|J|\) | configuration | \(N\) | triples per fixed \(A\) |
|---:|---|---:|---:|
| 1 | \(J\) is the repeated point | 60 | 105 |
| 1 | \(J\) is not the repeated point | 62 | 735 |
| 2 | repeated point in \(J\), the pair in \(L\) | 64 | 90 |
| 2 | repeated point in \(J\), the pair outside \(L\) | 65 | 120 |
| 2 | repeated point outside \(J\), the pair in \(L\) | 66 | 480 |
| 2 | repeated point outside \(J\), the pair outside \(L\) | 67 | 150 |

The first two rows sum to 840 for each of the 32 choices with \(|J|=1\); the
last four sum to 840 for each of the 224 choices with \(|J|=2\). Thus this
branch has 215,040 reduced labelled families in total, and none has the
scheduled \(f_8=94\).

This proves the \(a=1\) branch impossible.

## 6. Exact scope and reproducibility

Both Dehn--Sommerville branches surviving the summed-link bound are now
excluded independently. The proof only needs:

- decoding the five-edge graph;
- integer \(g_4\) and face-incidence arithmetic;
- the four minimal graph covers;
- the explicit finite system (4.2); and
- direct enumeration of five-element transversals.

In particular, it does not depend on the archived Z3 `unsat` result for
record 54.

The accompanying verifier reconstructs the source schedule, all 235 active
configurations, the 1,422 incidence states, the 159 canonical orbits, both
exact histograms, and the six \(a=1\) types. Reproduce with

```text
PYTHONPATH=src python3 scripts/certify_h9_14_record54.py --verify-only
PYTHONPATH=src pytest -q tests/test_h9_14_record54_certificate.py
```

The certificate establishes a combinatorial graph-shadow obstruction. It
makes no claim about the remaining 24 shadows and no general claim that all
ordinary simplicial 9-polytopes with fourteen vertices are excluded.
