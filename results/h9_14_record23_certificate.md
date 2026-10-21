# Solver-independent exclusion of H9/14 record 23

> **Current status.** Later record-74, record-46, record-37, and record-34 certificates,
> together
> with the multiwedge realizations of records 43, 44, 47, and 35, give the
> current partition \(121+6+15\); see `results/h9_14_facet_completion.md`.

## Result

**Record 23 is independently impossible.**  Its graph6 string is

```text
M???????????w?w??
```

and its missing-edge graph is

\[
H=K_{3,2}[A,B]\sqcup9K_1,
\qquad
A=\{0,1,2\},\quad B=\{12,13\},\quad
I=\{3,\ldots,11\}.
\]

The proof uses strict local `g`-vectors, Macaulay inequalities, blocker
duality, two exact bicoloured-incidence censuses, and a separate compressed
positive-branch enumerator.  It does not trust the archived Z3 `UNSAT`
disposition.

**Historical stage snapshot.** At the time this certificate was added, it
changed the theorem-facing partition of the 142 graph shadows to

\[
117\text{ independently impossible}
\; +\;2\text{ explicit ordinary-polytopal selections}
\; +\;23\text{ theorem-undecided}.
\]

This excludes one graph shadow.  It is not a complete compact H9/14
non-existence theorem.

## 1. Strict vertex-link `g4` leaves four schedules

Let `a` be the number of missing four-faces.  There are 1,128 independent
five-sets in `H`.  The six Dehn--Sommerville schedules are

| `a` | `f3` | `f4` | missing independent five-sets | `f8` | sum over vertices of link `g4` |
|---:|---:|---:|---:|---:|---:|
| 0 | 699 | 1,121 | 7  | 98 | 36 |
| 1 | 698 | 1,115 | 13 | 96 | 30 |
| 2 | 697 | 1,109 | 19 | 94 | 24 |
| 3 | 696 | 1,103 | 25 | 92 | 18 |
| 4 | 695 | 1,097 | 31 | 90 | 12 |
| 5 | 694 | 1,091 | 37 | 88 | 6 |

For an 8-polytope boundary,

\[
g_4=f_3-6f_2+21f_1-56f_0+126.
\]

Double-counting faces through the fourteen vertices gives

\[
\sum_v g_4(\operatorname{lk}v)=36-6a.
\]

Every link is not a simplex, and the equality case in the generalized lower
bound theorem would create an inherited missing face larger than the Lanner
bound.  Hence every summand is at least one.  The schedules `a=4,5` have sum
below 14 and are impossible.

## 2. The missing clutter is forced onto twelve vertices

If every graph-independent four- and five-set is provisionally a face, the
baseline vertex-link values are

\[
\beta_v=
\begin{cases}
5,&v\in A,\\
6,&v\in I,\\
1,&v\in B.
\end{cases}
\tag{2.1}
\]

The exact local missing-four intervals give zero A4 incidence at each vertex
of `B`.  An isolate-only A4 has an independent five-extension through a
vertex of `B`, which subtracts one from its baseline one.  A minimal A5
meeting `B` does the same.  Strict positivity therefore implies:

1. every extra minimal nonface lies in `W=A union I`; and
2. every missing four-face meets `A`.

The two inclusion-minimal covers of `H` are exactly `A` and `B`.

## 3. Blocker duality gives the exact facet formula

Let `G` be the clutter on `W` consisting of all A4s and the additional
minimal A5s.  Let `Z` be the subfamily of A5s contained in `I`.

Every transversal of `H` contains either `A` or `B`.  Thus the absence of a
transversal of order at most four is equivalent to

\[
\bigcap_{F\in Z}F=\varnothing,
\qquad
\tau(G)\geq3.                                               \tag{3.1}
\]

The first condition handles extensions of the three-cover `A`; the second
handles extensions of the two-cover `B`.  Since two five-subsets of the
nine-set `I` always intersect, (3.1) also implies `|Z|>=3`.

Once (3.1) holds, every five-element transversal is minimal.  Such
transversals split into

- `B` plus a three-subset of `W` hitting `G`; and
- `A` plus a two-subset of `I` hitting `Z`.

The would-be overlap `A union B` misses `Z`.  Therefore the exact number of
facet complements is

\[
N=\tau_3(G)+\tau_2(Z).                                      \tag{3.2}
\]

If `N` is not the scheduled `f8`, purity or the facet count fails.

## 4. The `a=0` bicoloured-incidence census

Here `G` consists of seven distinct A5s.  Encode the complements of these
five-sets inside `W` as the rows of a `7 by 12` binary matrix.  Every row has
seven ones.  A pair of columns is disjoint precisely when the corresponding
two vertices hit every A5.  Consequently the no-two-transversal condition in
(3.1) is exactly the condition that every column pair has a common one.

The command

```sh
genbg -q -Y1 -d7:0 -D7:7 7 12 49:49 \
  | /tmp/record23_a0
```

generates bicoloured bipartite graphs with seven row vertices, twelve column
vertices, 49 edges, row degree exactly seven, and at least one common row
neighbor for each column pair.  `genbg` quotients by colour-preserving
isomorphisms, i.e. by row and column permutations.

The checker independently reconstructs and verifies the bipartition, graph6
padding, every row degree, the `-Y1` property, and distinctness of the seven
rows.  It then expands every choice of three columns to serve as `A`.  This
may retain automorphic duplicates, but it omits no `S_3 times S_9`-coloured
object.  Since only the absence of 98 is used, automorphism weights are not
needed.

For a complement column `C_v`, strict local positivity translates to the
following lower bounds:

| original subset type | complement-union lower bound |
|---|---:|
| isolate vertex | `|C_v|>=2` |
| active vertex | `|C_v|>=3` |
| isolate pair | `|C_u union C_v|>=3` |
| pair meeting `A` | `|C_u union C_v|>=4` |
| isolate triple | union size at least 4 |
| triple meeting `A` | union size at least 5 |

For a selected active triple `A`, the rows in the intersection of its three
complement columns are exactly `Z`.  Condition (3.1) says this intersection
is nonempty and every isolate complement column meets it.

An empty intersection of three complement columns counts one member of
`tau3(G)`.  An isolate-column pair whose intersection with the `Z` rows is
empty counts one member of `tau2(Z)`.  Thus the checker evaluates (3.2)
directly.

The frozen run gives

```text
generated_graph_count             2807714
distinct_row_graph_count          2765889
locally_valid_left_choice_count  11456748
target_facet_count_occurrences           0
maximum_facet_count                     94
maximum_three_transversal_count         73
maximum_two_transversal_count           21
```

The complete facet-count histogram is

```text
46:1, 48:5, 49:12, 50:48, 51:182, 52:584, 53:1571,
54:3790, 55:8573, 56:19210, 57:36278, 58:68074,
59:113644, 60:187321, 61:274566, 62:394375, 63:516204,
64:673387, 65:779147, 66:921627, 67:947799, 68:1027350,
69:942371, 70:950004, 71:767830, 72:741564, 73:519701,
74:491360, 75:298732, 76:280029, 77:142750, 78:143051,
79:64749, 80:58532, 81:28013, 82:26302, 83:9434,
84:8271, 85:4674, 86:3137, 87:627, 88:691, 89:319,
90:649, 91:66, 92:54, 93:33, 94:57
```

The required count 98 is absent.

## 5. The `a=1` independent full-local audit

There are three A4 orbits, indexed by `k=|Q intersection A|=1,2,3`.  The A4
shadow contains eight five-sets, so the schedule forces five additional
minimal A5s.  The exact vertex capacities for these five A5s are

| vertex type | capacity |
|---|---:|
| active in `Q` | 2 |
| active outside `Q` | 3 |
| isolate in `Q` | 3 |
| isolate outside `Q` | 4 |

For defense in depth, this branch was checked in two different coordinate
systems.

### Primary full-local census

One degree-eight complement row represents `Q`; five degree-seven rows
represent the additional A5s.  The independent command is

```sh
genbg -q -Y1 -d7:0 -D8:6 6 12 43:43 \
  | /tmp/record23_a1_audit
```

The checker verifies the row profile, distinct A5s, that no A5 contains
`Q`, the complete A4 shadow, every vertex/edge/triple strict local inequality,
the applicable Macaulay inequalities, (3.1), and (3.2).  It then expands all
three-column choices for `A`.

```text
generated bicoloured graphs       75890
structurally admissible graphs    73919
locally valid active choices      86241
target 96 occurrences                 0
maximum N                            90
maximum tau3                         69
maximum tau2                         21
```

Its complete histogram is

```text
54:2, 55:6, 56:11, 57:38, 58:82, 59:225, 60:473,
61:800, 62:1453, 63:2099, 64:3618, 65:4128, 66:5878,
67:6816, 68:7620, 69:7513, 70:7503, 71:7766, 72:6287,
73:5516, 74:5107, 75:3804, 76:2815, 77:2060, 78:1851,
79:789, 80:907, 81:416, 82:367, 83:86, 84:63, 85:63,
86:51, 87:7, 88:3, 89:7, 90:11
```

### Relaxed internal cross-check

The separate self-contained positive-branch enumerator keeps only the
vertex-link `g4` capacities on the chosen A5s, so it enumerates a strict
superset of the primary states.  It canonicalizes the isolate-only subfamily
`Z` and then expands every remaining active A5 before applying the final
blocker test.  For `k=1,2,3` it finds respectively

```text
133726, 508884, 9515916 blocker-valid representative states,
support 56..92, 56..92, 54..92,
target 96 count 0 in all three cases.
```

The relaxed maximum 92 is a useful replay fixture rather than a discrepancy
in the theorem.  For example,

```text
Q  = {0,3,4,5}
A5 = {3,4,5,6,7}, {3,4,5,6,8}, {7,8,9,10,11},
     {0,1,2,3,4}, {0,1,2,5,6}
tau3=71, tau2=21, N=92.
```

It satisfies every relaxed vertex/blocker condition but fails the stronger
local gates at the isolate pair `{3,4}` and at triples containing it: the
corresponding higher local `g` entry is zero.  Thus the primary maximum 90 and
relaxed maximum 92 reconcile exactly.

The required count 96 is absent in both censuses.

## 6. The `a=2,3` compressed census

The self-contained positive-branch program writes an A4 as a nonempty active
part in `A` together with isolate membership columns.  It solves

\[
\sum_Sn_S=9,
\qquad
\sum_{S\ni i}n_S=4-|R_i|,
\]

rejects repeated rows, applies the exact A4 pair/triple local bounds, and
canonicalizes under `S_3(A) times S_a(rows) times S_9(I)`.

Let `p` count pairs of A4s sharing an isolate triple.  The five-shadow size is
`8a-p`; hence the number of additional minimal A5s is

\[
b=(7+6a)-(8a-p)=7-2a+p.                                    \tag{6.1}
\]

### `a=2`

There are 17 exact A4 orbits: 16 with `p=0,b=3`, and one with `p=1,b=4`.
After every eligible A5 completion, strict vertex `g4`, and the final blocker
test, 16,778 representative states remain:

```text
63:120, 65:180, 67:108, 68:540, 69:1260, 70:660,
71:1739, 72:1610, 73:1950, 74:1660, 75:3520,
77:2053, 78:450, 79:15, 81:178, 82:10, 83:115,
84:555, 86:45, 87:10
```

The required count 94 is absent.

### `a=3`

There are 144 exact A4 orbits: 132 with `p=0`, 11 with `p=1`, and one with
`p=3`.  Formula (6.1) gives `b=1,2,4`, respectively.  Since (3.1) requires
at least three isolate-only A5s, only the unique `p=3` orbit survives.  It has
the form

\[
Q\cup\{0\},\quad Q\cup\{1\},\quad Q\cup\{2\}
\]

for an isolate triple `Q`.  Its 1,935 blocker-valid A5 states have histogram

```text
73:270, 74:360, 75:540, 78:180, 84:540, 86:45.
```

The required count 92 is absent.

The counts in this section are symmetry-representative states, not labelled
automorphism-weighted families.  Duplicates are harmless and no orbit is
omitted; target absence does not use orbit weights.

## 7. Reproduction and trust ledger

Compile the three independent checkers with

```sh
clang++ -O3 -std=c++20 scripts/certify_h9_14_record23_a0.cpp \
  -o /tmp/record23_a0
clang++ -O3 -std=c++20 scripts/independent_record23_a1_audit.cpp \
  -o /tmp/record23_a1_audit
clang++ -O3 -std=c++20 scripts/independent_record23_positive_a.cpp \
  -o /tmp/record23_positive_a
```

Run the first two with the `genbg` commands above and run the third directly.

The trust boundary is explicit:

- The `a=0` and primary `a=1` coverage claims rely on the established
  `nauty/genbg` bicoloured graph generator.  The C++ checkers validate every
  property of the emitted stream but do not independently prove that
  `genbg` emitted every isomorphism class.
- No automorphism weight is trusted or required: every column colouring is
  expanded, duplicates are allowed, and only target absence is claimed.
- The positive-branch cross-check is self-contained finite C++ and has no
  graph-generator or SAT/SMT dependency.
- All calculations trust the upstream strict-link theorem, local Macaulay
  necessity, the Lanner order-five bound, and blocker duality.
- The archived Z3 result for record 23 is corroboration only and contributes
  no logical step to this certificate.

## 8. Frozen hashes and verifier commands

The portable archive check does not require `genbg`:

```sh
PYTHONPATH=src:. python3 scripts/certify_h9_14_record23.py --verify-only
```

The full heavy replay additionally checks the frozen `genbg` executable,
compiles all C++ sources, regenerates every incidence graph, and compares the
four exact transcript digests:

```sh
PYTHONPATH=src:. python3 scripts/certify_h9_14_record23.py \
  --verify-only --recompute
```

The canonical certificate digest is

```text
7b13a6b215c4cf6c2c871361c984f150f044a7107bce5b5fb9f484096360da5c
```

The frozen enumerator hashes are

```text
genbg
638224266e86d68d1cd720f18e8b7904d3571849f3da8c6130b1d0e58f515f83

scripts/certify_h9_14_record23_a0.cpp
b376f5d00e70d6c7e11138324e05b109f1d89ffb095fdfd92e8d1bffb7840a50

scripts/certify_h9_14_record23_positive.cpp
ba3ef1a4bbaadb680f6de96e98e610d9fa3cb312804fce04f51131b9d3ab36f5

scripts/independent_record23_positive_a.cpp
fb674b757d97c68d2522f30b3c93845312450e1126635d36bc8952331ad8aac1

scripts/independent_record23_a1_audit.cpp
0fb98ae8fba3bf495507f5e5006ed7dcdf42af35d837ef661fd97a14e01fe55a
```

The deterministic raw-output hashes are

```text
a=0 exact genbg census
bfd872ff8d7edd07271eabfd8575b035793d393ba4ae3843c96a776ecdad0d68

a=1 primary exact genbg census
8f1e607254fe6b03ee46769ac85534d13f0580bec0602a4494161d86bae8b08d

a=1 independently implemented replay
99b7df9981c7c16cb441c1eac1fb94137925117c0c661fa1f9c6c7410d2c795e

positive-a relaxed self-contained census
9203373d6be07b7b10689600b4452592a3850685af974faeba43b887d2202806
```
