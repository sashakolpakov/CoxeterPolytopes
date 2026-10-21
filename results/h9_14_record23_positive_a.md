# Record 23: exact positive-`a` blocker census

This note records the independent finite calculation in
`scripts/independent_record23_positive_a.cpp`.  It covers the
Dehn--Sommerville branches with one, two, or three missing four-faces.  It
uses only the vertex-link `g4` capacities, the minimal-nonface shadow, and
blocker duality.  In particular, it is a deliberately *weaker* filter than
the full vertex/edge/triple local-link census: excluding a branch here also
excludes it under the stronger gates.

## 1. Graph and Dehn--Sommerville reduction

The record-23 missing-edge graph is

\[
H=K_{3,2}[A,B]\sqcup9K_1,
\qquad
A=\{0,1,2\},\quad B=\{12,13\},\quad
I=\{3,\ldots,11\}.
\]

Its full automorphism group is

\[
\operatorname{Aut}(H)=S_3\times S_2\times S_9,
\]

and its only inclusion-minimal vertex covers are `A` and `B`.  The schedules
and summed vertex-link `g4` values are

| missing fours `a` | missing independent fives | minimal A5 count before the A4 shadow | facets | sum of vertex-link `g4` |
|---:|---:|---:|---:|---:|
| 0 | 7  | 7  | 98 | 36 |
| 1 | 13 | determined below | 96 | 30 |
| 2 | 19 | determined below | 94 | 24 |
| 3 | 25 | determined below | 92 | 18 |
| 4 | 31 | -- | 90 | 12 |
| 5 | 37 | -- | 88 | 6 |

Strict positivity in all fourteen vertex links gives a lower bound of 14 on
the last column, so `a=4,5` are impossible.

When all graph-independent four- and five-sets are provisionally faces, the
baseline vertex-link values are

\[
\beta_v=
\begin{cases}
5,&v\in A,\\
6,&v\in I,\\
1,&v\in B.
\end{cases}
\]

The local missing-four intervals force every A4 to avoid `B`.  An A4 made
only of isolates has an independent five-extension through either vertex of
`B`, which lowers that vertex's baseline from one to zero.  Likewise a
minimal A5 meeting `B` lowers its local value from one to zero.  Consequently
every extra minimal nonface lies in

\[
W=A\sqcup I,
\]

and every A4 meets `A`.

## 2. Complete blocker decomposition

Let `G` be the clutter on `W` consisting of the missing four-faces and the
additional minimal missing five-faces.  Let `Z` be the subfamily of A5s that
are contained in `I`.

Every transversal of `H` contains either `A` or `B`.  Therefore there is no
transversal of order at most four precisely when

1. `Z` has empty total intersection; and
2. `G` has no transversal of order at most two on `W`.

Indeed, `A` already hits every member of `G` meeting `A`, and it can be
extended by at most one point.  The members it misses are exactly `Z`.
Similarly, `B` misses every member of `G` and can be extended by at most two
points of `W`.

Once these conditions hold, every five-element transversal is minimal.  The
five-element graph covers split without overlap into

- `B` plus a three-subset of `W` hitting `G`; and
- `A` plus a two-subset of `I` hitting `Z`.

(The possible overlap `A union B` does not hit `Z`.)  Thus the exact number
of facet complements is

\[
N=\tau_3(G)+\tau_2(Z).                                      \tag{2.1}
\]

If this differs from the scheduled facet count, then either the complex has
a larger minimal transversal, contradicting purity, or all minimal
transversals have order five and the scheduled facet count is wrong.

Because two five-subsets of the nine-element set `I` always meet, the first
blocker condition implies

\[
|Z|\geq3.                                                    \tag{2.2}
\]

## 3. Missing-four orbit compiler

For `a<=3`, write a missing four-face as an active part `R_i`, a nonempty
subset of `A`, and an isolate part.  The nine isolates are encoded by their
row-membership columns `S subset [a]`, with multiplicities `n_S`.  The
compiler enumerates the exact integer system

\[
\sum_S n_S=9,
\qquad
\sum_{S\ni i}n_S=4-|R_i|.
\]

It rejects repeated rows and applies the exact A4 local predicates: a pair
containing an active vertex occurs in at most two A4s, and a triple containing
an active vertex occurs in at most one A4.  The remaining pair/triple bounds
are automatic for `a<=3`.  Canonical keys under

\[
S_3(A)\times S_a(\text{rows})\times S_9(I)
\]

are obtained by minimizing the ordered active parts and the isolate-column
multiplicities.  This gives

| `a` | exact A4 orbits | overlap statistic `p` |
|---:|---:|---|
| 1 | 3 | active-part sizes 1, 2, 3 |
| 2 | 17 | 16 with `p=0`, one with `p=1` |
| 3 | 144 | 132 with `p=0`, 11 with `p=1`, one with `p=3` |

Here `p` counts pairs of A4s sharing an isolate triple.  Each A4 has eight
independent five-extensions in `W`.  The local triple predicate ensures that
the only pairwise overlap between these shadows comes from such an isolate
triple, and no five-set lies in three shadows.  Hence the shadow has size

\[
8a-p,
\]

and the number of additional minimal A5s is

\[
b=(7+6a)-(8a-p)=7-2a+p.                                    \tag{3.1}
\]

For `a=3`, (2.2) and (3.1) remove the 143 orbits with `p=0,1`.  The unique
remaining orbit has `p=3` and, up to symmetry, is

\[
Q\cup\{0\},\quad Q\cup\{1\},\quad Q\cup\{2\},
\]

where `Q` is one isolate triple.

## 4. Exact positive-branch census

For a fixed A4 family, let `Sh` be its union of independent five-shadows.
For each vertex `v in W`, the enumerator permits at most

\[
c_v=\beta_v+6r_v-|\{F\in Sh:v\in F\}|-1                 \tag{4.1}
\]

chosen minimal A5s through `v`, where `r_v` is its A4 incidence.  This is
exactly strict vertex-link `g4>=1`.

The program exhausts the isolate-only subfamily `Z` with empty intersection,
adds every eligible remaining A5 not containing an A4, checks (4.1), applies
the final two-transversal test to the *complete* family, and evaluates (2.1).
No partial-family blocker condition is used.

### `a=1`

For one A4, its active-part size `k=1,2,3` gives the three A4 orbits.  Its
shadow has size eight, so `b=5`.  Formula (4.1) becomes

| vertex type | A5 incidence capacity |
|---|---:|
| active vertex in A4 | 2 |
| active vertex outside A4 | 3 |
| isolate in A4 | 3 |
| isolate outside A4 | 4 |

For `|Z|=3,4`, `Z` is canonicalized under

\[
S_{|Z|}\times S_{4-k}\times S_{5+k};
\]

the orbit counts for `k=1,2,3` are respectively `(56,1976)`, `(33,1009)`,
and `(14,334)`.  All remaining active A5 choices are then expanded.  For
`|Z|=5`, direct capacity-pruned enumeration is faster.

The exact representative-state supports are:

| `k` | blocker-valid representative states | support of `N` | target 96 |
|---:|---:|---:|---:|
| 1 | 133,726 | 56--92 | 0 |
| 2 | 508,884 | 56--92 | 0 |
| 3 | 9,515,916 | 54--92 | 0 |

These are representative-state counts, not automorphism-weighted labelled
family counts: different symmetry stages deliberately retain some duplicate
completions.  That does not affect the absence of 96.

As a replay fixture, the relaxed maximum 92 for `k=1` is attained by

\[
Q=\{0,3,4,5\}
\]

and A5s

\[
\begin{split}
&\{3,4,5,6,7\},\quad\{3,4,5,6,8\},\quad
  \{7,8,9,10,11\},\\
&\{0,1,2,3,4\},\quad\{0,1,2,5,6\}.
\end{split}
\]

It has `tau3=71`, `tau2=21`, and no transversal of order at most four.

### `a=2`

All 17 A4 orbits are exhausted.  The 16 `p=0` types have `b=3`; the unique
`p=1` type has `b=4`.  There are 16,778 blocker-valid representative states,
with histogram

```text
63:120, 65:180, 67:108, 68:540, 69:1260, 70:660,
71:1739, 72:1610, 73:1950, 74:1660, 75:3520,
77:2053, 78:450, 79:15, 81:178, 82:10, 83:115,
84:555, 86:45, 87:10
```

The required count 94 is absent.

### `a=3`

Only the unique `p=3` orbit reaches `b=4`.  Its exact A5 census has 1,935
blocker-valid representative states:

```text
73:270, 74:360, 75:540, 78:180, 84:540, 86:45.
```

The required count 92 is absent.

## 5. Scope and trust ledger

- The positive-`a` executable is self-contained ISO C++20 and uses only
  integer masks, permutations, and finite loops.  It has no external solver
  or graph-generator dependency.
- Its quotient counts are used only for coverage and target absence; no
  automorphism weights are asserted.
- It trusts the upstream geometric input that every relevant link has strict
  `g4>=1`, and the standard blocker correspondence between maximal faces and
  minimal transversals.
- It deliberately omits the stronger edge-link, triple-link, and local
  Macaulay gates on the chosen A5s.  This enlarges its state space.  Since the
  target facet counts are already absent, those stronger gates are not needed
  for the positive-`a` exclusion.
