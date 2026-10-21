# Record 31: the entire compact `a=2` schedule is excluded

## Certified statement

Status: **certified compact-only no-go**. This is a theorem-track advance,
not a new hyperbolic Coxeter polytope.

Record 31 has missing-edge graph `M?????????O?O?o??`, the path

```text
11--0--13--1--12
```

together with nine isolated vertices. Its full automorphism group on the
fourteen facets is `C2 x S9`. In the compact Dehn--Sommerville schedule
`a=2`, every completion must have 21 missing five-subsets and face vector

```text
(14, 87, 319, 769, 1275, 1463, 1126, 531, 118).
```

The exact conclusion is:

> No completion of Record 31 in the compact `a=2` schedule is the nerve of a
> compact hyperbolic Coxeter 9-polytope.

This closes one schedule, not the whole graph shadow. The compact schedules
`a=0,1` remain open. No finite-volume noncompact chamber is excluded.

## Exhaustive gate chain

The positive-A4 archive contains 49 full graph-action orbits with two
missing four-faces. Their complete disposition is

| exact gate | survivors |
|:---|---:|
| positive-A4 orbit representatives | 49 |
| strict face-link feasibility | 18 |
| pure scheduled blocker and ridge degree two | 12 |
| connected, every-link mod-two homology audit | 12 |
| unrestricted compact Coxeter-label gate | 0 |

The 12 ridge-valid targets arise from positive-A4 orbit indices

```text
0:1, 2:1, 3:1, 5:1, 34:6, 40:2.
```

Every one has minimal-transversal histogram `{5:118}` and ridge-degree
histogram `{2:531}`. Thus these are not casualties of a loose facet-count
filter: they reach the strongest combinatorial and topology gates currently
used here before failing the compact Coxeter-label condition.

## The analytic face-link obstruction

The first reduction is a theorem, not an enumerative heuristic. Let $K$ be
the boundary complex of the simplicial polar of a putative compact
9-polytope. If $F\in K$ has $s$ vertices, then its face figure has
dimension

\[
D=9-s.
\]

Every minimal nonface of every link has order at most five. Indeed, a
minimal nonface $T$ of \(\operatorname{lk}_K(F)\) lifts to a minimal
nonface $M\subseteq F\cup T$ of $K$, with $M\setminus F=T$; the Coxeter
minimal nonfaces are Lanner supports and have order at most five.

Now suppose $g_r=0$ in such a polytopal link. The Murai--Nevo equality
theorem gives an \((r-1)\)-stacked triangulation. An inclusion-minimal
interior simplex then has order at least

\[
D-r+2.
\]

It determines a missing face of the link. Therefore, whenever
$D-r+2>5$, equality is impossible and $g_r>0$. For the faces used in this
census, $s=1,2,3$ and $r=5-s$, so respectively

\[
(D,r)=(8,4),(7,3),(6,2),
\qquad g_{5-s}>0.
\]

This is the strict local lemma proved in
[`h9_14_local_link_constraints.md`](h9_14_local_link_constraints.md), applied
after the two A4s are fixed.

There is also an exact numerical form of the gate. For a fixed face $F$,
let

\[
p_F=g_{r-1},\qquad b_F=g_r
\]

be the link entries before choosing the remaining minimal A5s, and let
$t_F$ count the chosen A5s containing $F$. The final entry is
$g_r=b_F-t_F$. Strict positivity and the Macaulay inequality give

\[
\boxed{
\max\{0,b_F-p_F^{\langle r-1\rangle}\}
\le t_F\le b_F-1.
}
\]

For 31 of the 49 A4 orbits, an edge or triangle $F$ already has $b_F=0$.
Its required incidence interval is consequently `[0,-1]`, so no choice of
the remaining A5s can repair it. The excluded indices are

```text
4, 7, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
26, 27, 28, 29, 30, 31, 32, 36, 38, 39, 41, 43, 45, 46, 48.
```

The residual indices are

```text
0, 1, 2, 3, 5, 6, 8, 9, 10, 11, 33, 34, 35, 37, 40, 42, 44, 47.
```

For each A4 orbit the implementation derives all 420 link constraints through
faces of sizes one, two, and three. A separate standard-library Python replay
recomputes every $f$-, $h$-, and $g$-prefix and recovers exactly the same
31/18 partition and the same numerical contradictions.

## Lossless incidence census

The remaining 18 representatives require between three and six additional
A5s after their forced A5 shadows. The completion search imposes, exactly:

- every one of the 420 local incidence intervals;
- destruction of every forbidden transversal of order at most four;
- the scheduled 118 order-five facet complements;
- purity of the complete blocker; and
- degree two at every one of the 531 ridges.

Two different quotient implementations replay this search. Production uses a
colored incidence graph and nauty canonicalization. The second traversal uses
a direct binary-matrix invariant: it explicitly tests the path flip, keeps A4
and A5 row colors separate, and sorts the nine isolate-incidence columns. It
also reverses the pivot tie and candidate order. Their search transcripts are
different (17,015,439 versus 13,253,521 visited nodes), but their final
canonical ridge-target sets agree exactly: 12 orbits on each side, with no
orbit present on only one side.

The direct final comparison is performed under the complete `C2 x S9`
action. Thus equality is an equality of target orbits, not an equality of
accidental printed representatives.

## The 12 topology survivors

All 12 targets pass a literal downstream reconstruction. For every target it
recomputes the 118 facet complements and 531 ridge complements, then verifies

- the complete scheduled face vector;
- connectedness of the dual graph;
- spherical Euler characteristic in every link; and
- exact mod-two homology of every link.

This is deliberately recorded as a necessary topology audit. A connected
mod-two homology sphere is not, by itself, claimed to be a PL sphere, the
boundary of a convex polytope, or an ordinary polytopal realization.

## Why all 12 fail the compact Coxeter gate

For every target, the label join enumerates every classical Lanner placement
on its order-four and order-five minimal nonfaces: 108 labelled order-four
placements and 420 labelled order-five placements. At every stage it checks
ellipticity on the restrictions cut out by all 118 actual rank-nine vertices.
No state cap is used. Every target has terminal sparse state count zero.

The zero is unrestricted with respect to omitted finite non-right labels.
Write a finite Coxeter Gram restriction as $I-A$, where $A$ is the
nonnegative cosine-adjacency matrix. The sparse join assigns $m=2$, hence
zero adjacency, to every not-yet-forced intersecting pair. Replacing any such
entry by $m\ge3$ only increases $A$. Perron--Frobenius monotonicity cannot
turn a restriction with spectral radius at least one into an elliptic one.
Thus labels omitted by the sparse representation cannot repair any exhausted
state.

An independent downstream implementation regenerates the 108 and 420
placement tables, uses the opposite support-start and tie rules, rechecks the
literal blocker and ridge laws, and again obtains zero uncapped terminal
states for all 12 targets.

## Scope and reproducibility

The certificate makes only a compact `a=2` claim. In particular:

- it performs no parabolic-subdiagram or Vinberg endpoint census;
- it makes no inference about noncompact finite-volume chambers;
- it does not promote the topology audit to a convex-polytopality theorem;
- it does not close Record 31 while `a=0,1` remain open; and
- it constructs no new hyperbolic Coxeter polytope.

The frozen certificate is
`data/enumerations/h9_14_record31_a2_compact_no_go.json`. Its canonical
certificate SHA-256 digest is

```text
b8d8d6e7bd3ea4f16425a8f13c2e5449957717c2e78fe5172db734b38e6608b5
```

The quick source-provenance, digest, and semantic verification is

```sh
PYTHONPATH=.:src python3 scripts/certify_h9_14_record31_a2_no_go.py
```

The full independent incidence, link, topology, and label replay is

```sh
PYTHONPATH=.:src python3 scripts/certify_h9_14_record31_a2_no_go.py --fresh
```

The focused regression tests are in
`tests/test_h9_14_record31_a2_no_go.py`.
