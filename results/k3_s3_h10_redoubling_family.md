# An infinite K3 divisor-lattice redoubling family in `H10`

## Result

**TOOT.**  The 13-wall chamber `K3-S3-C2-H10-13` contains two parallel
doubling walls.  Consequently, for every integer $N\geq 1$, it produces a
finite-volume noncompact Coxeter polytope $P_N\subset H^{10}$ with


```text
facets                    N + 12
reflection-subgroup index N
volume                    N * volume(P_1).
```

For odd $N$, reversing the path gives the same lattice-decorated Gram type.
For every even $N$, the two possible path phases give two nonisomorphic
polytopes, already distinguished by the number of $m=4$ edges.  Different
values of $N$ have different facet counts and volumes.  Thus this is an
infinite positive family, not a finite list of speculative candidates.

Every wall is an explicit primitive class in



\[
 F=U+E_8(-2)+\langle-2\rangle=\operatorname{NS}(X)
\]

for a very general member of the previously constructed explicit
nine-parameter K3 family.  This is a divisor-lattice realization.  It does
not assert that the new chambers are nef cones, that every wall is an
irreducible curve, or that the reflections are surface automorphisms.

## The two walls and the normalization trap

In the standard basis of $F$, take one-based seed walls 10 and 13:


\[
 a=(-2,-2,0,0,0,0,0,-1,0,1,1),\qquad
 b=(0,0,0,0,0,0,0,0,0,0,-1).
\]

They satisfy


\[
 a^2=b^2=-2,\qquad (a,b)=2,
 \qquad (a,b)^2=4,
 \qquad \frac{(a,b)^2}{a^2b^2}=1.
\]

The raw square $4$ must not be mistaken for the normalized ratio.  The
ratio is exactly $1$, so the walls are **parallel/asymptotic**, not
ultraparallel.  Their sum


\[
 \delta=a+b=(-2,-2,0,0,0,0,0,-1,0,1,0)
\]

is isotropic and represents their common ideal direction.

The complete seed relation matrix has exactly these two doubling walls:
every finite wall they meet has an even Coxeter label.  Since the walls do
not meet in $H^{10}$, Allcock's finite-subtree theorem applies.  The group


\[
 \langle s_a,s_b\rangle\cong C_2*C_2
\]

has a bi-infinite Cayley line, and its finite subtrees are precisely finite
paths.  See Daniel Allcock,
[*Infinitely many hyperbolic Coxeter groups through dimension
19*](https://doi.org/10.2140/gt.2006.10.737), Theorem 4.1 and Corollary 4.2.

## Closed boundary proof

Nine seed walls are orthogonal to both $a$ and $b$, hence fixed by both
reflections.  Let $c$ and $d$ be seed walls 9 and 11.  The remaining
outer walls of the phase-`A` $N$-path are


\[
\begin{aligned}
 q_{2m}&=c+2m^2a+2m(m-1)b,\\
 q_{2m+1}&=d+2m(m+1)a+2m^2b,
\end{aligned}
\]

for $q_0,\ldots,q_N$.  The two endpoint walls are


\[
 b,\qquad Na+(N-1)b.
\]

Phase `B` swaps $(a,b)$ and $(c,d)$.  The formulas follow immediately by
induction from


\[
 s_a(xa+yb)=(-x+2y)a+yb,
 \qquad
 s_b(xa+yb)=xa+(2x-y)b.
\]

They also prove that the endpoint roots always have square $-2$,
divisibility $2$, and inner product $2$.  Therefore every descendant is
again redoublable.  Its boundary contains


```text
9 fixed walls + (N+1) moving walls + 2 endpoint walls = N+12 walls.
```

The certifier literally constructs the union of reflected seed chambers and
compares its oriented boundary with this closed list.  It verifies the
formulas through $N=12$, and the displayed reflection recurrence is the
symbolic proof for arbitrary $N$.  Every tested boundary spans $F$ with
index one; for arbitrary $N$, all roots remain primitive and all
reflections integral because they are images of seed roots under integral
lattice isometries.

## Exact fingerprints for every `N`

Here label `0` means parallel, label `1` means ultraparallel, and labels
`2,3,4` have their usual Coxeter meanings.

For $N=2k+1$, the relation-label counts are


\[
 (n_0,n_1,n_2,n_3,n_4)
 = (k+2,\;2k^2+3k,\;65+19k,\;9+2k,\;2),
\]

and the mirror-type counts are


\[
 \#(-4,2)=9+k,qquad \#(-2,1)=2+k,qquad \#(-2,2)=2.
\]

For $N=2k$, phase `A` has


\[
 (n_0,n_1,n_2,n_3,n_4)
 =(k+2,\;2k^2+k-1,\;56+19k,\;8+2k,\;1),
\]

while phase `B` has


\[
 (n_0,n_1,n_2,n_3,n_4)
 =(k+1,\;2k^2+k-1,\;55+19k,\;8+2k,\;3).
\]

Their mirror-type counts are respectively


\[
 (\#(-4,2),\#(-2,1),\#(-2,2))
 =(9+k,1+k,2)
\]

and


\[
 (8+k,2+k,2).
\]

The $m=4$ counts $1$ and $3$ prove that the two even phases are not
isomorphic.  For odd $N$, exact labelled-graph and lattice-decorated Gram
isomorphism implement path reversal.

## The first descendants

The first new-to-project-ledger examples are:

| polytope | facets | mirror types `(-2,1);(-2,2);(-4,2)` | labels `0;1;2;3;4` |
|---|---:|---:|---:|
| $P_{2,A}$ | 14 | `2;2;10` | `3;2;75;10;1` |
| $P_{2,B}$ | 14 | `3;2;9` | `2;2;74;10;3` |
| $P_3$ | 15 | `3;2;10` | `3;5;84;11;2` |
| $P_{4,A}$ | 16 | `3;2;11` | `4;9;94;12;1` |
| $P_{4,B}$ | 16 | `4;2;10` | `3;9;93;12;3` |

Allcock's theorem already proves finite volume for every row and every later
$N$.  As a method-independent calibration, the two 14-facet chambers were
also run through the exact Vinberg endpoint criterion and complete
extreme-ray enumeration:

| chamber | finite vertices | cusps | elliptic rank-9 edges | endpoint types `(finite,ideal)` | hyperideal rays |
|---|---:|---:|---:|---|---:|
| $P_{2,A}$ | 41 | 5 | 252 | `(2,0):164; (1,1):82; (0,2):6` | 0 |
| $P_{2,B}$ | 52 | 6 | 318 | `(2,0):212; (1,1):96; (0,2):10` | 0 |

Every elliptic edge has exactly two finite/ideal endpoints, every full
parabolic neighborhood passes Vinberg's test, every ideal link has rank 9,
and the equality extreme rays span each wall hyperplane.  Thus the 14 walls
in each case are irredundant.

## Index, volume, and noncompactness

Allcock's Corollary 4.2 identifies the boundary-reflection group of $P_N$
as an index-$N$ subgroup of the seed reflection group.  Hence


\[
 \chi(W_N)=-\frac{N\cdot527}{33443020800}
\]

and, in the standard $H^{10}$ Gauss--Bonnet normalization,


\[
 \operatorname{vol}(P_N)
 =\frac{N\cdot527\pi^5}{987614208000}.
\]

The seed group is nonuniform.  Every finite-index subgroup contains powers
of its parabolic elements, so the descendants remain nonuniform; therefore
all $P_N$ are noncompact as well as finite-volume.

## Lineage: what is and is not established

The exact bounded audit scans all 49 dimension-10 records in the pinned
CoxIter lanes: 27 ordinary records and 22 index-two-doubling records.  Their
facet counts range up to 44.  Every family member with a present header has
zero exact edge-labelled relation-graph matches.  For absent headers the
nonmatch is immediate.

One coarse collision is deliberately retained.  $P_{2,B}$ has the same
global label histogram as


```text
testing_index2/10-2,1,1,1,1,1,1,1,1,1,1-12-13_1.coxiter
```

but their node-profile multisets differ and a full labelled-graph
isomorphism test returns false.  Histogram comparison alone would have been
insufficient.

The known seed $P_1$ is itself absent from these pinned lanes.  That is a
useful negative control: the scan proves only a bounded corpus nonmatch, not
worldwide novelty.  Accordingly, the descendants are recorded as
**new-to-the-project-ledger, worldwide priority open**.  The separate
[hostile two-ray audit](h10_two_k3_redoubling_rays_hostile_audit.md) also
shows that the analogous ray on
$U+E_8(-2)+\langle-4\rangle$ is classical and exactly present in CoxIter;
the present $F$-ray is not that known ray in disguise.

## What this adds to the algebraic-surface program

This gives a precise reusable positive criterion:

> If a hyperbolic Néron--Severi lattice carries a complete finite-volume
> reflection chamber with two disjoint doubling walls, then one algebraic
> surface family weakly carries an infinite tower of finite-volume Coxeter
> polytopes obtained from every finite subtree of the doubling-wall Cayley
> graph.

For two parallel walls, the Cayley graph is a line and the tower is the
explicit $P_N$ family above.  This does not yet invert an arbitrary Coxeter
polytope into a surface.  It does establish a meaningful surface-to-polytope
functor: the softer algebraic datum is the Picard lattice plus a complete
reflective chamber and a redoubling configuration; the rigid outputs are
infinitely many exact hyperbolic Coxeter polytopes.

## Reproduction

- [certificate generator](../scripts/certify_k3_s3_h10_redoubling_family.py)
- [exact JSON certificate](../data/enumerations/k3_s3_h10_redoubling_family.json)
- [focused regression tests](../tests/test_k3_s3_h10_redoubling_family.py)
- [underlying full-Picard K3 carrier](k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.md)


```bash
PYTHONPATH=src:. python3 scripts/certify_k3_s3_h10_redoubling_family.py
PYTHONPATH=src:. pytest -q tests/test_k3_s3_h10_redoubling_family.py
```

Frozen payload SHA-256:

`11da6289847d134ffe56081604fbf6df8779f48f3984fd2a4a0e2aee4b8309f9`
