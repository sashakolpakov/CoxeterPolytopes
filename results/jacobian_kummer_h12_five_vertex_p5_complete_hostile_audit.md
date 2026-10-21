# Hostile audit: complete Jacobian--Kummer `H12` five-path census

Status: **PASS / FREEZE** (2026-08-10).  No mathematical correction or
lineage downgrade is required.

This audit used the compact frozen transcripts and selected exact rebuilds;
it did not rerun the 83-minute `BBBB` production job.  Every replay was run
under a 0.5 GiB watchdog.  The largest sampled RSS was 0.308 GiB.

## Denominator and histogram

The exact 384-element action was rebuilt and Burnside was recomputed for the
ten binary words modulo path reversal:

| word | orbit count |
|---|---:|
| `AAAA` | 142 |
| `AAAB` | 1,007 |
| `AABA` | 1,084 |
| `AABB` | 5,819 |
| `ABAB` | 6,304 |
| `ABBA` | 3,256 |
| `BAAB` | 3,081 |
| `ABBB` | 36,247 |
| `BABB` | 36,648 |
| `BBBB` | 108,159 |

Their sum is exactly `201,747`.  Grouping by the number `B` of B-edges and
using the independently checked formula `f=752+9B` gives

```text
752^142, 761^2091, 770^18460, 779^72895, 788^108159.
```

The compact `BABB` and `BBBB` transcripts contain respectively 36,648 and
108,159 consecutive rows and cover respectively 14,008,320 and 41,154,000
raw colorings by their recorded orbit sizes.  Their degree-digest statistics
replay as

| shard | degree classes | collision groups | collision members |
|---|---:|---:|---:|
| `BABB` | 34,583 | 1,950 | 4,015 |
| `BBBB` | 100,935 | 6,366 | 13,590 |

Every collision member has an exact color-refinement invariant and those
invariants are distinct inside its degree class.  The one extra refined
`BBBB` row is the cross-star candidate below.  The frozen `B<=2` and `ABBB`
archives independently have 20,693 and 36,247 distinct exact ordinary-Coxeter
type invariants.

## Cross-stratum and cross-shape separation

Facet count separates all `P5` strata except `ABBB` from `BABB`.  An
independent scan found 35,202 `ABBB` and 34,583 `BABB` degree-profile digest
classes with empty intersection.  Since unequal digests imply unequal exact
degree bytes, this is a rigorous nonisomorphism filter, not a hash-equality
argument.  Thus the four path strata are pairwise separated.

Against the complete star census, the same scan found zero common degree
classes at 779 facets for star--`BABB`, and exactly one at 788 facets for
star--`BBBB`.  The unique pair is

```text
star row 8269:  W77, W78, W83, W154
P5 row 20043:   W77, W83, W84, W152
```

A fresh matrix rebuild recovered 788 facets on both sides and the same literal
degree-profile SHA-256

```text
2fda6a517a41c74209234f5b5cf08220b9df8882a10df6a162815f74062d0f94.
```

The literal degree bytes are equal, but their exact stable-refinement words
are unequal, with digests

```text
P5:   5c3b654612dd3289904b2a0dff47619b4009caaadf1529e1a09de8178c98e3a4
star: 71fb5048964565b706f1c726fcf52fed5e2375510ef215b1315c39e6e5a43c61.
```

Both refinements have 395 final cells.  Distinct isomorphism invariants prove
the matrices nonisomorphic even though neither refinement is discrete.  With
the separately frozen literal star--`ABBB` separation, this validates exactly
the claimed high-B union of 197,926 types.  Star versus `B<=2` remains open,
so no star-plus-complete-`P5` theorem is licensed.

## Hyperbolic and carrier gates

The inference to genuine hyperbolic Coxeter polytopes is sound.  The pinned
seed is independently certified finite-volume noncompact in `H12`, and the
112 chosen colors are pairwise disjoint doubling walls.  Allcock's finite-
subtree theorem makes each properly colored five-vertex path union a Coxeter
fundamental polytope, and its index corollary gives reflection index five.
The union consists of five interiors-disjoint seed chambers, so its volume is
exactly five times the finite seed volume.  Equivalently, its reflection group
has finite index five in the nonuniform seed lattice; hence it remains
nonuniform and the polytope is noncompact, not infinite-volume.

The K3 wording also stays within what is proved.  Every descendant boundary
is made of primitive integral roots in
`U(2)+D5(-1)^2+<-4>` and retains a Smith-`1^13` boundary core, so the same
seven-dimensional full-Picard family is a weak divisor-reflection carrier.
Nothing here proves nef-cone equality, effectivity of every root, or a
surface-automorphism action.

## Lineage verdict

The positive theorem is a complete census inside the frozen JK-H12 factory.
It does not claim a new dimension, commensurability class, or worldwide
priority.  Allcock's exponential abundance in `H12` is explicitly credited,
and the underlying 172-facet seed remains only apparently absent from the
bounded literature/corpus audit.  Therefore “new-to-project; worldwide
priority open” is the correct lineage formulation.

