# A rank-27 rule using normalized A1 cross weights

The new exact rule has maximum vertex weight **67009/1000 = 67.009** in
dimension 27, below the strict face-count threshold **143/2 = 71.5**.
It pays the current conservative face catalogue. Complete finite-attachment
coverage and independent review are separate obligations; this report does
not assert a dimension theorem.

The [producer](../scripts/fresh_sqrt5_all_a1_cross_20261015.py) and
[certificate](fresh_sqrt5_all_a1_cross_20261015.json) use 1,343 parameter-marked
graphs of corner rank at most 27 from the frozen 2,397-literal catalogue.
They check all 14,128 dimension-feasible exact spectator-A1 cases. The rule
has 6,696 nonnegative bonus variables, including zero values. Generation
and saved replay pass. A separate review may bind these frozen files;
its completion is not implied by this primary report.

## Rule and exact vertex bound

Retain weight one on connected pairs at distance one or two and the two
H4-only extra units in the preceding rule. Connected-pair bonuses now
depend on the exact number a of entire isolated A1 components in the full
vertex diagram. On a disconnected marked pair consisting of an A1 root
and a node in an entire non-A1 component, add the specified rooted-component
bonus u_a divided by a. All A1–A1 pairs and unlisted bonuses have weight zero.

For fixed a>0, each node in a non-A1 component makes exactly a pairs with
the A1 roots. Its total contribution is therefore u_a. The cross terms
remain additive over non-A1 components. For each exact a, an unbounded rank
dynamic program enumerates all non-A1 forests on the remaining 27−a
dimensions. The case 27−a=1 is impossible and omitted. All finite component
types allowed over Q(√5) are included; no component-size or forest-count
sample substitutes for this computation.

For a face with core corner rank r, the orthogonal spectators have rank
27−r. Their exact A1 count is common to all four corners; each corner's
own A1 count is added separately. A residual non-A1 rank of one is
impossible. Each literal face payment uses u_a/a, rather than the total
vertex contribution u_a. This distinction is essential.

The floating-point linear program only proposes weights and separating
forests. Every face inequality is repaired and checked in rational
arithmetic, then every bonus is rounded upward to a multiple of 1/1000.
Upward rounding preserves face payments. The complete exact forest
calculation is repeated afterward and gives 67009/1000. No optimality or
dual certificate is claimed for this feasible rule. All bonuses are
nonnegative, so the reviewed triangle and persistent short/H4 payments
remain valid.

## Remaining coverage

The input family census includes reviewed dotted obligations but does not
yet include every terminal state of the active finite attachment search.
The new rule must be evaluated directly on every such zero-touch state,
with its exact dimension-27 spectator counts. Earlier numerical payments
cannot be transferred automatically. A complete coverage composition and
separate implementation review are required before this numerical result
can exclude dimension 27.

```sh
python3 scripts/fresh_sqrt5_all_a1_cross_20261015.py --verify
```
