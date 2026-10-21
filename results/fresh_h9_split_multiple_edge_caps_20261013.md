# Multiple-edge additions and a common escape beyond the whole bank

The 3,857 numerical no-removal proposals can close all 32 originally open directions of the selected 72-wall chamber. However, **none of their subsets completes that chamber to finite volume**. A rational interval reconstruction of every proposal certifies one common future strict spacelike point satisfying all 72 old and all 3,857 new inequalities. This excludes this literal bank, including mutually incompatible selections; it does not exclude further generations or arbitrary new walls.

## Constructive numerical search

Evaluating the caps on the light-cone endpoints of the old open edges gives 3,715 caps closing one direction and 142 closing two. The latter use four pairs of old directions. A search for mutually compatible additions, retaining every old wall, finds 100-wall configurations with 28 additions closing all 32 old directions.

Numerical traversal of their existing walls reaches 36 finite vertices, 38 finite-to-finite edges, four independent graph cycles, and 248 new open directions in two components. This is local closure followed by exposure of a larger frontier, not a compactness result. These 100-wall configurations have not received exact all-pair Coxeter certificates.

The [selection scout](../data/enumerations/fresh_h9_split_multiple_edge_caps_20261013_scout.json) and [producer](../scripts/fresh_h9_split_multiple_edge_caps_20261013.py) preserve the candidate IDs, old-edge coverage, selected compatibility cliques and reached frontiers. The numerical compatibility comparison uses labels 2 through 120 and strict dotted pairs.

An additional numerical diagnostic finds 312 old spacelike vertices surviving all 3,857 proposals. It also finds 34 spherical eight-wall blocks of the old chamber, all already represented in the original reached edge components. A direct polar construction at the selected remote spacelike vertex tests 12,826 profiles with at most four non-right attachments labelled 3, 4 or 5; none passes even the sixteen-seed pairing screen. These are bounded diagnostics only, archived in [the remote-gap scout](../data/enumerations/fresh_h9_split_remote_diagnostics_20261013_scout.json), reproduced by [this script](../scripts/fresh_h9_split_remote_diagnostics_20261013.py).

## Rigorous exclusion of the reconstructed bank

The [certificate](fresh_h9_split_multiple_edge_escape_20261013.json) works in the independent seed basis `[0,1,2,3,4,5,6,7,8,12]`. It reconstructs the seed Gram from its formula and checks inertia `(9,1)`, as well as the common strict timelike center. All old normal coefficients are enclosed from their archived exact `Q(phi)` quadratic expressions.

For each proposal, the eight edge pairings and the outside-anchor pairing are exact algebraic numbers. Nine independent linear equations determine a line in ten covector coordinates. Outward-rounded integer interval Gaussian elimination encloses its constant and direction vectors. The norm-two equation on this line is an exact quadratic with nonzero leading coefficient and strictly positive discriminant.

The checker encloses both real branches. Exactly one branch lies wholly inside a rational coordinate box of radius `10^-6` around the numerical proposal; the other branch has an entire coordinate interval disjoint from that box. This defines the literal reconstructed root without relying on floating equality. Its norm is exactly two by the quadratic definition. Its other proposed Coxeter angles are not assumed or needed for this exclusion.

The common witness has rational covector coordinates. Its exact square is

\[
\frac{257343836782858514032459241}{16000000000000000000000000}
-\frac{5624985589642989781161529}{640000000000000000000000}\varphi>0.
\]

Its inner product with the strict center is negative. Every old and reconstructed new inequality is strictly satisfied. For the new walls, the smallest certified positive pairing exceeds `0.0005029` in the saved normalization. The center also satisfies every new inequality strictly.

The segment from that timelike center to the spacelike witness therefore reaches the future light cone with every wall inequality strict. The chamber contains an open patch of ideal boundary and has infinite volume. Any subcollection of these new inequalities still contains this patch. This argument does not require the full bank to be mutually Coxeter-compatible.

The calculation uses rational intervals scaled by `10^70`; multiplication, division and square roots round outward using integer arithmetic. All sign conclusions follow from interval separation or exact `Q(phi)` signs. There are 131 distinct linear systems and 3,857 reconstructed branches. The [checker](../scripts/fresh_h9_split_multiple_edge_escape_20261013.py) replays the saved witness independently of the floating proposal step:

```sh
PYTHONPATH=src:. python scripts/fresh_h9_split_multiple_edge_escape_20261013.py --verify
```

The saved replay passes. A separate source review prompted the explicit unused-branch separation check. A [separate execution at 90-digit precision](fresh_h9_split_multiple_edge_escape_independent_20261013.json) also passes: all 38,570 new covector intervals nest inside the 70-digit bounds, all 3,857 branch choices agree, and every old and new wall retains strict center and escape pairings. This is a different-precision replay of the reviewed interval implementation, not a second implementation of its arithmetic. The certificate is a new exclusion of a specified bank, not a full implementation of Allcock's algorithm or an obstruction for the retained seed at every possible root.
