# Mixed indicators and prism-base transfers checked together

The union of the currently exported charge families still has an exact dual
lower bound above the dimension13 target12:

\[
 \frac{2896299884219560}{231058250350149}
 =12.534933852526\ldots>12.
\]

This closes the scope gap between the earlier separate certificates. It
includes the inherited pair, simplex clique and disjoint-port charges,
**all** mixed three-port indicators, and both represented-base prism
transfer families. It remains a barrier for this partial model, not a
dimension bound or evidence that its face masses are geometrically realizable.

## Rechecking the union

Start from the exact repaired face multipliers in the
[shared prism-base transfer certificate](fresh_sqrt5_prism_shared_transfers_20261017.md).
Those multipliers activate410 simplex rows absent from the earlier mixed
indicator check. The [new checker](../scripts/fresh_sqrt5_mixed_shared_transfers_20261017.py)
therefore reconstructs all1,225 positive simplex rows and1,582 positive quad
rows from their source geometry. It checks all1,858 positive quad provenances
and reproduces their original pair coefficients exactly.

The2,807 positive rows use672 vertex forests. All286 literal triples at
each forest are grouped into52,736 marked-forest orbit columns. As before,
the integer cost of an orbit column is its size, and its face coefficient
counts actual aggregate corners occupying at least two pairs of the triple.
The complete mixed coefficient vectors agree across every positive quad
provenance; this is checked rather than assumed.

Only three mixed columns require an increase in the starting vertex masses.
Their total correction is `142693/14000000000`. The checker reconstructs all
inherited pair, clique and port loads and verifies that the starting vertex
masses are their exact maxima. It then checks every mixed column against
the enlarged vertex masses. Forests outside this positive support carry
zero face load.

The same face multipliers are also checked directly against all2,845 shared
and3,840 separate-dimension signed transfer columns. Every signed dual load
is nonpositive; its exact rational value is saved. These columns have zero
vertex cost, so increasing the vertex masses leaves their inequalities
unchanged. Dividing the face objective by the new total vertex mass gives
the displayed lower bound.

## Provenance and scope

Mixed feature equality across provenances does **not** establish equality
of their prism-base maps. The signed transfers retain coefficient-wise
minima over every provenance, as in the frozen exports. The certificate
therefore concerns that explicitly weakened merged model. It does not claim
that one simultaneous split of every row mass satisfies the stronger
individual full-star transfer constraints.

Forced missing base types, their new port charges, and other missing face
families remain outside this model. No optimizer was run. The calculation
uses the same exact feature reconstruction method as the preceding mixed
certificate and is not an independent arithmetic implementation of it.

Generation and a separate complete saved replay passed. The
[certificate](fresh_sqrt5_mixed_shared_transfers_20261017.json) binds all direct
sources and can be recomputed with:

```sh
python3 scripts/fresh_sqrt5_mixed_shared_transfers_20261017.py --verify
```

**FINAL/FROZEN.** The existing checkpoint and every predecessor remain unchanged.
