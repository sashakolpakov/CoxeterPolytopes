# Independent review of the nonabelian rotation obstruction

**PASS:** every coherent added Coxeter wall retaining the specified
133-wall H20 chamber has a dotted contact in each block of five consecutive
chain walls, and therefore at least five dotted contacts in the 26-wall
chain. The proof covers arbitrary real proposed normals and all finite
angle labels. It is a necessary contact-pattern restriction, not a
finite-volume construction or a dimension bound.

The [reviewer](../scripts/fresh_h20_chain_abelian_gate_20261016_review.py) and
[certificate](fresh_h20_chain_abelian_gate_20261016_review.json) import no
producer arithmetic. The primary implementation uses standalone integer
polynomial algorithms; the review reconstructs its polynomial with SymPy
and independently factors the finite-field reductions. The
[primary report](fresh_h20_chain_abelian_gate_20261016.md) supplies the full
geometric notation and source references. The frozen exact chain geometry
and isolated parameter remain inherited dependencies.

The review checks the raw radicands, denominators and residue roots of the
field tower over `Q(sqrt2,sqrt79,sqrt11,sqrt3)`. Nonsquare unit residues at
313 and 359 prove the two successive degree doublings, so the tower is a
field of degree 64. Simple-root Hensel lifting validates the later
degree-one local specializations; reducible formal tower components are
not being treated as field embeddings.

Starting with the chain's central endpoint, the reviewer applies twelve
symbolic rotations and independently obtains the two Chebyshev expressions.
Substitution into the saved endpoint equation and two squarings gives the
degree-50 polynomial in `C=a^2/(a^2+2)`, with coefficients quadratic in the
saved element `w`. Every integer coefficient matches the primary output.
Only the forward implication, that the actual isolated parameter is a
root, is needed; no claim that squaring introduces no extra roots is used.

At the two degree-one primes, independent symbolic factorization verifies:

| Prime | Factor degrees | Other checks |
|---:|---|---|
| 2017 | 50, irreducible | Unit leading coefficient, squarefree reduction |
| 359 | 13 and 37, each irreducible | Unit leading coefficient, squarefree reduction |

All tower coefficients are integral locally, with unit denominators. The
first reduction proves irreducibility over the degree-64 field `K`.
The second gives unequal residue degrees above one prime in `K(C)/K`;
squarefreeness supplies the required locally maximal unramified polynomial
order. A Galois extension has equal residue degrees there, so `K(C)/K`
is not Galois. If `C` lay in the maximal abelian extension of the rationals,
then `Q(C)/Q` would be finite abelian Galois and its compositum with any
`K` would be Galois over `K`, a contradiction. No assumption that `K/Q`
is abelian is required.

Pairings `f_j` with consecutive chain roots satisfy
`delta_(j+2)+delta_j=2C delta_(j+1)`, where `delta_j=f_(j+1)-f_j`.
If five consecutive contacts are finite or parallel, all five values lie
in the maximal abelian extension. Either an internal difference is nonzero,
expressing `C` as a ratio of such values, or both internal differences vanish.
In the latter case both outer differences vanish and the recurrence makes
all chain contacts equal. The first case contradicts nonabelianity; the
second is excluded by the [independently reviewed uniform-contact theorem](fresh_h20_uniform_chain_cap_gate_20261016_review.md).

The parallel endpoint is covered explicitly: `t=2/a` lies inside the
closed interval used by the uniform theorem's positive blockers, and its
support cost is still strictly between 2 and `141/70`. The reviewer checks
these rational interval dependencies and binds the independent uniform
review. Thus the finite-angle theorem is not silently extrapolated to a
missing endpoint. Finally the first 25 chain positions contain five
disjoint blocks of five, requiring at least five dotted relations.

Generation and fresh saved replay passed. This review establishes the
stated dependency and arithmetic scope; external mathematical validation
and any broader nonexistence or construction claim remain separate.

```sh
python3 scripts/fresh_h20_chain_abelian_gate_20261016_review.py
python3 scripts/fresh_h20_chain_abelian_gate_20261016_review.py --verify
```
