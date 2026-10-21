# Rank-162 obstruction to the fixed-omission local estimates

The [exact standard-library construction](../scripts/fresh_nikulin_rank162_local_obstruction_20261016.py) and [certificate](fresh_nikulin_rank162_local_obstruction_20261016.json) give two full Coxeter Gram configurations at affine rank **N=162**, hence ambient hyperbolic dimension 163. Both satisfy the existing endpoint budgets and all present normalized wall pairings are strictly below 14. They therefore obstruct reducing the reviewed 163 bound by rank-sensitive bookkeeping alone while preserving these local data and omission choices.

These are **not finite-volume polytopes**. Each displayed chamber has an explicitly certified null ray satisfying every wall inequality strictly, and consequently an open patch on the ideal boundary and infinite volume. This is part of the certificate, not an unresolved volume test.

The first configuration replaces one E6 filler in the frozen sharp rank-163 local pattern by D5:

- two ordinary E7 factors with both contact rows at node 7;
- two ordinary B2 factors with both rows at node 0;
- one avoid-node A1 with coincident omission 0 and forbidden node 1;
- 23 ordinary E6 factors, with rows at the respective omitted nodes 2 and 4;
- one ordinary D5 factor, with rows at the respective omitted nodes 0 and 5.

All affine kernel heights are exactly one, so the common height ratio is one. The endpoint costs are exactly `(2,2)`. The E6 energies are `4/3`; the D5 energy is `5/4`, compared with target `10/9`. Thus the rank is `14+4+1+138+5=162`, and

\[
\sum W=23\frac43+\frac54=\frac{383}{12},\qquad
\sum T=\frac{322}{9},\qquad
\sum(T-W)=\frac{139}{36}.
\]

The exact balanced Lorentz pairing is

\[
-(x,y)=\frac12\sum W-2=\frac{335}{24}
=14-\frac1{24}<14.
\]

A second configuration replaces the two E7 and two B2 factors by two ordinary E8 factors at node 8 and one B2 factor at node 0. It has the same rank, energy and pairing, but endpoint costs `(16/9,11/6)`. Its two frame Schur complements are strictly positive, `2/9` and `1/6`; both selected endpoint diagrams are elliptic. Thus the obstruction is not confined to saturated budgets or ideal endpoint degeneracies.

The matrices are rational. Simply laced roots have square two; the B2 long root is rescaled to square four, with the unchanged Coxeter label four. The checker reconstructs every affine block, its positive kernel, every grounded inverse and positive pivot, both endpoint costs, all contacts, and the full Gram. It checks the nonpositive off-diagonal entries and their Coxeter values, including the ultraparallel endpoint pair. All normalized negative pairings are strictly below 14, not merely the selected endpoint pairing.

A block Schur calculation proves signatures `(163,1,29)` and `(163,1,28)` respectively. Every affine kernel has the same endpoint pairings `(-1,-1)`, so the common-height linear relations hold in the full Gram. After grounding the affine blocks, the remaining nondegenerate block has determinant minus one and exactly one positive and one negative direction; all other directions are explicitly radical.

The certificate then constructs a second null vector in the same causal sheet as the cusp ray. A small rational perturbation puts it strictly inside every affine wall, and an exact cusp-ray correction restores null norm. Its pairing with every present wall is negative. The maximum such pairing is `-1/9000` in the E7 configuration and `-1/15000` in the E8 configuration. Strictness proves the infinite-volume conclusion.

Generation and saved replay passed:

```sh
python3 scripts/fresh_nikulin_rank162_local_obstruction_20261016.py --verify
```

The source and JSON are frozen for separate review. No claim is made that these configurations have the full incidence or additional framing walls required of a finite-volume polytope. Changing omissions, using further framing pairs, or enforcing completion of the other faces remains available. In particular, arbitrary root combinations need not be actual facets, so their pairings cannot be substituted into the narrow-facet hypothesis without an incidence argument.
