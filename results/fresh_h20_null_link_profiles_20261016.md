# Complete incident-contact domain at the new H20 cusp

Any norm-2 acute Coxeter wall through the persistent null ray nu is right
to the affine G2 triple. Its contacts with the other 13 old incident
walls lie in **125 explicit finite profiles or two continuous single-root
tails**. The [source](../scripts/fresh_h20_null_link_profiles_20261016.py)
and [data](fresh_h20_null_link_profiles_20261016.json) pass generation and
fresh saved replay. This is a contact-domain theorem, not an escape-cutting
wall or a completed exclusion of such walls.

The exact old geometry is bound through the
[five-wall cusp completion](fresh_h20_affine_a17_completion_20261016.md).
The 13 spherical incident roots form

```
A2[3,4] + A5[6,7,8,9,10] + A4[12,13,14,15] + A1[17] + A1[31].
```

Write d_i=-r.n_i>=0 for their norm-2 contact magnitudes and J for the
inverse of their rational Gram matrix. The quotient `nu^perp/R*nu` is
positive definite. Orthogonal projection to the span of these 13 roots
therefore gives the necessary inequality `d^t J d<=2`. Every entry of J
is nonnegative. For an Ar path its entries are
`min(i,j)*(r+1-max(i,j))/(r+1)`, using indices1 through r.

Every nonzero Coxeter contact has magnitude at least1. A dotted contact
has magnitude greater than2, and every inverse diagonal is at least1/2,
so no incident contact can be dotted. The four smallest inverse
diagonals are `1/2,1/2,2/3,2/3`, whose sum is7/3>2. Consequently at most
three of the 13 contacts are nonzero.

Every nonisolated inverse diagonal is at least2/3. Thus its magnitude is
at most sqrt(3), so its finite label is at most6. If at least two contacts
are nonzero, the other contact already costs at least1/2; an isolated
A1 contact then satisfies `d_i^2/2<=3/2`, again giving label at most6.
Parallel contacts of magnitude2 are possible only at a single isolated
A1 and are retained explicitly.

The checker evaluates all19,604 profiles with one, two or three nonzero
positions and labels3 through6 over the exact squarefree-radical basis
`1,sqrt(2),sqrt(3),sqrt(5),sqrt(6),sqrt(10),sqrt(15),sqrt(30)`.
Exactly124 pass the budget; adding the all-right profile gives125.
All remaining possibilities are covered by the two families

```
only root17 nonright, magnitude d in [sqrt(3),2];
only root31 nonright, magnitude d in [sqrt(3),2].
```

These continuous intervals deliberately include non-Coxeter values.
They cover all labels at least6 and the parallel endpoint, with label6
overlapping the finite census. No finite label cutoff is imposed on the
actual candidate walls.

The domain is necessary and generally not sufficient: all old
nonincident contacts, mutual compatibility with any other new wall,
coherence and an escape-cutting inequality still have to be checked.

```sh
python scripts/fresh_h20_null_link_profiles_20261016.py --verify
```
