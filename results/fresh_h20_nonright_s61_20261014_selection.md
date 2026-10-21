# Numerical selection from the nonright-S61 bank

The saved 847-proposal bank has **three robust finite-contact pairs**, all right-angled, and no numerical finite corner involving two proposed caps. A common strict spacelike point survives the intersection of **all 133 retained walls and all 847 proposals**, even before mutual compatibility is imposed. Thus this bounded bank does not appear to close the global gap. These are numerical diagnostics, not exact exclusion or construction certificates.

The [selection script](../scripts/fresh_h20_nonright_s61_20261014_selection.py) reads the [bank](fresh_h20_nonright_s61_20261014_scout.json) without modifying it. The [saved output](fresh_h20_nonright_s61_20261014_selection.json) binds the bank, frozen 133-wall geometry and relevant loader sources. No new contact events were generated.

## Compatibility and local geometry

Write each proposal as `r=−Σc_iω_i+tU+xV+zW`, with the residual frame of signature `(-,+,+)` and support inverse `J`. The mutual pairing is

`(r,r′)=cᵀJc′−tt′+xx′+zz′`.

The maximum reconstructed norm error is below `4×10⁻¹⁵`. At a pairing tolerance of `10⁻¹⁰`, the 358,281 unordered proposal pairs split into 8,195 dotted pairs, three right pairs, 209,877 obtuse pairs and 140,206 forbidden numerical angle pairs. Near fits to π/7, large dihedral labels, or right angles at a looser `10⁻⁶` tolerance are rejected. The tolerance is a screening rule, not a proof that an algebraic angle equals or differs from a Coxeter angle.

The right pairs are:

| Proposal IDs | First support contacts | Second support contacts | Shared anchors 65,95 |
|---|---|---|---|
| 103,811 | walls 1,7 at π/3 | walls 9,10 at π/3 | π/5, π/6 |
| 104,812 | walls 1,7 at π/3 | walls 9,10 at π/3 | π/6, π/5 |
| 105,813 | walls 1,7 at π/3 | walls 9,10 at π/3 | π/6, π/6 |

Within each pair the saved residual vectors coincide. There is an exact **conditional** cancellation behind the numerical right angle: both support projection squares are `46/11`, and the mutual projection product is `24/11`. Equal residual vectors of square `2−46/11=−24/11` therefore give mutual pairing zero. This rational calculation is checked directly from `J`; exact equality and global validity of the physical anchor recipes still require the constructor's separate verification.

Each pair's first member has no other compatible proposal in the numerical graph. Consequently none can be enlarged to a bank clique. Proposals 44, 45 and 63 have no compatible partner at all. The sole three-member maximal clique is `[100,805,810]`, with all three pairings dotted. The graph otherwise has 8,195 two-member maximal cliques and 171 isolated vertices.

Each right pair has 20 common old finite-angle neighbors. Of their 190 subsets of size 18, 37 avoid an old dotted or parallel pair; every resulting 20-root candidate Gram fails the numerical positive-definiteness test. Thus none of the three finite-contact pairs produces a new finite corner in this screen. This does **not** classify vertices incident to only one proposal. If the three right pairs are exactly accepted, the common-center acute-projection argument would give actual right ridges, but no such exact acceptance is asserted here.

There is also a conditional exact explanation that excludes positive-semidefinite 20-root supports. The common neighbors are S61 plus walls 65 and 95. Including either anchor gives a three-node star with two mutually right cap nodes, both meeting that anchor with the same label at least 5. Its Schur value is `2−4cos(π/m)^2 < 0`. Omitting both anchors forces all 18 S61 roots, but those together with either cap already have negative residual square `−24/11`. This also precludes a bounded ideal vertex containing both caps, **provided** the physical recipes and complete neighbor lists are exactly certified. It is not an exact verification of the numerical bank.

## A gap for the entire relaxed bank

The old normalized normals are reconstructed from 100-digit physical input. Cholesky coordinates on the positive support span, together with the Lorentz residual frame, give a well-conditioned chart with time coordinate one. The numerical polyhedron has a strict timelike interior point of square approximately `−0.5444743`, with minimum normalized inequality slack approximately `0.04270755`.

Across 209 deterministic linear-programming objectives, a strict spacelike proposal is found with

`square ≈ 5.9617546381`,

and minimum normalized slack approximately `4.2707552×10⁻⁷` against all 980 inequalities. The JSON saves both its chart coordinates and a physical-coordinate proposal. This point was made strict by moving a small distance toward the common interior point. All caps were included, ignoring their mutual incompatibility, so it is a necessary diagnostic for every compatible subset of the saved bank. No exact sign certificate for these 980 inequalities is claimed.

The bank's original limitations remain: only selected one- and two-node support profiles, eight nearest anchors at each positive relaxation, and anchor labels 2–6 were tested. In particular, 59 profiles with spacelike residuals were omitted from that causal optimization and are outside this selection result. Other anchors, residual profiles, new walls or changed retained cores remain possible routes.

Generation and numerical replay pass. Replay compares stable discrete outcomes and rechecks the strict numerical gap; it does not demand bitwise equality of floating-point LP coordinates.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/fresh_h20_nonright_s61_20261014_selection.py --verify
```

No exact cap, finite-volume example or all-parameter impossibility theorem is produced by this diagnostic.
