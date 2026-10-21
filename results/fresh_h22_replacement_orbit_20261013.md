# The complete deleted-wall replacement orbit and a full-plane obstruction

The one-wall replacement extends coherently to its **entire 5,760-element stabilizer orbit**, giving **12,331 actual facets and 5,760 certified finite vertices in H22**. It cuts 5,760 specified inherited cusps. Nevertheless it leaves the same strict escape, and an exact positive-cone identity proves a stronger negative result: **no finite collection of acute-compatible walls selected from any of the 5,760 transported three-dimensional residual planes can close that escape**, even with arbitrary real coefficients, heights, norms and mutual angles.

This closes a substantial but explicitly restricted construction family. It is not a general obstruction to changing the old core, and the 12,331-facet chamber still has infinite volume.

## Preserved core and complete orbit

Start from the [813-facet replacement](fresh_h22_core_replacement_20261013.md). Its retained core consists of 811 of the original 812 incidence-branch facets: old wall 1 is deleted, while all 602 signed parent caps remain. This branch does not retain the later 860-wall construction.

The independently certified old symmetry group has order 241,920. Its full stabilizer of old wall 1 has order **5,760**. Since this wall belongs to one short-wall incidence part, the stabilizer preserves both incidence parts. The verifier explicitly checks that it preserves the entire parent bank: **1,733,760 transported mask memberships**, covering all 280 `A=3` and 21 `A=4` recipes under every stabilizer element.

The old center and the second incidence-center vector are fixed, as follows directly from

```
3h = sum(point walls) + sum(line walls),
21vprime = 5 sum(point walls) + 2 sum(line walls).
```

Consequently this is a symmetry of the retained 811-facet core, not an assumption about the full old group preserving a chosen half of the parent cusps.

For each stabilizer element `g`, transport the exact vectors `p,δ` from the single replacement. They satisfy

```
(gp)²=-2484,  (gδ)²=0,  gp·gδ=-46.
```

At the common parameter `t=350`, put

```
u_g=-gp+350gδ,       U=29716,
r_g,±=(u_g/sqrt(U), ±1/sqrt(2)).
```

The checker certifies 5,760 distinct old cusp rays, cap projections and finite-vertex vectors. The resulting 11,520 signed caps have **236,160 right and 75,456,000 dotted new unordered relations**, including relations with the retained core. Every wall has square 2 and remains strict at the common timelike center. Thus all 12,331 inequalities define actual facets.

## Exact compatibility for a whole parameter interval

The entire orbit is coherent for **every common real `t>=350`**, not just one parameter value. For the representative and each transported direction define

```
a_g = p·gp,
b_g = p·gδ + δ·gp,
c_g = δ·gδ,
U(t)=92t-2484.
```

For a nonidentity cusp image, `c_g<0`. The larger, same-height cap pairing is dotted exactly when

```
f_g(t)=a_g-b_g*t+c_g*t²+3U(t)<0.
```

At 350 the largest such polynomial value is **−1,033,910**, and its largest derivative is **−7,664**. Negative leading coefficients prove the entire tail. Opposite-height pairings are smaller; the two caps over the same cusp have mutual pairing zero. Stabilizer invariance reduces all orbit pairs to these representative checks.

The representative cap/core relations also remain valid on that tail. Every nonzero numerator is linear in `t`; the checker proves that its negative value divided by `sqrt(U(t))` increases from an already strictly dotted value. The 20 prescribed right relations persist identically.

Each orbit member has an actual finite vertex

```
(gp+296gδ,0)
```

at `t=350`, with exactly its 20 retained incident old walls and its two new caps active. All other orbit caps are strict there; their largest pairing numerator is **−1,441,562**. The incident Gram is positive definite. The same argument and a derivative check prove these vertices throughout the common parameter tail. No complete vertex or ideal-face census of the resulting chamber is claimed.

## A strict escape for arbitrary parameters

The saved integer vector `E=(100δ₀+h,40)` still has square **1586**, center pairing **−814**, and strict pairing with the entire retained core. At `t=350`, its largest old-coordinate cap numerator over the whole orbit is **−193,616**.

The underlying ideal point `(δ₀,0)` also remains feasible and has an explicitly unbounded horospherical link. Deleted wall 1 was not incident there. Its 24 unchanged old incident walls have mutually orthogonal affine components of ranks **5+7+7+1=20**, with checked positive null relations. All 602 parent caps are strict there (largest pairing **−2**), and all 11,520 replacement caps are strict (largest old numerator **−1800**). The link has dimension 21 but its incident normals span only 20 quotient directions; the pure vertical vector, of square 2, is orthogonal to every incident wall. It supplies an exact one-dimensional lineality direction. Thus this particular old end is not repaired by the full orbit; this is one checked link, not a new global cusp census.

Even before the stronger plane argument below, there is an exact all-parameter exclusion for the fixed-right pair family. In the representative plane, the retained wall 188 has numerator `192-2t`, so acuteness requires `t>=96`. Reversing the old projection is impossible: retained wall 5 would require `t<=4`, while positive old projection square requires `t>27`. At 96, every transported escape numerator is negative and

```
N_g(96)² - 3200U(96) >= 851825424 > 0.
```

The smallest derivative of this squared margin is **37,860,944**, and its quadratic coefficient is positive. Therefore the same escape is strict for all independently chosen admissible real parameters, in either vertical sign. No mutual compatibility assumption is needed for this negative conclusion.

## Stronger obstruction: the entire residual planes

For each `g`, consider arbitrary real normals in

```
L_g=span{gp,gδ,vertical}.
```

Three retained walls have the following pairing functionals in these coordinates:

```
g*old5:       (-4,   -1, 0),
g*old188:     (-192, -2, 0),
g*parent229:  (-423, -8, 2).
```

Here `parent229` means parent recipe 19 at positive height, in the original 812 indexing. Every transported parent is explicitly located in the retained bank. The functional matrix has determinant **−368**, so these three functionals are independent.

Let `eP=(E,gp)` and `eD=(E,gδ)`. Then on the whole plane `L_g`, exactly

```
(E,r) = α_g (g*old5,r) + β_g (g*old188,r) + 40 (g*parent229,r),

β_g=(-eP+4eD-15640)/184,
α_g=-eD-320-2β_g.
```

All 5,760 identities are checked as integer equalities. Their coefficients satisfy

```
184α_g >= 31464 > 0,
184β_g >= 14260 > 0.
```

Every normal acute to the retained core therefore pairs nonpositively with `E`. Equality would force all three independent functionals to vanish, hence would force the normal itself to be zero. Every nonzero proposed wall in any `L_g` is consequently **strictly negative** at `E`.

This proves that arbitrary finite additions from the union of these planes cannot produce finite volume while retaining the 811-facet core. It covers changed cap norms, free real heights, changed mutual angles and unbounded finite-angle labels; it is stronger than a parameter sweep. It does **not** cover linear combinations mixing different planes or normals released from additional old-wall orthogonality constraints. Those are the concrete ways to leave the proved obstruction.

## Reproduction and dependencies

- [Exact orbit and plane checker](../scripts/fresh_h22_replacement_orbit_20261013.py)
- [Frozen certificate, including every positive-cone identity](fresh_h22_replacement_orbit_20261013.json)
- [Underlying one-wall replacement and independent arithmetic review](fresh_h22_core_replacement_20261013.md)
- [Complete old corner/group census](fresh_h21_source_corner_walk_20261013.md)
- [Independent source review and fresh replay of the orbit, plane, and cusp proofs](fresh_h22_replacement_global_20261013_review.md)

```
python scripts/fresh_h22_replacement_orbit_20261013.py --verify
```

The producer and fresh replay pass. Integer array operations have pre-evaluation overflow bounds, and no floating-point predicates are used. The complete old group and its geometric action are a hash-bound independently verified dependency; new parent-bank invariance, cap transport, all orbit pair reductions, actual finite vertices and every escape identity are checked here. Source-level arithmetic reconstruction is shared with the earlier standalone replacement checker, so this is not being labeled a second independent arithmetic implementation.
