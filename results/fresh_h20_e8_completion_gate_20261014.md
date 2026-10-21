# All-residual-vector restrictions on E8 contacts

**Several complete contact families cannot cut the frozen 79-wall chamber's strict escape, even when the residual H12 vector is arbitrary. Other families remain viable.** In particular, the two-contact patterns `{3,6}` and `{4,5}` used by the separate constructive search are not excluded. This certificate supplies pruning constraints rather than a noncompletion theorem for the whole retained core.

The [exact checker](../scripts/fresh_h20_e8_completion_gate_20261014.py), [certificate](fresh_h20_e8_completion_gate_20261014.json), and [rational dual recipes](fresh_h20_e8_completion_gate_20261014_recipes.json) pass generation and saved replay. The checker imports no optimizer. A [numerical SOCP scout](fresh_h20_e8_completion_gate_20261014_scout.json) only proposed nonnegative rational wall multipliers; every accepted bound is subsequently proved by exact algebraic inequalities.

The retained input is the [79-wall chamber](fresh_h20_e8_cusp_completion_20261014.json). Its earlier [independent review](fresh_h20_e8_cusp_completion_20261014_review.md) is linked here without changing the previous round's frozen report.

## Contact conventions and conclusions

Number the eight E8 positions from zero, corresponding to retained wall indices

```text
position:  0   1   2   3   4   5   6   7
wall:     12  13  14  15  18  26  33  67.
```

Normalize a proposed wall to norm two and write its E8 pairings as `-c_i`. A nonzero Coxeter contact has `c_i=2cos(pi/m)` for a finite label `m>=3`, or `c_i>=2` for parallel/dotted hyperplanes. In particular `c_i>=1`. A zero coefficient means a right contact. The following statements require all unlisted E8 coefficients to be zero and all retained-wall pairings to be nonpositive.

| E8 contact support | Exact restriction for a cap cutting the saved escape |
|---|---|
| One position: 0, 1, or 7 | Impossible for every finite label and every parallel/dotted magnitude |
| One position: 2 | Only the pi/3 magnitude `c_2=1` remains possible |
| Two positions: `{0,1}`, `{0,2}`, `{1,2}`, `{1,3}`, `{1,7}`, `{2,3}` | Impossible for all independently chosen magnitudes at least one |
| Two positions: `{0,7}` | Impossible when their magnitudes are equal and at least one |

These restrictions cover all real normal vectors with those contacts, without a lattice, coordinate-field, height, or H12 face-support restriction. They are not limited to residual vectors on source vertices or edges. They do not classify the other two-contact patterns, the unequal `{0,7}` case, or higher-support contacts.

Every excluded cap pairs **strictly negatively** with the saved spacelike witness. Thus any finite completion using only these excluded families retains that escape. This is a restriction on that union of contact families, not on all possible completions.

## Why the completed H12 face gives a global bound

Let `omega_i` be the E8 fundamental weights, `I` the inverse Cartan matrix, and

```text
w = sum_i c_i omega_i,    r = -w + y,    y perpendicular to E8.
D = w^2 = c^T I c,       y^2 = 2-D.
```

The fourteen retained facets orthogonal to E8 are exactly the complete finite-volume [H12 source](fresh_h12_a12_source_20261014.md). Their cone has 35 timelike extreme rays and exactly two null rays. The acuteness conditions on these facets place `y` in its closed future causal cone. This conclusion does not assume that the proposed cap preserves the predecessor's chosen strict center.

For a nonzero contact, positivity of the E8 inverse and its minimum diagonal two imply `D>=2`. An all-right E8 contact would instead have `w=0` and `y^2=2`, which is impossible inside this causal cone. The escape's E8-orthogonal projection is itself future timelike and is checked exactly.

For any nonnegative rational multipliers `lambda_j` on retained walls `n_j`, set

```text
X = E - sum_j lambda_j n_j,
Y = projection of X to E8-perpendicular,
kappa = -Y^2,
a_i = -(omega_i,X).
```

Each saved recipe verifies `Y` is future timelike and `kappa>0`. Since `(r,n_j)<=0`, Lorentz Cauchy gives the bound

```text
(r,E) <= a^T c - sqrt(kappa*(c^T I c - 2)).
```

This is an all-parameter bound on every residual vector `y` of the required norm. The numerical optimization was performed over the larger convex set `y^2<=2-D`, but the exact proof above works directly with the norm equality.

## Eleven exact dual certificates

Write `M=kappa I-aa^T`. A positive value of

```text
c^T M c - 2 kappa
```

forces the displayed upper bound to be negative, regardless of the sign of `a^T c`. The checker establishes this over whole unbounded contact domains:

- For positions 1, 0, and 2, it checks a positive diagonal coefficient and the strict inequality at the respective thresholds `c^2=1,2,2`. Monotonicity then covers every larger magnitude.
- For position 7, the linear coefficient itself is strictly negative, so every positive feasible magnitude is excluded.
- For the six independently variable pairs, every entry of the corresponding two-by-two matrix is strictly positive, and its value at `(1,1)` exceeds `2 kappa`. The inequality therefore remains strict on the entire orthant `c_i,c_j>=1`.
- For equal contacts at `{0,7}`, the sum of the four matrix entries is positive and already exceeds `2 kappa` at common magnitude one. Scaling covers all larger common magnitudes. This argument deliberately does not assert the unequal-contact case.

All coefficients, future signs, and strict margins are checked in exact `Q(sqrt(2))(sqrt(79))` arithmetic. The norm equation and Lorentz inequality introduce no unverified radical comparisons.

## The exceptional minimum-weight contact

Position 0 at pi/3 has `D=2`, hence `y^2=0`. A null vector in the source cone is zero or a nonnegative multiple of one of its two null extreme rays: any sum of distinct future causal rays is timelike. Therefore only two one-parameter families remain.

For each ray `nu`, write `r=-omega_0+t nu`. The exact certificate selects an old wall with

```text
offset = (-omega_0,n)>0,   denominator = -(nu,n)>0.
```

Acuteness forces `t>=offset/denominator`, also excluding `t=0`. Since `(nu,E)<0`, the target pairing is at most its value at that lower bound. The checker proves the numerator of that value is strictly negative for both rays, without floating-point division. This closes pi/3 at position 0. All other nonright Coxeter contacts there have magnitude at least `sqrt(2)` and are covered by the continuous dual certificate. Thus no finite-label cutoff is involved.

## Constructive implication and limits

The full residual optimization leaves five possible one-node positions before further contact gates, with position 2 restricted to pi/3. It also leaves many multi-contact patterns. The separate scout found coupled candidates at `{3,6}` and `{4,5}`, consistent with these restrictions. Positive relaxed objective values alone are not certificates of valid Coxeter walls; full pair and geometry checks remain necessary.

This result does not force deletion of the retained core. It identifies contact patterns worth avoiding and shows why restricting a search to source vertex rays can miss useful residual vectors on lower-dimensional source faces.

```sh
python scripts/fresh_h20_e8_completion_gate_20261014.py --verify
```
