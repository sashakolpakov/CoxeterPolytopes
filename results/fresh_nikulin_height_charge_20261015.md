# Near-balanced unit charge with a shared defect budget

**Research proof candidate, pending independent review:** every finite-sided finite-volume real hyperbolic Coxeter polytope has dimension at most **165**. The new step counts the few factors where unit charging can fail. Each such factor consumes at least 3/2 of a shared budget of four, so at most two can occur. This improves the accounting in the frozen, reviewed [166 candidate](fresh_nikulin_coupled_charge_20261015_dimension166.md); it does not alter that package or assert a new uniform local charge coefficient.

The [checker](../scripts/fresh_nikulin_height_charge_20261015.py) and [certificate](fresh_nikulin_height_charge_20261015.json) pass generation and saved replay. They import no producer arithmetic. This is producer verification; independent review remains necessary. No publication-priority, current-best, or formal-verification claim is made.

## Local variables and the aggregate statement

Use exactly the preset omissions and geometric hypotheses of the [reviewed four-thirds proof](fresh_nikulin_unit_charge_20261015_review.md) and [common-height proof](fresh_nikulin_coupled_charge_20261015_review.md). An affine factor has rank r, Gram G with root squares two, positive kernel δ normalized by min δ=1, and omissions i,j. Its two contact rows p,q are nonnegative and nonzero; each nonzero entry is at least one. Each endpoint's retained-support cost, summed across all factors, is at most two. Put

\[
\sigma=\delta\cdot p,\quad\tau=\delta\cdot q,\quad
M=\max(\sigma/\tau,\tau/\sigma),\qquad
C=p_{S_i}^{T}G(S_i)^{-1}p_{S_i}+q_{S_j}^{T}G(S_j)^{-1}q_{S_j},
\]
\[
W=\min(\sigma,\tau)^2E,
\qquad T=2r/9\quad\text{or}\quad2(r-1)/9
\]

in the unique factor avoiding the narrow-facet node. The ratio M is common to all factors. The same omission choices are used in both height ranges below.

For **M≤31/30**, every factor has W+C≥T, with the following controlled exception: a failing factor must be one of explicitly listed E7/E8 support patterns, must have C≥3/2, and always satisfies

\[
W+C\ge T-2/9. \tag{1}
\]

Because the shared total cost is at most four, three failing factors are impossible. Thus, writing N=Σr=n−1 and summing all targets,

\[
\sum W\ge\frac{2N}{9}-\frac29-4-2\frac29
=\frac{2N-42}{9}. \tag{2}
\]

This does not contradict the exact balanced E8 counterexample to unit charge. That example is one of the expensive exceptions and remains admissible in the local domain.

## The height estimate

If σ≥h and τ≥k, then the ratio bound B=31/30 gives

\[
\min(\sigma,\tau)^2\ge
H(h,k):=\max\{\min(h,k)^2,\max(h,k)^2/B^2\}. \tag{3}
\]

All arguments below permit arbitrary real contact magnitudes with the zero-or-at-least-one gap. No finite-label cutoff is used.

## Complete exceptional support coverage

The checker reconstructs the literal affine A1, E6, E7, E8, F4 and G2 Grams, their positive kernels, every retained inverse and all 40 ordinary/avoid-node records. These records exactly match the frozen preset choices.

If C≥T, unit charging is immediate from W≥0. Otherwise every retained nonzero coordinate has inverse diagonal below T. All inverse entries are nonnegative, so setting each retained nonzero coordinate to one gives a rigorous support-cost lower bound, including every same-row cross term. The checker enumerates **every subset** of the cheap-coordinate list. No support-size cutoff is imposed. Precisely 71 singleton and 16 multiple-contact occurrences survive C<T; none has more than two retained coordinates.

For a singleton, orient the rows as p=a e_i+b e_l and q=c e_j, with b,c≥1 and a=0 or a≥1. Let R_uv be the normalized point-potential resistance and put

\[
E(u)=(1-u)R_{ij}+uR_{lj}-u(1-u)R_{il},\qquad0\le u\le1.
\]

When a=0, the lower bound is H(δ_l,δ_j)R_lj+c_il. When a≥1, it is H(δ_i+δ_l,δ_j) min E+c_il. The convex quadratic minimum is evaluated exactly at its endpoints and any interior stationary point. Both lower bounds meet T in **all 71 occurrences**, proving unit charging for every singleton.

One multiple-contact pattern is also harmless: ordinary E7 with omissions (3,6), retained nodes 4 and 7 on the respective rows. The cost is at least one and

\[
E(u,v)=\tfrac12(1+u^2+v^2).
\]

If both omitted coordinates vanish, u=v=1, so W+C≥5/2. If either omitted coordinate is nonzero, one height is at least three, and (3) gives W+C≥9/(2B²)+1. Both bounds exceed T=14/9.

The remaining 15 occurrences need no energy estimate beyond W≥0:

| Record | Support occurrences | Minimum cost | Maximum possible defect T−W−C |
|---|---:|---:|---:|
| Ordinary E7 | 4 | 3/2 | 1/18 |
| Ordinary E8 | 4 | 14/9 | 2/9 |
| Avoid-node E8 with omissions (3,7) | 7 | 3/2 | 1/18 |

The seven avoid-node occurrences are the separate prescribed forbidden-node records, not seven simultaneously used narrow factors. The larger uniform allowance 2/9 is safe for all records. Every possible defect consumes at least 3/2 of actual cost, proving the assertion used in (1) and (2). Zero retained contacts satisfy the inherited point-energy bound, checked again in every record. The coincident-omission A1 case has target zero.

## All classical A ranks

For affine A_r set m=r+1 and choose opposite cycle nodes 0 and h=⌊m/2⌋. The kernel entries are one. The point resistance and retained inverse entries are

\[
R(k)=\frac{k(m-k)}m,\qquad
J^{(i)}_{kl}=\frac{R(k-i)+R(l-i)-R(k-l)}2,
\]

with arguments reduced modulo m. These are the frozen all-rank cycle identities. An explicit Euclidean realization is v_k=e_1+⋯+e_k−(k/m)(e_1+⋯+e_m), with v_0=0. Its inner products are min(k,l)−kl/m=(R(k)+R(l)−R(k−l))/2; the checker verifies this identity for every used point pair.

For m=2,…,17 the checker lists every retained coordinate with diagonal at most two, enumerates every subset, and retains precisely the supports whose unit-coordinate cost is at most two. Entrywise nonnegative inverses make this a complete necessary list for actual endpoint rows. For each pair with total minimum cost below T, it checks both possibilities—zero or nonzero—for each omitted coordinate.

Let U,V be the resulting full nonempty supports. The normalized endpoint potentials lie in conv{v_k:k∈U} and conv{v_k:k∈V}. Their heights are at least |U| and |V|. For each of **580 cases**, the certificate supplies an exact rational closest-pair dual witness proving E≥e. Combining it with (3) and the complete minimum support cost gives

\[
W+C\ge H(|U|,|V|)e+C_{\min}\ge T. \tag{4}
\]

The dual witnesses are independently checkable without trusting an optimizer. For the signed Gram H of the two point lists and the two simplex-sum constraints A, the saved nonnegative barycentric vector x and multipliers λ satisfy Ax=(1,1), Hx+Aᵀλ≥0, and complementary slackness. Positive semidefiniteness of H is supplied by the explicit Euclidean point-potential realization. These conditions prove the global quadratic minimum. Singular minimizer sets do not cause a completeness assumption: one checked dual witness suffices for each complete support case. No floating-point solver is used.

For all larger cycles a single analytic inequality replaces the finite calculation. The individual cost budget two excludes every retained node at cyclic distance d≥3, because 3(m−3)/m>2 for m≥10 and the inverse diagonal increases to the antipode. For d≤2, with ε=m mod 2, the frozen linear-loss identity gives

\[
c=\frac{d(m-d)}m,\quad
\ell\le\frac{(m+\varepsilon)d}m,
\quad\frac\ell c\le\frac{m+\varepsilon}{m-2}.
\]

Since C≤4 and W≥E, it follows that

\[
W+C\ge R(h)-\frac{4(2+\varepsilon)}{m-2}. \tag{5}
\]

For even m≥18 the excess of the right side over T is

\[
\frac{m^2+6m-304}{36(m-2)}>0.
\]

For odd m≥19 it is

\[
\frac{m^3+6m^2-457m+18}{36m(m-2)}>0.
\]

Substitution m=18+2t or 19+2t yields numerators with strictly positive constant terms and nonnegative remaining coefficients for every t≥0. Both symbolic substitutions are checked exactly. This proves all ranks; it is not an extrapolation from the 580 cases. Rotating the same opposite omission pair avoids a prescribed node and proves the stronger ordinary target whenever m≥3. The remaining affine A1 avoid-node case is already immediate.

## Classical B, C and D

Use the complete 24 prefix-coordinate classes from the independently reviewed four-thirds proof. Re-evaluate their thresholds at C<T, rather than T/(4/3). All nondistant contacts have ℓ≤c, so the linear-loss identity directly gives unit charging. The other contacts have cost r/4, r/2, (r−1)/4 or (r−1)/2, which exceeds T=2r/9 or 2(r−1)/9 and is impossible in this range. Every index-class substitution is checked symbolically with independent nonnegative rank and index parameters.

The same point resistances exceed the respective targets. Modified B2 has minimum positive retained cost 1/2>T=2/9, and zero-cost point energy one. Thus every B/C/D rank satisfies W+C≥T without using height balance. The proof inherits the previously reviewed completeness of these coordinate classes; the new strict threshold comparisons are recomputed.

## Global split and dimension 165

The reviewed geometric chain supplies actual framing walls for a Nikulin narrow facet, one common M across the affine product, total cost at most four, the unique avoid-node loss 2/9, and the strict norm-two pairing bound −(x,y)<14. These are inherited hypotheses, not consequences of the local support calculations.

In the range M≤31/30, equation (2) and the unchanged Lorentz identity give

\[
-(x,y)\ge M\frac{N-39}{9}. \tag{6}
\]

For M≥31/30, retain the frozen four-thirds estimate with exactly the same preset omissions:

\[
-(x,y)\ge M\frac{N-43}{9}. \tag{7}
\]

At N=165+t, t≥0, (6) is at least 14+t/9, while (7) is at least

\[
\frac{31}{30}\frac{122+t}{9}
=14+\frac1{135}+\frac{31t}{270}>14.
\]

Both contradict the strict narrow-facet bound. Therefore the ideal-vertex branch excludes n≥166 and gives n≤165. The separately reviewed compact-facet branch still gives n≤143. Their combination is the stated candidate n≤165, with no arithmetic assumption and no restriction to simple cusps.

The improvement comes from an integer count of expensive exceptional factors. A smaller uniform local coefficient alone could not give it: the exact balanced E8 obstruction requires at least 8/7, whose scalar threshold still rounds to 166. The present argument neither assumes that 8/7 is a valid uniform upper coefficient nor changes any frozen omission rule.

```sh
python scripts/fresh_nikulin_height_charge_20261015.py
python scripts/fresh_nikulin_height_charge_20261015.py --verify
```

Generation, saved replay, direct source manifests, syntax and local links pass. The full unbounded geometry remains the explicitly inherited reviewed chain. This new aggregate proof and its complete finite exceptional checks are ready for independent scrutiny; prior 166 and earlier artifacts remain unchanged.
