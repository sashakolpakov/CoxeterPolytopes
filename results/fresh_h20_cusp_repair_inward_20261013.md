# Independent inward-cut audit and a golden-field block repair

The inward replacement method correctly preserves a nonempty finite-volume region. Independent exact replay passed for the saved 12-step construction. Its remaining 7,536 pair conflicts do not yield a Coxeter example. A new simultaneous two-wall repair gives a quartic-field plane with target angle π/5, but an untouched genuine bad ridge proves this particular repaired polytope is also non-Coxeter.

All wall indices are zero-based in the original 580-wall six-tile archive.

## Containment and redundancy

If `v=r+t s` with `t<0`, and `s` is retained, then `r=v+(−t)s`. Thus `v≤0,s≤0` imply the deleted old inequality `r≤0`. If a strict timelike center survives, the new region is full-dimensional and contained in the previous finite-volume region. A later inequality can be removed whenever its normal has an exact nonnegative expansion in the currently retained normals. Sequential removal remains valid by induction; later removals may substitute for terms used in earlier certificates.

The independent replay imports no geometry or arithmetic code from the construction. It anchors the first 580 normals and the center to their authoritative archives, then checks each field identity, coefficient sign, new positive norm, center sign, section equation, fresh index, and declared parent angle. It recomputes the final pair-conflict count.

The archived path has 12 cuts and no successful redundancy removals. It retains 580 inequalities and a rank-21 set of rational normals. Its first cut is `3r0−√2 r74`, of squared norm `48+24√2` and center pairing `−48+10√2`, making angle π/8 with `r74`. Later cuts alternate orthogonal projections toward `r2` and `r3`; because those walls are nonorthogonal, repeating these projections does not achieve simultaneous orthogonality in finitely many steps. This is a limitation of that greedy path, not an obstruction to other inward constructions.

## Simultaneous block repair

Put `a=r0`, `b=r74`, `p=r2`, `q=r3`. Their Gram is

\[
\begin{pmatrix}
4&-4&0&0\\
-4&6&-2&0\\
0&-2&2&-1\\
0&0&-1&2
\end{pmatrix}.
\]

Its determinant is −8. Only the `{p,q}` block is positive definite; the full four-root block must not be treated as elliptic.

The exact orthogonal projection of `b` away from `p,q` is

\[
b'=b+\tfrac43p+\tfrac23q,
\quad (b',b')=\tfrac{10}{3},\quad(a,b')=-4.
\]

Simply projecting the first π/8 cut gives `3a−√2 b'`. It preserves the center and containment but loses its target angle with `b`. Instead consider `v=a+t b'`, where `t<0`. This is orthogonal to both `p,q`, and

\[
\cos^2\angle(v,b)
=\frac{(-4+10t/3)^2}{6(4-8t+10t^2/3)}
=\frac59+\frac{4}{9(4-8t+10t^2/3)}
\in(5/9,2/3).
\]

The only possible Coxeter angle in this interval is π/5. Its unique negative solution is the real root

\[
1475t^4-7080t^3+10476t^2-4752t-324=0,
\qquad t\approx-0.0599372348773694.
\]

The exact replay proves irreducibility, isolates the root first in `(−3/50,−59/1000)`, then refines that interval by 160 rational bisections. With the integral coefficient convention

\[
V=3a+t(3b+4p+2q),
\]

its squared norm is `36−72t+30t²>0`, and its center pairing is `−48−42t<0`. The exact containment identity is

\[
a=\tfrac13V+(-t)b+(-4t/3)p+(-2t/3)q.
\]

All coefficients are positive. Replacing wall 0 with `V`, while keeping the other 579 old inequalities, therefore gives another nonempty finite-volume polytope. The new hyperplane is a genuine facet: at a relative-interior point of the old wall-0 facet, the other three terms are strictly negative, hence `V>0`, whereas the center has `V<0`.

The new-versus-old exact pair screen finds 361 compatible pairs, 18 positive pairings, and 200 other intersecting angles in the declared finite alphabet `{2,3,4,5,6,10}`. No redundancy pruning was attempted for this quartic candidate, and these raw hyperplane-pair counts alone are not treated as actual ridge counts.

There is a stronger, independent rejection. The old genuine ridge on walls `23,100` has Gram `[[6,−2],[−2,6]]`, hence squared cosine `1/9`. Its saved relative-interior point has square `−762753024`, satisfies every other old inequality strictly, and pairs with `V` as `−209088(1+t)<0`. Thus that actual bad ridge survives the new cut. The resulting polytope is definitely not Coxeter, irrespective of redundant-inequality questions.

## Arithmetic scope

Coordinates in `Q(√2)` do not by themselves imply arithmeticity. A useful conditional test for a future successful candidate is: if it retains a π/8 pair and rational facet normals spanning the full Lorentzian section, the cyclic product `4cos²(π/8)=2+√2` forces a nontrivial field embedding, while the rational Lorentzian subspace remains indefinite under that embedding. This contradicts the conjugate positive-semidefiniteness condition in [Vinberg's criterion, as stated by Bogachev–Douba, Theorem 2.1](https://msp.org/agt/2025/25-4/agt-v25-n4-p14-s.pdf). This is a conditional nonarithmetic test, not an arithmetic classification of the current non-Coxeter regions.

## Replays

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_cusp_repair_inward_audit_20261013.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_cusp_repair_golden_block_20261013.py --verify
```

Both passed. The corresponding exact records are `data/enumerations/fresh_h20_cusp_repair_inward_audit_20261013.json` and `data/enumerations/fresh_h20_cusp_repair_golden_block_20261013.json`.
