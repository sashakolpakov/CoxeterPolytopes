# A changed conductor-16 host and a certified first H9 edge

**Superseded host status:** the [complete edge-walking continuation](fresh_h9_new_host_next_20261013.md) now proves this full H9 host nonreflective by an exact infinite chamber stabilizer. The bounded computations below remain valid historical artifacts.

The new host has an independently certified first mirror on a previously open edge. Its normal has coefficients at most 108, replacing a much later compatible proposal with coefficients in the trillions. Four symmetry copies give an 18-wall H9 chamber prefix with four actual finite vertices. Exact spacelike escapes remain in both this prefix and its 10-wall H5 face. **No finite-volume polytope or global nonreflectivity theorem has been obtained for this host.**

The authoritative current chamber is [the corrected construction](fresh_h9_new_host_20261013_corrected.json), with [separate geometry replay](fresh_h9_new_host_20261013_corrected_independent.json). The [complete first-edge certificate](fresh_h9_new_host_20261013_first.json) and [independent full census](fresh_h9_new_host_20261013_first_independent.json) justify the first-mirror claim.

## The changed arithmetic form

Let

\[
K=\mathbb Q(T),\quad T=2\cos(\pi/8),\quad T^4-4T^2+2=0,
\qquad \pi=2+T,\qquad \alpha=(T-1)\pi=T^2+T-2.
\]

The lattice is

\[
L=\langle-\alpha\rangle\oplus F_4\oplus F_4\oplus\langle1\rangle,
\]

where each displayed F4 block has diagonal 2 and consecutive off-diagonal entries \(-1,-\sqrt2,-1\). Its determinant is 1. The time coefficient has signatures \((+,-,-,-)\), and \(N_{K/\mathbb Q}(\alpha)=-2\); thus the form has signature \((9,1)\) at the chosen place and is positive definite at all other real places. The one-block version has dimension 5.

This does not repackage the [globally excluded old host](fresh_h9_quartic_next_frontier_20261013.md), whose coefficient was \(T-1\). The determinant ratio is \(\pi\), whose rational norm is 2 and hence cannot be a field square. In even ranks 10 and 6, multiplying a quadratic form by a scalar does not change its determinant squareclass. The respective new and old forms are therefore not scalar-similar over K. This is a comparison with these specified repository forms, not a worldwide novelty claim or a classification of every commensurable lattice.

The previous [exact field structure certificate](fresh_h9_glued_host_20261013_structure.json) supplies the integral basis, class number one, unique dyadic prime, and the fact that totally positive units are squares. The new dual exponent is \(\pi\). Consequently a primitive reflection norm divides \(2\pi\), whose dyadic valuation is 5. The six safe valuation representatives are

\[
1,\quad \pi,\quad T^2,\quad 2/\pi,\quad 2,\quad 2\pi.
\]

The last is impossible primitively: reflection integrality for norm \(2\pi\) forces every spatial coordinate into \(\pi O_K\), because the spatial Gram is unimodular. If the time coordinate were a dyadic unit, its square term would have valuation 1, while the spatial square has valuation at least 2 and the proposed norm has valuation 5. Thus the time coordinate is also divisible by \(\pi\). The complete primitive norm pool is the remaining five representatives, modulo square units.

The full zero-height census has five dominant root orbits in H9 and three in H5. These are precisely the two F4 root orbits per block and the scalar reflection; there is no missing zero-height mixed-block mirror. The [independent audit](fresh_h9_new_host_20261013_independent.json) verifies the structural identities and repeats the complete height-100 census: 11 time controls, 14 H9 walls and 8 H5 walls. Both initial prefixes have exact escapes.

## Finding and certifying the first mirror

At the coordinate origin, three of the five H5 edge directions already have finite endpoints. The two middle-F4 directions remain open and are related by diagram reversal. The chosen edge fixes initial wall indices `[0,2,3,4]`; its direction is the second F4 fundamental weight, of norm 6.

A complete targeted screen through root height 10,000 considered all 30 finite-angle coupling profiles, all five primitive norms and all bounded time coordinates. None of its 391 quadratics had a field square root. The [independent screen](fresh_h9_new_host_20261013_edge_independent.json) uses different trace balancing and coefficient enumeration: 374 discriminants have nonsquare rational norm, and exact number-field polynomial factorization excludes the remaining 17. This bounded failure was not interpreted as a global obstruction.

The subsequent relative-norm calculation is

\[
z^2-6\alpha x^2=6R,
\]

where R is the residual norm of a coupling profile. There are 16 distinct targets. [The PARI transcript](fresh_h9_new_host_20261013_relative_raw.txt) certifies the degree-8 field

\[
X^8-288X^4-1728X^2-2592
\]

with `bnfcertify = 1`, class number 1 and signature `(2,3)`. Complete absolute integral norm representatives, rather than a particular rational norm solution, were obtained for the relevant norms up to 20,736. The first three fundamental units lie in K; the fourth has relative norm one. It therefore generates the relative unit kernel, up to sign. Totally positive unit norm discrepancies are corrected by exact base-field square roots.

[The relative-orbit artifact](fresh_h9_new_host_20261013_relative.json) checks 52 orbit points up to the earlier proposed endpoint, of which 18 are primitive reflective candidates. Four reach the first endpoint. This is a certified-CAS orbit calculation, with that dependency explicit. The final first-edge theorem does **not** depend on its completeness: the endpoint it found is close enough for a separate direct interval enumeration.

That complete interval census includes arbitrary bounded integral pairings to the four edge walls, every necessary norm class, nonsquare residual classes, both signs of all time coordinates, both roots of every quadratic, and exact reflection-integrality and primitivity tests. It checked:

| Quantity | Count |
|---|---:|
| Integral pairing profiles | 35 |
| Quadratic equations | 15,177 |
| Field-square discriminants | 47 |
| Primitive signed solutions | 42 |
| Folded segment mirrors | 25 |
| Mirrors after finite-edge unfolding, excluding the separately counted line-containing mirrors | 35 |
| H5 edge / initial / endpoint stabilizer mirrors | 5 / 25 / 20 |
| H9 edge / initial / endpoint stabilizer mirrors | 29 / 49 / 44 |

Every positive segment intersection occurs at parameter 1. The complete endpoint stabilizer selects the norm \(2/\pi\) normal

\[
\begin{aligned}
r={}&(-14-9T+24T^2+14T^3;\\
&-31-18T+54T^2+30T^3,\\
&-63-36T+108T^2+60T^3,\\
&-60-32T+102T^2+56T^3,\\
&-30-16T+51T^2+28T^3;\ 0).
\end{aligned}
\]

The displayed vector is in the H5 coordinates; the H9 lift inserts a zero second F4 block. Its root height is about 122,357.343. The endpoint has spherical type \(I_2(16)+A_2+A_1\) in H5, with a further F4 in H9. All other candidate simple normals at the same point generate smaller stabilizers, so merely selecting an arbitrary compatible cap there would miss the complete local chamber.

The H9 lift is proved by a complete second-F4 block census: for each normalized primitive norm q, every dominant spatial vector with square at most q at all embeddings has square 0 or q. A genuinely mixed crossing root would require its nonzero F4 square to be strictly between these at every embedding, which is impossible. Pure second-F4 mirrors contain the whole segment and are included in the stabilizer counts.

The [independent checker](../scripts/fresh_h9_new_host_20261013_first_independent.py) uses separate Fraction field arithmetic, a different exact unit-balancing procedure, direct five-equation elimination, all quadratic-tower roots, and independent finite Weyl unfolding. It passed the complete census, including both endpoint stabilizers. This establishes a first edge of the full reflection chamber, not only a subgroup corner.

## Corrected chamber and next frontier

The earlier four far caps are retained only in [their historical proposal artifact](fresh_h9_new_host_20261013_completion.json). The current [corrected chamber](fresh_h9_new_host_20261013_corrected.json) starts from the original complete prefix and adds the four symmetry copies of the certified first mirror. It has 18 H9 walls and a 10-wall H5 face. All 153 H9 pairs and 45 H5 pairs pass exact Coxeter-angle gates, including the new label 16. Four rank-nine positive definite incident Grams certify actual H9 finite vertices; the two corresponding H5 vertices are also checked against every inequality.

Exact escapes remain. The H5 escape has norm

\[
43-47T-12T^2+14T^3>0,
\]

and the H9 escape has norm

\[
\frac{125999-82423T-36730T^2+24194T^3}{289}>0.
\]

Both have saved strict future spacelike perturbations and strict timelike interior points. Thus these finite prefixes have infinite volume. The full arithmetic chamber remains unresolved.

From the new actual H5 corner `[0,1,3,4,8]`, the four next open supports are

- `[1,3,4,8]`;
- `[0,3,4,8]`;
- `[0,1,4,8]`;
- `[0,1,3,8]`.

Exact lookahead reaches the light cone before any existing wall in each of these directions. The remaining edge `[0,1,3,4]` returns to the original finite corner. This is a live frontier against the current walls; it does not claim those directions lack undiscovered arithmetic mirrors.

The natural chamber-symmetry test did not reproduce the old host's global obstruction. Initial F4 reversal is an integral isometry. The new endpoint's A2 swap is reflection in a norm-6 difference of incident roots and has two nonintegral entries \(\pm\sqrt2/3\); its inverse is the same matrix. Both transformations fix the initial timelike axis, and their product has trace \((4+T^2)/3\). Exchanging the two I2(16) rays would require \(\sqrt\pi\notin K\). [The exact matrix audit](fresh_h9_new_host_20261013_symmetry.json) records this scoped failure, not a theorem excluding other chamber symmetries.

## Reproduction and scope

From the repository root:

```sh
PYTHONPATH=src:. python scripts/fresh_h9_new_host_20261013_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_new_host_20261013_edge_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_new_host_20261013_first_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_new_host_20261013_corrected_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_new_host_20261013_symmetry.py --verify
```

The first-edge verification is a complete finite enumeration and requires no PARI installation. The optional discovery transcript can be reproduced with `gp -q scripts/fresh_h9_new_host_20261013_relative.gp`; the exact subsequent orbit processing is [relative.py](../scripts/fresh_h9_new_host_20261013_relative.py). The bounded edge screen's separate audit uses SymPy for 17 exact number-field irreducibility checks. The corrected geometry checker also rejects a forged normal and a timelike replacement for the escape.

All searches and audits for this checkpoint are complete. No larger uniform height bound was run, and the globally excluded old quartic host was not re-enumerated.
