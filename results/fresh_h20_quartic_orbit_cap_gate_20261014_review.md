# Independent review of the fixed quartic cusp endpoint gate

**PASS.** The [separate checker](../scripts/fresh_h20_quartic_orbit_cap_gate_20261014_review.py) reconstructs the fixed section and re-enumerates the necessary endpoint conditions using SymPy algebraic fields and rational root isolation. It imports neither the producer nor its arithmetic. Its [certificate](fresh_h20_quartic_orbit_cap_gate_20261014_review.json) binds the [primary certificate](fresh_h20_quartic_orbit_cap_gate_20261014.json), source data, complete corner walk and both scripts. Generation and saved replay pass.

The result concerns local record 0 at `u=1/2`, retaining its 22 cusp facets. It excludes finite-volume completion by projected, spacelike normals represented by **integral ambient vectors of norm 2 or 4**. It does not exclude arbitrary real normals, other sections, other ambient norms or deletion of seed facets. The later [rational-normal gate](fresh_h20_quartic_rational_cap_gate_20261014.md) is a separate strengthening.

| Independently checked item | Count |
|---|---:|
| Integral old positive-semidefinite triple cases | 120 |
| Right-contact trials | 240 |
| Right-contact angle passes | 12 |
| Common linear roots for two label-3 equations | 6 |
| Unresolved quadratic common factors | 0 |
| Accepted necessary angle records, including duplicates | 15 |
| Distinct necessary survivor families | 3 |
| Cusp edge germs | 252 |
| Edge germs with no finite or ideal opposite endpoint | 168 |

The retained triangle walls are `[23,50,200]`. Their projected pair angles are all `pi/3`. For a possible endpoint cap, the normalized three-wall determinant is `3/4-x²-y²-xy`, where `x,y` are its nonnegative normalized contact magnitudes. If both contacts are nonright Coxeter angles, `x,y>=1/2`, so positivity forces both labels to be 3 and the determinant to vanish. If one is right, the other label is at most 6. Label 5 is included explicitly by its quadratic squared-cosine identity and the correct root-sign test. Thus the finite alphabet is justified by the endpoint Gram itself; no assumption about the absence of a quadratic subfield is needed here.

An endpoint's old ambient triple Gram equals its projected Gram plus a positive rank-one matrix. It must therefore be positive semidefinite. The checker independently tests every principal minor of every integer triple within the exact Cauchy bounds, and reproduces all 120 cases. It recomputes every polynomial gcd over the quartic field and compares every accepted pairing, determinant and identity solution with the archive.

The final rational identity uses the missing integer pairing and the cusp pairing as two unknowns. Full-rank cases are solved over the rationals. Rank-one cases are retained, not discarded: the resulting Diophantine equation is solvable precisely when the appropriate gcd divides its constant. A solvable nontrivial congruence has negative cusp pairings in its arithmetic progression. After deduplication, only the pair `[23,50]` survives this necessary test: two all-right norm families and one norm-four, two-label-3 family. These are not constructed caps. Neither `[23,200]` nor `[50,200]` has a possible endpoint cap.

The geometric passage to infinite volume is also sound. The checker reconstructs three orthogonal positive affine components of quotient ranks 11, 2 and 6. Every one-wall deletion has a positive-definite Gram; each component's positive weighted normal sum is a positive multiple of the same future null vector `D`. Consequently their horospherical link is a bounded product of three simplices, with `12*3*7=252` genuine edge germs.

Any additional acute facet has cusp pairing at most zero, by the positive relations. If its cusp pairing were zero, all 22 retained pairings would vanish. Their projected normal span is `D-perp` within the section, of dimension 20, so the extra normal would be proportional to `D` and null. This contradicts its required positive square. Positive-proportional copies of retained facets are simply redundant copies, not new walls. Every genuine extra facet therefore has strictly negative cusp pairing and leaves the original link and edge germs intact sufficiently near the cusp.

For each of the 168 blocked germs, a first finite intersection with a new wall, or an ideal endpoint on a new wall, would require precisely the forbidden positive-semidefinite triple. A finite collection of additions thus remains strict at that germ's opposite null endpoint. Continuing slightly beyond the light cone and then perturbing toward the proposed completion's own strict interior point gives a strict spacelike escape. This proves the stated finite-volume obstruction without assuming the original saved center survives arbitrary additions.

```sh
python3 scripts/fresh_h20_quartic_orbit_cap_gate_20261014_review.py
python3 scripts/fresh_h20_quartic_orbit_cap_gate_20261014_review.py --verify
```

The arithmetic is independent of the producer's Fraction polynomial implementation, but both checks use the same frozen geometric inputs. This review does not repeat the full 91-direction classification or the complete source corner walk.
