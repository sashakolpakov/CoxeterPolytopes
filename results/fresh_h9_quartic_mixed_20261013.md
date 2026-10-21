# Mixed F4 roots and a finite endpoint in the conductor-sixteen H9 host

No compact or finite-volume H9 example was found. This round constructs a coherent **174-wall H9 Coxeter cone**, including four explicit Pell mirrors that close a previously escaping edge at an actual finite vertex. The complementary H5 face has 70 walls. Both displayed cones still have exact strict spacelike escapes, so both have infinite volume. The full arithmetic reflection group remains unresolved.

The host is the [previously introduced conductor-sixteen lattice](fresh_h9_glued_host_20261013.md),

```
T = 2 cos(pi/8),       T^4 - 4T^2 + 2 = 0,
alpha = T - 1,        L = <-alpha> + F4 + F4 + <1>.
```

Here each F4 block has diagonal 2 and successive off-diagonal entries `-1,-sqrt(2),-1`, with `sqrt(2)=T^2-2`. This round does not claim a new lattice or group. It searches directions omitted by the earlier one-block construction and then uses finite-vertex conditions to construct a large-height mirror directly. The complete necessary primitive norm pool remains `1, 2+T, T^2, 2T^2-T^3, 2`. Every retained root is integral, primitive, reflective, totally positive in norm, and passes every exact Coxeter pair gate. Intersecting pairs use labels `2,3,4,6,8,16`; an acute pair alone is insufficient.

| Construction | Explicit search scope | H9 walls | H5 walls | Exact geometry |
|---|---|---:|---:|---|
| [Full two-block census](fresh_h9_quartic_mixed_20261013_h1000.json) | All dominant roots through height 1,000 | 62 | 26 | Strict spacelike escapes |
| [Targeted mixed slice](fresh_h9_quartic_mixed_20261013_target.json) | Time, omega0 in one block, omega1 in the other; height 100,000 | 86 | — | Strict spacelike escape |
| [Neighboring-weight slices](fresh_h9_quartic_mixed_20261013_face_target.json) | Three two-covector H5 slices and their symmetries; height 100,000 | 170 | 68 | Strict spacelike escapes |
| [Relative-Pell endpoint](fresh_h9_quartic_mixed_20261013_endpoint.json) | Four explicit symmetry-related roots | **174** | **70** | Actual finite endpoint; other strict escapes remain |

The complete height-1,000 census has 61 time/norm controls and 58 primitive dominant candidates. Twelve have both F4 blocks active. Forty new walls enlarge the previous 22-wall cone to 62 walls. The [height-300 calibration](fresh_h9_quartic_mixed_20261013_h300.json) recovers all previous 22 walls and finds no primitive two-block candidate there.

The speed improvement is an exact unit change of coordinates before applying the trace ellipsoid bound. A numerical choice selects an integral unit and its integral inverse; multiplication by that unit is then an exact bijection. Every embedding bound is checked exactly. The [separate Fraction census](fresh_h9_quartic_mixed_20261013_independent.json) uses a different unit-selection procedure and spatial enumeration, and reconstructs all 61 controls and all 58 candidates, rather than checking only saved walls.

The targeted mixed slice has norm `-alpha*x^2 + 2*u^2 + 6*v^2`. It requires both spatial covectors to be strictly positive and the scalar to vanish. Its complete 626-control census yields three primitive roots, whose eight-position symmetry orbits give 24 mutually compatible new mirrors. Their heights are approximately 17,162, 21,552, and 43,105. All field coefficients have absolute value at most 87.

The face search treats covector pairs `(0,2),(1,2),(2,3)` in one F4 block, together with diagram reversal and block swap. Both covectors are strictly positive; boundary one-weight planes and the endpoint-only pair `(0,3)` are outside this census. Completing a square in the positive inverse-Gram minor reduces each search to a one-variable alphabet and an exact field-square test. The 1,878 slice controls yield 36 primitive base roots and 84 new compatible symmetry images. These 84 and the preceding 24 mirrors are all mutually compatible. No wall replacement is needed. The [independent replay](fresh_h9_quartic_mixed_20261013_face_target_independent.json) reverses the completed-square variable and tests squares through the quadratic tower over `Q(sqrt(2))`; it verifies the complete census and all 14,365 H9 pair gates.

There is also a narrow [global local obstruction](fresh_h9_quartic_mixed_20261013_local_obstruction.json). A normal confined to time and a single inner F4 fundamental weight would satisfy `6*u^2-alpha*x^2=q`. The two odd dyadic-valuation norm classes are impossible at the prime above 2. For the remaining three classes, all 81 residues modulo 3 fail `-alpha*x^2=q`. Thus this binary plane contains no primitive reflective normal at any height. It does not exclude neighboring covectors, scalar directions, or the full host.

The productive alternative is to fix the finite couplings at an edge of type `A2+A1+A1`. For a root norm `q`, let `(a,b,c,z)` be its covectors against the two A2 walls, the other F4 wall, and the scalar coordinate. Finite-angle integrality and positive definiteness require

```
a,b,c in {integral sqrt(q*r/2)},
z     in {integral sqrt(q*r/4)},
r in {0,1,2,3,2+sqrt(2),2+T},

q_perp = q - (2/3)(a^2+a*b+b^2) - c^2/2 - z^2 > 0
```

at every real embedding. The complete small profile census has 30 possibilities: respectively 4, 1, 5, 2, and 18 for the five norm classes. These are necessary endpoint profiles, not a classification of their integral roots.

One norm-two profile, `(a,b,c,z)=(0,0,1,0)`, is constructive. With free inner covector `t`, the norm equation becomes `6t(t+1)=alpha*x^2`. Set `x=3s` and `t=(y-1)/2`. The relative Pell equation is

```
y^2 - 6*alpha*s^2 = 1,

y = -17023554335 - 9213081216*T
    + 29061025044*T^2 + 15727713420*T^3,
s = -7548113028 - 4085009334*T
    + 12885434936*T^2 + 6973547136*T^3.
```

Both have positive distinguished embeddings, and `t` is integral. In time/F4/F4/scalar coordinates the root is `(3s, F4_inverse*(0,0,t,1), 0,0,0,0, 0)`. Its norm is exactly 2 and its reflection is integral. PARI supplied the relative unit; [the discovery script](../scripts/fresh_h9_quartic_mixed_20261013_endpoint_pell.gp) is archived. The exact certificate and independent verifier require no PARI. Its height is approximately `2.0269e22`; this is an explicitly solved endpoint equation, not an enlarged uniform height census.

All four block-swap and diagram-reversal images retain all 170 old walls and one another. For the displayed root, the point

```
v = (1, [alpha*x/(3*y)]*omega2, 0,0,0,0, 0)
```

has square `-alpha/y^2 < 0`, lies on the new mirror, and satisfies all 174 inequalities. Exactly nine independent walls meet there. Their Gram matrix is positive definite, with diagram `F4+A2+A2+A1`; the forced H5 endpoint has diagram `A2+A2+A1`. This is an actual finite vertex, rather than merely a positive definite subdiagram whose intersection might lie outside the chamber.

The new mirrors cut the former H5 escape, but a different one survives. The final H5 escape square is `6+7T-3T^3 > 0`, approximately 0.0083813. The H9 escape square remains

```
(-113752 + 68526*T + 201861*T^2 - 111278*T^3)/578 > 0,
```

approximately 0.087903. Both archives include strict spacelike perturbations and a common strict timelike interior. Since the host is definite at the other real embeddings, it is anisotropic over its field; a finite-volume completion would be compact. Nothing here proves that further compatible mirrors cannot complete it.

The final [Fraction-only replay](fresh_h9_quartic_mixed_20261013_endpoint_independent.json) passes the Pell identity, all four normal formulas, all 15,051 H9 pairs and 2,415 H5 pairs, both escapes, every finite-endpoint inequality, and the incident rank-nine positive definite Gram. A forged root is rejected. The construction agent additionally checked the large Pell identity with independent integer polynomial multiplication and reviewed the endpoint formula; this was a separate narrow review, not a repeated full pair census.

```sh
python scripts/fresh_h9_quartic_mixed_20261013_local_obstruction.py --verify
python scripts/fresh_h9_quartic_mixed_20261013_target_independent.py --verify
python scripts/fresh_h9_quartic_mixed_20261013_target_independent.py --face --verify
python scripts/fresh_h9_quartic_mixed_20261013_endpoint_independent.py --verify
```

The first target replay takes about 30 seconds, the full face census about 214 seconds, and the final endpoint replay about 1.3 seconds. `--geometry-only` on the target checker verifies saved geometry without repeating the census. Producers require `PYTHONPATH=src:scripts`. The independent arithmetic verifiers import no producer or repository number-field arithmetic.
