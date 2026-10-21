# H20 sections built around an entire Coxeter cusp

This round constructs **35 exact noncompact finite-volume H20 sections with a complete Euclidean Coxeter cusp link**. Eight also have a certified, deliberately prescribed global ridge angle. Every section has a different actual forbidden ridge, so none is a Coxeter polytope.

The mechanism is different from the earlier isolated angle repairs. It preserves the complete rank-19 cusp diagram while moving the section. Moreover, **27 irrational section normals are proved to require at least five original wall normals**, placing them outside the previous searches with at most four-wall support. No general exclusion of higher-support sections is claimed.

## Avoiding a completed route

The suggested chamber-symmetry fold was already covered by the repository's [complete cyclic symmetry census](borcherds_char2_all_cyclic_element_completion.md): all 28 element conjugacy classes of the 241,920-element chamber group were examined, and the largest nonidentity fixed hyperbolic space is H14. A nontrivial subgroup has fixed space contained in that of any nonidentity element. Thus there is no H20 fixed-space fold to test anew.

The construction below instead uses cusps of the original 210-wall H21 chamber, the classical example described in [Borcherds, *Automorphism groups of Lorentzian lattices*, Example 5](https://math.berkeley.edu/~reb/papers/al/al.pdf). It uses the frozen coordinate packet rather than importing a new chamber from the literature.

## A whole cusp survives the deformation

Let `u` be an ambient ideal vertex, so `u²=0` and every incident wall is orthogonal to `u`. Its incident-wall Gram decomposes into affine Coxeter components. Choose a bridge edge `a—b` in one component, whose positive kernel marks are `k_i`, and solve

```
(e₀,r_a)=k_b,    (e₀,r_b)=−k_a,
(e₀,r_i)=0 for the other roots of that component,
e₀ lies in the component's root span,    (e₀,h)=0.
```

The kernel relation makes the equations consistent. The last condition chooses a unique representative modulo `u`. Now move the hyperplane using

```
e(t)=e₀+t u,        e(t)²=e₀²=q>0,        (e(t),u)=0.
```

For every incident wall, `(e(t),r_i)` is constant. Therefore its projected Gram

```
H_ij=(r_i,r_j)−(e(t),r_i)(e(t),r_j)/q
```

is **independent of t**. Some bridge cuts split an affine component into two affine Coxeter components. For example, one selected cut changes an affine F4 component to affine A1 and G2; the two affine E8 components stay unchanged. Another changes an affine D7 component to two affine B3 components. The search retains a cut only after exact verification of all local angles, rank 19, positive affine kernel vectors and the corresponding positive multiples of `u`.

The coefficient of `u` in each positive kernel sum varies with t. Write the scaled projected normals as `s_i=q r_i−(e,r_i)e`; then each component satisfies

```
Σ k_i s_i = μ(t)u,       k_i>0.
```

If `μ(t)≤0`, no strict interior point can exist: the left side pairs negatively with such a point, while the right side pairs nonnegatively because `u` is future null. If all these coefficients are positive, the affine components give a bounded Euclidean cusp link. A separately certified strict timelike point verifies that the global H20 section exists. Thus the full cusp, not merely its formal Gram, is checked.

## Finite experiment and exact outcomes

Five previously archived ideal vectors also lie in the original chamber. Their affine diagrams supply 87 bridge cuts; 26 distinct rational normal planes pass the local Coxeter-link gate. This is a census for those five literal cusps, not a classification of all ambient cusps.

A numerical scout chooses eight of these pencils and one global target pair per pencil. The producer then isolates **every real root** of the exact target-angle equations for labels `2,3,4,5,6`. Numerical suggestions never establish geometry or exclude candidates.

| Exact gate | Count |
|---|---:|
| Selected cusp pencils | 8 |
| Real roots of the declared angle equations | 46 |
| Finite-volume sections with a full Coxeter cusp | 35 |
| Of these, with the targeted actual ridge certified | 8 |
| Sections rejected by another actual forbidden ridge | 35 |
| No interior, proved by a cusp-cone identity | 5 |
| Target pair lacks positive norms or an acute orientation | 5 |
| Wrong label-five conjugate | 1 |

The 35 section parameters comprise eight rational, thirteen quadratic and fourteen quartic cases. The other 27 algebraic target pairs are not certified as actual ridges, so no global angle-repair claim is made for them. Their separate actual forbidden-ridge witnesses still decide that the sections are non-Coxeter. There is no surviving or unresolved candidate in this declared list.

## A concrete irrational example

Use archived cusp index 1 and the bridge at ambient walls 5 and 6. The central normal has square 4, is orthogonal to the cusp vector, and satisfies `(h,u)=−12`. Set

```
t = sqrt(10)/5 − 1/2,
t²+t−3/20 = 0,
e=e₀+t u.
```

The complete cusp link is the product of affine `A11`, `B3`, `B3` and `A2` factors, of total rank 19. Actual ridge 4/6 has angle `π/3`. But actual ridge 1/23 has scaled Gram

```
[[8, −4], [−4, 12t+32/5]],
```

so its squared cosine is

```
5/(30t+16) = 5/(1+6sqrt(10)),
```

strictly between `1/4` and `1/2`. There is no Coxeter angle between `π/4` and `π/3`. The archive contains exact timelike points certifying both ridges and every other wall inequality. This section is a verified finite-volume non-Coxeter example with a genuine Coxeter cusp; its normal cannot be expressed using any four original walls.

## Why this exceeds four-wall support

Let `E=span_Q(e₀,u)`. If a rational span of four original walls contains `e₀+t u` for irrational t, it contains both rational vectors `e₀,u`. Consequently the four wall images in the quotient by E have rank at most two.

The independent checker tests this necessary condition modulo 1009, explicitly verifying that E retains rank two and that its rational coordinates have denominators invertible modulo this prime. If `E⊂span_Q(W₄)`, the combined rational matrix `[E;W₄]` has rank at most four. Reduction modulo the prime cannot increase this rank, so the four wall images modulo E have rank at most two. This argument does not assume that the rank of W₄ is preserved, or that containment itself survives reduction. The checker groups all projective quotient rays, both signs included, into their two-dimensional spans. Zero images, rank-one cases and repeated rays are also included. Thus every possible rational containment appears in the finite candidate list.

For seven planes no four-wall set survives the modular filter. For the eighth exactly one set survives, and an exact rational rank calculation rejects containment. Thus none of the eight normal planes lies in any four-wall span. The rational-root theorem proves the 27 relevant parameters irrational, completing their minimum-support bound of five. **The eight rational specializations are not assigned this stronger bound.**

## Geometry and independent replay

A strict timelike point in `e⊥` proves that each accepted section has full dimension. Finite volume is inherited from the ambient finite-vertex Klein polytope: its hyperplane section has only finite vertices and inherited ideal vertices, and a finite simplex triangulation has finite hyperbolic volume. The positive affine-kernel certificates show that `u` is an actual ideal vertex with compact horospherical link, so these sections are noncompact.

Every forbidden ridge has its own timelike point, a positive definite projected two-normal Gram, and all remaining wall inequalities verified. Additional zero inequalities, where present, are certified nonnegative combinations of the two ridge normals. Raw pair-Gram conflicts are not substituted for actual incidence.

- [Producer](../scripts/fresh_h20_new_route_20261013.py) and [exact archive](../data/enumerations/fresh_h20_new_route_20261013.json).
- [Independent Fraction verifier](../scripts/fresh_h20_new_route_20261013_verify.py) and [PASS record](../data/enumerations/fresh_h20_new_route_20261013_verify.json).

```sh
python scripts/fresh_h20_new_route_20261013.py --limit 8
python scripts/fresh_h20_new_route_20261013.py --replay-pencils
python scripts/fresh_h20_new_route_20261013_verify.py
```

The second command retains the eight literal pencils and repeats their exact root and geometry calculations. The independent replay uses only the standard library and earlier standalone verification helpers. It imports no producer, SymPy, NumPy or LP solver. It checks all 46 roots by Sturm sequences, all cusp identities, all 35 geometric rejections, the support exclusions, and hostile mutations of a ridge point and cusp kernel. Both scripts compile. Other agents reviewed the modular support argument, including preservation of E's rank, possible drops in the wall matrix's rank, and the irrational-parameter condition.

The unresolved broader construction problem is whether another cusp-preserving normal plane, another target system, or a different assembly yields a Coxeter H20 polytope. This batch supplies a way to enforce an entire valid cusp during that search; it does not exclude those possibilities.
