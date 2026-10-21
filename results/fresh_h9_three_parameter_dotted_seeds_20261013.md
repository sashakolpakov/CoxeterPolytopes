# Three metric parameters and caps disjoint from seed walls

These searches find no finite-volume candidate. They change both the seed
metric and the allowed cap pairings; all search conclusions are numerical
and restricted to the recorded proposals.

For two groups of eight matched wall pairs in total, let Z be the group
incidence matrix and let M be a symmetric two-by-two matrix. The seed is

```
K = Z M Zᵀ,
G = [[2I, -2I-K], [-2I-K, 2I]].
```

Unlike the earlier equal-diagonal family, M has three independent entries.
Every tested entry exceeds 2 and is at most 60. The negative eigenvalue of
diag(√n) M diag(√n) lies in (−4,0), giving seed rank 10 and signature (9,1).
The two group sizes are inherited from the starting records, 4+4 and 5+3.

A cap has prescribed negative pairings −x and −y with the two sets of seed
walls. Kernel consistency requires x+y=Zt. Put d=x−y, d̄=Zᵀd and D=diag(n).
Its square is

\[
-\tfrac12 t^{\mathsf T}M^{-1}t
+\tfrac18d^{\mathsf T}d
-\tfrac1{32}\bar d^{\mathsf T}(M^{-1}+D/4)^{-1}\bar d.
\]

Polarization gives pairings between two caps. The implementation solves
three equations in the entries of M⁻¹, using their analytic Jacobian, then
reconstructs every prescribed profile orbit at the resulting metric.

| Search | Root attempts | Accepted metric proposals | Whole banks numerically bounded | Compatible selections |
|---|---:|---:|---:|---|
| Two cap norms and their mutual angle | 43,670 | 4,671 | 19 | 12,011 maximal cliques, all with escapes |
| Three cap norms | 149,491 | 36,559 | 55 | 50,061 maximal cliques, all with escapes |
| Three cap norms, including dotted seed pairings | 85,916 | 43,598 | 0 | Every whole bank has an escape |
| Symmetric dotted-pairing curve intersections | 5,815 parameter proposals | 5,815 | 7 | All seven numerical coverage systems unsatisfiable |

The first two rows use magnitudes 0,1,√2,φ, where φ=(1+√5)/2. Their largest
compatible selections have six caps, and the smallest selected maximum
Klein square is approximately 5.828427. The first row contains 4,058
unequal-diagonal metrics and the second 33,238; none of their whole banks
is numerically bounded. These observations are not a theorem about all
unequal-diagonal metrics.

There is an exact motivation for the dotted-pairing experiment. Within one
group, transpose a pair of assignments 0,D in a norm-two cap. The difference
of the two cap normals has square D², so the two caps have pairing
2−D²/2. For a finite seed-pairing alphabet with maximum difference below 2,
this transposition gives an incompatible positive pairing. Taking
D=√(4+2√2) instead makes the cap pair have angle π/4, while D>2 makes its
specified seed pair disjoint. Both dotted searches use magnitudes 0,1,φ,D.

The three-equation dotted search samples 100,000 triples. The symmetric
search uses 300 samples of the earlier norm curves, with 0<b<1 and
b+2<a<20. It can miss tangencies and unsampled crossings. For each of its
seven whole-bank-bounded proposals, a Boolean solver enforces pairwise
compatibility and iteratively adds clauses requiring some selected cap to
cut a numerical spacelike witness. It reaches unsatisfiability after
2–11 witnesses. This excludes those numerical selection systems; it is
not an exact algebraic exclusion.

The root searches keep one starting point per selected equation system,
skip singular initial Jacobians in the three-norm searches, and do not
enumerate every real root. Metric deduplication, pair labels, ranks and
Klein geometry use tolerances. Full banks include incompatible caps;
their geometry is only a necessary screen. No positive exact construction
or complete parameter-space census is claimed.

Artifacts:

- [Two norms and an angle](../data/enumerations/fresh_h9_three_parameter_seed_20261013_scout.json).
- [Three norms](../data/enumerations/fresh_h9_three_norm_seed_20261013_scout.json).
- [Dotted seed pairings](../data/enumerations/fresh_h9_dotted_seed_caps_m4_20261013_scout.json).
- [Symmetric dotted intersections](../data/enumerations/fresh_h9_dotted_symmetric_m4_20261013_scout.json) and [coverage clauses](../data/enumerations/fresh_h9_dotted_symmetric_m4_20261013_cover.json).

Reproduce from the repository root:

```sh
python scripts/fresh_h9_three_parameter_seed_20261013.py
python scripts/fresh_h9_three_parameter_seed_20261013.py --screen
python scripts/fresh_h9_three_norm_seed_20261013.py
python scripts/fresh_h9_three_norm_seed_20261013.py --screen
python scripts/fresh_h9_dotted_seed_caps_20261013.py --m 4
python scripts/fresh_h9_dotted_seed_caps_20261013.py --m 4 --screen
python scripts/fresh_h9_dotted_seed_caps_20261013.py --m 4 --symmetric
python scripts/fresh_h9_dotted_seed_caps_20261013.py --m 4 --symmetric --screen
```
