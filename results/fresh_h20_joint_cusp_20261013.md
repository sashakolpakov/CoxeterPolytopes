# Joint angle repairs with an entire Coxeter cusp

This round gives two exact results. **Every interior section in one full affine-A2 normal span is non-Coxeter, for all real coefficients.** A separate construction produces **12 finite-volume H20 sections with irrational cusp directions and complete Coxeter cusp links**; eight also realize a prescribed actual global right angle. All twelve have another certified forbidden ridge. No Coxeter polytope was found.

The two proofs have independent standard-library/Fraction replays. The local direction scouts are exploratory and have narrower claims, described below.

## A two-parameter construction, and its persistent obstruction

The original Borcherds chamber has a recorded cusp whose diagram is affine `A11 × D7 × A2`. The triangular A2 factor has walls 60, 99 and 125, with Gram

```
[[ 4, -2, -2],
 [-2,  4, -2],
 [-2, -2,  4]].
```

Their sum is twice the archived future null cusp vector `u`. Thus an arbitrary interior line section of this triangular horospherical factor leaves a segment. The other two factors remain unchanged, giving a complete rank-19 Coxeter cusp `A11 × D7 × A1`.

In particular, both coefficients in

```
e = r60 + s r99 − t r125,       s,t > 0
```

can vary independently. The omitted wall is globally redundant in the section, since

```
r125 = (r60+s r99)/t    modulo e.
```

The remaining two projected walls have a positive affine-kernel relation

```
(1+t) Pi(r60)+(s+t) Pi(r99) = 2t u,
```

where `Pi` is orthogonal projection to `e⊥`. Together with the unchanged affine factors, this certifies an actual bounded horospherical link, not merely formal local angle values.

However, ambient walls **1 and 22 define an actual section ridge for every positive s,t**. This is certified by two strict ambient face points: `z+` vanishes exactly on `{1,22,60,99}`, and `z−` exactly on `{1,22,125}`. Their squares are −80 and −91. Moreover,

```
(e,z+) = 26t > 0,
(e,z−) = −23(1+s) < 0.
```

Consequently

```
p = 26t z− + 23(1+s) z+
```

is future timelike, belongs to `e⊥`, and vanishes on exactly walls 1 and 22. The projected normals are independent and lie in positive-definite `p⊥`, so this is an actual dihedral ridge.

Its squared cosine is

```
A/(2B),
A = s²+st−s+t²+t+1,
B = 2s²+2st−2s+t²+2t+2.
```

The exact identities

```
2A−B = t² > 0,
B−A = (s−1/2)²+3/4+st+t > 0
```

place this value strictly between `1/4` and `1/2`. Hence the angle lies strictly between `π/4` and `π/3`, and cannot be a Coxeter angle.

The same formulas and separate actual-face certificates cover both other choices of negative coefficient:

| Negative wall | Positive walls | Persistent forbidden ridge |
|---|---|---|
| 60 | 99, 125 | 3 / 27 |
| 99 | 60, 125 | 9 / 8 |
| 125 | 60, 99 | 1 / 22 |

Every normal in this three-dimensional real wall span that cuts the strict interior has mixed coefficients. After an overall sign and scaling, nonzero coefficients have one of these three forms. Zero coefficients reduce to the [previous universal two-wall result](fresh_h20_global_20261013.md). This therefore excludes the **entire displayed affine-A2 normal span**, including noncentral and algebraic normals. It does not exclude other affine spans or general H20 sections.

Strict section interiors are also constructed by interpolating faces `{i,j}` and `{k}` whose zero sets are disjoint. Finite volume follows from the finite-vertex Klein section argument used in the previous section reports. The positive cusp identities prove noncompactness.

## Irrational quotient directions give twelve new exact sections

The direction scout also found two cuts of the rank-three affine factor at recorded cusp 4, using walls `{34,46,204,210}`. Their quotient directions involve `sqrt(2)`. After dropping a certified redundant facet, the factor becomes an affine G2 triangle, with angles `π/2, π/3, π/6`; the other cusp factors have ranks 10 and 7.

Translate either normal by `tau u`, where `u` is this cusp vector. All local projected cusp angles stay constant. For a global target pair with exactly one wall incident to the cusp, its right-angle equation is linear in tau over `Q(sqrt(2))`. Six targets were selected for each direction and solved exactly. This gives twelve literal sections:

| Exact result | Count |
|---|---:|
| Irrational quotient directions | 2 |
| Finite-volume noncompact sections | 12 |
| Complete rank-19 Coxeter cusp links | 12 |
| Prescribed right angle certified at an actual ridge | 8 |
| Rejections by a different actual forbidden ridge | 12 |

The other four prescribed pairs satisfy the algebraic orthogonality equation, but their incidence was not certified. No actual target-angle claim is made for those four. Their separate actual forbidden-ridge witnesses still reject their sections.

A concise exact example is

```
E = 2sqrt(2) r34 + (3−sqrt(2)) r46
    − sqrt(2) r204 + 2sqrt(2) r210,
E² = 54.
```

Wall 204 is redundant in `E⊥` by a strictly positive combination of the other three displayed inequalities. Actual ridge 2/46 is right-angled; the complete cusp link is Coxeter. But actual ridge 1/169 has squared cosine **27/52**, strictly between the values for labels 4 and 5. The archive includes exact timelike ridge points and every wall inequality.

For every section, the replay checks the redundancy identity, strict interior, every local angle, connected positive affine kernels, their positive multiples of `u`, and the rejecting actual ridge. A connected Gram with nonpositive off-diagonal entries and a positive kernel vector is congruent to a weighted graph Laplacian; this certifies its rank and positive semidefiniteness. The cusp proof therefore includes the full geometry of the omitted facet and retained factors.

All twelve normals require **exactly four original walls**: their displayed expressions use four, and an exact quotient-span test excludes every span of at most three original walls. The argument uses the irrationality of `sqrt(2)`: any rational wall span containing the normal must contain both its rational coefficient vectors.

This batch covers only the twelve declared linear target equations. It does not exclude other translations, other global angle targets, or other irrational local directions.

## How the local directions were selected

Two scouts supplied proposals:

- Mixing two two-node directions from distinct affine factors: 44,032 oriented choices and 7,216 positive rational mixing-ratio candidates; no retained-facet pattern passed the local affine Coxeter gate. Exactly two prescribed endpoint facets were omitted. Exceptional cases without an endpoint-to-zero edge were restricted to rational roots. This is not an exclusion of all mixed-factor sections or all facet subsets.
- Mixing two directions on a length-two path inside one factor: 1,049 local facet patterns. This found the continuously variable A2 family and the two usable irrational directions above. Other local algebraic outputs can have redundant or inactive facets that were not geometrically certified.

These are proposal-generating searches. The universal A2 theorem and twelve irrational geometric certificates are independently verified; the report does not promote the entire exploratory local census to a universal nonexistence theorem.

Revalidation also corrected a diagram name in the [preceding cusp report](fresh_h20_new_route_20261013.md): its rank-seven affine factor is **D7**, with two trivalent nodes. Its earlier exact Gram and kernel certificates were already correct.

## Artifacts and replay

- [A2 producer](../scripts/fresh_h20_joint_cusp_a2_20261013.py), [exact certificates](../data/enumerations/fresh_h20_joint_cusp_a2_20261013.json), [independent verifier](../scripts/fresh_h20_joint_cusp_a2_verify_20261013.py), [PASS record](../data/enumerations/fresh_h20_joint_cusp_a2_verify_20261013.json).
- [Irrational-cut producer](../scripts/fresh_h20_joint_cusp_irrational_20261013.py), [exact certificates](../data/enumerations/fresh_h20_joint_cusp_irrational_20261013.json), [independent verifier](../scripts/fresh_h20_joint_cusp_irrational_verify_20261013.py), [PASS record](../data/enumerations/fresh_h20_joint_cusp_irrational_verify_20261013.json).
- Local scouts: [cross-factor script](../scripts/fresh_h20_joint_cusp_20261013.py) and [data](../data/enumerations/fresh_h20_joint_cusp_20261013.json); [two-direction script](../scripts/fresh_h20_joint_cusp_directions_20261013.py) and [data](../data/enumerations/fresh_h20_joint_cusp_directions_20261013.json).

```sh
python scripts/fresh_h20_joint_cusp_a2_20261013.py
python scripts/fresh_h20_joint_cusp_a2_verify_20261013.py
python scripts/fresh_h20_joint_cusp_irrational_20261013.py
python scripts/fresh_h20_joint_cusp_irrational_verify_20261013.py
```

Both independent replays use only the standard library and earlier standalone verification helpers, importing no producer, SymPy, NumPy, or numerical solver. A forged ambient-face point and a forged redundancy identity are rejected. A second agent reviewed the universal A2 actual-ridge and cusp arguments at source level and found no issue; this was not a second producer run. All six new scripts compile.
