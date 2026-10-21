# Intrinsic cusp caps: endpoint review and global closure test

The [endpoint profile archive](fresh_h20_intrinsic_caps_20261014_endpoint.md) passes source/proof review and a [separate accepted-record arithmetic check](../scripts/fresh_h20_intrinsic_caps_20261014_review.py), with [saved certificate](fresh_h20_intrinsic_caps_20261014_review.json). The checker uses SymPy radical expressions and rational matrix arithmetic, independently verifies all 260 accepted profiles and 252 transports, and replays the intrinsic geometric identities below. It does not repeat the producer's rejected-profile enumeration. Generation and saved replay pass.

The literal cusp diagram is **affine A11 × affine A2 × affine E6**, not D6. After scaling the two unchanged factors to inradius1, the triangle has an independent positive inradius `s`. This supplies an intrinsic construction space without restricting added normals to rational vectors in the old ambient realization.

| E6 deletion | Finite edge diagram | Literal edges | Finite profiles | PSD-singular profiles |
|---|---|---:|---:|---:|
| Outer node | A11 × A2 × E6 | 108 | 19 | 9 |
| Middle node | A11 × A2 × A5 × A1 | 108 | 54 | 14 |
| Center node | A11 × A2 × A2³ | 36 | 71 | 93 |

The middle type additionally has the isolated `I2(m)` family for every integer `m>=7`, together with its parallel limit. A singular profile is only a necessary ideal-endpoint condition; it is not a bounded-cusp certificate.

The reduction to three edge types is geometric, not an assumption that the original arithmetic realization has every symmetry. Rotations of the affine A11 cycle and the affine A2 triangle, together with permutations of the three E6 arms, give explicit Gram- and mark-preserving permutations covering every omitted triple. They preserve each factor separately and hence preserve its positive null relation and its metric scale. They extend to isometries of the intrinsic Lorentz space: their action on `D-perp` preserves its degenerate form and the common future null ray, and the complementary null vector is determined by its required pairings and zero square. The archive checks all 252 transports explicitly.

## Completeness of the local contact gate

Normalize all seed roots and a proposed cap to square2. Let B be the positive-definite Gram of the nineteen walls on a cusp edge. Write the cap's nonpositive contact vector as `-c`, with nonzero entries `c_i=2*cos(pi/m_i)`. A finite or ideal endpoint requires

```text
cost = c^T B^(-1) c <= 2.
```

Every checked inverse has nonnegative entries and diagonal at least1/2. Every nonright finite contact has magnitude at least1, so there can be at most four nonzero contacts. Setting their magnitudes to1 gives a rigorous lower bound and safely rejects supports before the label search.

A label at least7 is impossible if that old wall has a neighbor or the cap has another nonright contact. In either case a three-root principal minor is at most `1-cos²(pi/m)-1/4<0`. The only exception is a sole contact to an isolated A1; it gives the explicitly retained arbitrary-label tail. All remaining labels belong to `{3,4,5,6}` and are exhausted by the producer. This is a geometric PSD bound, not a field-degree or arithmetic cutoff. The fixed radical sign enclosures are rigorous, and unresolved signs cause assertion failure rather than acceptance.

Label5 is available intrinsically. There are two H3 attachments in the outer type, two H3 attachments and two H4 attachments in the middle type, and eight H3 attachments in the center type. The middle type also has one isolated `I2(5)` profile. An H3 attachment has

```text
cost = (3+sqrt5)/3,
2-cost = (3-sqrt5)/3 > 0.
```

An H4 attachment joins the cap by label3 to an old A2 endpoint and by label5 to the isolated A1. Its complete finite vertex diagram is `A11 × A5 × H4`, and

```text
cost = (17+3sqrt5)/12,
2-cost = (7-3sqrt5)/12 > 0.
```

These are the middle-type records43 and51. There is no H4 attachment on the other two edge types.

## A complete one-parameter cap representation

Let D be the original future null cusp and V the other future null direction perpendicular to the nineteen edge roots, normalized by `D.V=-1`. For a saved contact profile put

```text
p = -sum_i (B^(-1)c)_i * root_i,
rho2 = 2 - p^2.
```

Every cap with those contacts, positive norm2 and negative cusp pairing is

```text
r(h) = p + h*V - rho2*D/(2h),     h>0.
```

Indeed `r(h)^2=2`, `r(h).D=-h`, and its nineteen contacts are fixed. When `rho2>0`, the opposite endpoint is

```text
X(h) = V + rho2*D/(2h^2),
X(h)^2 = -rho2/h^2 < 0.
```

It is a genuine finite vertex once the remaining three seed walls pass their pair gates. Those three walls are strict at X independently of h, by their positive affine relations. When `rho2=0`, the proposed endpoint is V itself and further incident facets may be needed to bound its link.

For a factor with positive relation `sum m_i*root_i=mu*D`, its omitted wall has cap pairing

```text
(sum_retained m_i*c_i - mu*h)/m_omitted.
```

This proves the archive's three omitted-wall scalar formulas. All sufficiently large h make all three relations strictly dotted. Thus local existence is automatic for every positive Schur profile; the hard problem is choosing caps that are mutually Coxeter-compatible and cover the global gaps.

In upper half-space, set `R=sqrt2/h`. The ball center is the old cusp-link vertex plus `p/h`. If a factor's inradius is `a`, its omitted-wall distance divided by R is

```text
(a * sum_marks / m_omitted)/R
    - sum_retained(m_i*cos(pi/m_i))/m_omitted.
```

It must be a finite Coxeter cosine, or at least1 for a parallel/dotted relation. This treats the entire parameter interval; it does not require a finite nonedge anchor. Two norm2 caps with centers z_i,z_j and radii R_i,R_j have pairing

```text
(R_i^2 + R_j^2 - |z_i-z_j|^2)/(R_i*R_j).
```

Thus acuteness requires `|z_i-z_j|^2 >= R_i^2+R_j^2`; disjoint balls give dotted pairs once the distance exceeds `R_i+R_j`. The intervening finite-angle equalities remain exact requirements.

## A finite global certificate, rather than a vertex-coverage claim

Let C be the compact 19-dimensional product of the three Euclidean alcoves. For a finite compatible cap selection define

```text
L_i(z) = 2*z_i.z + R_i^2 - |z_i|^2.
```

The chamber has finite volume exactly when the closed balls cover all of C, equivalently `max_i L_i(z) >= |z|^2` for every z in C. This assertion assumes a full-dimensional chamber and retains the bounded original cusp link.

Here is a direct Lorentz proof. In the affine chart `D.X=-1`, write `X=V0+z+tD`; then `X^2=|z|^2-2t`. Seed inequalities give precisely `z in C`, while cap inequalities give `2t>=max_i L_i(z)`. The chart's only feasible direction with `D.X=0` is D itself, by the full-support positive affine relations. Therefore the ball cover is equivalent to the entire feasible polyhedral cone being future causal. Its finitely many extreme rays are then finite or ideal hyperbolic vertices; triangulation gives finite volume. Conversely a failed cover gives a spacelike feasible point, and a small perturbation gives a strict escape.

This criterion has a finite exact certificate. Intersect C with every power cell where one `L_i` is maximal. On such a cell, `L_i(z)-|z|^2` is concave, so its minimum is attained at a cell vertex. Check every vertex of every nonempty cell, including lower-dimensional cells; C is compact, so there is no recession or unbounded-cell case. Nonnegative values certify finite volume. A negative value certifies failure. Equality locations are handled by the same polyhedral-cone argument, rather than being assumed to have bounded links from rank alone.

Checking only the 252 original vertices, or closing all 252 original edge germs, is necessary but insufficient: the power cells introduce additional vertices. The [golden construction](fresh_h20_intrinsic_caps_20261014_golden_seed.md) is an explicit example of that distinction.

```sh
python3 scripts/fresh_h20_intrinsic_caps_20261014_review.py --verify
```
