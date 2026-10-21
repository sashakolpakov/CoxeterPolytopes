# Complete finite-endpoint profiles for the mixed D4 chamber

The [exact census](fresh_h9_mixed_d4_endpoint_profiles_20261014.json) gives **1,825 finite profiles and 72 unbounded dihedral tails** on all 52 open rank-eight edges of the frozen [twelve-wall chamber](fresh_h9_mixed_d4_coupling_20261014.md). These are complete necessary local possibilities for a new finite endpoint, including H3, H4 and arbitrary I2(m) components. They are not globally accepted caps or a finite-volume construction.

The [producer](../scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014.py) uses the frozen independent quartic arithmetic to recheck the physical old-wall labels, and exact arithmetic in `Q(sqrt5)` for normalized support inverses and Schur complements. It binds the [complete finite-corner frontier](fresh_h9_mixed_d4_coupling_20261014_frontier.md). No old ray census or numerical cap search is repeated.

The eight priority releases at the new polar vertices have the following local possibilities. The paired rows differ by replacing source wall 5 with wall 6; this statement concerns the displayed support diagrams, not a claimed global chamber symmetry.

| Edge IDs | Old rank-eight components | Finite profiles | Tails |
|---|---|---:|---:|
| 18, 45 | A3 + A3 + 2A1 | 58 each | 2 each |
| 19, 46 | A3 + A2 + 3A1 | 108 each | 3 each |
| 20, 47 | A3 + I2(5) + 3A1 | 62 each | 3 each |
| 21, 48 | H3 + A3 + 2A1 | 32 each | 2 each |

Every profile lists its nonright contacts, exact projection square and positive residual norm, and all spherical component types. The explicit finite list uses labels 3, 4 and 5; each tail consists of a single contact with **every integer label `m >= 6`** to an isolated old A1, all other support contacts being right. Thus labels 3–5 are not counted twice.

Concrete mixed contacts are available. On edge 18, contacts to walls `(1,5,11)` all at π/3 give endpoint type `D6 + A3` and residual norm `1/4`. This includes a nonright contact to the new polar wall 11. On edge 19, contacts `(0,2,9)` all at π/3 give `E7 + 2A1`, residual `1/12`. These are local recipes only; the remaining four old walls still impose angle and geometric inequalities.

## Why the label and support census is complete

Normalize the eight old support roots and the proposed root to square 2. Nonzero Coxeter Gram entries have magnitude at least 1. The old support Gram `G` is positive definite.

1. A cycle has nonpositive quadratic form on its all-ones vector. Therefore the new node can meet at most one vertex in each old connected component.
2. Four contacts contain a four-leaf star whose quadratic form is nonpositive on the vector with center coefficient 2 and leaf coefficients 1. Thus the new node has at most three contacts.
3. With three contacts, a label at least 4 gives squared edge magnitudes at least `2,1,1`; the star Schur complement is at most `2−(2+1+1)/2=0`. All three labels must therefore be 3.
4. A label at least 6 together with any second adjacent edge gives a three-node path with squared edge magnitudes at least `3,1`, again with Schur complement at most zero. Hence labels at least 6 occur only when the new node has one contact, to an isolated old A1. This is precisely the separate I2(m) tail.

These arguments use principal quadratic forms, rather than a bounded label search or an unverified catalog cutoff. Possible additional edges only decrease the displayed nonnegative-vector quadratic forms.

For all remaining contacts, write `c_i=2cos(π/m_i)`. Since the contacts lie in different old connected components, all their inverse cross terms vanish. The new nine-root Gram is positive definite exactly when

`rho = 2 − Σ c_i² (G⁻¹)_ii > 0`.

The exact squared magnitudes are `1`, `2`, and `(3+sqrt5)/2`. The checker exhausts 10,753 finite trials, accepts the 1,825 positive Schur complements, and separately records all 72 tails. Each tail has `rho=2sin(π/m)^2>0` for every finite integer `m >= 6`. Affine endpoints, where the Schur complement is zero, are not part of this census.

## Exact oriented one-parameter recipes

Each edge also contains a physical, unnormalized Lorentz frame. Let `X` be its finite start, `q` the saved spacelike other cone ray, `A=−X²`, and

`Y = A q + (X,q)X`.

The checker verifies `A>0`, `Y²>0`, `(X,Y)=0`, and perpendicularity of both vectors to all eight support roots. It also verifies that `Y` pairs strictly negatively with the released old wall, fixing the outgoing direction. Set `U=X/sqrt(A)` and `V=Y/sqrt(Y²)`.

For any accepted profile, form its support projection

`p = −Σ (G⁻¹c)_i s_i`.

Every root of square 2 with those support pairings has residual form

`r = p + sqrt(rho)(sinh(t)U ± cosh(t)V)`.

Retaining the finite starting vertex requires `t>0`. To meet the outgoing edge, the sign must be positive; the negative sign meets the backward half-line and violates the released old wall. The corresponding future finite endpoint is

`v(t) = cosh(t)U + sinh(t)V`.

With `T=exp(t)>1`, acuteness to any remaining old root `w` becomes the quadratic inequality

`(sqrt(rho)/2)(U+V,w)T² + (p,w)T + (sqrt(rho)/2)(V−U,w) <= 0`.

This supplies an exact continuous feasibility test before imposing the remaining discrete angle conditions. No such global parameter feasibility or first-intersection assertion is made by this artifact. The normalized profile scalars are in `Q(sqrt5)`; the physical frame vectors use the inherited quartic field specified separately in the JSON. Further square roots in the normalized recipe are kept explicit.

Generation and saved replay pass:

```sh
python scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014.py --verify
```

The coverage is exactly the 52 saved open directions incident to existing finite vertices. It does not classify vertex-free edges, arbitrary new facets elsewhere, affine endpoints, or compatible simultaneous completions. All prior artifacts remain unchanged.
