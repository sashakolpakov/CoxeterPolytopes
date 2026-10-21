No cap in the **308 single-weight/finite-vertex-ray families** can cut the saved global escape while remaining Coxeter-compatible with the retained 87 walls. The corresponding **14 single-weight/ideal-vertex-ray families** also supply no new compatible cutter. This exact all-parameter result closes these particular transverse families; mixed E7 contacts and residual directions along H13 edges or interiors remain outside it.

The base is the [87-wall E7 cusp completion](fresh_h20_e7_intrinsic_caps_20261014_cusps.json), whose completed H13 face has 44 finite vertices and two ideal vertices in the [85-wall parent](fresh_h20_e7_intrinsic_caps_20261014.json). Its E7 walls, in current numbering, are `[12,13,14,15,18,31,65]`. Their fundamental-weight squares are `[7/2,12,6,2,15/2,3/2,4]`. In particular, the previously named old wall 67 is current wall 65.

The [exact checker](../scripts/fresh_h20_e7_transverse_20261014.py) and [certificate](fresh_h20_e7_transverse_20261014.json) bind the frozen inputs. The checker recomputes all seven weights and all 44 vertex norms, incident sets, positive-definite rank-20 frames and signs against all 87 walls. It then certifies:

| Exact finite-ray exclusion | Families |
| --- | ---: |
| Escape bound at the minimum possible Coxeter coefficient | 174 |
| One retained-wall acuteness bound already prevents cutting | 131 |
| Only the forbidden coefficient interval between 1 and √2 could cut | 3 |
| Total | 308 |

For a finite face vertex V and an E7 fundamental weight ω, put `q=−V²>0` and `d=ω²`. The norm-two family is

`r(c)=−cω+sqrt((dc²−2)/q)V`.

It is orthogonal to the 19 retained walls obtained by releasing the relevant E7 wall from that vertex. Its pairing with the released norm-two wall is `−c`. Acuteness forces `c>0`, and a Coxeter relation therefore requires `c=1` or `c≥√2`; the latter includes every higher finite label, parallel contact and dotted separation. No label cutoff is used.

This also exhausts the possible sheets for a valid addition with that profile. Since every coefficient of an E7 fundamental weight in the simple roots is positive, every interior point x of the retained chamber satisfies `ω·x<0`. The vertex V is future timelike. The negative square-root sheet therefore has positive pairing with every such x, and the zero-root boundary does too. Only the displayed positive sheet with `dc²>2` can bound a chamber interior.

Write `r(c)/c=−ω+uV`, where `u=sqrt((d−2/c²)/q)` increases with c. For the saved escape E, define `A=−ω·E` and `β=−V·E`; both are strictly positive in all 308 cases. A cutter requires `u<A/β`. For an old normal N, write `a=−ω·N`, `b=−V·N`, with `b≥0`. Acuteness requires `a−ub≤0`. Thus a positive a with zero b immediately excludes the family; otherwise `aβ−Ab≥0` implies `u≥a/b≥A/β`, preventing cutting. These are the 131 one-wall certificates. In the other 174 immediate cases, the escape expression `A−βu` is already nonpositive at `c=1` and decreases thereafter.

Only vertices 40, 41 and 42 at current wall 65 require the last step. Each has a strictly negative escape pairing at `c=√2`, excluding every larger c. At `c=1`, a saved retained wall has strictly positive pairing, excluding that value. The exact checker proves both inequalities by field arithmetic and squaring; rounded numerical windows are not used in the proof.

For a null face ray ν, a norm-two cap `r=−cω+tν` must have `c²=2/d`. Ten of the 14 profiles give `0<c<1`, and two give `1<c<√2`, so their E7 contact is forbidden. The remaining two have `d=2,c=1`. Their pairing with the corresponding retained filler `−ω+2ν` is identically `+2`, independently of t. They therefore cannot be new acute-compatible facets. At `t=2` the normal merely repeats the existing facet, which is already strictly negative at E.

The same escape has square `1336028658157091827/250000000000000000>0` and is strict against all original 87 walls. Every admissible new cap in the specified families has pairing at most zero with it. For any finite coherent addition with a common strict timelike point H, the point `E+εH` is strictly feasible and remains spacelike for sufficiently small positive ε. Hence no finite completion using only these families has finite volume. This argument does not assert uniform strictness for an infinite collection and does not apply after removing retained walls.

The preceding [numerical scout](../scripts/fresh_h20_e7_transverse_20261014_scout.py) and [output](fresh_h20_e7_transverse_20261014_scout.json) tested the five coefficients `1,√2,φ,√3,2` and selected contact equations. It found no cutter. Its three continuous windows were approximately `(1.11354,1.18829)`, `(1.01722,1.05671)` and `(1.00310,1.02751)`. Those exploratory values led to the exact coefficient-gap proof above; the certificate does not import or rely on the scout.

Generation and saved-result replay:

```sh
python scripts/fresh_h20_e7_transverse_20261014.py
python scripts/fresh_h20_e7_transverse_20261014.py --verify
```

Both commands pass. No previous source, certificate or report was modified, and no new finite-volume polytope is claimed.
