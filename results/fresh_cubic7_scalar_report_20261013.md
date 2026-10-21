# Equal-block and scalar completion in the cubic-seven host

At the unchanged height bound 1000, allowing one equal-coordinate dihedral block together with the scalar coordinate found 25 new mirrors compatible with the original 25-wall cone. Seven cut the escape shared by the previous 158- and 122-wall constructions. The resulting cones have 182 and 146 walls, but still have exact spacelike escapes. No finite-volume example is claimed.

Write τ=2cos(2π/7), β=2−τ and D=2−τ². The slice has an integral saturated basis with coordinates

```
(x, a/D, a/D, z),
q = 2a²+z²−2τx².
```

D is a unit. The complete primitive norm pool is {1,2}. The ambient primitive classes β and 2β cannot occur here: reflection integrality at β forces x,a,z all divisible by β, contradicting norm valuation one at β. The previously proved ambient exclusion of primitive norm4 applies as well.

All three real embeddings bound x, a and z. For height H=1000 the chosen-embedding bounds are x²≤Hq/(2τ), a²≤q(H+1)/2 and z²≤q(H+1). At each other embedding they are x'²≤-q'/(2τ'), a'²≤q'/2 and z'²≤q'. Exact energy joins enumerate the norm equation. All four block positions, including a=0 and z=0, are covered; roots are deduplicated by their exact rays.

The original-25-compatible dictionary has four norm-one and 21 norm-two roots. Every saved root passes integrality, primitive reflection, old-wall and exact Coxeter-angle gates. One simple norm-one family has

```
x=2+3τ+τ², a=1+2τ+τ², z=1+4τ+2τ²,
height=20+58τ+26τ².
```

All four block embeddings of this family are accepted. Twenty-four of the 25 new mirrors are compatible with each previous completed cone. Adding them gives 182 walls from the 158-wall base and 146 from the 122-wall base.

The excluded twenty-fifth mirror is a useful pure-scalar norm-two cap:

```
x=1+2τ+τ², a=0, z=4τ+2τ²,
height=5+14τ+6τ².
```

It conflicts with the previous high Pell cap. Rebuilding from the original 25 walls using the complete two-block dictionary plus all scalar-slice roots selects this earlier cap and rejects the old Pell cap by the exact cross-angle gate. This alternative also has 146 walls.

All three final cones have the same certified spacelike escape, of norm

```
(-2885-26821τ+26557τ²)/5041 > 0.
```

The vector is supported on one equal block, one asymmetric block and the scalar line. Thus the new scalar walls remove the earlier witness, but a different escape remains. These are finite-cone infinite-volume certificates, not full-host nonreflectivity proofs.

The [complete artifact](fresh_cubic7_scalar_h1000_20261013.json) stores the original-25 dictionary, compatibility with both prior bases, the union dictionary used for replacement, all three selected wall sets, every finite-angle count, strict interiors and exact escapes. The [producer](../scripts/fresh_cubic7_scalar_search_20261013.py) ran in 14.34 seconds. The coordinating agent maintains the independent cubic replay separately.

```sh
PYTHONPATH=src:. python3 scripts/fresh_cubic7_scalar_search_20261013.py --height 1000 --output results/fresh_cubic7_scalar_h1000_20261013.json
```

Allowing an asymmetric block simultaneously with the equal block and scalar is a larger slice. It is pursued separately at the same height bound, rather than enlarging this slice's cutoff.
