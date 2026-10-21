# A complementary135th H20 wall cuts E0

The frozen134-wall chamber admits an exact further Coxeter wall cutting
the common escape E0 from the consecutive-arc obstruction. It is dotted
to all26 rotational-chain walls and to the preceding arc cap. The new
chamber has135 actual facets, retains all known finite vertices and the
three tracked proper cusps, and still has infinite volume: a different
explicit strict spacelike escape survives.

The [checker](../scripts/fresh_h20_complementary_cap_20261016.py) and
[certificate](fresh_h20_complementary_cap_20261016.json) pass generation
and fresh saved replay. Separate independent review is pending. No frozen
predecessor is changed.

## The recipe and its provenance

This is proposal44 of the old
[bounded nonright-S61 bank](fresh_h20_nonright_s61_20261014.md).
Its exact quadratic recipe was already reconstructed for that bank's
common-escape proof. That earlier proof did not certify every remaining
Coxeter contact or compatibility with the later134th wall. Those are new
obligations, checked here. This is reuse of an existing proposal, not a
claim of discovering an additional numerical family.

Use the frozen unit Lorentz frame x,e,f perpendicular to S61. In the
notation of the [arc-cap construction](fresh_h20_dotted_arc_cap_20261016.md),
set

```
r = y-omega5,                    y=t*x+z*f,
omega5²=30/11,                   y²=-8/11,
r²=2.
```

The old wall95 is √2 e, so r is exactly perpendicular to it. Write
p=x.n0 and q=f.n0 for the exact normalized old-wall coefficients. Since
omega5.n0=−6/11, the desired angleπ/6 with wall0 is

```
p*t+q*z = h,       h=-6/11-√3.
```

Put alpha=h/p and beta=−q/p. Then t=alpha+beta*z, and z is the positive
root of

```
(1-beta²)*z²-2*alpha*beta*z+(8/11-alpha²)=0.
```

The certificate proves the leading coefficient and discriminant positive,
isolates this branch in

```
0.9491<t<0.9493,       0.4167<z<0.4169,
```

and excludes the other branch, whose z is negative. The norm and the
prescribed contacts follow algebraically from these defining equations.
Outward448-bit rational intervals verify all remaining strict signs;
interval overlap alone is not used to assert the defining equalities.

The complete134 contact table is:

| Relation | Count |
|---|---:|
| right |18|
| π/3, with support wall5 |1|
| π/6, with old wall0 |1|
| dotted |114|

The right walls are the other17 support roots and wall95. Every other
retained wall is dotted, including all26 chain walls. The weakest dotted
old contact is wall65, with pairing approximately−2.028739677.
The mutual product with cap133 is approximately−4.689590041, strictly
dotted. Thus the cutter lies outside the three-nondotted-chain family
excluded by the new field argument.

## Coherence, facets and finite geometry

The same saved strict future timelike center works for all135 walls;
its new pairing is approximately−0.949193826. For each norm2 wall n_i,

```
H_i=H-(H.n_i)*n_i/2
```

is future timelike, lies on n_i, and is strictly negative on every other
wall by acuteness. This supplies135 actual facet-interior witnesses.
The old finite vertices survive by the acute-extension preservation
lemma; the three tracked proper cusps have strictly negative new
pairings, approximately−15.63175,−24.67293 and−17.39273, so their old
local links are unchanged.

There are **no new finite vertices** on cap134, and this statement is
complete for the new facet. Its20 finite old neighbors are S61 together
with0 and95. A finite20-root corner would use19 of these neighbors.
The old pair0/87 is dotted, so it must omit0 or87.

- Omitting0 retains all18 support roots and r. Their Schur value is
  `2−30/11=−8/11`, so the Gram cannot be positive definite.
- Omitting87 retains the triple(r,0,1), whose Gram is
  `[[2,−√3,0],[−√3,2,−1],[0,−1,2]]` and has positive null vector
  `(√3,2,1)`. This affine subdiagram also prevents positive definiteness.

Cap133 is dotted to r and cannot be an incident finite neighbor. This
exhausts the possible new finite corners. The chamber has at least170
finite vertices inherited from the previous package; no complete ambient
vertex or exterior-ray census is claimed. Absence of a new finite vertex
does not affect validity of the new facet or its escape-cutting role.

## The old gap is cut; another exact gap remains

For

```
E0=x-e/100-omega5,
```

the new pairing is approximately+1.778078910, rigorously positive. The
previous strict epsilon perturbation of E0 is cut too. The preceding
cap133 cuts the earlier847-bank escape, while this cap cuts the escape
that survived cap133. The two cuts therefore address different regions.

The following simple vector survives all135 walls:

```
E1=x-e/100-(49/50)*f-omega17
   -(1/100000)*sum_{j!=17}(omega_j),
E1²=1842000049/4000000000 >0.
```

Every pairing is certified strictly negative. Its pairings with caps133
and134 are approximately−0.453488869 and−1.357477175, and it pairs
negatively with the common future center. The ray from the center toward
E1 crosses the future null cone within strict wall inequalities, leaving
an open ideal-boundary patch. Hence the135-wall chamber has infinite
hyperbolic volume.

The result is an exact coherent complementary-wall extension, not a
finite-volume example or a claim that the remaining cap families have
been exhausted. Further walls must address the new escape while meeting
both existing caps with allowed Coxeter relations.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_complementary_cap_20261016.py --verify
```
