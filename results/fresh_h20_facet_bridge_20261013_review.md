# Independent source review of the H20 facet-bridge obstructions

**PASS.** I read the [exact checker](../scripts/fresh_h20_facet_bridge_20261013_exact.py), the relevant inherited pairing/point helpers and the [lane report](fresh_h20_facet_bridge_20261013.md), then reran:

```sh
python scripts/fresh_h20_facet_bridge_20261013_exact.py --verify
```

The replay reproduced [the frozen certificate](../data/enumerations/fresh_h20_facet_bridge_20261013_exact.json): all14 direction rejections and all175 reflected-normal inequalities passed. This is an independent mathematical/source review and a fresh replay using the existing arithmetic, not a second field implementation.

The perpendicular-direction proof covers its stated continuous family. The18 neighbor equations have rank18 in the21-dimensional ambient space; the independent facet175 equation leaves a positive-definite plane. Normalizing a candidate and old wall18 reduces their pairing to a unit-circle dot product. The exact interval bound on the old wall's normalized projection is strictly below `cos²(pi/9)` and below one, so neither a dotted/parallel pair nor a finite label of at least9 is possible. The projection bound is also above `cos²(pi/8)`, ensuring that each label2 through8 has exactly two real branches. The formula in the checker encloses both branches, including both right-angle directions. Each has a rigorously positive old pairing or a cosine strictly between consecutive allowed values. No center or actual-ridge assumption enters this necessary pair-compatibility rejection.

The reflection-bank argument has the correct signs. With `P²>0`, every predecessor normal satisfies `(P,r)<=0`, while the rechecked strict escape has `(r,X)<0` and `(P,X)<0`. Thus

```text
(s_P(r),X) = (r,X) - 2(P,r)(P,X)/P² < 0.
```

The same future spacelike point survives every reflected inequality and every current inequality. Therefore any full-dimensional hyperbolic chamber defined by a selection from that union has infinite volume, including selections deleting current walls. Pair compatibility of the reflected bank with the current walls is unnecessary for this negative conclusion.

These are two specific construction obstructions. They do not rule out nonorthogonal attachments, changes to neighbor relations, other reflections, or general finite-volume completion in H20.
