# Independent review of the coordinated residual caps

**PASS.** The [reviewer](../scripts/fresh_h20_dominant_axis_cluster_20261016_review.py)
and [certificate](fresh_h20_dominant_axis_cluster_20261016_review.json)
reconstruct the [156-wall cluster](fresh_h20_dominant_axis_cluster_20261016.md).
Generation and a fresh saved replay pass. No producer arithmetic code is
imported; numerical scores are not proof inputs.

Starting with the independently reconstructed 154-wall chamber, the reviewer
derives s_star=1/(2*s153) and c_star=sqrt(1+s_star^2). The orthogonal residual
frame proves exact norm2, the right contacts and the two pi/3 contacts.
All309 new pairings are replayed:38 right,2 pi/3 and269 strictly dotted.
The mutual product is exactly -2*s_star^2<-2, so the new walls are compatible.

Four separate literal four-root Gram inverses construct the new finite
vertices, rather than the producer's two-root projection formula. The
sixteen-root support is projected independently. All156 pairings at each
point are checked, including the other new wall. Every other possible
nineteen-neighbor support contains the dotted pair142/153. No finite vertex
can contain both new walls because they are dotted. This proves exactly
four new finite vertices, giving at least180 in the ambient chamber.

The strict common future interior and facet projections certify156 actual
facets. All ten new-wall/proper-cusp pairings are negative. The three
vertices of the compact (2,4,6) triangle remain strictly inside both new
halfspaces, preserving the whole triangle by convexity.

The exact residual cuts are z141,z143 <= s_star/c_star < 1. The same strict
spacelike E6 survives all156 inequalities, proving infinite volume. The
independently checked geometry does not establish completeness of the
ambient finite-vertex census, a quantitative volume reduction, or a new
finite-volume construction dimension.

```sh
python3 scripts/fresh_h20_dominant_axis_cluster_20261016_review.py --verify
```
