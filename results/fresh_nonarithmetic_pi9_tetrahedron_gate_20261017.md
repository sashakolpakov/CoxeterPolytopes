# The new π9 face cannot close with one tetrahedron facet

The frozen60-wall chamber has an actual H3 face supported by `[0,1,2,3,58,59]`, of type `H4+I2(9)`. Old5,6,55 are genuine sides meeting at a finite corner; their diagram is `A2+A1`. No single added coherent Coxeter wall can close this face as a compact tetrahedron.

The [checker](../scripts/fresh_nonarithmetic_pi9_tetrahedron_gate_20261017.py) and [certificate](fresh_nonarithmetic_pi9_tetrahedron_gate_20261017.json) reconstruct that finite corner from the [60-wall geometry](fresh_nonarithmetic_pi8_quad_completion_20261017.md). The support is positive definite. Every additional actual facet of this face must be right to the whole support: a finite spherical extension cannot attach nontrivially to H4, and the inverse diagonal of I2(9) exceeds2. Therefore all its intrinsic side angles are native Coxeter angles.

For a compact tetrahedron all proper principal Gram blocks are positive definite and the full Gram is Lorentzian. A connected minimal indefinite four-node diagram cannot have a finite label at least6: any nonright edge adjoining that edge would already make a three-node principal block nonpositive. Thus the64 contact triples with labels2..5 exhaust all possibilities for the fourth side against the retained three sides.

Exactly two triples survive, in side order `[5,6,55]`: `[5,2,3]` and `[2,5,3]`. Both are marked versions of the compact `[3,5,3]` chain. Their normals are uniquely determined on the allowable orientation of the [complete π9 endpoint conic](fresh_nonarithmetic_pi9_endpoint_domain_20261017.md). In that parameterization, the new55 contact is

`−s(t+1/t)=−1`, with `s=sqrt((2−p²)/2)` and `t>1`.

The surviving root is `t=(1+sqrt(1−4s²))/(2s)`, approximately2.3828518817. The other root is below1 and fails the old55 endpoint inequality. Both candidate normals satisfy their prescribed contacts and have all other old contacts strictly dotted except49. Their49 product is approximately−0.9686402337, strictly between−1 and0. This excludes every finite Coxeter label, without a field or angle cutoff.

Generation and a fresh saved replay pass, together with syntax, source-hash and report-link checks. This is a primary exact certificate sharing the physical interval engine of its frozen inputs; no separate arithmetic review is claimed. The conclusion only excludes a one-facet compact tetrahedron closure. A triangular prism or another compact polyhedron using more walls remains open, as does ambient finite-volume completion.
