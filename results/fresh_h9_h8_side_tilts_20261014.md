# The fixed H8 source has no nonright side lifts

Every Coxeter sidewall inducing a facet of the literal compact 11-facet Bugaenko H8 source must meet its containing hyperplane **at a right angle**. Thus tilting those eleven sidewalls does not provide a new attachment mechanism. Additional ambient caps disjoint from the H8 facet remain possible and can change the inward-end coverage.

The [checker](../scripts/fresh_h9_h8_side_tilts_20261014.py) and [certificate](fresh_h9_h8_side_tilts_20261014.json) reconstruct the source Gram from its physical normals, check the literal graph and all exceptional angle gaps, and pass generation and saved replay. The [independent Fraction review](fresh_h9_h8_side_tilts_20261014_review.md) also passes, including all 121 Gram entries, five exact gaps and all 2,048 right/nonright masks. This package has no dependency on the previous attachment optimizers or cap-direction scouts.

Let r be the containing facet normal and u_i a projected source normal, both of square two. At a generic finite point of a source facet, the ambient local two-dimensional Coxeter wedge supplies its unique adjacent sidewall. Its contact angle with r is some finite `theta_i=pi/m_i`, and its positively normalized normal is

`s_i=sin(theta_i)u_i-cos(theta_i)r`.

Consequently

`(s_i,s_j)=sin(theta_i)sin(theta_j)(u_i,u_j)+2cos(theta_i)cos(theta_j)`.

On a source pi/3 edge, if exactly one side lift is nonright, the ambient dihedral cosine is `sin(theta)/2`, strictly between zero and one half. No Coxeter angle has that cosine. Right versus nonright status must therefore be constant on every component of the source's pi/3 graph.

That literal graph has components `{0}`, `{1,2,3,4,5,6,7,9,10}`, `{8}`. The large component contains source-right pair 1,3. If the component were nonright, this pair's ambient scalar product would be positive, contradicting acuteness. Every node of that component is therefore right.

The remaining nodes 0 and 8 have source pi/5 edges to nodes 1 and 7 respectively. A nonright lift at either node would give ambient cosine `cos(pi/5)sin(pi/m)`. For m=3 or 4 its square lies strictly between 1/4 and 1/2, hence its cosine lies between the consecutive allowed values for pi/3 and pi/4. For m>=5 its cosine is positive and strictly less than one half. All finite nonright labels are excluded.

The theorem concerns the fixed source facet itself. A wall disjoint from that facet need not be perpendicular and may create new truncating balls. Such walls, changed old walls and different compact H8 source geometry are not excluded. This is not a nonexistence theorem for compact H9 polytopes.
