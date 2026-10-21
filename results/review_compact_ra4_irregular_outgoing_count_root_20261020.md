# Root review: irregular outgoing edges and ridge excess

20 October 2026. Complete written review of
`compact_ra4_irregular_outgoing_facet_excess_count_20261020.md`, SHA256
`f57760ea1bf0e023857ee8798bbdfba68ae3c31a2c3dcfccd13ee6589f297795`.
PASS. No scientific execution or unrestricted volume minimum.

For a pentagonal perpendicular ridge, the two actual outgoing endpoint
walls are adjacent, giving exactly coth(s_v)*coth(s_w)=phi. The defect
therefore changes sign, and a pentagonal ridge cannot separate a zero
defect from a nonzero one. No one-sided bound on outgoing lengths is
assumed.

The proper marked-set cut has at least three edges: its corresponding
dual icosahedral subgraph is nonempty and has even degrees, hence a cycle
of length at least three. In the all-marked case each intrinsic pentagon
needs a missing pentagonal-ridge edge; twelve pentagons and two incidences
per edge give at least six. A nonpentagonal ridge cannot intersect one
regular D in two edges.

I checked the inherited three-consecutive-L-edge closure argument in
the primary compact RA4 source. It concerns the actual closing walls of
the right-angled polygon. Thus a nonpentagonal k-ridge has at most
floor(2k/3) dodecahedral edge carriers. This is at most4(k-5), including
the extremal k=6 case. Counting the actual perpendicular incidences
proves3*N_irregular<=4W.

The primary source's small-degree lemma proves each positive facet has
excess at least two. If its excess were two, its two weighted graph
edges would each have weight one, and the reviewed degree-two theorem
would force W>=20. Therefore W<20 gives excess at least three for every
positive facet. The exact sum of all facet excesses is2W, yielding the
stated bounds N_uniform>=f-7W/3 and N_uniform>=f-2W respectively.
The universal capacity f>=45 makes the latter at least seven for W<=19.

The final center-containment consequence correctly invokes the separately
reviewed uniform-edge theorem, with no claim that the centers differ or
that the matched whole120-cells are contained. Low volume alone does not
bound W by19, so this does not establish a universal minimum or eliminate
arbitrary outgoing lengths in every possible smaller-volume example.
