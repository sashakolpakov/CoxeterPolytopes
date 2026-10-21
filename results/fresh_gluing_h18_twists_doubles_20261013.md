# Fresh H18 gluing and parent-double experiment

No new dimension or new nonarithmetic lower-dimensional example was found.
This experiment changes the construction, rather than increasing an earlier
Vinberg cutoff.

The common 19-wall H17 gluing facet has weighted Gram automorphism group
of order two. The two parent chamber groups have orders 24 and 6. Their
gluing-facet stabilizers both have order two and restrict onto the entire
facet group. Therefore the relative gluing double-coset space has exactly
one element: every direct twist is equivalent to the already reconstructed
Vinberg mixture.

There are 14 even left-parent walls away from the gluing facet, and two
even right-parent walls. Modulo the original seam stabilizers, these give
10 left-parent and one right-parent double constructions. Each was formed
exactly by reflecting the relevant parent, retaining the other parent, and
adding the reflected image of the original gluing wall as an exterior wall.
The resulting H18 chambers have 55 or 56 facets for a doubled left parent,
and 43 for a doubled right parent. Every pair passes the exact acute and
Coxeter-angle gates. Finite volume follows from their decomposition into
three copies of the established parent chambers. These are derived gluing
examples, not claims of worldwide novelty.

For the original seam-crossing H17, H16, H15, and H13 faces, an existing
forbidden-angle witness remains unless one of its two facet normals was
removed by the double. Exactly one witness is removed in the complete
first layer: doubling left wall `L18` affects the H13 face
`{L2,L3,L5,L19,L21}`. Its new 26-facet section still fails the Coxeter gate:
the retained pair `L1,L20` has exact normalized pairing square `1/36`.
All other 12,786 witnesses survive unchanged. Thus none of these 11 new
construction classes repairs the targeted lower-dimensional face angles.

An additional structural observation suggests why further ordinary doubles
are unlikely to help here: all 22 original left norm-two walls have a
label-three neighbor among those same 22 walls. Such a core cannot be the
first wall removed by an even doubling. After the one replacement above,
every archived rejection uses only these norm-two walls. This observation
has not been packaged as an independently audited arbitrary-depth theorem.

The scope is all relative direct twists and all single-parent even doubles
preserving the original seam, followed by the stated original seam-face
tests. It does not classify all Coxeter faces of the new chambers or all
possible gluing/unfolding constructions.

Replay:

```sh
PYTHONPATH=src:. python3 scripts/fresh_gluing_h18_twists_doubles_20261013.py
```

The script is
[`fresh_gluing_h18_twists_doubles_20261013.py`](../scripts/fresh_gluing_h18_twists_doubles_20261013.py);
the exact result archive is
[`fresh_gluing_h18_twists_doubles_20261013.json`](../data/enumerations/fresh_gluing_h18_twists_doubles_20261013.json).
The original parent sources and arithmetic conventions are recorded in
[`frontier_literature_2026_10_13.md`](frontier_literature_2026_10_13.md).
