# Independent review of the 35-wall cap and compact three-face

**PASS:** the new cap gives 35 actual facets, a complete census of 98 finite
vertices and 263 finite edges, and a compact tetrahedral three-face. A new
strict spacelike escape proves that the chamber still has infinite volume.
The same certificate excludes every coherent escape-cutting cap perpendicular
to the six-root support of that compact face. This is a fixed-support
obstruction, not an exclusion of arbitrary completions.

The [reviewer](../scripts/fresh_nonarithmetic_escape34_20261016_review.py) and
[certificate](fresh_nonarithmetic_escape34_20261016_review.json) read the
[primary certificate](fresh_nonarithmetic_escape34_20261016.json) without
importing its producer. They use the frozen independent octic arithmetic and
the preceding independent radical-extension implementation, with a separate
outer radical class and sign comparison. Exact reverse-order matrix
elimination and positive-definiteness checks come from the
[preceding review](fresh_nonarithmetic_wall3_cap_20261016_review.md).

The review checks the literal normal, its norm two, all 34 old-wall contacts,
and all 35 relative-interior facet witnesses. There are eight right contacts,
one label-three contact, and 25 dotted contacts. The old escape is cut, while
all 95 old finite vertices are strictly preserved. The complete nine possible
new-facet vertex supports yield exactly three new finite vertices. The
independently checked predecessor census therefore gives the stated complete
finite census, with 356 open finite-corner germs.

For the white support `[0,1,2,3,5,12]` and boundary `[6,7,8,34]`, the reviewer
reconstructs the projected boundary Gram. All its three-root principal
submatrices are positive definite and every entry of its full inverse is
negative. Its four negative Gram-dual rays are future timelike and satisfy
all retained wall inequalities, with precisely the required three active
boundary walls. Their positive cone is the whole tetrahedral cone; hence all
remaining walls are redundant on this compact three-face. No assertion that
its intrinsic dihedral angles are Coxeter angles is needed or made.

The replacement escape is independently reconstructed from the old-wall
support `[0,1,2,4,5,6,7,9,17]`, normalized against the center, and perturbed
by `H/8192`. It is strictly spacelike, in the future affine chart, and strictly
negative against all 35 normals. Its projection perpendicular to the six
white roots is a strictly positive combination of the four finite dual rays;
the reviewer checks this identity and every coefficient exactly.

To justify the all-real obstruction, let `u` be any coherent added acute wall
and choose a common strict future timelike interior point `h` of the retained
walls and `u`. At an old finite vertex, its positive-definite incident Gram
has entrywise nonnegative inverse. Projecting `h` to that vertex ray writes
a positive multiple of the vertex as `h` plus a nonnegative combination of
the incident outward normals. Pairing with `u` is strictly negative because
`(u,h)<0` and every old-new product is nonpositive. Thus every coherent acute
addition strictly preserves these finite vertices, even if the old selected
center is removed. If `u` is perpendicular to the six-root support, its
pairing with the new escape equals its pairing with the positive finite-ray
combination, and is strictly negative. It cannot cut that escape.

Generation and fresh saved replay passed. Source and certificate hashes are
bound by the certificate. The inherited chamber geometry and predecessor
finite census are explicit dependencies; no complete census of exterior
rays, finite-volume result, or new dimension is asserted.

```sh
python3 scripts/fresh_nonarithmetic_escape34_20261016_review.py
python3 scripts/fresh_nonarithmetic_escape34_20261016_review.py --verify
```
