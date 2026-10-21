# Review of the retained-cusp field obstruction

**PASS.** The [standard-library Fraction checker](../scripts/fresh_h22_incidence_field_obstruction_20261013.py) replays the [frozen certificate](fresh_h22_incidence_field_obstruction_20261013.json) successfully. Source and proof review found no scope or algebra issue. This review is a fresh invocation and mathematical review of that implementation, not a separately coded certificate producer.

The physical null ray is `δ=(a,2√11)` in the form `diag(Q,1)`. Its old part a is a strictly positive combination of **all** 210 old normals: uniform coefficient `1/10000` plus 22 positive excess coefficients. Those 22 normals are independently checked to span the full old 22-dimensional space. Thus a new normal nonpositive against every old normal can have `(a,r_old)=0` only when `r_old=0`.

If a projective representative of r has coordinates in a real field F not containing √11, orthogonality to δ implies both `r_last=0` and `(a,r_old)=0`. Otherwise their ratio would put √11 in F. The positive-combination identity then forces r to vanish. Rescaling does not change this argument or the zero pairing. The field may differ between proposed walls; a common field for the entire extension is not required.

Crucially, the retained norm-two parallel pair satisfies `δ=c4+old26`. Every acute-compatible new wall therefore has nonpositive pairing with δ. A genuine added wall over one of the allowed fields must have **strictly** negative pairing; it cannot cut the cusp away or add an incident facet. Redundant copies of retained walls change nothing and are not needed in this argument.

The inherited complete-incidence certificate is hash-bound. The checker reconstructs its 18 incident normals, checks normal rank 18 and Gram rank 17, and verifies a positive-norm vector orthogonal to δ and all incident normals. Hence the horospherical quotient, of dimension 21, has nonzero recession. All 564 incidence-core/binary walls pair strictly negatively with δ by the exact bound `-25+4√11<0`, so they leave this link unchanged. Any further **finite** collection of permitted-field walls is likewise absent sufficiently deep in the cusp and cannot give finite volume.

The theorem retains the old 210 facets and c4 and imposes acute compatibility with them. It does not forbid walls whose projective field contains √11, removal of those retained facets, or arbitrary H22 reflection polytopes.

```sh
python scripts/fresh_h22_incidence_field_obstruction_20261013.py --verify
```
