# Root review: local ideal-edge corners in mixed three-faces

19 October2026. **Full written mathematical review PASS**, under the
explicit inherited actual-wall and attachment inputs. I read the full
new report and the relevant earlier actual-move, symmetric-quad and
weighted-triangle arguments. No census, optimizer, checker, Lean or
TeX process ran.

Source `right_angled_mixed_nu11_all_ideal_corner_exclusion_research_20261019.md`,
SHA256 `fc047b78e9a8fc9eaf283bd5d689b69d72a8bbf0976a4765600ffbe8625955b0`.

At a finite three-face corner with three ideal incident edges, the
three ideal neighbors are distinct. Every incident polygon that is
not a bad triangle contributes at least one unit of polygon excess.
Each bad triangle has its unique opposite wall `q-ni-nj`. Its
ideal–ideal edge partner cannot also be bad unless the entire
three-face is P3. The partner walls for different pairs are distinct
and have strict product -1 with q, so are distinct from all incident
facets. Thus the three units are counted on distinct facets. No
assumption about other finite–finite edges is used.

I checked the Euler elimination. At `nu<=11`, excess at least three
and `p>=1` force the sole non-P3 candidate `(F,c,p)=(7,4,2)`.
Neither finite vertex can have a finite neighbor, so the already
reviewed seven-facet/four-cusp rigidity applies to the whole
three-face and gives D3. The former `(8,3,6)` allowance has excess
two and is excluded directly, without classifying its other corners.

For the unrestricted error bounds, different induced-path corners
use distinct finite–ideal edges and distinct bad triangles. Sharing
an ideal–ideal edge would invoke the explicit P3 attachment theorem.
This gives `4Z<=nu-2`. All unbacked paths now have `nu>=12`,
and the integer endpoint at twelve gives `sum M<=2E12`.
The companion bound follows from `3A3<=4c` and `p>=A3`, yielding
`sum A3<=5E12`. Removing the earlier finite–finite-edge terms is
therefore valid; the old terms were conservative allowances rather
than known examples.

For the strengthened conditional corollary, every triple at the
initial all-ideal-edge vertex is P3 or D3, so its bad-graph complement
is a matching. The exact P3 identities preserve cross edges, and the
quad identities preserve the missing-pair cross edges. In either
move every untouched and new direction has a surviving bad neighbor,
so its incident edge is ideal at the opposite vertex. This proves
closure of the actual finite move component inside this class. It
does not assume global absence of finite–finite edges.

Weights one on P3 moves and one half on each missing-pair quad are
undirected. Their actual displacement coefficients are the printed
`a,b_un,b_pair`. I recomputed them from triangle multiplicities and
the half-weight quad displacement. The printed causal margins and
strictly positive time coefficient apply to every matching size in
dimensions9–12, with the nonzero null endpoint at n9/r0 handled
explicitly. Finite undirected cancellation contradicts the sum of
nonzero future causal vectors on the closed component.

Consequently, under global `nu3<=11`, every finite vertex in
dimensions9–12 has a finite–finite incident edge. The unrestricted
dimension ceiling remains12: high-nu carriers and corners with fewer
ideal directions remain unresolved. The new implications are reviewed
written mathematics, not new script results, revised TeX/PDF claims
or geometric Lean theorems.
