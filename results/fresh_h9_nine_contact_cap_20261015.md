# A nine-contact cap in the original quartic field

An exact wall cuts both the 29-wall chamber's escape and the escape that survives the two separate octic alternatives. It has norm two, seven right contacts, two label-three contacts and twenty dotted pairings. All 30 facets remain actual facets, and the wall creates one finite endpoint.

The [producer](../scripts/fresh_h9_nine_contact_cap_20261015.py) and [certificate](fresh_h9_nine_contact_cap_20261015.json) solve nine contacts on
`[0,1,2,5,6,8,10,28,3]`, with pairing minus one at walls 8 and 3 and zero at the other seven. Write the solution as \(p-tv\), with \((p,v)=0\) and the last coordinate of \(v\) equal to one. The norm equation gives

\[
t^2=\frac{2740+1284\sqrt5}{121},\qquad
t=\frac{5+7\sqrt5}{11}\sqrt{6+\frac{16\sqrt5}{5}}.
\]

Thus nine contacts suffice to reconstruct a wall in the existing quartic field. Requiring ten linear contacts would have missed this recipe.

The [independent checker](../scripts/fresh_h9_nine_contact_cap_20261015_review.py) and [review certificate](fresh_h9_nine_contact_cap_20261015_review.json) use separate SymPy field arithmetic and outward sign bounds. They reverse the contact solve and verify the quadratic identity, all pair relations, 900 facet signs, all endpoint signs and both escape cuts. Generation and saved replay pass.

The [complete extension census](fresh_h9_nine_contact_extension_20261015.md) finds a new escape, so this addition does not produce finite volume.

```sh
python scripts/fresh_h9_nine_contact_cap_20261015.py --verify
python scripts/fresh_h9_nine_contact_cap_20261015_review.py --verify
```
