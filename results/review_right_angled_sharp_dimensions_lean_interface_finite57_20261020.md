# Independent source review: sharp right-angled conditional Lean logic

20 October 2026. Full source read only; no Lean compilation or printed
axiom audit was executed in this lane.

Reviewed `formal/CoxeterProofs/RightAngledSharpDimensions.lean`, SHA256
`e1844d73a3f07de6cc512e514b948dafef4253e272b3f98de36148543554e323`.

**PASS for the explicitly conditional interface and count logic.** Actual
kernel acceptance remains a separate scheduled verification. No source
edit is requested.

Only dimensions, facet counts and the descent index use natural numbers.
The integer lists describe exact polynomial coefficients and the doubled
four-vector Gram matrix. `heightBound k F q` is an opaque real-geometric
predicate, not an integer substitute for a height. Its meaning and the
unformalized strict bound are stated directly in the module header.

`ProjectionGeometry` explicitly supplies actual finite interior points,
their finite starting bounds, actual facet projection, inherited right
angles/finite volume, and the next real bound. The three-step theorem
applies this field at `(k,n)=(0,8),(1,7),(2,6)`, composing dimensions
9→8→7→6 and bounds 0→1→2→3. It does not assume a final six-face or
nonexistence conclusion. The actual-facet witnesses are obtained at each
step even though only the terminal data are needed by its conclusion.

The integer polynomial multiplication checks the factored cubic and
its normalized-kernel coefficient expansion correctly. The diagonal
value 27 and positive coefficients are exact integer checks. The concrete
doubled four-vector Gram sum is minus four, so its asserted all-ones
nonnegativity is inconsistent by integer arithmetic.

The important real boundary is explicit: `WrittenCubicGeometry` supplies
the tensor-count inequality `N*N<=27*N` and the Euclidean Gram-positivity
consequence at equality. Neither the full tensor contraction nor the
real equality argument is claimed to be formalized by those small
integer checks. `StrictCodeSixConsequences` does not assume the final
bound 26; the code derives it by cancellation, the exact four-vector
Gram contradiction and integer separation of 27 from 26.

The universal-eight theorem takes only the published Dufour facet input;
it has no Musin or Kolpakov assumption. Conversely the ideal-four theorem
uses the separate Musin/Kolpakov structure and opposite-pair center
geometry, without the Dufour or cubic structures. The shared radial
structure contains the two actual code reductions, not external
cardinality bounds. The abstract center/strict-height incompatibility
fields remain transparent local real-geometric inputs.

Actual nine-face and entirely ideal five-face inheritance are exposed
at the higher-dimensional endpoints. No field assumes no nine-polytope,
no ideal five-polytope, or either final dimension bound. The source imports
only `Std`, has no new `axiom` or `sorry`, and requests printed audits of
its twelve public conclusions. This is not a claim that the actual
geometry or published theorems have been kernel formalized.
