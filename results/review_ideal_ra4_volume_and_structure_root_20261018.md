# Root mathematical review of the consolidated ideal-four manuscript

18 October 2026. **Pass for mathematical consolidation, relative to the
stated published inputs and the previously reviewed component proofs.**
Source: `manuscripts/ideal_ra4_volume_and_structure.tex`, SHA256
`e503ad2cb8edafd827f735a62568a3372407e9dfdaf3f50d09b05fb3bc75f025`.
The complete source, including all five TikZ figures, was read. This is
a written mathematical and translation review, not a new script run,
PDF inspection, historical classification proof or Lean compilation.

The paper now has one primary chain of reasoning. The published minimum
and volume formula are identified as such. Defect arithmetic and the
opposite-cusp injection exclude all intermediate volumes. The equality
section independently forces exactly eight D9 facets and octahedral
remaining facets. Only after that reduction does the general Oct/D
theorem identify a two-cell chain and the reflection double.

The full minimum-defect case split survives consolidation. In particular,
the square-antiprism positive-neighbor count, the minimum-four matching
obstruction, the large-minimum finite-vertex contradiction, both N=7
branches and the cuboctahedral medial reconstruction have not been
replaced by an equality assertion or an unbounded census. The cusp-width
calculation supplies the nineteen-facet injection contradiction in the
last N=7 branch. The nine- and ten-vertex classification entries and
three-dimensional rigidity remain explicit external inputs.

The general argument retains the actual closed triangle surface, the
primal/dual-tree proof of chi≤2, the seven-facet lower injection, the
simple six-cusp graph and the unit-normal cube equations. Convexity
establishes the actual octahedral section. The cut proof addresses
outside facets, nonempty halves, absence of new finite vertices and
right angles. The termination count and geometric disjointness are
retained, as are direct all-octahedral recognition, the D4 common-neighbor
obstruction, path structure and the no-twist placement argument.

I separately checked the count formulas: m cells give
(18m+6,72m+24,68m+28,14m+10), with8(m−1) doubled facets and
6m+18 ordinary facets. At m=2 these give (42,168,164,38), thirty ordinary
and eight doubled facets. The double's halfspace-union argument now
explicitly proves containment using the reflected normal product; its
converse restores the gluing inequality on each side. Enumeration of
the finite model vertices remains an explicitly computational input.

The consolidation removes the former equality-specific second copy of
the cube-normal and cutting argument. It does not remove the independent
equality-spectrum reduction needed to apply the general theorem.
Duplicate equation tags were removed. No unrestricted decomposition,
universal24-cell-volume divisibility, compact120-cell minimum or new
dimension bound is asserted.

The new source is unfrozen and its first PDF is pending the resource
gate. The separate old frozen chain paper is preserved. Conditional
Lean sources and their assumption register are prepared but uncompiled;
the written geometric interfaces have not been instantiated as concrete
hyperbolic objects in Lean. An independent consolidation review is a
separate record, not attributed to this root read.

Subsequent execution: the same source compiled to16 pages/five figures
with zero warnings, and the core script/conditional-Lean package passed
in `ideal_ra4_unified_core_verification_final_20261018.json`. The separate
independent consolidation review also passes. These later executions do
not expand the geometric scope of this written review.
