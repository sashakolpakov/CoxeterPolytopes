# Exact union count for the frozen OG10 projected-wall families

The zero-`E8`, scalar-unit, and mixed-through-`b^2<=6` certificates contain
`4,178,946` family-indexed wall evaluations, but that sum is not a count of
distinct geometric covectors.  The scalar and mixed families overlap.

For `w=phi^n*v`, the coefficient of `phi` is `b=F_n*v`.  Thus `b^2<=6`
overlaps the frozen scalar family only in four nonempty signature-admissible
rows: `(n,v^2)=(1,2),(2,2),(2,4),(2,6)`.  Their wall counts are respectively
`9,600`, `5,760`, `86,400`, and `107,520`, totaling `209,280`.

The exact geometric-union count is therefore

\[
66+2{,}706{,}240+1{,}472{,}640-209{,}280
=3{,}969{,}666.
\]

The zero-`E8` family is disjoint from both nonzero-`E8` families.  This audit
changes no no-go theorem or sign census; it only prevents family-indexed
counts from being misreported as distinct covectors.  Reproduce it with
`scripts/audit_og10_golden_projected_wall_family_overlap.py` and
`tests/test_og10_golden_projected_wall_family_overlap.py`.
