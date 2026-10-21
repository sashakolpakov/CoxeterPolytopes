# Clarification of the B-row bound in the A-endpoint reduction

19 October 2026. Paper-only correction; no scientific job was run.
The two source-bound A-type reports are unchanged.

In Section 4 of
`arithmetic_an_same_endpoint_quad_reduction_20261019.md`, substituting
the B-row parameters `t=1/N`, `h=sqrt(2)/(2N)`, and `x=1/2` into
equation (1) gives

    y=(-N+6+2sqrt(2))/(2(N-3)).

The table prints the strictly larger quantity
`(-N+7+2sqrt(2))/(2(N-3))`. It is a valid upper bound, rather than
the exact value at the minimum. Both are negative for N=17,18, so
the branch exclusion and the full A-type report's use of that
exclusion are unaffected. The sentence introducing the table should
be read as “upper bounds,” with “exact” removed for this row.

This discrepancy was independently identified by resource_cleanup
while reviewing the full A-type consumer proof. It does not settle
the remainder of that independent review.
