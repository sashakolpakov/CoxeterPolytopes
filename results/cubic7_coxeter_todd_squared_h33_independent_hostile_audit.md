# Independent hostile audit: cubic-seven K12^2 H33 cutter

Terminal verdict: **PASS**.

This clean-room script imports no production search code.  It reconstructs
the cubic arithmetic and local norms `2,6`, the raw K12 Gram and ternary
radical gates, all 11 controls, and all five surviving necessary Gram types.
It independently reproduces the three earlier cross maxima.  At H33 it
checks both types: the relaxed cross maximum against the other type is
`-2-8theta-4theta^2`, while the cutter-type enumeration reproduces 8,173
Schur blocks, 6,960 minimum-four blocks, and maximum
`-4theta-2theta^2`.

The audit found a real proof gap in the first production certificate: its
`36/25` directional ledger omitted the other H33 type, whose relaxed bound
would be positive.  Production is repaired to scale `7/5`.  Its norm is
`196/75-2theta>0`, and exact necessary-Gram maxima give strict negative
pairings at all five types, including both H33 types.  This terminal PASS is
for the repaired bounded certificate only; it makes no compactness,
finite-volume, or Vinberg-termination claim.
