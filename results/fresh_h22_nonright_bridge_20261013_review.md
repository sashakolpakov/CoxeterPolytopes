# Review of the nonright H22 endpoint search

The exact spherical-profile replay passed on 2026-10-13: 1,951 one-contact,
3,706 two-contact, and 823 three-contact profiles, with 217 unbounded
dihedral families kept explicit. The command was

```sh
PYTHONPATH=src:. python scripts/fresh_h22_nonright_bridge_20261013_profiles.py --verify
```

Source review confirms that the edge inverse is the Schur complement of
the released diagonal entry in the full vertex inverse. Prescribed pairings
use the actual old root norms, including norm four. The extension test is
exact positivity of `2 - p^T G^-1 p`; nonnegative inverse entries justify
single-contact pruning. The finite spherical diagram classification justifies
at most three nonright contacts, only label-three trivalent nodes, and an
isolated old A1 for each unbounded dihedral tail. This is an exhaustive
abstract profile list for these 88 fixed supports.

The geometric scout uses the correct separate residual normalization:
`r=s+xv+ye`, `-x²+y²=R`, with
`x=sqrt(R)(t-1/t)/2` and `y=-sqrt(R)(t+1/t)/2`.
The anchor equation uses `-sqrt(2*q_anchor)*cos(pi/m)`, respecting the
anchor's norm. Both roots of each ordinary quadratic are considered.

The geometric failures remain bounded numerical observations. Fixed
tolerances, near-degenerate quadratic handling, the finite anchor alphabet,
and the sampled dihedral labels prevent an exact exclusion theorem. The
report states those limitations and claims no new wall or finite-volume
chamber. This review is a source audit and producer replay; it is not a
second implementation of the profile arithmetic or a rigorous replay of
the 85,889,220 numerical proposals.

Verdict: **PASS for the exact abstract profile census and stated scout
scope; no general nonright-bridge obstruction established.**
