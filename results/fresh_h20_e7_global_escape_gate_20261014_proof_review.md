# E7-face escape gate: independent proof review

**2026-10-14 — PASS.** No issue was found in the [primary proof](fresh_h20_e7_global_escape_gate_20261014.md). The primary saved replay and the [separate symbolic checker](../scripts/fresh_h20_e7_global_escape_gate_20261014_proof_review.py) both pass. The [review certificate](fresh_h20_e7_global_escape_gate_20261014_proof_review.json) binds the frozen source and data.

The two lifts exhaust the stated scope. The previously independently checked E7 Schur gate permits only the all-right profile or a single label 3 at wall 33 when the projection is spacelike. The twenty corner roots and their orthogonal timelike vector form a complete basis. Consequently a norm-two wall inducing an actual facet of this E7 face has exactly the displayed form, with projected norm 2 or 1/2. For its thirteen remaining contacts, the relaxation `{0,1} union [sqrt2,infinity)` contains every Coxeter label, including parallel and dotted pairs.

The positive-time conclusion does not require an additional apex-preservation hypothesis. The full `E7+A11+A2` corner has an entrywise nonnegative Gram inverse, independently reconstructed here. For a proposed common strict future center `p=uT-sum b_i Omega_i`, one has `u>0` and `b_i>0`. An acute proposed wall has `r=tT-sum c_i Omega_i`, with `c_i>=0`. Thus `r.p=-2tu+c^T G_corner^-1 b` cannot be negative when `t<=0`. This proves `t>0` and hence strict retention of the old finite apex.

The right-lift positive dual has the correct inequality direction: positive multiples of retained-wall inequalities plus fourteen strictly negative residual coefficients give a strictly negative escape pairing. For the nonright lift, the independently verified substitution identity is

`B*(r.oldwall) = (boundary_row.d-rhs) + 2*z*(r.E)`.

The two relevant values of `z` are positive in the frozen physical reconstruction. Under `r.E>=0`, acuteness therefore implies both saved linear inequalities. Their positive combination forces `d12<1`, hence `d12=0`; the wall-16 inequality then forces `d0>1`, hence `d0>=sqrt2`. Both discrete gaps are used in the correct direction, with no upper label cutoff.

The rank-one positive-definiteness condition makes the projected-norm lower bound strictly convex. The saved positive denominator and eleven positive constrained-gradient numerators prove the proposed point is its global minimum on the relaxed domain. Separate symbolic differentiation reproduces the free-coordinate minimizer, the value `11/6-2C^2/D`, and the strict margin `(4/3)D-2C^2>0`. Thus the minimum exceeds 1/2 and contradicts the required nonright projected norm. The argument works for `r.E>=0`, so it excludes zero pairing as well as positive pairing.

This review independently reconstructs the rational corner and symbolic optimization identities, and checks the saved strict coefficient signs using a separate recursive Fraction implementation. It does **not** independently reconstruct the ambient field coefficients; the [separate arithmetic review](fresh_h20_e7_global_escape_gate_20261014_review.md) does that with sparse Fraction arithmetic and outward rational intervals. The E7 profile census is a source-bound dependency on the [earlier independent face review](fresh_h20_e7_face_20261014_review.md).

The conclusion is limited to additions producing actual facets of the selected E7 face while the 87 old walls are retained. Every such wall strictly retains the same future spacelike escape. A finite, nonempty acute completion consisting only of these additions therefore has infinite volume. Walls outside this face-facet class, and replacements of old walls, remain outside the exclusion; this is not a global noncompletion theorem.

```sh
PYTHONPATH=src:. python scripts/fresh_h20_e7_global_escape_gate_20261014.py --verify
PYTHONPATH=src:. python scripts/fresh_h20_e7_global_escape_gate_20261014_proof_review.py --verify
```
