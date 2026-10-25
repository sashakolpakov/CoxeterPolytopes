# Case IDs for the A1/D5 witness in the selected-affine audit

25 October2026. Author: `/root/audit_finite57_algorithms`.
This corrects one sentence in my own source review, not the underlying
mathematical result. The frozen report and its pin manifest are unchanged.

In Section6 of
`audit_finite57_selected_affine_completion_chain_fresh_20261024.md`
(34,339 bytes; SHA256
`e740b5ccf7ed66546ceca59507a3c3fdedd6179069732d12f73ca6dbb37c3555`),
the sentence “The last two rows contain whole unmarked A11” must read:

> The two A1/D5 cases, **2371 and2539**, contain a whole unmarked A11
> at their alternate corner(3,13).

The table immediately before the sentence already assigns the alternate
A1/D5 payment to2371/2539 and the D5/A12-end payment to1978/2535.
The mistaken positional reference does not change that table.

I read the original fixed certificate's exact components at(3,13):

| Case | Marked components | Unmarked witness |
|---|---|---|
|2371|A1 on[3]; D5 on[10,11,12,13,14], marked long endpoint13|A11 on[0,1,4,5,6,7,8,9,16,17,18]|
|2539|D5 on[0,1,2,3,4], marked long endpoint3; A1 on[13]|A11 on[6,7,8,9,10,11,14,15,16,17,18]|

Thus both marked roots are outside the A11 witness. The final Boolean
`A11 OR (A9 AND A2)` guard supplies the required3/8 at these two corners.
No A11 witness is needed for the distinct D5/A12-end payments in1978/2535.
The common full-core seal supplies whole-component status, including
arbitrary additional supporting whites, as stated in the frozen review.

Evidence: `fresh_d5_d12_alternate_coupon_rebalance_20261018.json`, SHA256
`f49cc9c383e700818ca09d21f1e2eae596cf119b3ce126d172554d9eb39549ef`,
lines2013–2078 and4074–4145, under the fixed2371/2539 proof records.
The separate later source
`fresh_order17_d5_alternate_half_zero_transfer_20261018.py`, SHA256
`122507a456cd86743f6f90d2daf656f0cc8d9dcb0b23fdbf1f390a9114e126fe`,
also tests exactly these case IDs and the component multiset A1+A11+D5
in its `local` function. Its stronger D113 rule is not imported here.

This was a source/certificate read and streaming identity check only.
No worker, classification, enumeration, test, Lean or TeX process ran.
Root identified the imprecise sentence; this author supplies the correction.
