# Root review: isolated finite A-type continuations

20 October 2026. Complete independent written review of
`arithmetic_full_an_isolated_finite_continuation_packing_20261020.md`,
SHA256 `37fdc15a140a6d1649f940830c27be6c7108f09636faacffcc3185f437424f11`.
PASS for the specified isolated finite continuations. No execution or
dimension ceiling is claimed.

The revised ground-field normalization is valid: rescale roots within
K along the simply laced chain, then scale the quadratic form by a
totally positive scalar. The opposite root normalization has norm2x
and product-x, without adjoining sqrt(x). Its full-rank Gram gives
Schur complement x*(2-x*d_i). A proper extension of Q(x) would have
a nonphysical real embedding preserving this Lorentzian Gram, so the
actual ground field is exactly Q(x).

The elementary totative argument gives all possible bulk labels. I
checked the three cases of m, including the alternate exponent for
m=6k with k=2 mod3. The indicated exponent is coprime to m and lies
strictly between one and m/3. It yields a nonphysical conjugate greater
than one and contradicts conjugate positivity. The positive affine null
combinations then correctly exclude labels3,4,6 in the stated rank
range. The three remaining conjugates give exactly the listed windows.

I independently derived the shared timelike-residual pair formula and
the fact that its twice-unit product lies in the ground field. The
quadratic field restricts a finite pair to the displayed labels. The
AM--GM bound excludes adjacent indices. For distance two I independently
obtained equation (6) and every coefficient comparison. The exceptional
formal label12 solution N=8 is inconsistent even with its required
inverse-Gram entry. No pair at distance two survives.

The counting interval means the full contiguous set d_i>2, rather than
the potentially disconnected set left after a label window. The weighted
recurrence retains the original path indices and handles unavailable
ones by skipping them. These corrected interfaces are sound.

The theorem does not bound all short edges: isolated dotted, one-sided,
bridge and small-support continuations remain outside it. Shared closing
walls cannot be counted as multiple isolated walls. The source keeps
these limitations explicit and changes no uniform dimension ceiling.

Final scope correction reviewed on20 October: current source SHA256
`9ab34f073ea95bed2040bafe6b8573331721f9b1e75d65172f814637da046810`.
An isolated contact is counted at its old CONTACT INDEX only. The same
wall can also close another outgoing edge, where its continuation need
not be isolated. The revised opening and closing paragraphs say this
explicitly. The spacing theorem and proof are unchanged and retain
PASS; they do not count all directions closed by these walls.
