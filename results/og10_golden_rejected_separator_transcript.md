# Rejected OG10 projected-separator transcript

This note records a negative regression, not a candidate.

The proposed data claimed

\[
v=(-2,-4,-6,-5,-4,-3,-2,-3),\qquad
w=(4+7\phi)v,
\]

but supplied the sole nonzero \(E_8\)-covector entry as \(8+14\phi\).
In the fixed simple-root convention, direct multiplication gives

\[
G_{E_8}v=-e_7.
\]

Consequently the claimed \(w\) requires the covector entry
\(-4-7\phi\), not \(8+14\phi\).

For the supplied covector literally, exact inversion gives

\[
q=2p^2=782+1259\phi,
\]

and the identity-place Schur complement is negative.  The extended form has
signature \((9,2)\), not \((10,1)\).  Its apparent opposite evaluation signs
therefore do not define a valid carrier.

With the sole \(E_8\) entry corrected to \(-4-7\phi\), exact reconstruction
does give

\[
q=2-\phi,qquad
\operatorname{sig}(H)=(10,1),qquad
\operatorname{sig}(H')=(11,0).
\]

Moreover \(2p\) is integral, the resulting involution squares to one, and
\(J^tHJ=H\).  But its exact evaluations are

\[
400\ell(x_{\rm int})=-10659-17409\phi<0,
\qquad
13\ell(v_{61})=-259-415\phi<0.
\]

Thus the corrected carrier does **not** cut the escape ray.  There is no
positive wall, chamber, or hyperbolic example here.

The standalone regression is
`scripts/audit_og10_golden_rejected_separator.py`, exercised by
`tests/test_og10_golden_rejected_separator.py`.  No candidate JSON is frozen
for this rejected transcript.
