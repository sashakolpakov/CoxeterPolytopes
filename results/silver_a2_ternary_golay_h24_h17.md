# Silver `N(A2^12)` H24 continuation through real height 17

## Outcome

Neither of the two silver-field hosts

\[
 \langle-\sqrt2\rangle\oplus N(A_2^{12}),\qquad
 \langle-2\sqrt2\rangle\oplus N(A_2^{12})
\]

is finite-volume at the exact height-17 frontier.  Beyond the twelve affine
caps found in the height-10 screen, both hosts acquire the same complete
packet of 3,960 dominant norm-two walls at

\[
 H=8+6\sqrt2,\qquad y^2=18+12\sqrt2.
\]

All 3,960 walls are mutually Vinberg-acute: among the 15,677,640 ordered
distinct pairs, 15,677,508 pair negatively, 132 are orthogonal, and none
pair positively.  Nevertheless an explicit strict hyperideal point survives
the affine and 3,960-wall packets.  This is a deeper new host-family screen,
not a new compact example and not a global nonreflectivity theorem.

## Exact glue census

Write the two coefficients of every `A2` Dynkin-weight pair modulo three.
The rational and radical residue words must both lie in the extended ternary
Golay code.  Its 729 squared ordered word pairs collapse to only 288 exact
3-by-3 contingency tables; their multiplicities sum to

\[
 729^2=3^{12}.
\]

For each contingency table an exact dynamic program records the two
coefficients of the accumulated spatial square.  Six contingency types reach
the target after the affine-cap inequalities are imposed, giving exactly
3,960 labelled dominant vectors.  The same program proves that the preceding
control at `7+5*sqrt(2)` is empty.  Every other post-affine control below the
cutoff is empty or was already exactly rejected in the height-10 audit.

The same-height acuteness audit represents the 3,960 vectors by their 24
integral rational and radical Dynkin-weight coefficients.  It evaluates the
two coefficients of all pairings with integer matrix products.  Exact signs
are decided by comparing `a^2` with `2b^2`; floating-point arithmetic is not
used for the sign certificate.  The maximum distinct pairing is zero.

## Strict surviving points

The largest sum of all 24 Dynkin weights among the second-packet roots is

\[
 9+9\sqrt2.
\]

For the coefficient `a=sqrt(2)`, put the common Dynkin weight

\[
 c=\frac{\sqrt2}{5}
\]

in both coordinates of every `A2` block.  At time coordinate one its exact
Lorentz norm is

\[
 24c^2-\sqrt2=\frac{48}{25}-\sqrt2>0.
\]

Its maximum pairings with an affine cap and a second-packet wall are,
respectively,

\[
 -\frac65-\frac35\sqrt2,qquad
 -\frac25-\frac15\sqrt2,
\]

so it lies strictly inside every displayed halfspace.

For `a=2sqrt(2)`, take instead `c=2/5`.  Its norm and the two maximum
pairings are

\[
 \frac{96}{25}-2\sqrt2>0,qquad
 -\frac{12}{5}-\frac65\sqrt2,qquad
 -\frac25-\frac25\sqrt2.
\]

Thus both complete prefixes contain an open hyperideal patch and have
infinite volume at this cutoff.  Later walls may still cut these points.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_a2_ternary_golay_h24_h17.py
PYTHONPATH=src pytest -q tests/test_silver_a2_ternary_golay.py
```

The frozen archive is
`data/enumerations/silver_a2_ternary_golay_h24_h17.json`.
