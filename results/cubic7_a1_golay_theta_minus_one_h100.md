# Cubic-seven `A1^24` Golay unit seed: exact H=100 prefix

For

\[
L=\langle-(\theta-1)\rangle\oplus N(A_1^{24}),
\qquad \theta=2\cos(2\pi/7),
\]

the exact prefix through Vinberg height 100 has two and only two nonempty
positive-height packets.  It does **not** close: after both packets the cone
has 218,103,483 hyperideal extreme rays and no cusps.

The complete primitive mirror norms are `1,2`.  Indeed, the timelike
coefficient is a unit and the Niemeier block is unimodular, so a primitive
pairing ideal is the unit ideal and the norm ideal divides `(2)`.  The prime
`(2)` is inert.  The Minkowski bound is `14/9`, and the displayed full unit
signature rank makes totally positive units squares.

The complete control census contains 21 levels and 1,329 unordered square
partitions.  Testing all seven nonzero XORs of the three coefficient-parity
words against the Golay weights `0,8,12,16,24` rejects every pattern except:

- norm 1, control `2+5theta+2theta^2`, height
  `5+14theta+6theta^2`: one spatial coordinate `2theta+theta^2`, giving 24
  coordinate-cap walls;
- norm 2, control `2+8theta+4theta^2`, height
  `6+18theta+8theta^2`: an ordered pair of distinct coordinates
  `1+2theta+theta^2` and `2theta+theta^2`, giving `24*23=552` walls.

Both survivors have zero coefficient residues, hence explicit root-lattice
labelings.  The 552-root packet is mutually acute and acute against the 24
caps.  Its active spatial rows have exact rank 24.  The symmetric active ray

\[
(1;(\theta-1)^{24})
\]

has norm `49-97theta+48theta^2>0`, so it is an exact hyperideal extreme ray.
An exact vertex classification gives 218,103,784 total rays: 301 finite,
zero ideal, and 218,103,483 hyperideal.

This is an exhaustive height-prefix result, not a global nonreflectivity
claim.

- Backend: `src/coxeter_polytopes/search/cubic7_a1_golay.py`
- Audit: `scripts/audit_cubic7_a1_golay.py`
- Tests: `tests/test_cubic7_a1_golay.py`
- Archive: `data/enumerations/cubic7_a1_golay_theta_minus_one_h100.json`
