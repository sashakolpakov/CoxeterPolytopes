# Outside-author mathematical audit: literal H3+A6 end-edge period

21 October 2026. **Written-proof PASS at the explicitly inherited
finite-census and arithmetic-input boundary** for
`compact_h9_h3_a6_end_edge_period_obstruction_candidate_20261021.md`,
SHA-256 `b8985c4b7445f8e21a27775795491a740ba126850f7c7675233e938e104621d3`.
I independently reconstructed the profile, residue, near/far orientation,
integrality and full-chamber period arguments. No mathematical job,
scientific import, enumeration, Lean or TeX ran. No new example, global
dimension bound or broader lattice-class exclusion is asserted here.

The earlier 8839ead4 source was read first. The current source changes
only the ambiguous q1/q4 notation to R=1 and R=4 and supplies the
elementary discreteness paragraph; I read and checked both changes.
No mathematical blocker remains.

## Exact inherited boundary

All six source pins in the candidate match the bytes read:

| Input under results/, except scripts as marked | SHA-256 |
| --- | --- |
| fresh_h3_blocks_h9_report_20261013.md | `4545875776e23edf97c0a4a690d8eb0589bf76d3600fe6b7f5daf65081910b19` |
| fresh_h3_blocks_h9_a6_phi_h30_20261013.json | `ee0370bc5306d01d0f6954ae0c4759e1ead4194c34190e189ec06d88428b8d5c` |
| fresh_h3_blocks_h9_independent_audit_20261013.json | `d4acf35f13083ee44419c9e0bc7edd6cda2503d53ca93a7fbc07a38a5ff83863` |
| scripts/fresh_h3_blocks_h9_a6_search.py | `fcfe23a3beadf58242df05f0c91c0f1c93462b7247bedf2606f717609b4898b7` |
| scripts/fresh_h3_blocks_h9_independent_verify.py | `06124f9147cb6c70c43928033bc8c57227d80865f1deb3711098976fe2ddc8fd` |
| compact_h9_construction_reentry_frontier_20261020.md | `668a075c6001f0860d05fb848a8c2357de59ac133d286f8dbbf599c603a49df6` |

The 45-wall object is a bounded prefix, not a finite-volume chamber
certificate or an independently repeated positive-height census. The
saved independent audit repeats the zero census, controls and root/angle
checks; its own scope explicitly leaves full positive-height spatial
exhaustion source-reviewed. The new proof preserves that limitation.

I reviewed the producer's completeness bounds: both-embedding Cauchy
bounds enumerate every dominant covector; selected-prefix inverse-Gram
costs are valid Schur lower bounds; the additional full-inverse pruning
is used only in the physical embedding, where all entries and remaining
dominant covectors are nonnegative. Exact integral reconstruction and
the two block-energy join lose no candidate. The two-embedding integer
rectangle helper is conservatively bounded. These are source checks,
not fresh execution of the positive census or its transitive arithmetic.

For the first H3 edge, a wall meeting 0<=t<=1/2 satisfies
`phi^2*x^2 <= (phi^2/2)*(q+phi*x^2)`, hence height <=phi^3.
The stored levels below this have no candidates; at phi^3 the two
primitive q1 walls are exactly s and r, while the q4 copies are
imprimitive. Directly s is negative throughout the segment and r first
vanishes at P. The full initial zero-height germ plus this census gives
the actual complete P germ A6+3A1. The argument genuinely needs this
inherited complete-census input, rather than only 45 compatible roots.

## Norms, alphabet and finite profiles

Direct multiplication gives Gw=(0,1,0), w^2=2phi^2,
P^2=-1/(2phi), r^2=1, and pairwise orthogonality of a,b,r.
The time/H3 decomposition and its norm identity follow immediately.
The norm pool follows from primitive q|28, inert2 and7, and unit
rescaling. The PID and unit facts are applicable: coefficient rounding
has norm error <=5/16, and a totally positive unit in (1,phi^2) would
have difference from its conjugate strictly between0 and sqrt5,
impossible for its nonzero integral phi coefficient. Thus phi^2 is the
least totally positive unit >1 and all such units are squares.

I checked each finite contact entry using
`product^2=q*q_old*cos^2(pi/m)`, the field Q(phi), and integral
reflection divisibility. Labels6/10 and prime7 contacts do not produce
the required field squares. The apparent q4-to-norm1 products1,phi
fail divisibility by2. Multiple contacts to A5 would give a cycle in
the positive actual replacement and are forbidden.

Without an A5 contact the finite star condition is
`l^2+n^2+2p^2<2q`. At q2 this gives exactly the seven R rows printed,
including norm1-only4. At q1/q4 there are only the isolated and one4
rows; two4 contacts are affine. With an A5 contact, a rank>=6 finite
connected diagram is A/B/D/E. I reconstructed the entire table:
q1/q4 have only terminal B6; q2 has A6/D6/E6 without another contact,
A7/E7 with one3 on a norm2 A1, B7 with the norm1 terminal4, and
D8 with two3 A1 leaves. Central one-A1 attachment is affine E6;
internal two-A1 attachment has the affine two-fork. There is no
additional E8 or weighted-branch role. Compactness ensures the actual
first endpoint has a complete finite replacement, so this list applies
to the arriving wall, not an abstract extension.

## Residue filters and the sign-sensitive surviving branch

At7, every nonzero rational residue is a square in F49, while phi
is nonsquare because its norm is -1 in F7. This checks the rational-R
and rational-bracket exclusions. R=2-phi also has nonsquare norm.
For R=3-phi, reduction at sqrt5 gives T^2-2z^2=0 in F5; both vanish,
contradicting the exact valuation1 of -sqrt5. For prime7 norms the
orthogonal H3 coefficients force T=2x. The q7 parity, q28 primitive
content and q14 F9 norm arguments each check, including the stronger
divisibilities q|2x and q|14z.

I derived `Omega.s_j=-j`, `s_j^2=j(j+1)` and equation(3).
Its only remaining profiles have q2,d1,sigma0,j1 or5, T=2x and7|x.
For j1 the meeting condition is 6z-1>0; together with its norm
equation this forces z>1/3. Conjugation gives0<z'<1/3, so z is
totally positive. Since (z,3z-1)=1, every non7 prime valuation of z
is even and the7 valuation has the indicated parity. The PID and
totally-positive unit-square statement give z=s^2 or7s^2. In either
case s|x, and cancellation gives (x/s)^2=-phi modulo3. This works
even when s has arbitrary3-adic content; no coprimality assumption
at3 is hidden.

For j5, writing u=z-1 gives a DIFFERENT physical sign. Its meeting
condition is1+6u>0, its norm equation forces u>0, and conjugation
gives -1/3<u'<0. Thus the preceding total-positive-unit proof cannot
be recycled to exclude it. The current candidate correctly retains
this branch and does not claim that the edge has no finite endpoint.

## Integral full-chamber period and conclusion

Equation(4) gives Y^2-12x^2/(7phi)=1. With P^2=-1/(2phi) and
Omega^2=6/7 this proves both boost norms and its mixed product;
Y>1 gives the future sign. The two-sided lattice calculation is valid:
B(E1)-E1=7uOmega+2xP, and for an original time/H3 vector v,
`-2phi(v,P)=2phi^2*v_time-phi*v_H3,second` is integral.
Also6uP and (x/phi)Omega are integral because2P is integral and7|x.
The inverse changes only the boost sign x. Thus B and B^-1 preserve
the literal lattice, rather than just its quadratic K-space.

The integral A6 rotation sends E6 to-H. Hence S=BR fixes a,b,r and
sends the entire oriented old germ to `-E2,...,-E6,m`.
I checked m.B(P)=0 and B(P)=Y*(P+x/(phi*Y)*Omega), on the positive
edge ray. Since m is the actual first arriving facet, this is the
actual endpoint with the same complete outward germ. An integral
isometry preserves the full integral mirror arrangement; local
finiteness and the common complete germ therefore imply S(C)=C.
This step is not justified by the finite prefix alone: actual firstness
under the compactness assumption is essential and is supplied here.

The end-weight quotient `(Omega,R(Omega))/Omega^2=-1/6` and tr(R)=3
on the full ten-dimensional space give tr(S)=3+5u. The signs of u,u'
force a positive integral phi coefficient, hence u>sqrt5-1/3>5/3
and tr(S)>10. Thus S has infinite order. The new discreteness paragraph
is correct: the conjugate orthogonal group is compact, and bounding
both embeddings leaves finitely many O-valued matrices. A discrete
isometry subgroup stabilizing a compact set with nonempty interior
is finite, giving the required compactness contradiction.

Noncofiniteness additionally uses the explicitly named admissible
anisotropic arithmetic-lattice theorem: a cofinite subgroup of the
uniform integral group would have finite index and compact chamber.
I did not replace that published arithmetic input with a new proof
or claim a fresh independent scientific verification. Subject to these
stated inputs, the literal full integral host is excluded. Its unknown
first far mirror is handled conditionally for every surviving solution;
neither a Pell generator nor another lattice or commensurability class
is classified.
