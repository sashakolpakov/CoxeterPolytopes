# Root review: golden isolated-five stars and invariant lattices

20 October 2026. Complete written review; no mathematical worker,
enumeration, Lean or TeX execution.

Source: `results/golden_b_isolated_five_invariant_lattice_interface_20261020.md`,
SHA256 `e0e347a7642ca196ec416cc8889adf666bf4903a3ebfe865e9d65549257b4e10`.

**PASS, with the published literal cutoff G as an explicit input.**
This removes two specified vertex-star positions, not a dimension range.

I checked both caps separately. Their norms are two; their time and
nonzero spatial coordinates are units. Signed permutations and either
reflection therefore generate the algebra containing every `2E_ij`,
including the time row and column. Coordinate-content scaling over the
PID really does reduce every invariant lattice to `2L <= Lambda <= L`.

Modulo two, both cap formulas have the stated rank-one image. For the
five-coordinate cap, differences of five-subsets give every spatial
coordinate difference when n>5. Their common kernel is exactly D and
their image span is exactly E. An invariant subspace outside D contains
one entire image line and all its permutation translates, hence E.
This proves the four-item list even when D lies in E; no extra lifts
are lost because each lattice is the full inverse image of its subspace.

The norm-hull argument is valid at the inert prime (2):
`f(d)=n-1-2phi` has valuation zero or one, so a nonintegral element of
`Lambda_D/2` has nonintegral norm. Conversely the characteristic vector
is fixed by every isometry of L, since in F4
`B(d,y)=phi^2*t+sum x_i=sqrt(f(y))`.
Thus the two orthogonal groups are equal, and taking lattice duals
establishes the E case. This is actual group equality, not a transfer
of reflectivity along commensurability.

Finally, the full-rank star normalization gives time coefficients
`phi^-1` and `phi^2` at positions one and five. A larger totally real
ground field would duplicate its indefinite embedding. Arithmeticity
supplies an invariant lattice for the entire reflection group, so the
group containment needed to apply G follows in the same coordinates.
The full integral chamber refining a compact chamber is compact; G
then contradicts its existence.

I read Vinberg, *Hyperbolic reflection groups* (1985), Introduction,
printed p.38, and the displayed golden family in section5.6, p.64:
https://www.mathnet.ru/php/getFT.phtml?jrnid=rm&option_lang=eng&paperid=2140&what=fullteng .
I also checked Kellerhals (2014), p.474, equation(4.3) and its following
sentence, which explicitly gives the same literal form and cutoff:
https://homeweb.unifr.ch/kellerha/pub/CMFT.pdf .
The original Bugaenko1984 proof was not retrieved. G is a named
published input; this review does not reprove it. Positions two through
four, other finite contacts and global compact bounds remain open.
