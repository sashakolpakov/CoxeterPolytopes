# Outside-author review: the ramified silver facet bridge

**PASS, source-only and with the stated inherited silver-family input.**
I read the full source
`silver_ramified_diagonal_odd_facet_bridge_20261020.md`, SHA256
`ae39ae8dae80e7bd69c333e7e8b4a9f992a48c1396c06e182abc8e75814b23a1`,
and independently rederived all its lattice, crossing and containment
claims. No change is requested.

The primitive-pool argument is complete. The dual exponent p gives
`q|2p=p^3`. Valuation three forces every spatial coordinate into pO;
primitivity then makes the time coordinate a unit. Its contribution
has valuation one and cannot cancel the spatial terms of valuation
at least two. The surviving normalized positive norms are exactly
the stated necessary pool `1,2,q0`; every vector of such a norm is
reflective because `2/q` is integral.

The perpendicular module has the claimed basis: its defining equation
is `x1+x2=pa`, so `v0,pu,e3,...,en` generate it exactly. Adding u has
index ideal p. Direct products give

`v0^2=q0`, `u^2=v0.u=1`,

and hence `h=alpha(v0-u)` has norm `-alpha`, is orthogonal to u, and
gives an integral unit-determinant basis change of the over-lattice.
This proves the literal odd-silver facet model without asserting an
ambient integral equivalence.

I checked the complete crossing table from both real inequalities.
The trace bounds restrict B to zero or plus/minus one; the individual
physical inequalities remove every remaining option except those in
the table. In particular the q2 product `|t|=alpha^-1` really survives
and must not be dropped for having a nonnative acute angle `3pi/8`.

Since `p/q0=alpha`, the three nonzero-product projections are precisely
`pv-alpha c` and `pv-c` with the displayed norms `q0,X,q0`. All belong
to N. Each norm divides two, so the integral over-lattice is preserved
by the corresponding reflection. This proves the needed mirror-SUBSET
assertion; no converse extension assertion is being used.

The actual-facet step is valid for the full locally finite arrangement:
take any generic finite point on c-perp. The adjacent full chamber has
c as a facet. Chamber transitivity transfers hypothetical compactness
to it. Every intrinsic boundary mirror of that compact facet is among
the standard odd-silver mirrors just identified. A generic interior
point therefore selects a full standard chamber contained in the facet.
The infinite-volume input rules this out when n-1 is at least nine.
Disjoint ambient mirrors have constant sign on c-perp and create no
missing intrinsic boundary. Intrinsic Coxeter angles are unnecessary.

The exact inherited uniform input is
`odd_silver_stable_integral_isometry_and_uniform_obstruction_20261020.md`,
SHA256 `2ad82b9a52c7b81364a070b985932932eab72438922e7b0d717b392e5e4c5ee3`.
Its strong-approximation, integral identification and closed E8-host
arguments remain imported written theorems, not newly replayed here.
The optional nonreflectivity conclusion additionally uses the explicit
standard cocompact arithmetic-lattice theorem: conjugate positivity
makes the form K-anisotropic, and a finite-covolume subgroup of the
arithmetic lattice has finite index.

The proved scope is the literal lattice `M_n=<-sqrt(2)>+I_n`, n at
least ten. This does not classify arbitrary invariant lattices of the
same space, close the actual native-eight node-two star, or improve a
uniform dimension bound. No executable import, scientific run, build,
formalization or candidate-source mutation occurred.
