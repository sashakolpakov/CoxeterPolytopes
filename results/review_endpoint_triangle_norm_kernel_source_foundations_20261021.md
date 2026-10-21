# EndpointTriangleNormKernel: outside-author source and assumption review

21 October2026. **PASS at the complete static source/logic-map scope.**
No correctness blocker found. This is not compilation, an actual axiom
audit, or acceptance of either geometric application. No Lean/compiler,
test, worker import, scientific or TeX process was run. Sources and map
were read in full; no reviewed file was edited.

## Exact reviewed bytes and independence

| File | SHA256 |
|---|---|
| `formal/CoxeterProofs/EndpointTriangleNormKernel.lean` | `e187298414441c7d65013928dd06e4c00c9cc72fb58474962c8c6dbbe4ef8a95` |
| `docs/endpoint_triangle_norm_kernel_logic_20261021.md` | `c9beffd99c3d4b963573dad65fcfb28637b4d85de215c187a2e6b6886fef4a75` |

I did not author this kernel or map. I contributed earlier actual-point,
component and arithmetic-normalization foundations in the repository;
that shared background is not an independent new review of those premises.
For application-interface comparison I also read the two linked written
sources, without declaring a new full geometric review of them:

* short-I source `arithmetic17_a16_isolated_i_endpoint_triangle_exclusion_20261021.md`,
  SHA `917cdf51bdab0a12592aae1886abdbd94eb61cf548c6e61f1583db3a26e7b78c`;
* ordinary extension `arithmetic17_a16_every_rank_one_endpoint_triangle_exclusion_20261021.md`,
  SHA `17d9adde741465ac787c235daf3ed8869839413a0776da426c2036309cc309a5`.

The module imports only `Std`. Its ten theorem statements and their proofs
were checked; the generic ordered commutative ring is not an integer proxy
for the real geometric variables. The final Nat theorem is explicitly the
separate rational-degree case. Actual elaboration remains untested here.

## Independent derivation of every logical group

**Cycle expansion, lines26–36.** Put lambda=A*a-b*r. Expanding gives

    rho/4 = A²a² - b²r² + b² = h(A²-b²).

This displayed shorthand is only the hand derivation: the Lean statement
uses cleared polynomial equalities throughout. Multiplying by17 and using
the supplied time relation yields
`17*rho=4*h*(8*U-17)=(8*U-17)*v`.
Both mixed terms cancel, including the foreign-contact b terms. The three
explicit `congrArg` multiplications provide the nonlinear consequences
used by `grind`; no cancellation conclusion is an input.

**Integral witness, lines39–43.** The left side of the claimed identity
is `(15*(8*U-17)-17*(7*U-15))*v=U*v`, since120-119=1 and
-255+255=0. Thus the exact polynomial witness
`theta=15*rho-(7*U-15)*v` is derived. Integrality of that expression
when its inputs lie in O_K, and its identification with theta in a field,
are explicitly outside this generic-ring theorem.

**All three scalar bounds, lines45–69.** In the short-I physical range,
`0<U*v<2*U<17`; hence `0<theta<1`. At a nonphysical embedding,
`-4*U<U*v<0` and `8*U<17` imply `-1<2*theta` and theta<0.
In the widened ordinary-A1 range, `5*U*v<14*U<119`; therefore
`0<theta` and `5*theta<7`. Positivity and each multiplication by
positive U are supplied or explicitly derived. There is no division by
a possibly nonpositive expression, implicit numerical approximation,
or reuse of the short-I physical v<2 in the extended theorem.

**Products and norm contradiction, lines71–103.** The recursive product
starts at1, including n=0. Its induction retains every factor and proves
`0<=productN<=1` from factorwise bounds. Thus a positive a<1 gives
`norm=a*productN<=a<1`, contradicting the explicit input1<=norm.
There is no assumption that the product has already been bounded, or
that a field norm is small; those conclusions are derived from factors.

**Extended pair and full norm, lines106–124.** The separate pair theorem
gives `0<a*b<7/10<1` from a<7/5 and b<1/2. The last norm theorem
then calls the previous product contradiction with a*b as its small
factor; every remaining embedding still occurs in `productN`. In the
application d>=2 is needed to select this nonphysical b. The theorem
does not manufacture such an embedding or apply this pair argument to
degree one.

**Rational case, lines127–133.** Positive Nat a,b with a<=8,b<=2 have
`1<=a*b<=16`. A positive multiple of17 cannot lie in that interval.
The proof first derives product positivity/cap, then exposes the
divisibility witness and uses Presburger arithmetic; it assumes neither
Meyer nor a quadraticity statement. Relating rational algebraic integers
U,v to these positive naturals is an application step, not hidden in the
Nat statement.

The proof strategies are consistent with the statements: elementary
ordered-ring multiplication and polynomial arithmetic for the scalar
lemmas, induction for the finite product, and a direct Nat divisibility
contradiction. No `sorry`, newly declared axiom, unsafe computation,
imported geometric record or desired ambient bound appears in the source.
This static observation does not predict an actual printed-axiom list.

## Application boundaries and ordinary-A1 separation

The map accurately distinguishes the short-I physical v<2 from the
extended physical5v<14; both use the SAME nonphysical -4<v_sigma<0
and8U_sigma<17. For the short-I norm argument take the physical theta
as a and all absolute nonphysical theta values as the finite tail.
For the extension in degree>=2 take one absolute nonphysical theta as
b, and retain all remaining absolute conjugates in the tail. The map's
external norm identity/lower-bound premises are exactly what these uses
need. The geometric source's larger physical range is not being silently
substituted into `physical_I_cap`.

Map lines44–56 explicitly leave actual endpoints, raw normalizations,
cyclic integrality, the totally real field/embeddings, Schur bounds,
nonzero algebraic-integral norm interpretation and its lower bound to the
application. It also leaves the same-vector reference payments, all-wall
feasibility, complete finite-white classification and final acute-angle
contradiction in written mathematics. None of those is instantiated by
the generic norm lemmas. The ordinary-A1 extension requires its own
written review/root acceptance independently of the kernel's compilation.
This is faithful to the distinct linked sources and their exact scopes.

The order-ring witnesses and lower bounds are ordinary theorem parameters,
not named Lean axioms; nevertheless they remain substantive mathematical
inputs. In particular an eventual standard axiom printout would not prove
cyclic integrality, number-field norms, endpoint existence or the A16+A1
geometric exclusion end to end. No ambient bound is formalized here.

## Pending execution status

There are exactly TEN requested `#print axioms` commands, at lines137–146,
one for each substantive theorem reviewed above. At this review time they
have ZERO actual output credit. Map lines4–6 expressly say compilation
and the separate actual axiom audit are pending, and lines60–71 give two
serial resource-runner commands rather than a completed receipt. The
source header's reference to the logic map is consistent with that pending
status. No successful old module or dependency needs rerunning for this
static verdict. Root may separately admit one guarded module compilation
and one actual audit; any elaboration failure or actual axiom output must
be preserved/reviewed on its own evidence.

**Conclusion:** source/assumption-map PASS for all ten declarations, with
the ordinary-A1 extension kept separate and the pending actual verification
status stated honestly. No source correction requested.
