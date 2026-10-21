# Source preparation for the bounded omega2 finite-profile search

19 October 2026. **SOURCE PREPARED, NOT EXECUTED BY THIS LANE.** Only
source reading/writing, metadata/hash inspection and an AST syntax parse
were performed. The coordinator owns review and actual gated execution.

Source `scripts/search_h19_omega2_finite_profiles_20261019.py` has331 lines,
15,540 bytes, SHA256
`02f82e5ddfb72486593a761d8a1ca86fef8e5ae3ffaaa72e901dccf129967b3e`.
The new certificate is `results/h19_omega2_finite_profiles_20261019.json`.
No frozen source, certificate, report, catalogue or manuscript was edited.

## Independent mathematical/source review

The complete contact-profile input is
`results/compact_h9_discriminant19_edge_pencil_obstructions_20261019.md`,
SHA256 `41e78966ecb8137714347e79ddc5e434044409a85d5d98bcecd437350e8ace0a`.
Its congruence refinement is
`results/compact_h9_discriminant19_edge_profile_congruences_20261019.md`,
SHA256 `3fac696b113f386e60ae15debeb2ff3f643b8681034086cd306023c6b4f6cdff`.
I read both completely and checked the contact signs, block-inverse
identity, original norm equation and every class of the40-profile table.
The34 remaining necessary labeled profiles and exact K values agree.
This is a written PASS of those reductions, not of root existence or a
complete mirror census. Both reports are direct runtime hash pins.

For the native negative simple roots old0=-e_binary and
old(j+1)=-e_j, contacts beta,eta_j are **nonnegative magnitudes**.
Writing a new short root as r=(x,y,z) gives

    2y=beta-phi*x,
    z=C^(-1)(T*e2+eta), eta_2=0,
    R=(C^(-1))_2*eta, Q0=eta^T*C^(-1)*eta,
    A=30T+R,
    A^2-15D*x^2=K=60-15beta^2-30Q0+R^2.

The binary norm is (beta^2-D*x^2)/2. The spatial norm is
30T^2+2RT+Q0, so the displayed Pell equation is equivalent to norm2
once the coordinate definitions are imposed. The source verifies the
norm and every signed white contact again for every integral candidate.
It tests both A signs, T=(A-R)/30 in O and y=(beta-phi*x)/2 in O;
it does not treat the scalar Pell equation alone as sufficient.

The finite white star has type
A1(binary)+A2(nodes0,1)+A4(nodes3..6)+A1(node7).
The inverse Grams are entrywise nonnegative, and two contacted nodes in
one component already contribute at least2 to the Schur energy. Thus
there is at most one contact per component, of magnitude1 orphi.
The source enumerates all3*5*9*3=405 labeled component choices, removes
the zero choice and q>=2, and asserts exactly51 nonzero finite profiles.
The block-inverse identity is checked on each profile as
K=30(2-q); K is therefore strictly positive at the principal embedding.

**The modulo3 test uses the original norm equation.** Testing only
whether K is square modulo3 would be vacuous because K=R^2 there.
In F9 the reduced equation is

    beta^2+2Q0+T*R-D*x^2=1.

If R is nonzero it is soluble in T. If R is zero it requires
(beta^2+2Q0-4)/D to be a square, with D inverse2phi in F9.
The exact square set is {0,1,2,1+phi,2+2phi}. This leaves40 profiles.
Modulo(D), phi=15 in F19 and the Pell equation requires K square;
the exact F19 square set removes six further profiles. The asserted
counts are therefore405 ->51 ->40 ->34, retaining distinct node labels.
These are necessary filters only. All later integrality and wall tests
remain in the search.

## Bounded-domain completeness and exact square test

The requested range is

    x=a+b*phi in O, 0<x<=600, |x'|<=2.

Then b=(x-x')/sqrt5 lies strictly between-1 and301, so b=0..300.
Put m=(isqrt(5*b*b)-b)//2=floor(b*(sqrt5-1)/2).
Since x'=a-b*(sqrt5-1)/2 lies in[-2,2], the six integers
a=m-2,...,m+3 contain every possibility. All1,806 pairs in this
conservative box receive the exact field-sign inequalities. No floating
floor, approximate height or guessed Pell-unit recurrence is used.
The conjugate bound is also conservative for a canonical norm-two
lattice root: conjugating the completed binary square gives
2>=(−D')x'^2/2, and −D'>1, so |x'|<2 automatically.

For each retained profile and each admissible x, the source asks whether
K+15D*x^2 has a square root in Q(sqrt5). The rational-norm test and
both rational possibilities C^2=(A0 plus or minus sqrt(N))/2 exhaust
field square roots C+E*sqrt5; each rational square test uses integer
numerator/denominator square roots. The zero-rational-part case is
handled explicitly. Both signs of the resulting field square root are
checked. There is no separate T or A bound.

Every integral norm-two root then receives ALL26 exact native wall
tests, even when one wall or the ray already fails. Exact squared
normalized products distinguish right/native3,4,5,6 angles and strict
ultraparallel pairs. Parallel and positive products are rejected.
The complete field-angle argument in the pinned reports supplies why
no other finite label can occur for a canonical short candidate against
the current norm2/2b alphabet. The alphabet itself is inherited.

The old centre c is never a rejection filter. For x>0,
(p,r)=-D*x<0. Halving a positive rational epsilon constructs
p+epsilon*c strictly inside all27 wall halfspaces; exact norm and future
orientation are checked. The old-centre product/status is only recorded.

## Actual initial-edge endpoint of every survivor

The fixed old3 direction is p+tau*omega2. Its saved26-wall first-hit
parameter tau_old has a spacelike endpoint. A candidate must satisfy

    (p+tau_old*omega2,r)=-D*x+tau_old*A>0.

For every survivor the source constructs tau=D*x/A and checks
0<tau<tau_old. The exact endpoint identity is

    (p+tau*omega2)^2=-2D*K/A^2<0.

It checks the future orientation, all26 old-wall inequalities, the new
wall equality, and EXACTLY the eight expected old white walls active.
Thus no surviving scalar solution is silently called a finite endpoint
without checking the actual retained chamber. The Schur q<2 condition
gives the positive-definite eight-white-plus-new Gram. This local check
does not prove that the resulting chamber has finite volume.

## Inputs, execution and output scope

The source reconstructs the native10x10 form and inverse; it reads only
the12,259-byte native prefix and10,261-byte saved26-wall certificate.
Nineteen old roots are retained by excluding indices10,20,21. The seven
new roots are recovered from the saved `new_covectors` using G inverse,
with all norms/integrality/reflection criteria checked. The pure
standard-library field/inverse helper is the sole imported local module.
The26-wall checker, helper, native data and both profile reports are
direct pins; all seven dependencies recorded by the26-wall certificate
are also rehashed. No scientific producer or giant JSON is loaded.
The two complete compact H8 facets are inherited and are not enumerated.

The source has **no inner flock**. The reviewed wrapper must own it:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- \
  scripts/search_h19_omega2_finite_profiles_20261019.py --write
```

Separate replay replaces `--write` by `--verify`. The wrapper retains
the6GiB/40GiB start gates,3GiB/8GiB abort gates,59-second deadline,
1GiB live process-tree RSS cap and parent-death protection. The child
also has a45-second alarm. Generation uses exclusive creation; replay
recomputes and compares canonical JSON. Certificate size is capped
below1MiB. Timing/RSS have not been measured by this lane: the source
has only1,806 initial coefficient pairs and34 times the retained x count
in small exact square tests, with no solver or unbounded enumeration.

Any PASS or empty outcome is a result for the stated height600 finite
profile search only. No root completeness outside that bound, other
star, universal finite-volume obstruction, vertex preservation, compact
example, nonreflectivity theorem or dimension bound follows.
