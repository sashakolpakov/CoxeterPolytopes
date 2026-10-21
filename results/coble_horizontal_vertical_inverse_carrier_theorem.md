# Horizontal versus vertical Coble carriers: the corrected inverse target

**Status (2026-08-09): proved carrier-support theorem and a strict correction
to the proposed `C_free` strategy; not a completed surface-free `L4`
theorem.**

The low-memory executable consequence is frozen by
[`certify_coble_horizontal_vertical_support.py`](../scripts/certify_coble_horizontal_vertical_support.py),
the
[`JSON`](../data/enumerations/coble_horizontal_vertical_support.json), and
the
[`tests`](../tests/test_coble_horizontal_vertical_support.py).

The important conclusion is not another surface realization.  It is a
restriction on what kind of `POLYTOPE -> SURFACE` theorem can possibly reach
the high-dimensional Coble examples:

```text
characteristic-uniform/free incidence branch  =>  dimension at most 11;
H14 or H17 exact Coble branch                  =>  vertical in characteristic 3.
```

Thus the earlier plan to use one characteristic-uniform
free/evaluation-surjective class and then recover the full set
`{10,11,14,17}` cannot work.  A free horizontal class would give the
strictly sharper conclusion `{10,11}`.  The `H14` and `H17` search must use
characteristic-rigid determinantal or finite-geometry incidence, not a
generic Bertini tower.

## 1. Input already proved elsewhere

For a marked `(P,h) in D_CobVin`, let

\[
             \mathcal C_P=\coprod_\beta\mathcal C_{P,\beta}
\]

be the finite disjoint union of all incidence-stack branches produced by the
all-boundary/all-contraction compiler.  It is of finite type over

\[
                         Z=\operatorname{Spec}\mathbf Z[1/2].
\]

The uniform incidence inverse theorem proves, for every algebraically closed
field `k` of characteristic different from two,

\[
 \mathcal C_P(k)\ne\varnothing
 \quad\Longleftrightarrow\quad
 P\text{ has a marked strong-Coble full-nef carrier over }k,       \tag{1.1}
\]

and every point on the left returns literally the original oriented `P`.
This includes all walls, the integral marking, the complete face cone,
finite volume, and the hidden-root/hidden-`(-1)`-curve exclusions.

If the raw boundary frame has `b` components, then every resulting surface
has

\[
 K^2=-b,\qquad \rho=10+b,\qquad \dim P=9+b.             \tag{1.2}
\]

Exact return and the finite-index root gate make the surface automorphism
group finite.  Kondō's independent nine-class theorem is therefore
applicable after (1.1).

## 2. The characteristic-support theorem

Write `|im(C_P)|` for the set-theoretic image of the finite-type morphism
`C_P -> Z`.

> **Theorem 2.1 (horizontal/vertical support).**  For every
> `(P,h) in D_CobVin`:
>
> 1. if `b` is not in `{1,2,5,8}`, then `C_P` is empty;
> 2. if `b=5` or `b=8`, then
>    `|im(C_P)|` is contained in the closed point `(3)` of `Z`;
> 3. if an irreducible component of `C_P` dominates `Z`, then
>    `b in {1,2}` and `dim P in {10,11}`.

**Proof.**  A geometric point gives, by (1.1), an exact strong-Coble
full-nef carrier with finite automorphism group.  Kondō's rows have support

| characteristic | possible boundary counts |
|---|---|
| `3` | `1,2,5,8` |
| `5` | `1,2` |
| every other characteristic different from `2` | `1,2` |

This proves (1) and (2).  A component dominating `Z` has nonempty geometric
generic fiber, hence a point over `Qbar`.  The characteristic-zero row of
the same classification forces `b in {1,2}`.  Equation (1.2) gives
`dim P in {10,11}`.  `square`

There is a slightly stronger constructibility observation.  By Chevalley's
theorem the image of each branch is constructible.  A constructible subset
of the one-dimensional irreducible scheme `Z` which contains its generic
point contains a nonempty open.  A branch is therefore either horizontal,
with carriers in characteristic zero and all but finitely many odd
characteristics, or vertical, supported at finitely many odd primes.  For
`b=5,8`, Theorem 2.1 reduces that finite support to the single prime `3`.

This is an AG-native restriction on the compiled carrier stack.  It is not
a Coxeter-volume calculation and is not encoded in the definition of
`D_CobVin`.

## 3. Why this corrects, rather than revives, the old `C_free`

The earlier proposed free predicate required base-point-freeness, `H1`
vanishing, and evaluation-surjectivity for every successive curve system.
That predicate is false on the two high-dimensional positive controls:

```text
H14/60: the six relevant 3x3 line matrices have rank 2;
H17/48: all twenty-four relevant 4x4 graph matrices have rank 3.
```

The associated systems are composed with a pencil at the failed step.  The
old induction therefore cannot be used, and this report makes no such use
of it.

There is also a conceptual obstruction independent of those rank
calculations.  A genuinely characteristic-uniform sequence of free center
choices, projective-bundle curve choices, and dense Bertini opens produces
a parameter component dominating `Z`.  Theorem 2.1 then forces
`dim P<=11`.  Consequently no corrected horizontal free lemma can contain
either the `H14` or the `H17` Coble branch.  Their rank defects are not an
accidental implementation nuisance; they reflect necessary vertical
characteristic-three geometry.

The right split is therefore

```text
C_hor:  controlled horizontal determinantal/supersolvable grammar;
         desired AG pushback is dim(P) in {10,11};

C_vert,3: characteristic-three finite-geometry/determinantal grammar;
          possible high endpoints are H14 and H17.
```

The existing polar/chord, grid/projectivity, and Möbius cores belong on the
second line.  Their incidence equations force characteristic three.  They
must not be advertised as free extensions.

## 4. A uniform-open realization lemma, at its exact proved strength

The following lemma is elementary but useful because it isolates the exact
implication needed from a future P-only grammar.

> **Lemma 4.1 (uniform-open tower).**  Let a P-derived contraction branch be
> defined over `A=Z[1/2N]`.  Suppose its construction is a finite tower
> `T_r -> ... -> T_0 -> Spec A` such that:
>
> 1. `T_0` is geometrically integral and dominant over `Spec A`;
> 2. a free-center step is a nonempty open in a smooth family with
>    geometrically integral generic fiber;
> 3. a moving-curve step is a nonempty open in a projective bundle formed
>    from a locally free, base-change-compatible section space, on which the
>    required smoothness, geometric integrality, incidence, and avoidance
>    conditions have fiberwise-proper complement;
> 4. a forced intersection step is a section or finite flat multisection
>    whose prescribed transverse locus is nonempty over the generic point;
> 5. every step is taken on the same dominating irreducible component; and
> 6. at the end the P-derived boundaries are reduced, smooth, disjoint and
>    satisfy `sum B_i=-2K`, while every P-derived root/crossing wall is the
>    constructed irreducible curve of its declared class.
>
> Then the branch has a component dominating `Spec A`.  Its generic
> geometric fiber is a strong Coble root-seed carrier.  If the raw input is
> in `D_CobVin`, the carrier returns the same `P`; hence `b in {1,2}` and
> `dim P in {10,11}`.

**Proof.**  Inductively, a nonempty open in a geometrically integral generic
fiber is geometrically integral and nonempty.  Projective bundles preserve
this property; removing finitely many fiberwise-proper bad loci leaves a
nonempty generic open.  The same holds for the stipulated transverse
section/multisection.  Thus the final generic fiber is nonempty.

The boundary condition alone supplies the two Coble Hilbert assertions.  If
`D in |-2K|`, then `D.B_i=-4`, so every `B_i` is a component of `D`; their
sum already has class `-2K`, hence `D=sum B_i`.  If `A in |-K|`, then
`2A=sum B_i`, impossible because the latter is reduced.  Therefore

\[
                    |-K|=\varnothing,\qquad
                    |-2K|=\{\textstyle\sum B_i\}.
\]

The constructed surface is strong Coble.  The uniform `D_CobVin` theorem
then gives exact same-`P` return.  Theorem 2.1 gives the final bound.
`square`

This lemma does **not** say that acyclicity or expected dimension alone
implies its hypotheses.  Mnëv-type incidence behavior and the failed
composed-pencil examples forbid that shortcut.

## 5. The smallest acceptable P-only target class

The notation `C_free` should now be retired unless it explicitly means a
horizontal class.  A better target is `C_hdet`, defined by a fixed
**controlled determinantal construction grammar** applied to raw P-derived
contraction data.

Membership must be forced by the following data and no others.

1. The input lies in `D_CobVin`; the boundary frames and every contraction
   branch are recovered from the integral wall lattice.
2. A deterministic canonicalizer produces the point/curve/proximity
   dependency program.  It may use a supersolvable join/meet module, a
   moving projective-bundle module, or a rigid determinantal module, but the
   module list is fixed before seeing the classified output.
3. All symbolic evaluation, Fitting, Jacobian, incidence, and avoidance
   matrices are generated from the P-derived divisor/proximity data over
   `Z[1/2][t_1,...,t_s]`.  A generic rank is accepted by an exact nonzero
   minor or leading-monomial identity, not by supplied surface coordinates.
4. A moving module must satisfy the full relative closed-immersion and
   geometric-integrality gates.  A rank-one/rigid module instead needs a
   canonical kernel plus exact symbolic Jacobian and irreducibility
   identities.  Base-point-free plus generically separable is forbidden.
5. All forbidden coincidences are represented by finitely many explicit
   nonzero determinant polynomials on the same integral parameter space.
6. The final coefficientwise divisor identities are exactly those of the raw
   wall branch, including reduced `B=-2K` and every root/crossing class.

The key quantifier is that the compiler itself derives and checks this
program from arbitrary raw `P`.  Membership cannot be “there exists a
coordinate/Hilbert witness.”  That witness-defined class is the already
proved `C_FG_HilbExt`, whose strict inverse score remains only `L3`.

If a uniform realization theorem is proved for this grammar, Lemma 4.1 and
Theorem 2.1 immediately give the desired nonzero pushback:

\[
 \boxed{
 P\in C_{\rm hdet}
 \Longrightarrow
 P=\mathbf P(\operatorname{Nef}S)
 \Longrightarrow
 b\in\{1,2\},\quad\dim P\in\{10,11\}.}               \tag{5.1}
\]

Thus there are no `C_hdet` polytopes in dimensions `12` or higher.  This is
strictly stronger than the unrestricted Coble-carrier ceiling `17` and is
the correct prospective AG leverage delta for a horizontal class.

## 6. Exact same-P theorem chain and open obligations

The complete chain is now short enough to audit line by line.

| arrow | status | exact obligation |
|---|---|---|
| raw `P ->` finite all-contraction branches | **proved** | uniform `D_CobVin` compiler |
| raw branch `-> C_hdet` dependency program | **open** | arbitrary-input canonicalizer for the fixed grammar |
| `C_hdet ->` dominating incidence component | **open** | prove every grammar module satisfies Lemma 4.1; no existential coordinate witness |
| generic point `->` strong Coble surface | **proved once prior arrow holds** | reduced `B=-2K` boundary lemma |
| surface `->` exact same marked `P` | **proved** | Vinberg/root/cusp hidden-wall theorem |
| exact surface `-> dim P<=11` | **proved** | Theorem 2.1 plus Kondō |

There are three mandatory calibration tests before `C_hdet` can be promoted.

1. Run the raw canonicalizer on the one-boundary `H10/22` and `H10/30`
   inputs and the two-boundary `H11/26` input.  At least two distinct raw
   types must pass without importing their surface labels, equations, or
   stored Grams.  The explicit rational `H10/30` plane construction is a
   positive control, not proof that the generic grammar accepts it.
2. Verify that `H14/60` and `H17/48` fail the horizontal grammar for the
   structural rank/composed-pencil reason and are routed to `C_vert,3`.
3. Census a broader P-only denominator before Kondō is invoked.  Record the
   raw survivors by `b`, then apply (5.1).  Without this before/after table,
   the theorem has a genuine bound but no counted search delta.

For the vertical lane, use the same exact-return theorem but replace Lemma
4.1 by a characteristic-specific determinantal realization theorem.  A
branch with `b=5` or `8` must have elimination support at `(3)`.  The
existing five finite-geometry cores prove isolated positive controls; a
variable nonthin grammar and independent terminal pushback remain open.

## 7. Honest five-gate score

| gate | verdict |
|---|---|
| independent raw `D_CobVin` denominator | **pass** |
| finite exhaustive incidence stack | **pass** |
| exact return for every geometric point | **pass** |
| horizontal/vertical AG support theorem | **pass** |
| raw-P-forced nonempty `C_hdet` branch for every member | **open** |
| theorem on a non-singleton P-only `C_hdet` class | **open** |

Accordingly, Theorem 2.1 is a proved and useful `P -> carrier support -> AG`
restriction, and (5.1) is the precise payoff of the next realization lemma.
It is not yet the requested flagship because `C_hdet` has not been frozen as
a non-singleton arbitrary-input class with forced branch existence.

## 8. Replay

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_horizontal_vertical_support.py \
  --check-output data/enumerations/coble_horizontal_vertical_support.json
PYTHONPATH=src:. pytest -q tests/test_coble_horizontal_vertical_support.py
```

These checks are constant-size.  They do not invoke `cdd`, enumerate rays,
or load the large H14/H17 archives.
