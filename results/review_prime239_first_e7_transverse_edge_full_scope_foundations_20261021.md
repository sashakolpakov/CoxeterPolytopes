# Full review: norm239 first E7 transverse edge

21 October2026. **PASS at the explicit full-arrangement and arithmetic
input boundary.** I read the entire259-line source and independently
reconstructed its new geometry, original-module equations, finite profile
classification and every hand row. No worker, scientific import,
enumeration program, numerical test, compiler, Lean or TeX process ran.
No reviewed source or frozen paper was edited.

Reviewed source:
`results/compact_h9_prime239_first_e7_golden_transverse_edge_candidate_20261021.md`,
SHA-256 `2a207fb82414203e34625d323238abc4c297541b6e34f23c4d28b16b7c3efa15`,
10677B. Both direct input bindings match:

| Input | SHA-256 | Bytes |
| --- | --- | ---: |
| `compact_h9_prime239_full_base_and_first_binary_edge_candidate_20261021.md` | `879a3feba630072a4b127e68a3729c4da3fca4f3c44630d9aae623ddd899de02` |9217|
| `golden_odd_glue_full_base_and_discriminant_reduction_20261021.md` | `b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce` |12913|

I authored the first binary-edge input and previously reviewed/shared
related odd-glue and norm79 foundations. I did **not** author this new
transverse theorem. This is an outside-author review of its complete
new argument, not an independent discovery or self-review of its starting
germ. The actual E8+A1 germ, primitive norm alphabet2/2b, integral
arrangement discreteness/local finiteness and ordinary finite Coxeter
classification remain inherited inputs at their stated scopes.

## Original mirror, finite upper and parity

The native Gram multiplies the two displayed dual columns to unit
covectors at7 and6. I checked their squares2/6 and mutual product3.
The ray q+tau*Omega7 is the actual outgoing ray after releasing-e7;
its square is-2D+2tau². Its initial segment belongs to the same chamber
by the full starting germ, without a compactness assumption.

The proposed upper
`R=phi(q-m)/2+phi³*Omega7` has original binary coordinates
`(13+22phi,-10-17phi)`. Its binary covectors are
`(-54-87phi,-67-114phi)` and its only nonzero E8 covector is
phi³ at7. These are integral; R²=2, R.m=-phi and R.(-e7)=-phi³.
The original-lattice reflection is therefore integral and R is
primitive. I rederived `q.R=-Dphi`, `Omega7.R=2phi³`,
`tau*=D/(2phi²)` and the strictly positive finite gap
`1-tau*²/D=(7-4phi)/(4phi⁴)`. No firstness is used here.

The index at2 is respected. For clarity, if an integral binary part
is aX+bY, the inherited inverse formulas give its q,m coefficients

```text
T=a*(3+5phi)+b*(3+7phi),
kappa=a*(17+25phi)+b*(21+33phi),
binary part=(Tq-kappa*m)/2.
```

Thus T and kappa are integral; the original q,m vectors agree modulo2
with unit coordinate1+phi, forcing T=kappa modulo2O. Conversely that
congruence makes the binary coordinates integral. This proves the source's
condition(2) without replacing the original module by an enlarged
orthogonal direct sum. The E8 covectors similarly force h in O.

## Exactly five first-simple profiles

The source correctly classifies the first **simple chamber wall**, not
every arbitrary arrangement mirror meeting the transverse ray. Retaining
E7+m, a finite connected completion containing E7 has rank8 and is E8,
attached by3 at node6. It cannot also connect m in rank9. Otherwise
the new node is independent of E7 and makes I2(k) with m. Native
finite contact and conjugate definiteness give kappa=0,1,phi: for
nonzero contact,1<=kappa<2 and |kappa'|<2 force the phi coefficient
to be0 or1, leaving exactly1 orphi.

For norm2b, any nonzero retained norm2 product is b*eta. Its angular
invariant has rational norm239*N(eta)²>=239, whereas finite physical
intersection and positive definite conjugate require norm<16. Thus
that norm class is orthogonal to every retained wall, and has no E8
attachment. These facts give exactly the three epsilon0 norm2 rows,
one E8-completion row and one norm2b row. The norm equation

```text
D*T²=4h²+12*epsilon*h+12*epsilon²+kappa²-4
```

and collision formula follow directly from q²=-2D, m²=2 and the
two E8 duals. Acuteness against the still-actual old-e7 wall gives
h>=0. A new collision has T>0 and, when epsilon0, h>0. No missing
native label, negative own parameter or second retained contact is
hidden in the profile list.

## Every finite interval row

I independently obtained the physical/conjugate bounds from this finite
upper, including the exact h<=phi³ bound in the golden-contact row.
The coefficient interval calculations yield the stated **7/5/3/1/0**
row counts. The rational norm columns, in source order, recompute as:

- epsilon0,kappa0:0,-1,-5,-16,-11,-9,-45 for h²-1;
- epsilon0,kappa1:-11,-59,-191,-71,-11 for4h²-3;
- epsilon0,kappaphi:-19,-79,-239 for4h²+phi-3;
- epsilon1,kappa0:the unique h=phi-1 gives phi³, of norm-1;
- norm2b:the positive integral control j has no possible coefficient row.

In the first row zero forces T0 and no positive collision. In all other
excluded norm2 rows,239 does not divide the rational norm, contrary to
being D times the square of an O-integer. Thus the exclusion is exact
and needs no omitted square-root membership test. In the golden row,
the last target is exactly Dphi², giving the sole positive T=phi and
correct parity. This is R, at tau*, with no earlier candidate or tie.

For epsilon1 I checked the distinct formula
`DU²=(h+1)(h+2)` and `tau²/D=1-1/(2h+3)²`. The upper gives
0<=h<2.15 and the conjugate gives-2<h'<-1. Hence its phi coefficient
is1 and constant coefficient-1, leaving exactly the excluded h=phi-1.
This includes h0 rather than silently presuming positive own contact.

For norm2b, reflection integrality forces h=bj and parity forces T=2U.
The equations yield `DU²=h²-b`,
`j²<=4phi⁴/(44-13phi)<1.1²`, and |j'|<.45. The two embedding
bounds force the phi coefficient0 and then contradict a positive integer
j. All bounds follow from the displayed rational phi interval. No
unbounded search has been replaced by an unexplained height cutoff.

## Firstness and entire actual endpoint

An original-integral upper with finite collision guarantees a first new
simple chamber wall by that parameter. Local finiteness and finite-corner
simplicity make its retained E7+m data one of the five profiles just
exhausted. Consequently it is R. This argument does **not** reuse the
binary-axis mixed-E8 exclusion on a transverse ray: mixed roots are
represented by the E8-completion possibility, and no claim that all
transverse intersections are binary is made.

The endpoint `Q=2phi²q+D*Omega7` is a positive multiple of the ray
at tau*. I checked `Q²=-2D(7-4phi)<0`, its future sign,
orthogonality to E7,m,R and `Q.(-e7)=-D<0`. It remains feasible
against every original chamber wall by the established first-exit
argument. The retained walls and R have positive definite Gram
E7+I2(5); ordinary finite spherical chamber recognition gives the full
actual simple facet germ. Further roots in that finite local system
are arrangement mirrors, not additional chamber facets or unresolved
first-simple ties. No formal Gram is substituted for a realized point.

No correctness blocker was found. The theorem supplies one complete
new actual edge of the literal norm239 original module. Finite-volume
completion, a chamber period and nonreflectivity remain open. It gives
no new hyperbolic polytope example, form-class statement, dimension
bound, executed census or Lean credit. The two input scopes and original
dyadic glue remain necessary for every further continuation.
