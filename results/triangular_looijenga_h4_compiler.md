# A surface-free triangular Looijenga compiler in `H4`

## Result

There is now a genuine variable, polytope-first theorem in the class

```text
C_triangle(4,8) = primitive arithmetic Coxeter wall systems in I_(1,4)
                  with at most eight (-1)/(-2) walls and full wall span,
                  for which an anticanonical triangle is not supplied.
```

The compiler takes only the raw integral wall vectors.  It reconstructs the
canonical class by adjunction, exhausts all lattice-combinatorial boundary
triangles, constructs their plane or infinitely-near blowup models, exhausts
the period-root branches, and independently checks every exceptional and
boundary inequality in `C_D++`.  It does not import a surface, a boundary, a
known Coxeter diagram, or a Vinberg finite-volume decision.

The exhaustive answer is unexpectedly small:

> **Triangular `H4/8` inverse-carrier theorem.**  A raw wall system `P` in
> `C_triangle(4,8)` is the complete nef chamber of a marked rational
> Looijenga pair with an anticanonical triangle if and only if the compiler's
> decoration/character/cone locus is nonempty.  Equivalently, after integral
> relabelling, its weighted Gram matrix is one of exactly four types
> `T6a,T6b,T7,T8` listed below.

Thus this is not merely a criterion with an undecided residual ledger: the
finite residual ledger for **carrier detection** is empty.  This is not the
stronger `U_total` assertion that every raw polytope in the structural class
has such a carrier; most raw inputs may be rejected.  It transfers no
polytope nonexistence theorem or dimension bound.

The transferred polytope statement is also positive:

> There exist four distinct integral crystallographic Coxeter polytopes in
> `H4`, with respectively `6,6,7,8` facets, obtained as complete nef
> chambers of triangular Looijenga pairs.  Exact extreme-ray enumeration
> proves that all four have finite volume and are noncompact.

`T8` is the previously frozen eight-wall example.  The two six-wall types
and the seven-wall type are new positive full-nef carriers in this project.
The completed lineage audit identifies their abstract polytopes as known:
`T6a/T6b` are Tumarkin pyramid families `07_02/07_01`, and `T7` is Im Hof's
doubly truncated orthoscheme item 13, also Ma--Zheng `P_(2,1)`.  It also
proves that all four groups are canonical mod-two congruence preimages in
`W(I_(1,4))`, with indices `10,20,30,60`.  `T8` worldwide priority remains
unresolved.  See
[`triangular_looijenga_h4_worldwide_lineage.md`](triangular_looijenga_h4_worldwide_lineage.md).

## 1. What the raw polytope determines

Let

```text
L = I_(1,4) = ZH + ZE1 + ... + ZE4,
H^2=1, Ei^2=-1,
```

and let the primitive wall vectors be `w_j`, with `w_j^2` in `{-1,-2}`.
If these vectors span `L_Q`, adjunction gives a full-rank linear system

```text
k.w_j = -2-w_j^2.
```

Consequently `k` is unique.  The first carrier gate requires it to be
integral, characteristic, and of square five.  In the surviving orbit it is

```text
k = -3H+E1+E2+E3+E4.
```

This is a real inverse step: `k` is recovered from the wall norms and is not
attached to the input as a surface decoration.

Failure of this gate is a triangular-Looijenga carrier obstruction.  It does
not say that the abstract hyperbolic polytope fails to exist.

## 2. Why every search terminates

The orthogonal complement `k^perp` is negative definite of type `A4(-1)`.
Write

```text
u=(a,b,c,d,-3a-b-c-d) in k^perp.
```

Then

```text
-u^2 = z^T Q z,

Q = [ 8 3 3 3 ]
    [ 3 2 1 1 ]
    [ 3 1 2 1 ]
    [ 3 1 1 2 ],

Q^(-1) = (1/5)[ 4 -3 -3 -3 ]
                 [-3  6  1  1 ]
                 [-3  1  6  1 ]
                 [-3  1  1  6 ].
```

For a boundary component `d` of square `s`, adjunction and the triangle
alphabet give

```text
-2 <= s <= 3,
k.d = -s-2,
u = 5d-(k.d)k,
u^2 = 5(-s^2+s-4).
```

Hence `-u^2<=50`.  Cauchy in the displayed positive form gives the exact
bounds `|u_0|<=6` and `|u_i|<=7` for `i>0`.  The compiler therefore examines
a proved finite box and finds exactly `70` possible component classes.

The same argument at norm two gives the exact `20` roots of `k^perp`.
For a numerical exceptional class

```text
e=aH+sum c_i Ei,  e^2=k.e=-1,
```

Cauchy gives `5a^2-6a-3<=0`, hence `a` is `0` or `1`; exactly ten classes
occur.  These three finite shells close the termination proof.

## 3. All anticanonical triangles, reconstructed rather than supplied

The equations

```text
d1+d2+d3=-k,
di.dj=1 for i!=j,
-2<=di^2<=3
```

give exactly `165` unordered marked triangles:

| boundary squares | triangles | roots in `D^perp` |
|---|---:|---:|
| `(-2,-2,3)` | 60 | 2 (`A1`) |
| `(-2,-1,2)` | 60 | 2 (`A1`) |
| `(-2,0,1)` | 20 | 6 (`A2`) |
| `(-1,-1,1)` | 15 | 4 (`A1+A1`) |
| `(-1,0,0)` | 10 | 6 (`A2`) |

The last three families include the fifty direct blowups of a `P2` triangle,
with blowup distributions `(3,1,0)`, `(2,2,0)`, and `(2,1,1)`.  Stopping
there would miss 120 valid decorations.

The other two families have explicit line-conic models:

- For `(-2,-1,2)`, begin with a transverse line-plus-conic cubic.  Blow up
  one of their two nodes; its exceptional curve becomes the third boundary
  component.  Make the remaining three blowups at smooth boundary points.
- For `(-2,-2,3)`, use the same node blowup, put two further points on the
  strict line, and blow up a smooth point of the exceptional boundary away
  from its two boundary nodes.  This final point is infinitely near.

In a blowdown basis `(h;e1,...,e4)`, with a row `(a;m1,...,m4)` denoting
`ah-sum m_i e_i`, the five certified templates are

```text
(-1, 0,0): (1;0001), (1;0010), (1;1100)
(-1,-1,1): (1;0000), (1;0011), (1;1100)
(-2, 0,1): (1;0000), (1;0001), (1;1110)
(-2,-1,2): (0;-1000), (1;1011), (2;1100)
(-2,-2,3): (0;-1001), (1;1110), (2;1000).
```

The lattice has five exceptional blowdown bases.  Template matching, modulo
boundary and exceptional relabelling, gives at least one construction
witness for every one of the 165 triangles.  Thus the numerical triangles
are not being mistaken for geometric surfaces.

## 4. Exact period and effectivity compiler

For each reconstructed triangle put

```text
Lambda=D^perp,
Phi={alpha in Lambda: alpha^2=-2 and alpha is a deformation root}.
```

The compiler produces a primitive integral basis of `Lambda` and verifies
that every square-minus-two vector in this rank-two lattice is a deformation
root: each appears as an irredundant simple-root facet on a branch whose
`C_D++` inequalities are independently regenerated.

It enumerates every oriented simple system `Delta` of every closed root
subsystem.  If `R0=Z Delta`, the exact period stratum is

```text
{chi in Hom(Lambda/R0,C*): chi(beta_bar)!=1 for beta in Phi\R0}.
```

All root spans in this census are primitive.  The possible nonempty strata
are therefore

```text
(C*)^2 minus one, two, or three root hypertori;
C*;
C* minus {1};
one point.
```

The lifted Weyl chamber selects the chosen oriented simple roots as the
actual irreducible internal `(-2)`-curves.

Period theory alone does not control exceptional curves.  For each branch
the compiler independently enumerates all ten numerical exceptional classes.
An exceptional class is an irreducible wall precisely when it is dominant
against the effective boundary roots and `Delta`; every nondominant class
decomposes.  Finally every one of the ten exceptional inequalities and all
three boundary inequalities is expressed as a nonnegative integral sum of
the proposed walls.  The largest minimal sum has weight six.  This proves

```text
C_P = C_D++ intersect the selected simple-root halfspaces,
```

so there are no hidden exceptional, boundary, or nodal facets.

## 5. Exhaustive branch and Gram census

Before the facet bound there are `885` decoration/chamber branches.  Exactly
`25` generic branches have ten walls and therefore lie outside
`C_triangle(4,8)`.  The remaining census is

```text
860 decoration branches
260 distinct marked integral wall systems
4 weighted Gram types.
```

| type | wall norms | marked systems | decoration branches | finite + ideal vertices |
|---|---|---:|---:|---:|
| `T6a` | `3*(-2)+3*(-1)` | 120 | 240 | `5+2` |
| `T6b` | `2*(-2)+4*(-1)` | 60 | 120 | `4+3` |
| `T7` | `2*(-2)+5*(-1)` | 60 | 300 | `5+3` |
| `T8` | `1*(-2)+7*(-1)` | 20 | 200 | `5+4` |

The exact canonical Gram matrices are frozen in the JSON certificate.  Their
surface strata show that the association is emphatically not one polytope
to one rigid surface:

- `T6a`: boundary `(-2,-2,3)` with an `A1` branch and stratum `C*`, or
  boundary `(-2,0,1)` with full `A2` and a point stratum;
- `T6b`: generic `(-2,-2,3)` with a two-dimensional deleted-root torus, or
  `(-1,0,0)` with full `A2` and a point;
- `T7`: three carrier species, with strata `C*`, `C*\{1}`, and a point;
- `T8`: four carrier species, including two two-dimensional deleted-root
  tori and the earlier `C*\{1}` worked carrier.

This is the useful direction suggested by the question: the arithmetic
polytope determines a finite union of algebraic moduli strata, and the
nonemptiness plus cone equations are necessary and sufficient.

## 6. Hyperbolic-polytope proof, not a volume assumption

For each of the four Gram types the certificate independently verifies:

1. every reflection is an integral involutive isometry;
2. every pair is one of `m=2,3,4`, parallel, or ultraparallel in the forced
   `(-1)/(-2)` alphabet;
3. a positive strict interior vector exists;
4. projection of that vector to each wall gives a positive facet witness;
5. every extreme ray is obtained by exact intersections of four independent
   walls; and
6. every extreme ray is timelike or isotropic, with at least one isotropic
   ray.

Thus no type has a hidden hyperideal direction.  The exact finite/ideal
profiles in the table prove finite volume and noncompactness.  No Vinberg
finite-volume conclusion is fed into the construction.

## 7. Logical scope

The detector is total and exhaustive for the triangular-Looijenga carrier
question **inside** `C_triangle(4,8)`.  It is not `U_total`, not a claim that
every member of the raw class has a rational-surface carrier, and not a claim
that every arithmetic Coxeter polytope belongs to this class.  In particular:

```text
compiler accepts P  =>  an actual full-nef Looijenga carrier exists;
compiler rejects P  =>  no carrier of this triangular species exists;
compiler rejects P  !=> the abstract Coxeter polytope does not exist.
```

This separation is what makes the result safe for a primal-dual program:
the four positive rows enlarge the example boundary, while every negative
row is a precise carrier obstruction that can later be combined with an
independent universality theorem, but is not prematurely promoted to a
polytope no-go theorem.

## 8. External geometric inputs

The finite lattice, character, monoid, reflection, facet, and ray statements
are replayed by the certificate.  Actual surface realization uses:

1. Gross--Hacking--Keel, [*Moduli of surfaces with an anti-canonical
   cycle*](https://arxiv.org/abs/1211.6367): toric models, the root criterion,
   nodal versus period-trivial roots, and marked moduli as the lifted period
   torus.
2. Friedman, [*On the ample cone of a rational surface with an
   anticanonical cycle*](https://arxiv.org/abs/1207.7012): period surjectivity,
   numerical-exceptional effectivity, the complete negative-curve wall list,
   and the Looijenga-root criterion.
3. Friedman, [*On the geometry of anticanonical
   pairs*](https://arxiv.org/abs/1502.02560): period surjectivity and the
   root/effectivity criterion.

## 9. Replay

From the repository root:

```bash
python3 scripts/certify_triangular_looijenga_h4_compiler.py --check-archive
python3 -m pytest -q tests/test_triangular_looijenga_h4_compiler.py
```

The frozen artifacts are:

- `scripts/certify_triangular_looijenga_h4_compiler.py`;
- `data/enumerations/triangular_looijenga_h4_compiler.json`;
- `tests/test_triangular_looijenga_h4_compiler.py`.
