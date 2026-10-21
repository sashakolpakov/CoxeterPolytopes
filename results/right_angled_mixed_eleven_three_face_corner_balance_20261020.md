# Balancing ideal-edge degrees in mixed eleven-dimensional corners

20 October 2026. **New paper-only necessary count inequality, pending
independent review.** No script, census, optimizer, Lean or build ran.
This note changes no dimension ceiling. It uses the new one-cusp lemma
as an explicit input and improves an ordinary incidence lower bound;
it does not infer an average from a maximum-height corner.

Let `P` be a hypothetical finite-volume mixed right-angled
eleven-polytope. Write `f,c,a2,a3,B` as in the reviewed opposition
and polygon/three-face packets. Thus `f` counts finite vertices and
`B` counts actual canonical `P3` three-faces. Set

    E=sum_(X not P3)(nu(X)-10)>=0,
    R=sum_(X not P3)(8-F(X)),

and let `Xtot` count exceptional finite-vertex/three-face incidences,
excluding `P3` and `Q8`. The signed quantity `R` is not nonnegative.

## 1. Exact local and inherited inputs

The one-cusp ten-facet exclusion is proved in
`right_angled_one_cusp_ten_facet_exclusion_research_20261020.md`,
SHA256 `19590da4197ee04940df141a492ff89d69aea44b49bf87fe3c63699c252b56d4`.
Its complete exceptional budget gives

    Xtot<=18E+4R.                              (1)

For a three-face `Y`, write `A_j(Y)` for the number of its finite
vertices with exactly `j` incident edges ending at ideal vertices.
The reviewed small-face list has

|Type|Nonzero ideal-edge degree counts|
|---|---|
|P3|A3=2|
|Q8|A1=8|
|R7|A0=1, A2=3|

These are actual edge endpoint counts, not graph degrees in an
arbitrary abstract model. In particular each of these types satisfies
`A2-3A0=0`.

We also use the reviewed local theorem from
`right_angled_mixed_nu11_all_ideal_corner_exclusion_research_20261019.md`
and its coefficient-two successor
`right_angled_mixed_nu12_carrier_error_sharpening_research_20261019.md`:
at an all-ideal-edge finite corner, a face of `nu<=11` is `P3` or
its actual reflected double `D3`; the latter has exactly two finite
vertices and `nu=11`. For every other face of `nu>=12`,

    A3(Y)<=2*(nu(Y)-11).                       (2)

These geometric inputs and the one-cusp lemma remain written premises
here; no saved-certificate hash is a substitute for them.

## 2. Two global errors with controlled size

Define sums over actual three-faces

    D=sum_Y [A2(Y)-3A0(Y)],
    H=sum_(Y not P3) A3(Y).

Then

    D>=-30E,        H<=2E.                     (3)

For the first inequality, the small three-face types contribute zero.
For `nu=11`, the new one-cusp lemma and `F+2s>=12` give `F<=9`,
hence `p=4F-2nu-4<=10=10(nu-10)`. At `nu>=12`,
`p<=2nu-4<=10(nu-10)`. Thus every remaining face contributes
at least `-3p>=-30(nu-10)` to `D`. The compact dodecahedron
attains this bound, so it is not appropriate to improve it without
a further geometric or type-dependent argument.

For `H`, each `D3` contributes two corners against its one unit of
`E`. Other faces with `nu<=11` contribute none. Equation (2) bounds
each `nu>=12` contribution by `2(nu-11)<=2(nu-10)`. Summing proves
the second part of (3). No global restriction on finite--finite edges
or on the three-face spectrum is used.

## 3. A balanced corner row

At an ambient finite vertex `q`, let `k` of its eleven edges end
ideally, put `m=11-k`, and let `t=t(G_q)` be its actual `P3`
opposition degree. Let `X(q)` count its exceptional three-faces.
Every triple with zero or two ideal-ending directions is exceptional,
so

    X(q)>=C(m,3)+C(k,2)*m.                     (4)

The local contributions to `D,H` are exactly

    D(q)=C(k,2)*m-3*C(m,3),
    H(q)=C(k,3)-t.                             (5)

The second equality uses the actual `P3` equivalence: its two finite
corners have all three edges ideal-ending, and all graph triangles
give exactly those faces.

For every integer `0<=k<=11`,

    40*C(m,3)+28*C(k,2)*m
          -16*(165-C(k,3))>=0.                 (6)

The complete successive integer values are

    3960,2160,972,288,0,0,180,432,648,720,540,0.

Combining (4)--(6), with no assumption that every ideal triple is
canonical, gives the pointwise inequality

    31X(q)>=16*(165-t)+3D(q)-16H(q).

After summing over ordinary finite vertices and using (3),

    31Xtot>=16*(165f-2B)-122E.                 (7)

The balance term is important: the one center and three leaves of
every `R7` are counted together globally. They were treated only
through a worst local ratio in the earlier `11/23` inequality.

## 4. Elimination with the same actual counts

Equations (1) and (7) give

    680E+124R>=2640f-32B.                      (8)

Retain the exact previously reviewed coface and Euler identities

    E=B-10a3+9a2+11520c,
    R=8a3-9a2-2B,
    a3=(9/2)a2-(165/4)f-5760c.

Consequently

    R=-(5/4)B-(165/8)f+5760c-(3/4)E.          (9)

The strict weighted row `a2>(25/2)f+1930c` also gives

    360c<B-(75/2)f-E.                         (10)

Substitute (9) into (8):

    587E-123B-(10395/2)f+714240c>=0.

Since `714240=1984*360`, equation (10) makes its left side strictly
smaller than

    1861B-(159195/2)f-1397E.

We obtain the stronger necessary ordinary lower bound

    3722B>159195f+2794E,
    2B>(159195/1861)f+(2794/1861)E.            (11)

Both the intercept and the excess coefficient improve the original
`22E+4R` row. The exact positive differences of the respective cross
products are

    159195*1379-117825*1861=257580,
    2794*1379-2070*1861=656.

The separately strengthened unbalanced `18E+4R` estimate remains
useful too:

    2B>(117825/1379)f+(2254/1379)E.

It has a larger excess coefficient than (11); neither new inequality
should be discarded by silently comparing only the intercepts.

## 5. Remaining obstruction

The intercept `159195/1861` is about 85.543. This is a genuine
ordinary finite-vertex incidence gain, conditional on the explicit
new local lemma, but it is not a proof that dimension eleven is
absent. The unrestricted all-formal-row scalar obstruction remains,
and the actual-graph upper bound needed to close it is not supplied
here. A possible useful target is an upper estimate of the form
`2B<=84f+K*E` with a sufficiently small `K`; its universal geometry
and analytic proof remain open.
