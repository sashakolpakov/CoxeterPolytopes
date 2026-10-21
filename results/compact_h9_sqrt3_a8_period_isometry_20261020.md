# An integral infinite-order A8 period: conditional global stopping test

20 October 2026. Written proof only. No scientific execution, import,
root census, Lean, or build. The scalar facet lattice is the exact lattice
M=E8(O)⊥Ov, O=Z[s], s=√3, v²=1−s. Put ε=2+s, η=1+s and m=4+s.
Number E8 by the path α₁,…,α₇ with β attached at α₃.

**Conditional result.** If the new norm-2 mirror

    r=(23+13s)v−(4+3s)ω₂−ω₁−ωβ                    (1)

is the FIRST endpoint on the outgoing α₂ edge of the actual E8 corner,
then M is nonreflective. Consequently the original full integral host
⟨−√3⟩⊥⟨1⟩⁹ is nonreflective. The first-endpoint premise in this sentence
is a necessary remaining census obligation; it is not asserted here.

The parent has recorded actual generation and replay of the earlier A7
edge certificate `compact_h9_sqrt3_a7_opposite_edge_20261020.json`.
That certificate establishes that the original A8 corner and this E8
corner belong to the SAME full M chamber. This note uses that geometric
interface; it does not re-execute the certificate.

## 1. Two equal-norm A8 corner vectors

The original A8 simple roots are nᵢ=αᵢ for i≤7 and

    n₈=ηεv−ω₇−ηωβ.

Their future orthogonal vector is

    V₀=(11+8s)v−2εωβ,       V₀²=9(1−s).

In the original nine-coordinate M model, V₀=(9;m⁹). The endpoint of
the edge cut by (1) is

    V₁=(73+42s)v−(16+10s)ω₂,       V₁²=9(1−s).

The ordered new A8 simple-root chain is

    (f₁,…,f₈)=(α₁,r,β,α₃,α₄,α₅,α₆,α₇).           (2)

All eight roots have norm 2. The vector V₁ is orthogonal to every fᵢ
and is future timelike. These identities follow directly from the
E8 inverse Cartan matrix. In particular r²=2 and its white bridge is
the A8 extension, not merely an abstract positive subdiagram.

There is therefore a unique K-isometry A satisfying

    A nᵢ=f₉₋ᵢ (i=1,…,8),       A V₀=V₁.             (3)

The reversal in (3) is important. The corresponding unreversed chain
map need not preserve the integral glue.

## 2. The integral glue check

In the original constrained coordinates x=(t;z₁,…,z₉), Σz=mt, one has

    x=(t/9)V₀+Σ_{i=1}^8[-Σ_{j≤i}zⱼ+imt/9]nᵢ.

Consequently (3) preserves M if and only if

    (V₁+mΣ_i i f₉₋ᵢ)/9 ∈ M.                       (4)

Set S=Σ_i i fᵢ. Since Σ_i i f₉₋ᵢ=9Σ_i fᵢ−S, it suffices to check
(V₁−mS)/9. In the integral E8,v basis this vector has v coefficient

    −21−12s

and E8 coefficients, ordered α₁,…,α₇,β,

    (34+19s, 68+38s, 96+54s, 76+43s,
     56+32s, 36+21s, 16+10s, 48+27s).

Thus (4) is integral. This proves A(M)⊂M. Since A is an isometry,
its determinant has square 1; its integral matrix therefore has unit
determinant and integral inverse. Hence A is an actual lattice
automorphism, not merely a K-isometry. It preserves the future cone
because it sends the future vector V₀ to V₁.

## 3. Its order is infinite

The first seven images in the E8,v basis are

    Aα₁=α₇, Aα₂=α₆, Aα₃=α₅, Aα₄=α₄,
    Aα₅=α₃, Aα₆=β, Aα₇=r.

The sum of their diagonal coefficients is −20−12s. For an independently
checkable calculation of the remaining two diagonal coefficients set

    Kβ=5α₇+10α₆+15α₅+12α₄+9α₃+6β+3r,
    K₇=2α₇+4α₆+6α₅+5α₄+4α₃+3β+2r,
    U=α₁+K₇+ηKβ,       W=V₁+2εKβ.

The equations A n₈=α₁ and A V₀=V₁ give

    Av=[(11+8s)W−16εU]/9,
    Aβ=[(5+3s)W−(11+8s)U]/9.

Here (11+8s)²−16ε(5+3s)=9. The v coefficient of Av is 69+40s;
the β coefficient of Aβ is −44−25s. Therefore

    tr A=(−20−12s)+(69+40s)+(−44−25s)=5+3s>9.

A finite-order real matrix of rank nine has all eigenvalues on the unit
circle and trace at most nine. Thus A has infinite order.

## 4. Why FIRST endpoints would make this a global obstruction

Assume the first-endpoint premise. Then V₀ and V₁ are actual A8 vertices
of the same full integral reflection chamber C, with their outward
simple-root germs given by n₁,…,n₈ and f₁,…,f₈. Since A preserves M,
it permutes the FULL integral reflection arrangement. Equation (3) sends
the inward chamber germ at V₀ to the inward chamber germ at V₁.
Uniqueness of the chamber containing a given germ therefore gives A(C)=C.

A compact full-dimensional polytope has finite isometry group: its
finite vertex set is permuted, and an isometry fixing all its spanning
vertices is the identity. The infinite-order A therefore rules out a
compact C. Since the quadratic lattice is anisotropic over its number
field, a cofinite arithmetic reflection chamber would be compact. Thus M
would be nonreflective. The already proved full scalar-facet transfer
then excludes reflectivity of the original H⁹ lattice.

Without first-endpoint certification, A remains a proved integral
infinite-order isometry of M, but it might move C to a different chamber.
In that case no nonreflectivity conclusion follows. The immediate next
step is the complete eighteen-profile bounded edge census, retaining all
earlier candidates and resolving simultaneous candidates by the full
simple-root cone rather than by an arbitrary ordering.
