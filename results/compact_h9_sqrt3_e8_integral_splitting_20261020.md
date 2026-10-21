# Exact E8 splitting of the scalar-facet lattice

20 October 2026. Written algebra only: no scientific execution, import,
root census, Lean, or build. This continues the frozen scalar-facet note
`compact_h9_sqrt3_scalar_facet_next_edge_20261020.md`, SHA256
`f6be0faf51bdc8188a3f8363fe2c5e3ee35e16c356ba49dd1546574b5195567d`.

Write s=√3, ε=2+s, η=1+s, O=Z[s]. The exact integral lattice is

    M={ (t;z₁,…,z₉)∈O¹⁰ : Σz=(4+s)t },
    Q=Σzᵢ²−εt².

**Proposition.** There is an integral orthogonal decomposition

    M = E8(O) ⊥ O v,       v²=1−s.

This is an explicit unimodular change of basis, not merely a similarity
over K. It does not assume that the auxiliary root below is the first
wall met along the previously considered edge.

## The nine-vector integral basis

For i=1,…,8 put nᵢ=(0;−eᵢ+eᵢ₊₁). Define

    r=(5+3s; (4+2s)³, (3+2s)⁵, 2+s),
    v=(11+8s; (8+5s)⁸, 4+3s).

Both satisfy the defining coordinate sum. For r the sum is 29+17s;
for v it is 68+43s. Their squares check directly:

    Σr_z²=196+112s,     εr_t²=194+112s,
    Σv_z²=1155+664s,    εv_t²=1154+665s.

Thus r²=2 and v²=1−s. Also r·n₃=−1, r·n₈=−η, and all its other nᵢ
products vanish. The seven roots n₁,…,n₇ together with r consequently
have the E8 Cartan Gram matrix: an A7 chain with an extra leaf at its
third node, or arms of lengths 2,4,1. They are orthogonal to v.
For the last orthogonality identity, the spatial product r_z·v_z and
εr_t v_t both equal 473+273s.

Here is an independent check that there is no remaining lattice glue.
Use the free coordinates (t,z₁,…,z₈) of M. The vectors n₁,…,n₇ are a
unimodular basis of the kernel of

    (t,z₁,…,z₈) ↦ (t,S₈),       S₈=Σ_{i=1}^8 zᵢ.

The two remaining basis columns have quotient coordinates

    r ↦ (5+3s, 27+16s),
    v ↦ (11+8s, 64+40s).

Their determinant is

    (5+3s)(64+40s)−(11+8s)(27+16s)=−1.

Thus n₁,…,n₇,r,v are an O-basis of the entire lattice. This explicitly
settles the integral identification, independently of a determinant-square
argument. The determinant is indeed 1−s, as in the earlier M basis.

## A genuine full-arrangement E8 vertex

The earlier whole-lattice argument restricts primitive M reflection norms,
modulo unit squares, to 2 and 2ε. In the new basis write any root as
u+a v, where u∈E8(O), a∈O. Because E8 is even, u²∈2O. The root norm
is also in 2O. Since (1−s) has valuation one at η and 2 has valuation
two, the equation

    u²+(1−s)a² ∈ 2O

forces a∈ηO. This is a necessary condition on every root, at every
height; it is not a bounded-search observation.

The roots through the timelike line O v are exactly the 240 ordinary E8
roots. Indeed write u=A+sB with A,B∈E8(Z). Then

    u²=A²+3B²+2s(A·B).

If u²=2, the minimum nonzero E8 square 2 forces B=0 and A²=2.
If u²=2ε=4+2s, B≠0 would already contribute 3B²≥6 to the constant
coefficient, whereas B=0 cannot give the nonzero s coefficient. Thus
there is no norm-2ε root through v. Unit-square rescaling of a normal
does not alter its mirror and is already covered by this normalization.

The physical v is future timelike (t>0 and v²=1−s<0). Local finiteness
of the full integral reflection arrangement therefore makes each E8
chamber germ at [v] a germ of a genuine full reflection chamber. In
particular one may restart Vinberg's procedure at an actual E8 vertex,
even if an earlier mirror cuts off r from the original A8 corner.
This asserts existence of a finite corner, not compactness of its chamber.

For later enumeration the conjugate norm equation also gives the exact
necessary bounds

    |a′|² ≤ q′/(1+s),       q∈{2,2ε},

in addition to a∈ηO. These constraints concern the new integral time
coefficient a; it must not be confused with the old coordinate t.

## Reference check and remaining global problem

A focused primary-source search did not retrieve a theorem deciding
reflectivity of this exact E8(O)⊥⟨1−√3⟩ lattice. This is an access/search
boundary, not an assertion that its status is unknown worldwide.

[Bugaenko's 1992 Bielefeld abstract](https://www.math.uni-bielefeld.de/chronik/booklets/1992.pdf)
records the eight-dimensional example E8⊥⟨−(1+√5)⟩. Its field is
different. The earlier checked
[Guglielmetti thesis](https://homeweb.unifr.ch/kellerha/pub/DissRGuglielmetti-final.pdf)
row for the diagonal family ⟨−3−2√3,1,…,1⟩ supplies no integral
identification with the present lattice. Neither reference provides the
needed stopping theorem here.

The exact splitting removes the former index-nine A8-gluing barrier and
provides a much simpler global root model. It supplies neither an infinite
chamber stabilizer nor a compact full chamber. A finite prefix, a finite
corner, or a finite-height failure to find a root still settles neither
reflectivity nor the desired compact H⁹ construction.
