# Root review: three finite projections and the dimension-eight ceiling

20 October 2026. Complete written mathematical review PASS, with the
published lower bound27 in dimension six explicit. No execution or
end-to-end geometric Lean formalization is claimed here.

I read the full authoritative source
`right_angled_dimension8_descending_height_proof_20261020.md`, SHA256
`9470fa1fa3a973a12c412fbc1257537655536dba4ffc4e9ac23c320615babdbe`,
and independently rederived all its algebra, dimension conventions,
actual-wall checks and the spherical inequality.

## Actual finite-point descent

For an interior future unit point q, every unit outward facet normal
m has positive t(m)=-q.m. At a maximizing normal n, t=t(n)>0, and

    q'=(q+t*n)/sqrt(1+t²)

has square-1, is orthogonal to n, and is future because its product
with q is negative. For every OTHER actual normal m, the numerator's
product with m is q.m+t(n.m)<0 by actual interior strictness and
all-wall acuteness. This proves an actual relative-interior point in
the selected facet, not only feasibility in a guessed local star.

Every intrinsic wall of that facet is an actual ridge neighbor m_i
perpendicular to n. Its restriction remains unit, and its new height
is t(m_i)/sqrt(1+t²)<=t/sqrt(1+t²). Repeating intrinsically is valid;
the maximum is over the finite actual wall set at EACH step. After
three projections all heights are strictly less than1/sqrt3. The
strictness follows from the first finite t, including ties at maxima.

Starting in dimension nine ends in an actual six-face with a finite
relative-interior point. All dimensions involved are at least six,
so ideal-polygon angle inheritance is irrelevant. Finite volume and
right angles are the standard actual-face inheritance facts.

## Six-dimensional spherical code and elementary square identity

The six-face's intrinsic Lorentz space has dimension seven. Relative
to its chosen point its normals are (t_i,sqrt(1+t_i²)u_i), with
u_i in S5 in Euclidean SIX-space. Their all-wall products and the
strict height bounds give u_i.u_j<1/4. Directions are distinct.

The polynomial is exactly

    F(x)=(x+1/2)²(x-1/4)=x³+(3/4)x²-1/16,
    F(1)=27/16.

For unit u,v in R6, put A(u)=u tensor u-I/6 and
T(u)=u tensor u tensor u-S(u)/8, where S(u) has the three symmetric
Kronecker-delta terms. Direct contraction gives

    A(u).A(v)=x²-1/6,
    u^3.S(v)=S(u).v^3=3x,
    S(u).S(v)=3(6+2)x=24x,
    T(u).T(v)=x³-3x/8.

Thus the manuscript's ordered-pair sum is the exact sum of
N²/16, (3/8)|sum u_i|², (3/4)|sum A(u_i)|² and |sum T(u_i)|².
This proves the needed positive-kernel inequality directly; no
external Delsarte or harmonic-addition theorem is being assumed.

Since every off-diagonal F value is nonpositive, the sum bounds give
N²<=27N and hence N<=27. If N=27, equality forces every off-diagonal
F value to vanish. The strict upper product bound excludes1/4,
leaving every pair equal to-1/2. Four such vectors would have sum
of squared norm4+12(-1/2)=-2, a contradiction. Therefore N<=26.

## External input, contradiction and sharpness

The exact published input is: EVERY finite-volume right-angled
six-polytope has at least27 facets. It is explicitly stated in
Alexandrov, arXiv2303.09533v1, Section4, citing Dufour2010 Lemma1
and Proposition4; Section6 gives the same27 as
1+nu5, with nu5>=26. I checked that primary mathematical text and
the analogous derivation in Egorov, arXiv2603.29925v1, Section3.2.
This is not merely the existence of a27-facet model, and no
entire-ideal or arithmetic restriction is present. The direct
Dufour PDF was not retrieved in this pass; this access limitation
is honestly recorded in the source.

The actual six-face simultaneously has at most26 and at least27
facets, excluding the nine-polytope. Actual nine-face inheritance
excludes every higher dimension. The primary Potyagailo--Vinberg
2005 article, DOI10.4171/CMH/4, supplies finite-volume right-angled
examples through dimension eight, making the resulting ceiling
sharp. That construction is a separate published sharpness input.

The result is independent of the old internal ceiling10 and the
new entire-ideal ceiling4. The proof has no optimizer, enumeration,
field, vertex-type or restricted-facet-spectrum hypothesis.
