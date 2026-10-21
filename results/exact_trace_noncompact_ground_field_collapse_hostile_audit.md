# Hostile audit: exact cusps collapse the returned field to Q

**Verdict (2026-08-11): PASS for exact finite-volume noncompact chambers
over an admissible totally real form.**  The result is stronger than an
\(H^2\) cusp-packing bound: the returned coefficient field itself is
\(\mathbf Q\), in every dimension.  The proof uses the whole geometric cusp
packet and does not infer a rational line from an abstract affine subdiagram.

## 1. Finite volume produces the right fixed line

Choose an ideal accumulation point of interior chamber points and send it to
infinity.  Above a sufficiently high horosphere, supports not incident to
that point are bounded hemispheres and impose no condition.  Incident
supports are vertical, so the high end is

\[
                         C\times(T,\infty)
\]

for a full-dimensional Euclidean polyhedron \(C\).  Its hyperbolic volume is

\[
                  \operatorname {vol}_{\mathrm E}(C)
                  \int_T^\infty y^{-n}\,dy.
\]

For \(n\ge2\) the integral in height is finite and positive.  If \(C\) were
unbounded, its Euclidean volume would be infinite, contradicting finite
hyperbolic volume.  Thus \(C\) is bounded.  This argument does not assume
that the cusp is simple.

Write an incident Lorentzian normal as \(e_i=a_i u_0+w_i\), with \(w_i\)
in the positive-definite horospherical direction space \(E\).  The recession
cone of the bounded section is zero.  Hence the oriented \(w_i\) span \(E\)
and admit a nonzero nonnegative dependence

\[
                         \sum_i\lambda_iw_i=0.
\]

At an interior point all incident inequalities are strict, so
\(\sum_i\lambda_i a_i>0\).  Therefore a positive combination of incident
normals is a nonzero multiple of \(u_0\).  They span
\(E+\mathbf R u_0=u_0^\perp\), and their common orthogonal is precisely
\(\mathbf R u_0\).  The executable replay includes the nonsimple square
cusp in \(\mathbf H^3\) and cube cusps through dimension ten.

## 2. Projective modules and nonprimitive normals cause no gap

Let \(M\) be finite projective over the returned order and
\(V=M\otimes k\).  Every returned facet reflection is an endomorphism of
this \(k\)-vector space.  For all walls incident to the cusp, set

\[
                         W=\bigcap_i\ker(s_i-1).
\]

Equivalently \(W\) is the kernel of the single map
\(x\mapsto((s_i-1)x)_i\).  Flat scalar extension to the distinguished real
place commutes with this kernel.  Section 1 therefore gives

\[
                 W\otimes_{k,\sigma_0}\mathbf R
                 =\mathbf R u_0,
                 \qquad \dim_kW=1.
\]

This is an operator proof.  It never chooses a generator of a rank-one
normal module, so freeness of \(M\), primitive normals, divisibility, and
normal scaling are irrelevant.

## 3. The signed-place contradiction is exact

Choose \(0\ne u\in W\).  Its distinguished image lies on the null cusp
line, giving \(\sigma_0(h(u,u))=0\).  Since a field embedding is injective,
this means \(h(u,u)=0\) in \(k\), not merely numerically at one place.

For every other embedding \(\sigma\), the vector \(u^\sigma\) remains
nonzero and

\[
                          h^\sigma(u^\sigma,u^\sigma)=0.
\]

Admissibility makes \(h^\sigma\) positive definite, a contradiction.
Thus there is no other real embedding.  Total reality forces
\([k:\mathbf Q]=1\), hence \(k=\mathbf Q\).

## 4. Consequences and exact boundary

For an exact primitive full-trace complex-K3 carrier, the noncompact field
set is therefore the singleton \(\{\mathbf Q\}\).  Only after that field
collapse does K3 geometry give

\[
                            n+1=\rho(X)\le20,
                            \qquad n\le19.
\]

This carrier-relative bound does not replace the rational weak-log
\(n\le25\) branch and is not a theorem that every raw chamber has a K3
carrier.  The theorem also does not cover compact chambers, partial affine
subdiagrams, nonadmissible forms, carrier construction, Property 4, chamber
finiteness, or essential surjectivity.  Reflection maximality and saturation
are unused.

