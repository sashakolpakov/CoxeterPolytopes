# Independent hostile audit: exact cusps force the rational field

**Verdict (2026-08-11): PASS in the exact finite-sided,
finite-volume, admissible returned-field scope.** A genuine cusp supplies a
\(k\)-rational null line, while an admissible nonidentity conjugate forbids
every nonzero null vector. Thus the returned coefficient field is
\(\mathbf Q\). This is a field collapse, not merely a height estimate.

## 1. A nonsimple cusp still has one common fixed line

Move an ideal accumulation point to infinity. Above a sufficiently high
horosphere, every nonincident support is a bounded hemisphere and imposes
no condition. The end of the chamber is

\[
                         C\times(T,\infty),
\]

where \(C\subset\mathbf R^{n-1}\) is the full-dimensional Euclidean
polyhedron cut out by all incident vertical supports. Since

\[
 \operatorname{vol}\bigl(C\times(T,\infty)\bigr)
 =\operatorname{vol}_{\mathrm E}(C)\int_T^\infty y^{-n}\,dy
\]

and \(n\ge2\), finite hyperbolic volume forces \(C\) to be bounded.

Write the incident normals as \(e_i=a_i u_0+w_i\), with \(w_i\) in the
positive-definite horospherical space \(E\). Boundedness gives the zero
recession cone

\[
 \{z:\langle w_i,z\rangle\le0\text{ for every }i\}
\]

zero. Hence the \(w_i\) span \(E\), and they admit a nonzero nonnegative
dependence \(\sum_i\lambda_iw_i=0\). Evaluation at an interior point of
\(C\) gives \(\sum_i\lambda_i a_i>0\). Therefore the \(e_i\) span
\(E+\mathbf R u_0=u_0^\perp\), so the intersection of all incident
reflection fixed hyperplanes is exactly \(\mathbf R u_0\). No simplicity
or choice of primitive normals is involved.

## 2. The line descends to \(k\)

For the returned \(k\)-linear reflection operators, define

\[
 W=\bigcap_{i\in I}\ker(s_i-1)
   =\ker\!\left(V\longrightarrow\bigoplus_{i\in I}V\right).
\]

Field extension is flat, so kernel formation commutes with
\(k\hookrightarrow\mathbf R\). The distinguished base change of \(W\) is
the one-dimensional cusp line; consequently \(\dim_kW=1\). This operator
argument works for a finite projective order-module and never requires a
free order basis or a normal generator.

Choose \(0\ne u\in W\). At the distinguished place,
\(\sigma_0(h(u,u))=0\). Injectivity of \(\sigma_0\) implies
\(h(u,u)=0\) in \(k\). At every other embedding, \(u^\sigma\ne0\) and

\[
                         h^\sigma(u^\sigma,u^\sigma)=0,
\]

contradicting positive definiteness. Since \(k\) is totally real, it has no
other embedding and therefore \(k=\mathbf Q\).

## 3. Exact boundary

The result needs the full incident packet of an actual finite-volume cusp.
A partial affine diagram does not suffice. It applies to simple and
nonsimple cusps, requires neither reflection maximality nor a K3 surface,
and says nothing about compact chambers. For an exact primitive full-trace
complex-K3 carrier, it makes the noncompact field set the singleton
\(\{\mathbf Q\}\); only then does \(n+1=\rho(X)\le20\) give \(n\le19\).

The independent replay is
[audit_exact_trace_noncompact_ground_field_collapse_independent_hostile.py](../scripts/audit_exact_trace_noncompact_ground_field_collapse_independent_hostile.py),
with tests in
[test_exact_trace_noncompact_ground_field_collapse_independent_hostile.py](../tests/test_exact_trace_noncompact_ground_field_collapse_independent_hostile.py).
