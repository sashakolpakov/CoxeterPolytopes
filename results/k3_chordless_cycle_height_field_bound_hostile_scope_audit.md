# Hostile audit: bounded chordless-cycle field gate

**Verdict (2026-08-11): PASS, conditionally and with the height premise kept
literal.**  The
[main theorem](k3_chordless_cycle_height_field_bound.md) correctly removes
the pair-generation hypothesis.  It does not remove the bounded mixed
correspondence-degree hypothesis.

## 1. Directed cycle algebra

Individual Cartan entries need not lie in the returned order when displayed
normals are nonprimitive.  On a finite projective returned order-module,
cyclic products are nevertheless integral.  If
\(s_i=I-P_i\) with \(P_i=e_i\alpha_i\), expansion of
\(\operatorname {tr}_k(s_{i_1}\cdots s_{i_\ell})\) expresses it as the
rank \(m\) plus the alternating sum of cyclic products over all ordered
subsequences.  The product of reflections is an endomorphism of the finite
projective order-module, so its module trace lies in the order.  Induction
on \(\ell\), beginning with
\(\alpha_i(e_i)=2\), makes the unique full-length cyclic product integral.
The induction runs over all ordered words, with repetitions allowed, so
every repeated or non-induced proper subsequence has already been covered;
length two explicitly yields \(a_{ij}a_{ji}\).  No primitive normal scaling
is assumed.  If the returned lattice is not
projective over its order, integrality in the actual acting order and the
resulting multiplication operators must be separately certified.

For a chord \(ij\) splitting a directed cycle \(C\) into \(C_1,C_2\), the
two chord directions occur once each, while every perimeter direction
occurs once.  Hence

\[
 c(C_1)c(C_2)=c(C)\,a_{ij}a_{ji}.
\]

The chord pair is nonzero, so division is legitimate in the field.  Closed
walks with repeated vertices split multiplicatively.  Thus pair cycles and
chordless simple cycles generate the same field as all cyclic Cartan
products.  This is a field-generation statement, not an equality of
integral orders.

## 2. Why cycle length is bounded without a facet bound

For \(d\ge2\), at a positive-definite conjugate place, the principal Gram matrix on an
induced \(\ell\)-cycle is positive semidefinite of rank at most
\(m=n+1\).  Its adjacent off-diagonal entries are nonzero.  For a kernel
vector, the cyclic three-term recurrence is uniquely determined by its
first two coordinates.  Therefore nullity is at most two and

\[
                         \ell-2\le m.
\]

The proof does not assert that the whole facet Gram matrix has full rank,
nor that the number of facets is bounded.  When \(d=1\), the ground field is
already rational and the conjugate-place cycle bound is unnecessary.

## 3. Conjugate and sign audit

The returned trace form is explicitly
\(q=\operatorname {Tr}_{k/\mathbf Q}(\alpha h)\) with \(\alpha\) totally
negative.  Hence the distinguished Lorentzian \(h\)-summand becomes the
unique \((1,m-1)\) summand of \(q\), and all definite conjugate summands
become negative definite.  This is exactly the signed-place decomposition
used by the mixed-height inequality.

For a length-\(\ell\) Cartan cycle, root-length factors cancel around the
cycle.  Positive-definite Cauchy--Schwarz gives

\[
                         |\sigma(c(C))|\le2^\ell
\]

at every non-distinguished embedding.  Acuteness at the distinguished
place gives \((-1)^\ell\sigma_0(c(C))\ge0\).

With

\[
 \theta_C=(-1)^\ell c(C)+2^\ell,
\]

every conjugate is strictly positive.  A zero conjugate would force
\(\theta_C=0\) by injectivity of a field embedding, contradicting its
positive distinguished value.  If the signed distinguished cycle exceeds
\(2^\ell\), then \(\theta_C\) is dominant; otherwise all of its conjugates
are at most \(2^{\ell+1}\).  The raw signed cycle is never silently treated
as totally positive.

## 4. K3 constants

For an exact full-Picard **complex projective** K3 trace return,

\[
 d(n+1)=\rho(X)\le20.
\]

For \(n=2\), this yields \(D=6\), chordless length \(L\le5\), at most two
selected field generators, and primitive-combination coefficient bound
\(M_6=16\).  If

\[
 R_0=\max\{64,C/s_0\},
\]

the primitive element has every conjugate in \((0,32R_0]\).  Its
Vandermonde discriminant therefore satisfies

\[
                  |\operatorname {Disc}k|
                  \le(32R_0)^{30}.
\]

The exponent is \(6\cdot5\), not the number of generator cycles.  The
bounded monic integer coefficient boxes prove a finite field set, not only
a bounded field degree.

## 5. Exact remaining hypothesis

For each selected large affine cycle, the theorem assumes

\[
 H\in N=\operatorname {NS}(X),\qquad
 H^2\ge s_0,\qquad
 H\cdot\tau_{\theta_C}(H)\le C.
\]

The constants \(s_0,C\) are uniform across the carrier class.  The
one-positive-place identity then bounds the distinguished conjugate.
No argument in the report derives \(C\) from reflection maximality,
finite volume, or bounded \(H^2\).  The rank-four family demonstrates why:
its polarization square stays \(14\), while its mixed height is
\(32a-10\).

The theorem also explicitly requires:

- a complex projective K3 carrier (the characteristic-\(p\) Picard ceiling
  can be \(22\), so the displayed constants would change);
- a finite projective returned order-module whose facet reflections are
  integral, or a separate cyclic-integrality/multiplication certificate;
- the field generated by the returned cyclic Cartan packet;
- exact full-Picard trace, or an independently returned positive class in
  the trace summand; and
- a big-and-nef class rather than an arbitrary class of positive square if
  the number is to be interpreted as a polarized correspondence degree.

Subject to these premises, the theorem is valid for orthogonal and
nonorthogonal couplings alike.  It rules out a bounded-height
infinite-field rank-six family, but it does not rule out one whose canonical
mixed degree escapes.
