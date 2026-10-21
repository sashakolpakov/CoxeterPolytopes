# Independent hostile audit: the compact maximal (2,5,6) K3 PHYS no-go

**Verdict (2026-08-11): PASS, with one source-presentation caution and a
classification-free repair.** The frozen packet really is a compact rank-six
arithmetic reflection packet, its full codifferent trace-\((-2)\) shell has
exactly five elements, and the prime above \(3\) kills every simultaneous
same-vector K3 normalization. Common form similarity, arbitrary facet-line
representatives, fractional module choices, coefficient-order choices, and
facet permutations do not open an escape.

The mathematical claim of absolute **reflection** maximality also passes.
This audit proves it directly by area and an exhaustive four-index triangle
census, so it does not depend on interpreting the cited classification. The
owner's two-source route is consistent with the published tables, but the
sentence attributing the general triangle-overgroup lemma specifically to
Singerman should not carry the proof by itself. A second clean reviewer
independently reached the same arithmetic verdict and flagged exactly this
source-attribution issue.

The audit did not import the owner certifier and did not edit any owner byte.
It uses the basis \((1,\phi)\), derives the codifferent from the inverse trace
Gram, and uses ideal-valuation parity rather than the owner's explicit
quadratic-field square routine.

## 1. Frozen-byte and raw-packet reconstruction

The audited owner hashes are

~~~text
9fa8afdb6f8e3f5225c494d48e3e9130a9b86cfe2aa8f58b38b177b9f972b83e  results/k3_compact_maximal_256_phys_empty.md
452f2946ec042dd7380e1772e3355d520753b54b909dd970c910857338527d40  scripts/certify_k3_compact_maximal_256_phys_empty.py
d5d7dd20cb241b9dfa21c88327c62521862cfc4007be0a5965268ced449c5b03  data/enumerations/k3_compact_maximal_256_phys_empty.json
fd6842bccc25a5dbac21cc9121d8fc31ae6d808121b653b596255bc2ec5ef11a  tests/test_k3_compact_maximal_256_phys_empty.py
~~~

Put \(\phi^2=\phi+1\). Direct determinant expansion of

\[
 G=\begin{pmatrix}-6&0&3\\0&-2&\phi\\3&\phi&-2\end{pmatrix}
\]

gives \(\det G=6\phi\). Its exact LDL pivot signs are

\[
 (-,-,+)\quad\text{at }\phi={1+\sqrt5\over2},
 \qquad
 (-,-,-)\quad\text{at }\phi'={1-\sqrt5\over2}.
\]

The Cartan entries are \(-1,-3,-\phi,-\phi\), and the nonzero
normalized products are \(\phi^2/4\) and \(3/4\). Thus the type is
\((2,5,6)\), all reflections preserve the displayed
\(\mathcal O_k\)-module, and

\[
 {1\over2}+{1\over5}+{1\over6}={13\over15}<1.
\]

The chamber is compact and the trace lattice has signature \((1,5)\) and
rank \(2\cdot3=6\). Finally

\[
 \mathcal O_0=\mathbf Z[3,\phi^2]
             =\mathbf Z[\phi]=\mathcal O_k.
\]

This closes the coefficient-overorder lane before any K3 normalization is
attempted.

## 2. Classification-free absolute reflection maximality

The triangle area is

\[
 \operatorname{area}(P)
 =\pi\left(1-{1\over2}-{1\over5}-{1\over6}\right)
 ={2\pi\over15}.
\]

Suppose its side group \(W\) lay properly in a discrete reflection group
\(W'\), with index \(N\ge2\). The mirrors of \(W'\) refine the compact
chamber \(P\); a chamber \(Q\) of \(W'\) contained in \(P\) is therefore
compact, and \(P\) is a union of \(N\) copies of \(Q\). Hence

\[
                    \operatorname{area}(Q)={2\pi\over15N}.
\]

A compact Coxeter quadrilateral has area at least \(\pi/6\): relative to
four right angles, positivity forces at least one label \(m\ge3\). A
compact Coxeter \(s\)-gon with \(s\ge5\) has area at least
\((s/2-2)\pi\ge\pi/2\). But \(Q\) has area at most \(\pi/15\), so \(Q\)
must be triangular.

The smallest compact Coxeter triangle is \((2,3,7)\), of area \(\pi/42\).
This is elementary: for ordered \(2\le a\le b\le c\), the branches
\(a\ge3\), \((a,b)=(2,3)\), \((a,b)=(2,4)\), and \(a=2,b\ge5\) have
respective minimum candidates \((3,3,4)\), \((2,3,7)\), \((2,4,5)\), and
\((2,5,5)\). Therefore \(2/(15N)\ge1/42\), whence \(N\le5\).

For \(N=2,3,4,5\), the exact equation is

\[
 {1\over a}+{1\over b}+{1\over c}=1-{2\over15N},
 \qquad 2\le a\le b\le c.
\]

There is no cutoff hidden in the executable census: the equation forces
\(a\le3\), then
\(b\le2/(1-2/(15N)-1/a)\), and \(c\) is determined. The full result is

\[
\begin{array}{c|c}
N& (a,b,c)\\ \hline
2&(2,3,10)\\
3&\varnothing\\
4&\varnothing\\
5&\varnothing.
\end{array}
\]

At index two, the two chambers inside connected \(P\) share a side, so
\(P\) is a reflection double of \((2,3,10)\). The three possible shared
sides give

\[
\begin{array}{c|c}
\text{shared endpoint orders}&\text{double}\\ \hline
(2,10)&(3,3,5)\\
(2,3)&(\pi/10,\pi/10,2\pi/3),\ \text{not Coxeter}\\
(3,10)&\text{a quadrilateral}.
\end{array}
\]

None is \((2,5,6)\). Thus no proper discrete reflection overgroup exists.
This proves exactly the raw notion needed by the compiler: absolute
maximality among hyperbolic reflection groups. It does not need the stronger
claim that \(\Delta(2,5,6)\) is maximal among every Fuchsian overgroup.

As a secondary source check, the target is absent from Singerman's complete
nonmaximal triangle-subgroup patterns, and Maclachlan's arithmetic table
puts \((2,5,6)\) as the maximal triangle group of its class. See
[Singerman](https://doi.org/10.1112/jlms/s2-6.1.29) and
[Maclachlan, Table 4](https://msp.org/pjm/1996/176-1/pjm-v176-n1-p10-p.pdf).
Those citations support the conclusion, but the direct proof above removes
the attribution ambiguity from the audit verdict.

## 3. The codifferent shell is complete

In the integral basis \((1,\phi)\), the trace Gram and its inverse are

\[
 T=\begin{pmatrix}2&1\\1&3\end{pmatrix},\qquad
 T^{-1}={1\over5}\begin{pmatrix}3&-1\\-1&2\end{pmatrix}.
\]

Thus a trace-dual basis of \(\mathcal O_k^\vee\) is

\[
 d_0={3-\phi\over5},\qquad d_1={-1+2\phi\over5}.
\]

Every codifferent element is uniquely \(A d_0+B d_1\) with
\(A,B\in\mathbf Z\), and its trace is \(A\). Trace \(-2\) therefore fixes
\(A=-2\), with no bounded search. Writing \(c=B+1\) gives

\[
 -2d_0+B d_1=-1+{c\over\sqrt5}.
\]

Its two real embeddings are \(-1\pm c/\sqrt5\), so total negativity is
equivalent to \(|c|<\sqrt5\). Since \(2^2<5<3^2\), the entire shell is

\[
 \left\{-1+{c\over\sqrt5}:c=-2,-1,0,1,2\right\}.
\]

The adjacent excluded values \(c=\pm3\) have opposite signs at the two
embeddings. This proves both completeness and the absence of a boundary
case.

## 4. Every scaling and representative freedom cancels or squares

Let \(h_i=g(r_i,r_i)\). For the two sides at angle \(\pi/6\),
\(h_1/h_3=3\). Every realization of the same pinned reflection lines and
compatible irreducible form has the shape

\[
                  g' = \lambda DGD,
 \qquad \lambda\in k^\times,
 \qquad D=\operatorname{diag}(a_1,a_2,a_3),
\]

after choosing line representatives. Hence literal facet roots have
energies \(u_i=\lambda a_i^2h_i\), and necessarily

\[
                         {u_1\over u_3}
                     =3\left({a_1\over a_3}\right)^2.       \tag{4.1}
\]

The common similarity \(\lambda\) cancels. Replacing any line generator,
including by a fractional-module generator, contributes only a square.
Raw similarities preserve the normalized \(m=6\) line-pair class \([3]\),
and the three distinct labels make that pair intrinsic under facet
permutation.

Because \(\mathcal O_0=\mathcal O_k\), there is no compatible overorder.
For any literal root \(\delta\) in an \(\mathcal O_k\)-stable integral trace
module,

\[
 \operatorname{Tr}_{k/\mathbf Q}(\alpha g'(\delta,\delta))\in\mathbf Z
 \quad\text{for every }\alpha\in\mathcal O_k,
\]

so its energy lies in \(\mathcal O_k^\vee\). Weyl framing, period choices,
and later lattice overlattices cannot alter its trace square \(-2\).
Consequently (4.1) is invariant under every freedom actually present in the
full-trace, same-vector compiler.

## 5. The inert-prime obstruction has no valuation leak

The polynomial \(x^2-x-1\) takes the values \(2,2,1\) modulo \(3\), so it is
irreducible and \(3\) is inert in \(\mathbf Q(\sqrt5)\). For the unique prime
\(\mathfrak p=(3)\), normalized by \(v_{\mathfrak p}(3)=1\), every shell
element satisfies

\[
 u_c={c-\sqrt5\over\sqrt5},\qquad
 N(c-\sqrt5)=c^2-5\in\{-5,-4,-1\}.
\]

Neither numerator nor denominator has positive
\(\mathfrak p\)-valuation; hence \(v_{\mathfrak p}(u_c)=0\) for all five
values. Equation (4.1) would give

\[
 0=v_{\mathfrak p}(u_1/u_3)
   =1+2v_{\mathfrak p}(a_1/a_3),
\]

an even-odd contradiction. The executable audit repeats this for all
\(5^2=25\) ordered energy pairs. Thus the simultaneous diagonal root
normalization is impossible and **PHYS is empty**.

## 6. Exact scope of the conclusion

The counterexample disproves

\[
 \text{raw arithmetic reflection maximality}
 \Longrightarrow
 \operatorname{K3Comp}(P)\ne\varnothing
\]

for the primitive **full-trace, full-Picard, same-vector** K3 compiler. The
failure is before EMB, ROOT, or ample framing. It is not an embedding
artifact: if an even rank-six PHYS host existed, then
\(\ell(A_N)\le6\) and its K3 complement would have rank \(16\ge6+2\), so
the stable primitive-embedding criterion would apply.

The packet does **not** exclude a larger Picard lattice whose ambient roots
project to the raw walls, and it does not exclude a different algebraic-
geometric carrier or the full Property-4 disjunction. Raw maximality can
remove hidden native root walls only after PHYS, primitive embedding,
zero-projection exclusion, and wall-selective native saturation have been
supplied independently. The owner report states these nonclaims correctly.

## 7. Reproduction and presentation audit

The separate replay is
[audit_k3_compact_maximal_256_phys_empty_independent_hostile.py](../scripts/audit_k3_compact_maximal_256_phys_empty_independent_hostile.py),
with an exact
[JSON archive](../data/enumerations/k3_compact_maximal_256_phys_empty_independent_hostile_audit.json)
and independent
[regression tests](../tests/test_k3_compact_maximal_256_phys_empty_independent_hostile_audit.py).

The owner Markdown has balanced display delimiters, no malformed control
sequences found by the audit, working relative artifact links, explicit
scope nonclaims, and a visible declaration of its published maximality
dependency. Its TeX converts successfully through Pandoc. The only
presentation caution is the source attribution described in Section 2; it
does not affect the theorem after the direct reflection-area proof.

Run:

~~~text
python3 -m scripts.audit_k3_compact_maximal_256_phys_empty_independent_hostile --check-output
pytest -q tests/test_k3_compact_maximal_256_phys_empty_independent_hostile_audit.py
~~~
