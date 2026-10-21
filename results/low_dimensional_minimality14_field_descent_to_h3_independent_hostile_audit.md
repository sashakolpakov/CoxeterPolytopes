# Independent hostile audit: low-dimensional minimality-14 field descent

**Date:** 2026-08-12  
**Verdict:** **PASS after two explicit repairs.**

The audited statement is the following field reduction, and no stronger
reflection-group finiteness statement:

\[
 \mathscr K_4,\mathscr K_5\subset
 \mathscr K_3\cup\mathscr K_2\cup
 \mathscr K_{\rm tr}\cup\mathscr K_{\rm small}.
\]

Here the polygon, transition, and small-diagram field sets are finite.  Thus
compact dimension three is the only low-dimensional ground-field input not
proved *inside the audited note*.  In the new global argument it is supplied
by the separate general-fork theorem.

The first draft compressed two genuine proof obligations.  The terminal
version now:

1. defines orthogonal complements with the \(k\)-rational minus-one
   eigenspace lines of the reflections, rather than assuming that a real
   square-\((-2)\) normalization is \(k\)-rational; and
2. states compact Coxeter simplicity and constructs the small branch from
   an inclusion-minimal positive principal packet inside the complete
   finite-edge wall packet.

No remaining dimension, incidence, field-inheritance, or transition-family
gap was found.

## 1. Compact face geometry is exhaustive

A compact Coxeter chamber is simple because every finite vertex link is a
spherical Coxeter simplex.  It follows that two meeting facets share a
ridge and every face is simple.

If the first descent fails, a wall \(f\) meeting the narrow facet \(E\) has
nonzero incidence with it, and

\[
                    F=E\cap H_f
\]

is a compact codimension-two face.  If the second descent fails, some wall
\(g\) cuts \(F\) in a facet and is nonorthogonal to \(E\) or \(f\).  Since
\(\dim F=n-2\ge2\), choose a vertex of this facet and the unique edge of
the simple face \(F\) leaving it.  Its other endpoint is off \(H_g\).  In
the ambient simple chamber, the two endpoints have unique walls \(g,h\)
which do not contain the whole edge.  Hence \(g,h\) really are the
endpoint-exclusive pair; no simplicity or dimension assumption is hidden.

The common face of \(e,f,g\) is finite, so its Gram is negative definite.
It is connected because \(e.f>0\) and \(g\) meets at least one of \(e,f\).
The same endpoint-link argument makes both three-node endpoint deletions of
\(\{g,h;e,f\}\) negative definite.

## 2. Both literal descents retain the full field

For a \(k\)-rational reflection, its minus-one eigenspace is a \(k\)-line.
Choose \(k\)-vectors on the lines \(\ell_e,\ell_f\), and reserve
square-\((-2)\) representatives for normalized real Gram calculations.
Therefore

\[
 U_1=\ell_e^\perp,\qquad
 U_2=(\ell_e+\ell_f)^\perp
\]

are \(k\)-defined.  This removes the normalization error which would occur
if a square root needed to scale a root line to norm \(-2\).

In the first case every wall cutting \(E\) is orthogonal to \(\ell_e\); in
the second every wall cutting \(F\) is orthogonal to both lines.  Intrinsic
projection is consequently literal.  Norms, angles, parallelism, dotted
distances, Cartan cycles, and the \(<14\) bound are unchanged.  The original
reflections preserve respectively \(L\cap U_1\) and \(L\cap U_2\), which
are full projective arithmetic modules.  The boundary-wall normals span the
appropriate Lorentz space because the descended face is compact.

Their Gram packet therefore has a positive direction.  The independently
proved cycle-field inheritance lemma applies to this principal packet and
forces its cycle field to equal the ambient minimal field \(k\).  Thus the
descendants are exact arithmetic Coxeter chambers over the *full* field,
not projected Coxeter faces over an undeclared subfield.

## 3. The finite-edge trichotomy has no missing branch

Let \(u=g.h\) in square-\((-2)\) normalization.

- If \(u=2\), the endpoint normals are independent and their two-node Gram
  is parabolic.  Parabolic full-field inheritance gives \(k=\mathbf Q\).
- If \(u<2\), take an inclusion-minimal principal packet with a positive
  direction in the complete finite-edge packet.  It has rank \(>2\).  If
  every proper packet is elliptic, it is a Coxeter triangle or a
  rank-at-least-four Lanner packet.  If a proper packet is parabolic, it is
  the quasi-Lanner noncompact-simplex case; the exact cusp line makes its
  admissible \(k\)-form isotropic, hence \(k=\mathbf Q\).
- Otherwise \(u>2\).  The marked pair itself supplies a positive direction.

In the last branch the selected four nodes form one of the seven literal
supports in the independently exhaustive four-node archive.  They cannot
lie in schematic family \(5\): after deleting \(h\), the component
containing \(g\) contains the connected triple \(e,f,g\), whereas each
family-5 endpoint-deletion component containing the other marked endpoint
has only two vertices.  Hence the packet lies in family \(1,2,3,\) or \(4\).
Every selected wall meets \(E\), so every relevant normalized entry is
below \(14\).  The dotted pair makes the packet hyperbolic and therefore
full-field; the audited tail/fork integral-locus theorem makes its fields
finite.

This also verifies that the transition packet is not merely a support
subgraph: it carries the full cyclic-product algebra required by the field
argument.

## 4. Iteration and noncompact inputs

For a noncompact finite-volume input, the full incident cusp packet has a
one-dimensional common fixed null line defined over \(k\).  A nonzero
isotropic \(k\)-vector remains isotropic at every embedding, contradicting
definiteness at a nonselected real place.  Thus \(k=\mathbf Q\), including
nonsimple cusps; the compact simplicity argument is not being applied to
the noncompact case.

For compact \(n=4\), the two literal descents land in dimensions \(3\) or
\(2\).  For compact \(n=5\), they first land in dimensions \(4\) or \(3\),
and applying the same theorem to a dimension-four descendant lands in
dimension \(3\) or \(2\).  All descended faces are compact and retain
minimality \(14\), so this iteration is legitimate.  The remaining exits
are exactly the finite transition and small-diagram sets.

## 5. Scope and frozen bytes

This PASS certifies only the low-dimensional **field descent**.  It does not
bound facets, dotted heights away from the narrow star, native lattices,
dimension globally, chambers, or reflection groups.  It invokes neither a
global arithmetic-reflection finiteness theorem nor fixed-field/fixed-rank
reflective-lattice finiteness.

The terminal audited owner file is

```text
8157c0779705f5f265e75126d10361441bf2e2b504e522d917cbda848b75f019  results/low_dimensional_minimality14_field_descent_to_h3.md
```

It parses successfully with Pandoc.
