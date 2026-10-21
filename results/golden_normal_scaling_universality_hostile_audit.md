# Hostile audit: golden normal-scaling universality

**Audit date:** 2026-08-09  
**Verdict:** **PASS, with one certificate-provenance caveat.**  The stated
theorem is correct with its written quantifiers.  It closes normal scaling
and finite saturation for each *fixed, full-rank realization* in
`C_35^conn`; it does not prove that a polytope has a surface carrier.

## Mathematical audit

1. **Unit-square propagation is valid.**  On an `m=3` or `m=5` edge,
   `a_ij a_ji` is respectively `1` or `phi^2`.  Integral factors of a unit
   are units, and

   \[
   d_j/d_i=a_{ij}/a_{ji}=a_{ij}^2/(a_{ij}a_{ji})
   \]

   is a unit square.  Connectivity propagates one common root length.
   Positive choices at the distinguished embedding remove the sign
   ambiguity.  After common form homothety, the normalized Gram entries are
   exactly `2`, `-1`, `-phi`, and `0` on the prescribed finite edges.

2. **There is no hidden cycle parameter.**  Comparing two systems with the
   same geometric products gives `c_ji=c_ij^-1`.  Symmetrizability around a
   cycle gives `(product c_ij)^2=1`; oriented signs make every `c_ij`
   positive at the hyperbolic embedding, so the product is `+1`.  The unit
   cochain is therefore a vertex coboundary.  This proves uniqueness of the
   normalized **Gram**, not uniqueness of literal coordinate vectors.

3. **The finite-overclosure argument is valid and correctly scoped.**  If
   the wall span `R` has full rank and an integral lattice `M` contains it,
   then pairing elements of `M` with `R` gives
   `R subset M subset R^vee`.  The discriminant module is finite.  The same
   argument applies to the fixed trace form.  It would fail without full
   span, integrality, or a locked common scale; all three exclusions are
   explicit in the report.

4. **Bugaenko's dotted entry is genuinely forced.**  An independent SymPy
   reconstruction from the archived Coxeter relations reproduced

   \[
   D_0=(-9+7\phi)\left(x-\frac{36+2\phi}{31}\right)(x+2\phi),
   \quad
   D_1=(-16+12\phi)\left(x-\frac{4+2\phi}{5}\right)(x+2\phi).
   \]

   Their sole common root is `x=-2phi`; the resulting Gram has rank `9`.
   Independently, the first nine wall coordinates have determinant norm
   `-1`, and the ambient Gram determinant has field norm `-4`.

5. **Bugaenko's saturation claim is valid.**  The field discriminant module
   is the simple module `O/(2)` because `2` is inert.  Its only nontrivial
   candidate is the full dual, whose line generator has square
   `-1/(2phi)`, so it is not integral.  For the trace lattice, the unimodular
   basis change gives
   `<2> + <-2> + E8(-1)^2`.  Its discriminant classes have squares
   `1/2`, `-1/2`, and `0`; multiplication by `phi` acts modulo two as
   `[[0,1],[1,1]]` and cycles all three nonzero classes.  Hence the unique
   ordinary isotropic line is not `O`-stable, and no nontrivial `O`-stable
   even trace overclosure exists.

6. **The hostile boundaries are real.**  `p=3+phi` has norm `11`, so it is
   not a square in `F`, while its two real embeddings give the required
   Lorentzian/positive-definite pair.  The family
   `p_k=3+phi^(2k)` gives infinitely many distinct actual realizations of
   the same angle-only infinity diagram.  Common form homothety is harmless
   only geometrically: with fixed trace scalar it multiplies a rank-`2r`
   trace determinant by `m^(2r)`.

## Exact residual boundary

The surface-free definition is essential.  The result proves

```text
fixed full-rank wall realization in C_35^conn
  -> one normalized wall Gram
  -> finitely many integral/stable saturations.
```

It does **not** prove any of the following:

- existence of a primitive K3 or other surface embedding;
- absence of additional effective roots cutting the proposed chamber;
- equality of the nef chamber with the Coxeter chamber;
- existence of a successful saturation for every polytope; or
- finiteness from an angle-only diagram containing infinity edges.

Thus this is useful polytope-to-lattice leverage, but not yet the requested
polytope-to-surface necessary-and-sufficient theorem.

## Certificate audit and replay

The archive replay passed with canonical certificate digest
`b1792a92b6073a9ccb1e56fc3c973b38a1e3213af7528fe15b125bfe5b87c5a9`.
The focused test suite passed `7/7`, and byte compilation passed.

One nonfatal provenance caveat remains: the main certificate hard-codes the
trace discriminant action and dual-line square rather than deriving them
from the archived Gram.  I independently derived both, and
`scripts/bugaenko_trace_lattice.py` exactly verifies the claimed standard
trace lattice.  A future hardening should import or hash that derivation in
the main certificate; this affects self-containment, not the verdict.

Audited SHA-256 digests:

```text
1d114df0683da3d6949ed053eea5fe3e22acb92c5289a6fca84c86bb259dc3e8  results/golden_normal_scaling_universality.md
ac873044b1d31e815e5095ae834955deb49f59bd7ebe89f977740b09ee11b610  scripts/certify_golden_normal_scaling_universality.py
b5c3adcefc089d0851578c22541e21cb7917439ab3cf3749b240edb0e9bbf04f  data/enumerations/golden_normal_scaling_universality.json
a579c2a3870e2330ccafb187d0071986c87a5175788acd782f467f4cf3366d5c  tests/test_golden_normal_scaling_universality.py
1d7900fd72f41a7c1ff346ba17aa8c6170d24db4b7b3504acb79274c4b4c85b3  scripts/bugaenko_trace_lattice.py
c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c  data/known_examples/bugaenko_compact_h8.json
86a7a0aab71b9fbc20c63cda75ad81a7a22f7b9745a8708ac81676f15007272c  data/enumerations/polytope_first_golden_k3_leverage.json
```
