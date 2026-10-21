# Root review: stable silver isometry and diagonal-family obstruction

20 October 2026. Full independent mathematical/source review PASS of
`odd_silver_stable_integral_isometry_and_uniform_obstruction_20261020.md`,
SHA256 `2ad82b9a52c7b81364a070b985932932eab72438922e7b0d717b392e5e4c5ee3`.
The proof is relative to its explicit Hasse--Minkowski and
Kneser--Platonov inputs and the inherited closed silver-host proof.
It is not end-to-end Lean verification or external refereeing.

I independently verified the two local bases, every nonzero Gram
product, F4 determinant and basis determinant argument. Their denominator
supports cover all finite primes, including the ramified dyadic prime.
The local isometries may be chosen proper by composing with an integral
norm-one reflection when necessary. The two rational quadratic spaces
are already explicitly isometric, so the genus comparison takes place
in one fixed global space.

The finite local proper-spinor image contains every unit class. The
odd-prime sum-of-two-squares argument has a nonzero derivative. At the
dyadic prime I independently checked the norm1/7 modulo8 dichotomy,
the adjustment by alpha, Hilbert90 and the remaining rational square
classes. The norms5,-1 and-alpha have the stated integral local
representatives; the -1 representative uses the explicit2-adic square
criterion for-7. Each norm is a local unit, so its reflection preserves
the local lattice, and the norm-one composition makes it proper.

The global spinor step also checks. A prescribed a positive at the
nonphysical embedding is represented by V via Hasse--Minkowski on
V+<-a> of rank seven: all finite completions are isotropic and both
real completions are indefinite. Its last coordinate cannot vanish
because V is globally anisotropic. Principal ideals and the independent
unit sign patterns therefore remove all remaining valuation/sign
obstructions. The resulting adelic maps lift to Spin.

I independently opened Cao--Huang arXiv1906.06967v4, Theorem1.1, and
checked the exact strong-approximation criterion. Rank-six Spin(V) is
simply connected and absolutely almost simple, and its physical real
group is noncompact. No K-isotropy assumption is required. This proves
the claimed integral stable isometry; multiplying by-I corrects its
future sheet if needed. Iteration and the inherited explicit F4-squared
identification transfer the full-arrangement obstruction to L9.

The separately reviewed full coordinate-facet identity then descends
any hypothetical finite-volume L_n chamber to L9 for every n>=9. Thus
the literal diagonal family is excluded. No conclusion about all
Q(sqrt2) forms or the entire rational commensurability class follows.

## Exact local-algebra execution

The standalone checker
`scripts/check_odd_silver_stable_local_bases_20261020.py`, SHA256
`8e5535561ad64053e1b19c48b324179c9fd401797e3d8b0e8ebf13d15c1492f2`,
was fully read before execution. Its Fraction-pair operations, literal
five-by-five bases, two-sided inverses, denominator support, witness
identities and before/after source pins match the proof. It imports no
project algebra and performs no enumeration.

Generation and separate saved replay both PASS under the shared-lock
59s/1GiB/one-thread wrapper. They took0.061s and0.059s respectively,
with peak process-tree RSS55,083,008 bytes. Initial gates passed with
177.25GiB disk and at least8.63GiB available RAM. Ledgers:

* `results/odd_silver_stable_local_bases_write_bounded_20261020.json`
* `results/odd_silver_stable_local_bases_replay_bounded_20261020.json`

The checker verifies local algebra only. The2-adic criterion, spinor
theory, Hasse principle, strong approximation, face geometry and old
nonreflectivity remain written inputs; old enumerations were not rerun.
No global integral matrix is claimed. No build intermediates were
created, and no new Lean or TeX work was started.
