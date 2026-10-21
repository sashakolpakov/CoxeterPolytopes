# Odd silver diagonal family: proof and reproduction

20 October 2026. For `O=Z[√2]`, `α=1+√2`, the full integral reflection
group of `L_n=<-α>+I_n` is nonreflective for every `n≥9`. This excludes
finite-volume Coxeter groups whose wall reflections preserve these
literal lattices. It is not a uniform dimension bound over the field,
an exclusion of the entire rational commensurability class, or a new
compact example. Stop further searches in this diagonal family.

The [two-edge and full-face proof](../results/compact_h9_odd_silver_diagonal_two_edges_20261020.md)
establishes the complete primitive norm pool, actual B9 germ, two exact
first edges, and the full coordinate-facet identity `L_n→L_(n−1)`.
The [stable-isometry proof](../results/odd_silver_stable_integral_isometry_and_uniform_obstruction_20261020.md)
gives `L5≅<-α>+F4+<1>` by explicit local bases and spinor/approximation
theory. Applying it twice identifies L9 with the previously excluded
`<-α>+E8(O)+<1>`. The facet identity transfers this obstruction to all
higher n. The root's complete written reviews of
[the geometric/edge inputs](../results/review_compact_h9_odd_silver_two_edges_root_20261020.md)
and [the stable isometry and execution](../results/review_odd_silver_stable_uniform_obstruction_root_20261020.md)
pass. This is internal mathematical review, not external refereeing.

## Small exact replay

From the repository root:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- \
  scripts/check_odd_silver_stable_local_bases_20261020.py --verify
```

The checker owns no lock; `--parent-lock` is required. The wrapper enforces
the shared resource gates, one process/thread, 59 seconds and 1GiB RSS;
the child additionally has a 45-second alarm. No old producer is imported
or executed. The historical generation used the same command with
`--write`; its exclusive creation refuses to overwrite the retained
certificate. Keep that certificate unchanged.

The actual [generation](../results/odd_silver_stable_local_bases_write_bounded_20261020.json)
and [separate fresh-process replay](../results/odd_silver_stable_local_bases_replay_bounded_20261020.json)
both PASS, in 0.061 and 0.059 seconds. Peak process-tree RSS was
55,083,008 bytes; start gates passed at 177.25GiB disk and at least
8.63GiB available RAM. No build intermediates were created.
The [saved certificate](../results/odd_silver_stable_local_bases_20261020.json)
binds the four written source inputs and the checker before/after the run.

It recomputes the two five-column local bases over `Q(√2)` using
Fractions, their exact Gram matrices, determinants, two-sided inverses
and denominator support. It also checks the three local unit-spinor
witness norms, with the norm −1 interpreted as a polynomial identity
under `r²=−7` and the congruence `−7≡1 (mod8)`. This replay is a second
execution of the same implementation, not an independent implementation.

## Written and inherited inputs

The 2-adic square criterion, local spinor theory, Hasse–Minkowski and
the rank-seven representation argument are written inputs. Strong
approximation is the explicit Kneser–Platonov theorem as stated in
[Cao–Huang, Theorem 1.1](https://arxiv.org/pdf/1906.06967v4); the proof
checks rank-six Spin is simply connected/absolutely almost simple and
has a noncompact physical real factor. The checker does not establish
that global theorem or compute the resulting global integral matrix.
Its existence, and the future-sheet correction by −I, are written
conclusions. There is no end-to-end Lean verification and no new Lean job.

The earlier [E8 full-chamber proof](../results/fresh_silver_e8_h9_continuation_20261016.md),
[scope audit](../results/fresh_silver_e8_chamber_scope_20261016.md), and
[explicit F4-squared/E8 identification](../results/fresh_silver_e8_f4_identification_20261016.md)
are inherited reviewed inputs. Their old checkers, enumerations,
certificates and integral matrices were **not rerun** in this batch.
The new packet hashes the critical retained evidence without executing it.
Likewise the new edge geometry and lower-face completeness are written
proofs, not consequences of a local matrix test alone.

The [immutable packet](../results/odd_silver_diagonal_family_20261020_freeze.json)
records exact file sizes and streamed SHA256 hashes. Central status
documents remain outside that packet. No further L_n search is authorized
by this reproduction note.
