# Dimension bounds, verification status, and the right-angled lane

Checked2026-10-18 against the current repository and the primary sources
linked below. Dimensions mean hyperbolic dimension, not lattice rank.
All polytopes here are finite-sided and have finite volume. A mixed
polytope has at least one finite and one ideal vertex; an entirely ideal
polytope has no finite vertices. Noncompact includes both kinds.

These are the best published bounds located in this literature check,
not a certificate that no newer or unpublished theorem exists.

|Class|Highest known example dimension|Published upper bound|Current ceiling used here|Status of current ceiling|Gap to examples|
|---|---:|---:|---:|---|---:|
|Arbitrary finite-volume Coxeter|21|995, Prokhorov|58|New internally complete computer-assisted proof|37|
|Arbitrary compact Coxeter|8|29, Vinberg|29|Published; no improvement here|21|
|Arithmetic noncompact Coxeter|21|21, Esselmann; dimension20 excluded|21|Published and sharp; no improvement here|0|
|Arithmetic compact Coxeter|8|29, Vinberg|21|New internally complete computer-assisted proof|13|
|Arithmetic finite-volume, compact or noncompact|21|29, Vinberg, with sharper noncompact21|21|New compact proof plus published noncompact theorem; sharp for the union|0|
|Right-angled compact|4|4|4|Published and sharp; no improvement here|0|
|Right-angled entirely ideal|4|6, Kolpakov|6|Published; no improvement here|2|
|Right-angled mixed finite/ideal vertices|8|12, Dufour|12|Published general finite-volume theorem; no improvement here|4|

The current58 result is the
[completed simultaneous rank15 proof](fresh_a15_d15_dimension58_20261018.md).
It composes written geometric arguments with exact local certificates,
complete finite upper bounds, cusp estimates and separately checked
computations. Generation, separate fresh replays and three internal
mathematical composition reviews passed. The current compact arithmetic21
comes from the [field composition](fresh_arithmetic_dimension21_composition_20261016.md)
and its [global-scope audit](fresh_arithmetic_dimension21_hostile_20261016.md).
No unresolved mathematical step is currently identified in these final
chains. Internal reviews took place within the agent workflow, rather
than by external mathematician referees.

Neither new theorem has been formalized in Lean or another proof assistant.
No Lean sources, lakefile or lean-toolchain are present in the repository.
Exact script checking is not a formal verification of the geometric prose.
The published results are conventional refereed mathematical theorems;
this check makes no claim about formalizations elsewhere.

The [earlier163 proof](dimension163_proof_status_20261017.md) is also
internally complete but is superseded by58. Current improvements to the
dimension58 local weights are research toward57, not a proof of57.
No new finite-volume example has been established in the repository.

## Published comparison and examples

[Vinberg1984](https://www.mathnet.ru/eng/mmo444) supplies the compact29
ceiling and the arithmetic restrictions used by the repository.
[Prokhorov1986/87](https://doi.org/10.1070/IM1987v028n02ABEH000889)
supplies the unrestricted finite-volume995 ceiling, in the standard
formulation excluding dimensions at least996.
[Esselmann1996](https://doi.org/10.1006/jnth.1996.0141) concerns reflective
integral Lorentzian lattices; its hyperbolic noncompact arithmetic
consequence is dimension at most21 with20 excluded. Borcherds' example
attains21 in the noncompact class. It does not show compact sharpness.
[Belolipetsky2016, §§1,2,7](https://arxiv.org/pdf/1506.03111) explicitly
states these distinctions and discusses Bugaenko's compact8 example.
Arithmetic reflection lattices are of simplest type, by Vinberg's lemma;
this is not a dimension theorem for arbitrary lattices of that type.

[Potyagailo–Vinberg2005, §2](https://ems.press/content/serial-article-files/43032)
reproves nonexistence of compact right-angled polytopes above4, and
records the right-angled compact120-cell attaining4. Their §3 constructs
finite-volume right-angled polytopes through8. Their older general14
bound was improved to12 by
[Dufour2010](https://doi.org/10.1007/s10711-009-9454-2); the publisher's
original abstract states that bound explicitly. Dufour's full text was
not accessible in this check.

[Kolpakov2012, Corollary7, p1959](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-p.pdf)
excludes entirely ideal right-angled polytopes in every dimension at
least7. The ideal right-angled24-cell gives an example in4. The first
open dimensions for this class are therefore5 and6 in the sources
checked here. The sometimes repeated ideal upper8 is obsolete.

[Italiano–Martelli–Migliorini, §2.1](https://www.cambridge.org/core/journals/journal-of-the-institute-of-mathematics-of-jussieu/article/hyperbolic-manifolds-that-fibre-algebraically-up-to-dimension-8/27F73EB6174580E7192BD79A0268F9B5)
explicitly states that its P3 through P8 have both finite and ideal
vertices. Thus dimension8 belongs in the mixed row, not the entirely
ideal row. [Alexandrov2023](https://arxiv.org/pdf/2303.09533) states the
example frontiers4 and8 and proves a stronger Nikulin–Khovanskii-type
inequality. [Egorov2026](https://arxiv.org/pdf/2603.29925), §1, repeats
the compact4, ideal6 and general12 ceilings and improves counts of finite
and ideal vertices. Its introductory body does not assert an example
in12; a loose phrase in its abstract must not be read as such an example.

## Added research lane

Treat the three right-angled vertex classes separately. The compact
dimension question is settled at4 and is a calibration class. The
entirely ideal5/6 and mixed9–12 classes have the shortest open dimension
gaps and deserve a parallel existence/nonexistence effort.

The right-angled local models are especially rigid: a finite vertex
has a spherical simplex link, and an ideal vertex has a Euclidean cube
link. Faces inherit right angles. Start from these exact link incidences,
the Alexandrov inequality and Egorov's updated vertex lower bounds;
combine them with exact Lorentz-Gram feasibility and whole-face incidence
constraints. A vertex-count obstruction without an upper count is not
by itself a dimension exclusion. An exhaustive nonexistence computation
must prove the completeness of its combinatorial domain and permit
arbitrary facet count, unless its theorem explicitly restricts that count.

An independently delegated lane addresses the user's separate volume
question in dimension4: is the facet-double of the ideal right-angled
24-cell second in the ideal right-angled volume spectrum? Establish its
exact volume and incidence data, check for published answers, and then
seek either an intermediate-volume example or a universal exclusion.
Keep second-smallest volume, equality classification and uniqueness as
separate claims. No answer to that conjecture is asserted here.
