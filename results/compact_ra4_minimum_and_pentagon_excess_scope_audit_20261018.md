# Compact right-angled H4: minimum and pentagon-excess scope audit

18 October 2026. Bounded primary-literature follow-up and independent
algebraic check. No proof computation, optimizer, compilation, or Lean was
run. Earlier reports and frozen evidence were not changed.

**Verified literature boundary.** The regular compact right-angled
120-cell has volume 34pi^2/3; see Ma--Zheng, Remark2.7, p9 of
[v4](https://arxiv.org/pdf/1801.08814v4#page=9).
Kolpakov's thesis, Section4.4, pp73--75, asks its universal volume and
facet minimality. Lemma6, p74, gives

    Vol(P)=(v-4f+16)pi^2/12,

and its corollary proves f>120 implies Vol(P)>34pi^2/3. On p75 the stronger
assertion that a non-dodecahedral facet forces f>120 and v>600 is explicitly
a question. These are partial results, not a universal minimum theorem.
[Author's thesis](https://homeweb.unifr.ch/kellerha/pub/KolpakovThesis_Final.pdf).

The follow-up search did not verify a later universal minimum-volume,
minimum-facet, or second-volume classification theorem. This is an audit
limitation, not a proof that the questions remain open in October2026.
Two recent extremal papers do not fill that gap: Vesnin--Egorov's
[2025 theorem](https://doi.org/10.1134/S0037446625060047) is in H3;
Egorov's [September2026 preprint](https://arxiv.org/html/2609.04258v1#S1)
restricts its H4 theorem to common angle alpha<pi/2 and explicitly excludes
the right-angled endpoint. Its Theorem1.1 cannot be applied by taking an
endpoint limit.

**Independent pentagon-excess calculation.** Write (v,e,r,f) for the
face vector of any compact right-angled H4 polytope and put

    W=sum_(ridges R)(number_of_sides(R)-5).

The polytope is simple. Thus e=2v, r=v+f, and counting ridge corners gives
sum_R sides(R)=6v. Consequently

    W=v-5f>=0,
    Vol(P)=pi^2(f+W+16)/12.                        (1)

For a three-dimensional facet F, trivalence and Euler give
number_of_faces(F)=v(F)/2+2. Therefore its excess is exactly

    d(F)=sum_(two-faces R of F)(sides(R)-5)
        =(v(F)-20)/2=number_of_faces(F)-12,
    sum_F d(F)=2W.                                (2)

The nonnegativity is elementary because compact right-angled hyperbolic
polygons have at least five sides. The stronger local fact d(F)=0 or
d(F)>=2 follows from v(F)=20 or v(F)>=24; the first case is precisely the
regular right-angled dodecahedron. This three-dimensional gap is recorded
in Alexandrov--Bogachev--Vesnin--Egorov, Section1.2 of the
[2023 paper](https://www.mathnet.ru/php/archive.phtml?jrnid=sm&option_lang=eng&paperid=9740&wshow=paper).
In particular, the number of non-dodecahedral facets is at most W.

By (1), universal 120-cell volume minimality is equivalent to the new
universal inequality f+W>=120. The available W>=0 proves it for f>=120;
it leaves the possible f<120 cases. Pentagon excess does not itself give
the ideal argument's opposite-cusp injection or a cutting section.

**The double is a construction, not a universal second-volume theorem.**
Reflect the 120-cell across one dodecahedral facet. Its double has

    v=2*600-2*20=1160,
    f=2*120-2-12=226,
    W=30,
    Vol=68pi^2/3.

The twelve neighboring pairs of facets merge into doubles of the
dodecahedron across a pentagon. Each such three-polyhedron has30 vertices,
17 faces and d=5; the other214 facets remain dodecahedra, giving
sum d=60=2W. These counts independently agree with (1)--(2).

In the explicitly defined class of convex unions of whole chambers of
the regular right-angled 120-cell tessellation, volumes are N*34pi^2/3.
For N=2 the two chambers share a facet, and facet transitivity makes their
union isometric to this double. This is a statement about that class.
Neither universal membership in the class nor a universal second-volume
gap or equality theorem was verified. Even universal minimality of the
120-cell, if established separately, would not imply those stronger claims.
