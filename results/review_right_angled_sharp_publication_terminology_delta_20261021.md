# SharpRA publication successor: exact delta review

21 October 2026. **PASS.** Independently read the complete `diff -u`
between the preserved and publication TeX sources and checked their
hashes and sizes. The sole difference is line 319:

`four-facet $Q$` → `four-dimensional facet $Q$`.

|Source|Bytes|SHA-256|
|---|---:|---|
|`manuscripts/right_angled_sharp_dimension_bounds.tex`|15449|`3c3ddd8540593004bcd5fb57832afc77ffc639268e66cb0ab29302b286654310`|
|`manuscripts/right_angled_sharp_dimension_bounds_publication.tex`|15461|`0d607b48616429ae8ff7cf93f8cd8ead4564d4a2515f74178848c84fc6f86bab`|

The correction removes an ambiguity identified in the [final read](review_right_angled_sharp_final_read_20261021.md).
It changes no definition, theorem, formula, hypothesis, reference or
proof dependency. That source read therefore transfers to the successor
without a new mathematical argument or a renewed worker/Lean claim.

The same `right_angled_sharp` catalogue id now points to the publication
TeX, `manuscripts/build/right_angled_sharp_dimension_bounds_publication.pdf`
and `docs/right_angled_sharp_dimensions_publication_proof_map_20261021.md`.
The last has SHA-256
`cdc764a6e94ce0ee3b24ef33b38152881d54d55d8a71017fbe7d0c97ca27157c`.
It correctly gives the current six-job selection (two exact workers,
three Lean roots with 12+10+8 printed declarations, PDF) and distinguishes
the historical combined-core executions from the successor's new build.
Its originally ambiguous core/PDF attribution was corrected before this
hash was pinned. There remain six `primary_paper` bundles: the new stem
is a version of the same sharpRA paper.

The original freeze still names the original source, PDF and proof map
at their original paths. I independently checked those three artifacts
against their frozen byte sizes/hashes: source as above, PDF 207718 bytes,
SHA-256 `3f769dc4f35404fcaf6aead34f9f6adde83fc6fdbd29a2e458f6642d2b2d1150`,
and map 9142 bytes, SHA-256
`d485f9c5bf6be253f4b585716c6f4b4ef6eaa3143193dfb5826e105078cbc273`.
The freeze's current hash is
`e106da0806edcd724d38b41946e053c5127d9cdd1df9033d59a39b5e86104ad6`.
The root separately reports streaming all 44 frozen paths and finding
44/44 byte-size/hash matches; I did not repeat that complete stream.
The catalogue explicitly labels this packet as immutable predecessor
evidence, rather than silently rebinding it to the successor. This is
an appropriate way to preserve the freeze.

The [new bounded build ledger](right_angled_sharp_publication_terminology_build_20261021.json)
binds the successor source to PDF SHA-256
`c4449e3d93eb569c2c56bcc8876adeadfbb9e9a5a993fca63cd3cb12312b97ca`;
my independent PDF hash agrees. It records five pages, zero warnings,
0.914 seconds, 83,640,320-byte peak RSS, admitted disk/RAM of
161.44/7.49 GiB and removed temporary build files. It claims typesetting,
not mathematical verification. No fresh combined-core run or new visual
audit follows. No scientific, Lean or TeX job was run for this review;
only source/metadata reads and hashes, and this report, were produced.
