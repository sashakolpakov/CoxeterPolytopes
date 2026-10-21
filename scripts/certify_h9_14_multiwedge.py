#!/usr/bin/env python3
"""Freeze the multiwedge polytopality upgrades for records 43, 44, and 47."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_multiwedge import audit_all


ROOT = Path(__file__).resolve().parents[1]
PROOF_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_multiwedge.py"
)
OUTPUT = ROOT / "data/enumerations/h9_14_multiwedge_polytopality.json"
SELECTION_SOURCES = {
    43: ROOT / "data/enumerations/h9_14_record43_selection.json",
    44: ROOT / "data/enumerations/h9_14_record44_selection.json",
    47: ROOT / "data/enumerations/h9_14_record47_selection.json",
}
EXPECTED_SELECTION_CERTIFICATES = {
    43: "af521ed4e7476c95945df9b770f36a343dcf3fc79dabed557aaab0fd97a671a3",
    44: "e66e493d0c68f7ae28ffe03e63f96209276a8d83a5ceb27c498f35b42c4a1867",
    47: "960ef6564c6b8013caf0989c95a46cc8401c846d32d47c66c057091340301546",
}
EXPECTED_CERTIFICATE_SHA256 = "d5efd409cf5d5a5a62e9784f8c9e71eebb5a08ee20e8ea7f1d05e894d305a183"


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _certificate_digest(manifest: Mapping[str, object]) -> str:
    critical = dict(manifest)
    critical.pop("certificate_sha256", None)
    return _canonical_sha256(critical)


def _selection_provenance() -> list[dict[str, object]]:
    result = []
    for index, path in sorted(SELECTION_SOURCES.items()):
        archive = json.loads(path.read_text(encoding="utf-8"))
        embedded = archive.get("certificate_sha256")
        if embedded != EXPECTED_SELECTION_CERTIFICATES[index]:
            raise ValueError(f"record-{index} selection archive is not frozen")
        result.append(
            {
                "record_index": index,
                "path": str(path.relative_to(ROOT)),
                "file_sha256": _file_sha256(path),
                "embedded_certificate_sha256": embedded,
            }
        )
    return result


def build() -> dict[str, object]:
    audits = [audit.to_record() for audit in audit_all()]
    if [audit["record_index"] for audit in audits] != [43, 44, 47] or not all(
        audit["ordinary_polytopal"] for audit in audits
    ):
        raise AssertionError("multiwedge audits did not certify all three records")
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_multiwedge_ordinary_polytopality",
        "claim": (
            "H9/14 records 43, 44, and 47 are nerves of ordinary convex "
            "simple 9-polytopes with 14 facets."
        ),
        "construction": (
            "For each record, the nerve is the join of the boundary of a "
            "4-simplex with a simplicial multiwedge of a small polytopal "
            "sphere. Dually, take the product of a 4-simplex with successive "
            "facet wedges of the dual base simple polytope."
        ),
        "proof_lemma": {
            "statement": (
                "A facet wedge of a simple convex polytope is again simple "
                "and convex, and its dual boundary is the simplicial wedge. "
                "Iteration realizes every positive integral multiwedge. "
                "Joins of simplicial polytope boundaries are polytopal."
            ),
            "primary_reference": (
                "S. Choi and H. Park, Wedge Operations and Torus Symmetries, "
                "Tohoku Math. J. 68 (2016), Proposition 2.2(3) and Lemma 2.3"
            ),
            "url": (
                "https://www.jstage.jst.go.jp/article/tmj/68/1/"
                "68_1458248864/_pdf"
            ),
        },
        "proof_source": {
            "path": str(PROOF_SOURCE.relative_to(ROOT)),
            "sha256": _file_sha256(PROOF_SOURCE),
        },
        "selection_provenance": _selection_provenance(),
        "audits": audits,
        "algebro_geometric_corollary": {
            "claim": (
                "For each record there is a smooth projective toric 9-fold X "
                "of Picard rank five whose reduced toric boundary D is SNC, "
                "satisfies K_X + D = 0, and has dual intersection complex "
                "equal to the target nerve."
            ),
            "base_models": {
                "44": (
                    "an integral Delzant hexagon, for example the moment "
                    "polygon of the blow-up of P^2 at its three torus fixed points"
                ),
                "43": (
                    "the dual stacking sequence realized by toric corner cuts "
                    "of a large integral Delzant tetrahedron, including the "
                    "final cut at the new vertex {u,11,12}"
                ),
                "47": (
                    "the dual stacking sequence realized by toric corner cuts "
                    "at {u,1,2}, {u,2,11}, and {u,1,11}"
                ),
            },
            "delzant_wedge": (
                "With a primitive integral slack ell_F, use 0 <= t <= "
                "ell_F(x). The new primitive inward normals are e_{d+1} and "
                "(a_F,-1); together with the old normals they form a lattice "
                "basis at every new vertex. Primitive coefficient one is "
                "essential: replacing ell_F by k ell_F can introduce index k."
            ),
            "product": (
                "After the facet wedges, take the product with a Delzant "
                "4-simplex. The result has dimension nine and fourteen facets."
            ),
            "divisor_encoding": (
                "The fourteen torus-invariant prime divisors correspond to "
                "facets. A set has nonempty intersection exactly when its rays "
                "span a fan cone, equivalently when the corresponding moment-"
                "polytope facets have a common face. Thus Dual(D) is the nerve."
            ),
            "chow_ring": (
                "A^*(X) is Z[x_0,...,x_13]/(I_lin + I_SR), where I_SR is "
                "generated exactly by the six quadratic and four quintic "
                "squarefree monomials of the archived minimal nonfaces."
            ),
            "carrier_caveat": (
                "The labelled incidence complex is canonical for the toric "
                "boundary pair (X,D), or for a distinguished toric presentation; "
                "the abstract variety or abstract graded Chow ring alone need "
                "not recover the boundary labels."
            ),
            "polarized_toric_dictionary": {
                "compactness": (
                    "An ample toric divisor H has a compact moment polytope Q; "
                    "ordinary polytope compactness is projectivity plus the "
                    "choice of polarization."
                ),
                "top_intersection": (
                    "H^9 = 9! times the lattice-normalized volume of Q."
                ),
                "facet_intersection": (
                    "H^8 D_i = 8! times the lattice-normalized volume of the "
                    "facet corresponding to D_i."
                ),
                "invariant_measure": (
                    "The moment-map pushforward of symplectic volume is the "
                    "Duistermaat--Heckman/Lebesgue measure on Q."
                ),
                "gale_rank": (
                    "Pic(X)_R = R^14/M_R has dimension five; quotienting by "
                    "the polarization/support class has dimension four, equal "
                    "to the Gale rank 14-9-1."
                ),
                "coxeter_caveat": (
                    "These ordinary volumes and relation spaces do not supply "
                    "the rank-ten Lorentzian Coxeter form, angle labels, or "
                    "hyperbolic chamber volume."
                ),
            },
            "picard_rank": 5,
        },
        "updated_theorem_partition": {
            "independently_certified_impossible": 119,
            "certified_explicit_ordinary_polytopal": 5,
            "theorem_undecided": 18,
            "total_graph_shadows": 142,
        },
        "nonclaims": [
            "No Coxeter labels are supplied.",
            "No Lorentzian Coxeter Gram matrix is supplied.",
            "Ordinary polytopality does not imply hyperbolic Coxeter realizability.",
            "The toric divisor intersection ring does not by itself encode Coxeter angles.",
        ],
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_multiwedge_ordinary_polytopality":
        raise ValueError("unexpected multiwedge certificate kind")
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("multiwedge certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("multiwedge certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("multiwedge certificate is not frozen")
    return {
        "record_indices": [audit["record_index"] for audit in rebuilt["audits"]],
        "ordinary_polytopal_count": len(rebuilt["audits"]),
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    del argv
    manifest = build()
    OUTPUT.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    summary = verify(manifest)
    print(json.dumps({"output": str(OUTPUT), **summary}, indent=2))


if __name__ == "__main__":
    main()
