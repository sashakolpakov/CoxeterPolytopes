"""Verification of archived Vinberg transcripts and their stopping criterion."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Protocol, Sequence

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.search.bugaenko_vinberg import (
    BugaenkoVinbergSearch,
    build_bugaenko_vinberg_transcript,
    run_bugaenko_vinberg_search,
)
from coxeter_polytopes.search.h11_vinberg import (
    H11VinbergSearch,
    build_h11_vinberg_transcript,
    run_h11_vinberg_search,
)
from coxeter_polytopes.search.h19_vinberg import (
    H19VinbergSearch,
    build_h19_vinberg_transcript,
    run_h19_vinberg_search,
)

from .chamber import ChamberVerification, verify_chamber_geometry
from .vinberg_volume import (
    VinbergVolumeVerification,
    verify_vinberg_volume_criterion,
)


class VinbergTranscriptError(ValueError):
    """Raised when a wall-enumeration transcript is not reproducible exactly."""


@dataclass(frozen=True, slots=True)
class VinbergCompletionVerification:
    search: BugaenkoVinbergSearch
    chamber: ChamberVerification
    diagram: VinbergVolumeVerification


class VerifiedVinbergPrefix(Protocol):
    """Structural type returned by an exact Vinberg transcript replay."""

    roots: Sequence[Sequence[Golden]]


@dataclass(frozen=True, slots=True)
class VinbergPrefixStoppingVerification:
    """Join of an exact accepted-root prefix and a finite-volume chamber."""

    chamber: ChamberVerification
    diagram: VinbergVolumeVerification


def _first_difference(left: Any, right: Any, path: str = "$.") -> str | None:
    if type(left) is not type(right):
        return f"{path} has type {type(left).__name__}, expected {type(right).__name__}"
    if isinstance(left, Mapping):
        if set(left) != set(right):
            missing = sorted(set(right) - set(left))
            extra = sorted(set(left) - set(right))
            return f"{path} key mismatch: missing={missing}, extra={extra}"
        for key in right:
            difference = _first_difference(left[key], right[key], f"{path}{key}.")
            if difference is not None:
                return difference
        return None
    if isinstance(left, list):
        if len(left) != len(right):
            return f"{path} has length {len(left)}, expected {len(right)}"
        for index, (left_value, right_value) in enumerate(zip(left, right, strict=True)):
            difference = _first_difference(left_value, right_value, f"{path}[{index}].")
            if difference is not None:
                return difference
        return None
    if left != right:
        return f"{path} is {left!r}, expected {right!r}"
    return None


def verify_bugaenko_vinberg_transcript(
    transcript: Mapping[str, Any],
) -> BugaenkoVinbergSearch:
    """Rerun the exact search and compare every archived enumeration counter."""

    scan = transcript.get("control_coordinate_scan")
    if not isinstance(scan, Mapping) or not isinstance(scan.get("identity_bound"), int):
        raise VinbergTranscriptError("transcript has no integer control-coordinate bound")
    identity_bound = scan["identity_bound"]
    expected = build_bugaenko_vinberg_transcript(identity_bound)
    difference = _first_difference(transcript, expected)
    if difference is not None:
        raise VinbergTranscriptError(difference)
    return run_bugaenko_vinberg_search(identity_bound)


def verify_h11_vinberg_transcript(
    transcript: Mapping[str, Any],
) -> H11VinbergSearch:
    """Rerun a bounded H11 norm-two search and compare its full transcript."""

    scan = transcript.get("control_coordinate_scan")
    if not isinstance(scan, Mapping) or not isinstance(scan.get("identity_bound"), int):
        raise VinbergTranscriptError("transcript has no integer control-coordinate bound")
    identity_bound = scan["identity_bound"]
    expected = build_h11_vinberg_transcript(identity_bound)
    difference = _first_difference(transcript, expected)
    if difference is not None:
        raise VinbergTranscriptError(difference)
    return run_h11_vinberg_search(identity_bound)


def verify_h19_vinberg_transcript(
    transcript: Mapping[str, Any],
) -> H19VinbergSearch:
    """Rerun a bounded H19 norm-two search and compare its full transcript."""

    scan = transcript.get("control_coordinate_scan")
    if not isinstance(scan, Mapping) or not isinstance(scan.get("identity_bound"), int):
        raise VinbergTranscriptError("transcript has no integer control-coordinate bound")
    identity_bound = scan["identity_bound"]
    expected = build_h19_vinberg_transcript(identity_bound)
    difference = _first_difference(transcript, expected)
    if difference is not None:
        raise VinbergTranscriptError(difference)
    return run_h19_vinberg_search(identity_bound)


def _decode_roots(value: Any) -> tuple[tuple[Golden, ...], ...]:
    if not isinstance(value, list):
        raise VinbergTranscriptError("certificate simple_roots must be an array")
    roots = []
    for i, vector in enumerate(value):
        if not isinstance(vector, list):
            raise VinbergTranscriptError(f"certificate simple_roots[{i}] must be an array")
        try:
            roots.append(tuple(Golden.from_json(entry) for entry in vector))
        except (TypeError, ValueError) as exc:
            raise VinbergTranscriptError(
                f"invalid certificate root coordinate in simple_roots[{i}]"
            ) from exc
    return tuple(roots)


def _decode_wall_gram(value: Any) -> tuple[tuple[Golden, ...], ...]:
    if not isinstance(value, list):
        raise VinbergTranscriptError("certificate wall_gram must be an array")
    rows: list[tuple[Golden, ...]] = []
    for row_index, row in enumerate(value):
        if not isinstance(row, list):
            raise VinbergTranscriptError(
                f"certificate wall_gram[{row_index}] must be an array"
            )
        try:
            rows.append(tuple(Golden.from_json(entry) for entry in row))
        except (TypeError, ValueError) as exc:
            raise VinbergTranscriptError(
                f"invalid certificate wall_gram entry in row {row_index}"
            ) from exc
    return tuple(rows)


def verify_vinberg_prefix_stopping(
    verified_prefix: VerifiedVinbergPrefix,
    certificate: Mapping[str, Any],
) -> VinbergPrefixStoppingVerification:
    """Use finite covolume to stop an exactly replayed Vinberg prefix.

    ``verified_prefix`` must be the result of a method-specific exact
    transcript verifier such as :func:`verify_bugaenko_vinberg_transcript`,
    not an arbitrary wall sample.  This function joins that accepted-root
    evidence to two independent finite-volume checks of the same outer
    chamber: exact extreme rays and Vinberg's diagram criterion.  The
    stopping theorem then proves that no later accepted simple root exists.

    Unlike :func:`verify_vinberg_volume_criterion`, this is the legitimate
    route by which a not-yet-complete accepted prefix *becomes* a complete
    simple-root list.  It supports both compact and cofinite noncompact
    certificates.
    """

    if not hasattr(verified_prefix, "roots"):
        raise TypeError("verified_prefix must expose replayed accepted roots")
    replayed_roots = tuple(tuple(root) for root in verified_prefix.roots)
    if len(replayed_roots) != len(set(replayed_roots)):
        raise VinbergTranscriptError("verified accepted prefix contains duplicate roots")
    certificate_roots = _decode_roots(certificate.get("simple_roots"))
    if len(certificate_roots) != len(set(certificate_roots)):
        raise VinbergTranscriptError("certificate contains duplicate roots")
    if set(certificate_roots) != set(replayed_roots):
        missing = set(replayed_roots) - set(certificate_roots)
        extra = set(certificate_roots) - set(replayed_roots)
        raise VinbergTranscriptError(
            "verified prefix and chamber certificate have different walls: "
            f"missing={len(missing)}, extra={len(extra)}"
        )

    chamber = verify_chamber_geometry(certificate)
    wall_gram = _decode_wall_gram(certificate.get("wall_gram"))
    dimension = certificate.get("dimension")
    if not isinstance(dimension, int):
        raise VinbergTranscriptError("certificate dimension must be an integer")
    diagram = verify_vinberg_volume_criterion(
        wall_gram,
        dimension,
        # The exact replay lists every facet of the *current outer
        # polyhedron*, so the diagram criterion may evaluate that polyhedron.
        # If it is finite-volume, Vinberg stopping then proves the distinct
        # conclusion that no later ambient simple root can be accepted.
        simple_roots_complete=True,
    )
    if diagram.finite_volume is not True:
        raise VinbergTranscriptError(
            "accepted prefix does not pass Vinberg's finite-volume diagram test"
        )
    if bool(chamber.cusps) != (diagram.decision == "finite_volume_noncompact"):
        raise VinbergTranscriptError(
            "extreme-ray and diagram checks disagree about compactness"
        )
    return VinbergPrefixStoppingVerification(chamber=chamber, diagram=diagram)


def verify_bugaenko_vinberg_completion(
    transcript: Mapping[str, Any], certificate: Mapping[str, Any]
) -> VinbergCompletionVerification:
    """Join exhaustive height ordering to the independent compactness proof.

    The search and geometry verifier remain separate.  This final check only
    establishes that they concern exactly the same oriented norm-two walls.
    Vinberg's stopping theorem then makes the list complete: a later root in
    height order cannot cut a chamber already proved compact.
    """

    search = verify_bugaenko_vinberg_transcript(transcript)
    stopping = verify_vinberg_prefix_stopping(search, certificate)
    if stopping.diagram.decision != "compact":
        raise VinbergTranscriptError("Bugaenko certificate must be compact")
    return VinbergCompletionVerification(
        search=search,
        chamber=stopping.chamber,
        diagram=stopping.diagram,
    )
