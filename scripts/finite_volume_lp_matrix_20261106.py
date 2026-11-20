#!/usr/bin/env python3
"""Exact sparse LP transport and candidate checking for the finite-volume work.

Source preparation, 6 November 2026. No complete geometric matrix or dimension
claim is supplied by this module. Producers must state their coverage premises
and distinguish the full corner model from a restricted family of proved rules.
Run numerical or certificate work through run_proof_batch_20261024.py.
"""
from __future__ import annotations

import argparse
from array import array
from dataclasses import dataclass, field
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import re
import warnings

SCHEMA = "finite-volume-lp-exact-v1"
MAX_BYTES = 256 * 1024**2
MAX_LINE = 16 * 1024**2
MAX_ROWS = 1_000_000
MAX_NONZEROS = 5_000_000
MAX_VARIABLES = 500_000
NAME = re.compile(r"[A-Za-z][A-Za-z0-9_.:/-]*\Z")
RATIONAL = re.compile(r"[+-]?[0-9]{1,512}(?:/[0-9]{1,512})?\Z")


def rational(value):
    """Reject binary floats in the authoritative coefficient representation."""
    if type(value) not in (int, str, Fraction):
        raise TypeError("LP coefficients must be integers, rational strings or Fractions")
    if isinstance(value, str) and (len(value) > 1026 or not RATIONAL.fullmatch(value)):
        raise ValueError("Use bounded integer or numerator/denominator notation")
    result = Fraction(value)
    if max(result.numerator.bit_length(), result.denominator.bit_length()) > 1700:
        raise ValueError("Rational coefficient exceeds the declared size limit")
    return result


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def named(value):
    if not isinstance(value, str) or not NAME.fullmatch(value):
        raise ValueError(f"Invalid LP name: {value!r}")
    return value


def file_digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Repeated JSON key: {key}")
        result[key] = value
    return result


def parse(raw):
    return json.loads(raw, object_pairs_hook=unique_object,
                      parse_constant=lambda value: (_ for _ in ()).throw(
                          ValueError("Nonfinite JSON value")))


@dataclass(frozen=True)
class LinearRow:
    name: str
    coefficients: dict[str, Fraction]
    sense: str
    rhs: Fraction
    provenance: dict


@dataclass
class LinearProgram:
    """Variable registry; rows are streamed rather than retained in memory."""

    name: str
    metadata: dict
    variables: dict = field(default_factory=dict)
    objective: dict = field(default_factory=dict)

    def add_variable(self, name, lower=0, upper=None):
        named(name)
        if name in self.variables:
            raise ValueError(f"Repeated variable: {name}")
        low = None if lower is None else rational(lower)
        high = None if upper is None else rational(upper)
        if low is not None and high is not None and low > high:
            raise ValueError("Inconsistent variable bounds")
        self.variables[name] = (low, high)

    def set_objective(self, coefficients):
        if not set(coefficients) <= self.variables.keys():
            raise ValueError("Unknown objective variable")
        self.objective = {name: rational(value) for name, value in coefficients.items()
                          if rational(value)}

    def export(self, path, rows):
        """Write a new exact JSONL model; never overwrite an existing artifact."""
        target = Path(path)
        temporary = target.with_name(target.name + ".partial")
        if target.exists() or temporary.exists():
            raise FileExistsError("Choose a fresh model path")
        if not 0 < len(self.variables) <= MAX_VARIABLES:
            raise ValueError("Variable count exceeds transport limits")
        named(self.name)
        if not self.metadata.get("model_scope") or not self.metadata.get("coverage_premises"):
            raise ValueError("State model scope and coverage premises explicitly")
        names = list(self.variables)
        indices = {name: index for index, name in enumerate(names)}
        header = dict(kind="header", schema=SCHEMA, name=self.name,
                      metadata=self.metadata,
                      variables=[dict(name=name,
                                      lower=None if lo is None else str(lo),
                                      upper=None if hi is None else str(hi))
                                 for name, (lo, hi) in self.variables.items()],
                      objective=[[indices[name], str(value)]
                                 for name, value in self.objective.items()])
        count = nonzeros = size = 0
        checksum = hashlib.sha256()
        target.parent.mkdir(parents=True, exist_ok=True)
        owned_partial = False
        try:
            with temporary.open("xb") as stream:
                owned_partial = True
                raw = encoded(header)
                if len(raw) > MAX_LINE:
                    raise ValueError("Model header exceeds line limit")
                stream.write(raw)
                size += len(raw)
                for row in rows:
                    named(row.name)
                    if row.sense not in ("<=", ">=", "=") or not row.provenance:
                        raise ValueError("Each row needs a sense and provenance")
                    if not set(row.coefficients) <= indices.keys():
                        raise ValueError(f"Unknown variable in {row.name}")
                    terms = sorted((indices[name], rational(value))
                                   for name, value in row.coefficients.items()
                                   if rational(value))
                    item = dict(kind="row", ordinal=count, name=row.name,
                                coefficients=[[i, str(value)] for i, value in terms],
                                sense=row.sense, rhs=str(rational(row.rhs)),
                                provenance=row.provenance)
                    raw = encoded(item)
                    count += 1
                    nonzeros += len(terms)
                    size += len(raw)
                    if (len(raw) > MAX_LINE or size > MAX_BYTES
                            or count > MAX_ROWS or nonzeros > MAX_NONZEROS):
                        raise ValueError("Model exceeds a declared transport limit; shard it")
                    stream.write(raw)
                    checksum.update(raw)
                footer = encoded(dict(kind="end", rows=count, nonzeros=nonzeros,
                                      row_sha256=checksum.hexdigest()))
                if size + len(footer) > MAX_BYTES:
                    raise ValueError("Model exceeds byte limit")
                stream.write(footer)
            os.link(temporary, target)
            temporary.unlink()
        except BaseException:
            if owned_partial:
                temporary.unlink(missing_ok=True)
            raise
        return dict(path=str(target), sha256=file_digest(target),
                    variables=len(names), rows=count, nonzeros=nonzeros)


def read_model(path, expected_sha256=None):
    """Yield a validated header, rows and footer; refuse truncated transports."""
    with Path(path).open("rb") as stream:
        total = count = nonzeros = 0
        checksum = hashlib.sha256()
        transport = hashlib.sha256()
        raw = stream.readline(MAX_LINE + 1)
        if len(raw) > MAX_LINE or not raw.endswith(b"\n"):
            raise ValueError("Missing or oversized header")
        total += len(raw)
        transport.update(raw)
        header = parse(raw)
        if header.get("kind") != "header" or header.get("schema") != SCHEMA:
            raise ValueError("Unsupported model")
        named(header["name"])
        variables = header["variables"]
        if not isinstance(variables, list) or not 0 < len(variables) <= MAX_VARIABLES:
            raise ValueError("Invalid variable domain")
        names = [named(v["name"]) for v in variables]
        if len(names) != len(set(names)):
            raise ValueError("Duplicate variables")
        for variable in variables:
            low, high = variable["lower"], variable["upper"]
            if low is not None:
                low = rational(low)
            if high is not None:
                high = rational(high)
            if low is not None and high is not None and low > high:
                raise ValueError("Invalid variable bounds")
        objective_indices = set()
        for index, value in header["objective"]:
            if (type(index) is not int or not 0 <= index < len(variables)
                    or index in objective_indices):
                raise ValueError("Invalid objective index")
            objective_indices.add(index)
            rational(value)
        if (not header["metadata"].get("model_scope")
                or not header["metadata"].get("coverage_premises")):
            raise ValueError("Missing mathematical scope")
        yield header
        while True:
            raw = stream.readline(MAX_LINE + 1)
            total += len(raw)
            if not raw or not raw.endswith(b"\n") or len(raw) > MAX_LINE or total > MAX_BYTES:
                raise ValueError("Truncated or oversized model")
            transport.update(raw)
            item = parse(raw)
            if item.get("kind") == "end":
                if (type(item["rows"]) is not int or type(item["nonzeros"]) is not int
                        or item["rows"] != count or item["nonzeros"] != nonzeros
                        or item["row_sha256"] != checksum.hexdigest()
                        or stream.read(1)):
                    raise ValueError("Model footer mismatch or trailing content")
                consumed_hash = transport.hexdigest()
                if expected_sha256 is not None and consumed_hash != expected_sha256:
                    raise ValueError("Consumed model bytes differ from the expected hash")
                item["_transport_sha256"] = consumed_hash
                yield item
                return
            if (item.get("kind") != "row" or type(item.get("ordinal")) is not int
                    or item["ordinal"] != count):
                raise ValueError("Missing or reordered row")
            named(item["name"])
            if item["sense"] not in ("<=", ">=", "=") or not item["provenance"]:
                raise ValueError("Invalid row scope")
            previous = -1
            for index, value in item["coefficients"]:
                if type(index) is not int or not previous < index < len(variables):
                    raise ValueError("Repeated, unsorted or invalid coefficient index")
                if not rational(value):
                    raise ValueError("Zero coefficient must be omitted")
                previous = index
            rational(item["rhs"])
            count += 1
            nonzeros += len(item["coefficients"])
            if count > MAX_ROWS or nonzeros > MAX_NONZEROS:
                raise ValueError("Model exceeds transport limits")
            checksum.update(raw)
            yield item


def check_candidate(path, values):
    """Check only this LP, exactly. Geometric completeness remains a premise."""
    initial_hash = file_digest(path)
    records = read_model(path, expected_sha256=initial_hash)
    header = next(records)
    expected = [v["name"] for v in header["variables"]]
    if set(values) != set(expected):
        raise ValueError("Candidate variable domain differs from the model")
    vector = [rational(values[name]) for name in expected]
    violations = []
    violation_count = rows = 0
    minimum_slack = None
    for variable, value in zip(header["variables"], vector, strict=True):
        for side in ("lower", "upper"):
            bound = variable[side]
            if bound is not None:
                slack = value - rational(bound) if side == "lower" else rational(bound) - value
                if slack < 0:
                    violation_count += 1
                    if len(violations) < 20:
                        violations.append(dict(variable=variable["name"], bound=side,
                                               slack=str(slack)))
    for row in records:
        if row["kind"] == "end":
            continue
        lhs = sum((rational(value) * vector[index]
                   for index, value in row["coefficients"]), Fraction())
        rhs = rational(row["rhs"])
        slack = rhs - lhs if row["sense"] == "<=" else lhs - rhs
        if row["sense"] == "=":
            slack = -abs(lhs - rhs)
        minimum_slack = slack if minimum_slack is None else min(minimum_slack, slack)
        rows += 1
        if slack < 0:
            violation_count += 1
            if len(violations) < 20:
                violations.append(dict(ordinal=row["ordinal"], row=row["name"],
                                       slack=str(slack)))
    if file_digest(path) != initial_hash:
        raise ValueError("Model changed during exact verification")
    objective = sum((rational(value) * vector[index]
                     for index, value in header["objective"]), Fraction())
    return dict(status="EXACT_LP_FEASIBLE" if not violation_count else "EXACT_LP_VIOLATION",
                model_sha256=initial_hash, rows_checked=rows,
                variables_checked=len(vector), violation_count=violation_count,
                first_violations=violations,
                minimum_row_slack=None if minimum_slack is None else str(minimum_slack),
                objective=str(objective), model_scope=header["metadata"]["model_scope"],
                geometric_completeness_verified=False, dimension_bound_claimed=False)


def numerical_candidate(path, seconds=35, denominator=1_000_000):
    """One-thread numerical proposal, followed by independent rational checks."""
    import numpy as np
    import scipy
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix

    if not 0 < seconds <= 40 or not 1 <= denominator <= 10**9:
        raise ValueError("Invalid numerical proposal limits")
    before = file_digest(path)
    records = read_model(path, expected_sha256=before)
    header = next(records)
    n = len(header["variables"])
    objective = np.zeros(n)
    for index, value in header["objective"]:
        objective[index] = float(rational(value))
    banks = {kind: dict(i=array("i"), j=array("i"), a=array("d"), rhs=array("d"))
             for kind in ("ub", "eq")}
    for row in records:
        if row["kind"] == "end":
            continue
        bank = banks["eq" if row["sense"] == "=" else "ub"]
        index = len(bank["rhs"])
        sign = -1 if row["sense"] == ">=" else 1
        for column, value in row["coefficients"]:
            bank["i"].append(index)
            bank["j"].append(column)
            bank["a"].append(sign * float(rational(value)))
        bank["rhs"].append(sign * float(rational(row["rhs"])))
    matrices = {}
    for kind, bank in banks.items():
        if bank["rhs"]:
            matrices["A_" + kind] = coo_matrix(
                (np.frombuffer(bank["a"], dtype=np.float64),
                 (np.frombuffer(bank["i"], dtype=np.int32),
                  np.frombuffer(bank["j"], dtype=np.int32))),
                shape=(len(bank["rhs"]), n)).tocsr()
            matrices["b_" + kind] = np.frombuffer(bank["rhs"], dtype=np.float64)
    bounds = [(None if v["lower"] is None else float(rational(v["lower"])),
               None if v["upper"] is None else float(rational(v["upper"])))
              for v in header["variables"]]
    if file_digest(path) != before:
        raise ValueError("Model changed before solve")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = linprog(objective, bounds=bounds, method="highs-ds", **matrices,
                         options={"threads": 1, "time_limit": seconds, "presolve": True})
    report = dict(status="NUMERICAL_RESULT_ONLY", model_sha256=before,
                  scipy_version=scipy.__version__, solver_status=int(result.status),
                  solver_success=bool(result.success), solver_message=str(result.message),
                  warnings=[str(item.message) for item in caught],
                  dimension_bound_claimed=False, exact_infeasibility_proved=False)
    if result.x is None or not all(math.isfinite(float(value)) for value in result.x):
        return report
    values = {variable["name"]: str(Fraction(str(float(value))).limit_denominator(denominator))
              for variable, value in zip(header["variables"], result.x, strict=True)}
    check = check_candidate(path, values)
    if check["model_sha256"] != before:
        raise ValueError("Model changed after solve")
    report.update(values=values, exact_check=check,
                  status="EXACT_LP_FEASIBLE" if check["violation_count"] == 0
                  else "NUMERICAL_CANDIDATE_FAILED_EXACT_CHECK")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("solve", "verify"))
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seconds", type=float, default=35)
    args = parser.parse_args()
    if not os.environ.get("COXETER_CHECK_OWNER"):
        raise SystemExit("Use the repository's bounded research runner")
    if args.output.exists():
        raise FileExistsError("Choose a new result path")
    if args.mode == "solve":
        result = numerical_candidate(args.model, args.seconds)
    else:
        if args.candidate is None:
            parser.error("--candidate is required for verify")
        with args.candidate.open("rb") as stream:
            raw = stream.read(MAX_LINE + 1)
        if len(raw) > MAX_LINE:
            raise ValueError("Candidate exceeds size limit")
        candidate = parse(raw)
        if candidate.get("model_sha256") != file_digest(args.model):
            raise ValueError("Candidate is bound to a different model")
        result = check_candidate(args.model, candidate["values"])
        if result["model_sha256"] != candidate["model_sha256"]:
            raise ValueError("Verified model differs from the candidate's bound model")
    result["checker_sha256"] = file_digest(__file__)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({key: value for key, value in result.items() if key != "values"},
                     sort_keys=True, allow_nan=False), flush=True)
    if result["status"] != "EXACT_LP_FEASIBLE":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
