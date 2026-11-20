#!/usr/bin/env python3
"""Assemble and certify the global LP from retained local energy theorems.

Run only through run_bounded_research_check_20261018.py --parent-lock.
The input modules are reviewed theorem/domain transcriptions, not geometric
generators or replacements for their proofs. No old global bound is an input.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import signal
import sys
import tempfile
import time

STARTED = time.monotonic()
ROOT = Path(__file__).resolve().parents[1]
HELPER = "scripts/finite_volume_lp_matrix_20261106.py"
HELPER_SHA = "50a772e45c009f4aa0b1611d3b36eb24c806b468ebdfb1d8585acfee6ff9b88c"
MAX_BRANCHES = 20_000
MAX_SLOTS = 29
MAX_CUTS = 256
MAX_DOCUMENT = 64 * 1024**2


def require(ok, message):
    if not ok:
        raise ValueError(message)
    if time.monotonic() - STARTED > 53:
        raise TimeoutError("53-second entry allowance exhausted")


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       allow_nan=False) + "\n").encode()


def write_new(path, value):
    raw = encoded(value)
    require(len(raw) <= MAX_DOCUMENT, "Document exceeds its 64-MiB cap")
    with path.open("xb") as stream:
        stream.write(raw)
    return {"path": str(path.relative_to(ROOT)), "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest()}


def load_reviewed_module(relative, expected, index):
    path = ROOT / relative
    require(path.resolve() == path and path.parent == ROOT / "scripts"
            and path.suffix == ".py", "Use an unaliased repository source")
    require(path.stat().st_size <= 1024**2 and digest(path) == expected,
            "Reviewed catalogue source changed")
    name = f"existing_energy_catalogue_{index}"
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, "Cannot load catalogue")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    result = module.catalogue()
    require(digest(path) == expected, "Catalogue changed while evaluated")
    return result


def validate_and_join(catalogues, module_pins, threshold):
    branches, sources, premises, coverage = [], dict(module_pins), {}, []
    dimensions = set()
    for catalogue in catalogues:
        require(type(catalogue) is dict, "Expected a complete catalogue")
        require(type(catalogue["family"]) is str and catalogue["coverage"],
                "Declare the family and complete coverage premise")
        require(type(catalogue["dimensions"]) is list
                and all(type(n) is int and 1 <= n <= 29
                        for n in catalogue["dimensions"]), "Invalid actual ranks")
        declared = set(catalogue["dimensions"])
        dimensions.update(declared)
        coverage.append({"family": catalogue["family"],
                         "coverage": catalogue["coverage"]})
        for path, sha in catalogue["sources"].items():
            require(path not in sources or sources[path] == sha,
                    "Conflicting proof source bindings")
            sources[path] = sha
        for key, theorem in catalogue["premises"].items():
            require(key not in premises or premises[key] == theorem,
                    "Conflicting theorem identifiers")
            require(type(key) is str and theorem, "Missing theorem statement/domain")
            premises[key] = theorem
        require(type(catalogue["branches"]) is list, "Expected all branch records")
        seen_dimensions = set()
        for original in catalogue["branches"]:
            require(type(original) is dict and original.get("domain"),
                    "Every branch needs its complete applicability domain")
            n = original["dimension"]
            require(type(n) is int and n in declared, "Undeclared actual rank")
            seen_dimensions.add(n)
            branches.append(original)
            require(len(branches) <= MAX_BRANCHES, "Branch inventory cap reached")
        require(seen_dimensions == declared, "An actual-rank cover is absent")
    require(dimensions == set(range(min(dimensions), 30)),
            "The fixed catalogue must cover every rank through 29")
    require(min(dimensions) <= threshold <= 30, "Threshold outside fixed support")
    keys = set()
    for branch in branches:
        key = branch["key"]
        require(type(key) is str and key not in keys, "Repeated branch identifier")
        keys.add(key)
        slots = branch["slots"]
        require(type(slots) is list and 1 <= len(slots) <= MAX_SLOTS,
                "Invalid complete energy account")
        slot_ids = [slot["key"] for slot in slots]
        require(all(type(key) is str for key in slot_ids)
                and len(set(slot_ids)) == len(slot_ids), "Repeated energy slot")
        ranks = [slot.get("rank") for slot in slots]
        if all(type(rank) is int for rank in ranks):
            require(all(rank >= 0 for rank in ranks)
                    and sum(ranks) == branch["dimension"], "Rank account differs")
        rows, row_keys = branch["rows"], set()
        require(type(rows) is list and len(rows) <= MAX_CUTS, "Cut cap reached")
        for row in rows:
            require(row["key"] not in row_keys, "Repeated local energy cut")
            row_keys.add(row["key"])
            support = row["support"]
            require(type(support) is list and support
                    and len(set(support)) == len(support)
                    and set(support) <= set(slot_ids), "Invalid funded support")
            lower = Fraction(row["lower"])
            require(type(row["lower"]) in (str, int) and lower >= 0
                    and max(lower.numerator.bit_length(),
                            lower.denominator.bit_length()) <= 512,
                    "Use bounded exact nonnegative local floors")
            require(row["theorem"] in premises, "Unbound local theorem")
        exit_record = branch.get("exit")
        if exit_record is not None:
            require(type(exit_record) is dict and exit_record.get("hypotheses")
                    and exit_record["theorem"] in premises,
                    "A local contradiction needs its full theorem/domain join")
    for relative, sha in sources.items():
        path = ROOT / relative
        require(path.resolve() == path and ROOT in path.parents
                and path.is_file() and path.stat().st_size <= 4 * 1024**2,
                "Missing, aliased or oversized proof source")
        require(digest(path) == sha, "Retained proof source changed: " + relative)
    selected = [b for b in branches if b["dimension"] >= threshold]
    totals = {"branches": len(selected),
              "slots": sum(len(b["slots"]) for b in selected),
              "cuts": sum(len(b["rows"]) for b in selected),
              "support_terms": sum(len(r["support"]) for b in selected
                                   for r in b["rows"])}
    for key, cap in (("branches", MAX_BRANCHES), ("slots", 60_000),
                     ("cuts", 80_000), ("support_terms", 160_000)):
        require(totals[key] <= cap, "Complete assembly exceeds " + key + " cap")
    return {"schema": "existing-proof-global-energy-catalogue-v1",
            "threshold": threshold, "fixed_dimensions": sorted(dimensions),
            "external_ceiling": 29, "coverage": coverage,
            "sources": sources, "premises": premises,
            "branches": selected, "fixed_branch_count": len(branches),
            "assembly_counts": totals}


def assemble(catalogue, helper, catalogue_sha):
    metadata = {
        "model_scope": "global LP from proved local energy inequalities and complete domains",
        "coverage_premises": catalogue["coverage"],
        "retained_local_theorems": catalogue["premises"],
        "catalogue_sha256": catalogue_sha,
        "threshold": catalogue["threshold"], "external_compact_ceiling": 29,
        "same_original_vector_and_orthogonal_slots": True,
        "local_theorem_validity_is_retained_input": True,
        "no_global_dimension_exclusion_used_as_a_local_row": True,
    }
    primal = helper.LinearProgram("existing_proof_energy_primal", dict(metadata))
    dual = helper.LinearProgram("existing_proof_energy_capacity_dual", dict(metadata))
    primal.add_variable("unit_deficit", lower=1, upper=1)
    dual.add_variable("eta", lower=1)
    primal_rows, dual_rows, deltas = [], [], []
    for index, branch in enumerate(catalogue["branches"]):
        prefix = f"b{index}"
        mass, delta = prefix + ".mass", prefix + ".delta"
        primal.add_variable(mass)
        primal.add_variable(delta)
        deltas.append(delta)
        energies = {}
        for j, slot in enumerate(branch["slots"]):
            variable = f"{prefix}.energy{j}"
            primal.add_variable(variable)
            energies[slot["key"]] = variable
        provenance = {"branch": branch["key"], "dimension": branch["dimension"]}
        primal_rows.append(helper.LinearRow(prefix + ".deficit",
            {delta: Fraction(1), mass: Fraction(-1),
             **{variable: Fraction(1) for variable in energies.values()}},
            "=", Fraction(0), provenance))
        packing = {key: {} for key in energies}
        budget = {"eta": Fraction(-1)}
        for j, row in enumerate(branch["rows"]):
            lower = Fraction(row["lower"])
            primal_rows.append(helper.LinearRow(f"{prefix}.local{j}",
                {mass: -lower, **{energies[key]: Fraction(1)
                                  for key in row["support"]}},
                ">=", Fraction(0), {**provenance, "local_row": row["key"],
                                     "theorem": row["theorem"]}))
            if branch.get("exit") is None:
                weight = f"{prefix}.weight{j}"
                dual.add_variable(weight)
                budget[weight] = lower
                for key in row["support"]:
                    packing[key][weight] = Fraction(1)
        if branch.get("exit") is not None:
            primal_rows.append(helper.LinearRow(prefix + ".local_exit",
                {mass: Fraction(1)}, "=", Fraction(0),
                {**provenance, "local_exit": branch["exit"]}))
            # Its free exit multiplier is constructively eta; no capacity is used.
        else:
            for j, coefficients in enumerate(packing.values()):
                dual_rows.append(helper.LinearRow(f"{prefix}.capacity{j}",
                    coefficients, "<=", Fraction(1), provenance))
            dual_rows.append(helper.LinearRow(prefix + ".payment", budget,
                ">=", Fraction(0), provenance))
    primal_rows.append(helper.LinearRow("positive_total_deficit",
        {"unit_deficit": Fraction(-1),
         **{delta: Fraction(1) for delta in deltas}}, "=", Fraction(0),
        {"meaning": "normalize the strict total energy deficit without epsilon"}))
    # A zero objective can choose a tight payment ratio with a very large
    # denominator. Maximizing the sum of available local payments preserves
    # feasibility and favors saturated capacity vertices. Exact checks remain
    # authoritative even if the numerical proposal is not an optimum.
    dual.set_objective({variable: -value for row in dual_rows
                        if row.name.endswith(".payment")
                        for variable, value in row.coefficients.items()
                        if variable != "eta" and value})
    # eta is lower-bounded by one: a feasible point is the complete certificate.
    return primal, primal_rows, dual, dual_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("export-solve", "verify"))
    parser.add_argument("--catalogue", action="append", nargs=2,
                        metavar=("SOURCE", "SHA256"), required=True)
    parser.add_argument("--threshold", type=int, required=True)
    parser.add_argument("--directory", type=Path, required=True)
    args = parser.parse_args()
    require(1 <= len(args.catalogue) <= 2, "Use compact alone or arithmetic plus compact")
    require(digest(ROOT / HELPER) == HELPER_SHA, "Exact LP helper changed")
    import finite_volume_lp_matrix_20261106 as helper
    module_pins = dict(args.catalogue)
    require(len(module_pins) == len(args.catalogue), "Repeated catalogue module")
    catalogues = [load_reviewed_module(path, sha, i)
                  for i, (path, sha) in enumerate(args.catalogue)]
    catalogue = validate_and_join(catalogues, module_pins, args.threshold)
    catalogue["sources"][HELPER] = HELPER_SHA
    catalogue_bytes = encoded(catalogue)
    require(len(catalogue_bytes) <= MAX_DOCUMENT, "Catalogue exceeds its cap")
    catalogue_sha = hashlib.sha256(catalogue_bytes).hexdigest()
    primal, primal_rows, dual, dual_rows = assemble(catalogue, helper, catalogue_sha)
    directory = (ROOT / args.directory).resolve()
    require(directory.parent == ROOT / "results", "Use a fresh results child")
    if args.stage == "export-solve":
        directory.mkdir(exist_ok=False)
        artifacts = {"catalogue": write_new(directory / "catalogue.json", catalogue)}
        artifacts["primal"] = primal.export(directory / "primal.jsonl", primal_rows)
        artifacts["dual"] = dual.export(directory / "dual.jsonl", dual_rows)
        remaining = 53 - (time.monotonic() - STARTED)
        require(remaining >= 8, "Insufficient allowance for solver plus exact checks")
        solved = helper.numerical_candidate(directory / "dual.jsonl",
                                            seconds=min(25, remaining - 5))
        require(solved["model_sha256"] == artifacts["dual"]["sha256"],
                "Solver report is not bound to the exported dual")
        if solved["status"] == "EXACT_LP_FEASIBLE":
            require(solved["exact_check"]["model_sha256"]
                    == artifacts["dual"]["sha256"], "Exact check has another dual")
        artifacts["candidate"] = write_new(directory / "candidate.json", solved)
        final_hashes = {"catalogue.json": catalogue_sha,
                        "primal.jsonl": artifacts["primal"]["sha256"],
                        "dual.jsonl": artifacts["dual"]["sha256"],
                        "candidate.json": artifacts["candidate"]["sha256"]}
        status = solved["status"]
    else:
        require((directory / "catalogue.json").read_bytes() == catalogue_bytes,
                "Complete regenerated catalogue differs")
        with tempfile.TemporaryDirectory(prefix="energy-lp-replay-") as temp:
            temp = Path(temp)
            primal.export(temp / "primal.jsonl", primal_rows)
            dual.export(temp / "dual.jsonl", dual_rows)
            regenerated = {name: digest(temp / name)
                           for name in ("primal.jsonl", "dual.jsonl")}
            for name in regenerated:
                require(regenerated[name] == digest(directory / name),
                        "Fresh complete LP replay differs: " + name)
        raw = (directory / "candidate.json").read_bytes()
        require(len(raw) <= MAX_DOCUMENT, "Candidate exceeds cap")
        solved = json.loads(raw)
        require(solved["status"] == "EXACT_LP_FEASIBLE", "No exact candidate to replay")
        require(solved["model_sha256"] == regenerated["dual.jsonl"]
                and solved["exact_check"]["model_sha256"] == regenerated["dual.jsonl"],
                "Candidate is not bound to the freshly regenerated dual")
        exact = helper.check_candidate(directory / "dual.jsonl", solved["values"])
        require(exact["violation_count"] == 0
                and exact["model_sha256"] == regenerated["dual.jsonl"],
                "Fresh exact capacity check failed or checked another dual")
        final_hashes = {**regenerated, "catalogue.json": catalogue_sha,
                        "candidate.json": hashlib.sha256(raw).hexdigest()}
        artifacts, status = {}, exact["status"]
    for relative, sha in catalogue["sources"].items():
        require(digest(ROOT / relative) == sha, "Source changed before final seal")
    for name, sha in final_hashes.items():
        require(digest(directory / name) == sha, "Output changed before final seal")
    report = {
        "schema": "existing-proof-global-energy-lp-execution-v1",
        "stage": args.stage, "status": status,
        "threshold": args.threshold, "fixed_dimensions": catalogue["fixed_dimensions"],
        "selected_branches": len(catalogue["branches"]),
        "local_exits": sum(b.get("exit") is not None for b in catalogue["branches"]),
        "primal_variables": len(primal.variables), "primal_rows": len(primal_rows),
        "dual_variables": len(dual.variables), "dual_rows": len(dual_rows),
        "catalogue_sha256": catalogue_sha, "source_bindings": catalogue["sources"],
        "complete_artifact_hashes": final_hashes,
        "retained_coverage_premises": catalogue["coverage"],
        "artifacts": artifacts, "elapsed_seconds": round(time.monotonic() - STARTED, 3),
        "normalized_primal_infeasibility_certified": status == "EXACT_LP_FEASIBLE",
        "conditional_on_retained_geometry_local_theorems_and_complete_domains": True,
        "new_primitive_geometric_relaxation_solved": False,
    }
    write_new(directory / ("completed.json" if args.stage == "export-solve"
                           else "verified.json"), report)
    print(json.dumps({key: report[key] for key in (
        "stage", "status", "threshold", "selected_branches", "local_exits",
        "primal_variables", "primal_rows", "dual_variables", "dual_rows",
        "normalized_primal_infeasibility_certified", "elapsed_seconds")}))


if __name__ == "__main__":
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("Entry deadline")))
    signal.setitimer(signal.ITIMER_REAL, 53)
    main()
