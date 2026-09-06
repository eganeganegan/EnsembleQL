"""Command-line interface for EnsembleQL."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from .trajectory import load
from .query import explain


def _add_query_source(command: argparse.ArgumentParser) -> None:
    source = command.add_mutually_exclusive_group(required=True)
    source.add_argument("--query", help="query text")
    source.add_argument("--file", type=Path, help="path to a .eql query")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ensembleql", description="Query molecular events over time")
    subcommands = parser.add_subparsers(dest="command", required=True)
    query = subcommands.add_parser("query", help="execute an EnsembleQL query")
    query.add_argument("--topology", required=True, help="PDB topology")
    query.add_argument("--trajectory", required=True, help="trajectory file")
    query.add_argument(
        "--default-timestep",
        default="1ps",
        help="fallback frame spacing when trajectory timestamps are absent (default: 1ps)",
    )
    _add_query_source(query)
    query.add_argument("--format", choices=("table", "csv", "json"), default="table")
    explain_command = subcommands.add_parser("explain", help="show a query plan without scanning frames")
    explain_command.add_argument("--topology", required=True, help="PDB topology")
    _add_query_source(explain_command)
    explain_command.add_argument("--format", choices=("text", "json"), default="text")
    return parser


def _query_text(args) -> str:
    return args.query if args.query is not None else args.file.read_text(encoding="utf-8")


def _print_explanation(plan: dict) -> None:
    print(f"Streaming: {'yes' if plan['streaming'] else 'no'}")
    for title, key in (("Selections", "selections"), ("Observables", "observables"),
                       ("Frame predicates", "frame_predicates"), ("Temporal operations", "temporal_operations")):
        values = plan[key]
        print(f"{title} ({len(values)}):")
        for value in values:
            print(f"  - {value}")
    print("Plan:")
    print(plan["plan_tree"], end="")


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        text = _query_text(args)
        if args.command == "explain":
            plan = explain(text, topology=args.topology)
            if args.format == "json":
                json.dump(plan, sys.stdout, indent=2)
                print()
            else:
                _print_explanation(plan)
            return 0
        results = load(
            args.trajectory,
            topology=args.topology,
            default_timestep=args.default_timestep,
        ).query(text)
        records = results.to_records()
        if args.format == "json":
            json.dump(records, sys.stdout, indent=2)
            print()
        elif args.format == "csv":
            fields = ("start_ps", "end_ps", "duration_ps", "type", "selections", "metadata")
            writer = csv.DictWriter(sys.stdout, fieldnames=fields)
            writer.writeheader()
            for record in records:
                record["selections"] = json.dumps(record["selections"])
                record["metadata"] = json.dumps(record["metadata"])
                writer.writerow(record)
        else:
            print(results)
        return 0
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"ensembleql: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
