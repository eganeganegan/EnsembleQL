"""Command-line interface for EnsembleQL."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from .trajectory import load


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ensembleql", description="Query molecular events over time")
    subcommands = parser.add_subparsers(dest="command", required=True)
    query = subcommands.add_parser("query", help="execute an EnsembleQL query")
    query.add_argument("--topology", required=True, help="PDB topology")
    query.add_argument("--trajectory", required=True, help="XYZ trajectory")
    source = query.add_mutually_exclusive_group(required=True)
    source.add_argument("--query", help="query text")
    source.add_argument("--file", type=Path, help="path to a .eql query")
    query.add_argument("--format", choices=("table", "csv", "json"), default="table")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        text = args.query if args.query is not None else args.file.read_text(encoding="utf-8")
        results = load(args.trajectory, topology=args.topology).query(text)
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
