"""Create an auditable no-benchmark workflow for restricted PID environments."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re


def prepare(source, destination):
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError("Source and output must be different files")
    if source.resolve().parent != destination.resolve().parent:
        raise ValueError("Keep the generated workflow beside the source for relative paths")
    original = source.read_bytes()
    kept, removed = [], []
    for number, line in enumerate(original.splitlines(keepends=True), 1):
        if re.match(rb"^\s*benchmark:", line):
            expression = line.split(b"benchmark:", 1)[1].strip()
            if not expression or not re.fullmatch(rb'(?:f)?[\"\'][^\r\n]*[\"\']', expression):
                raise ValueError("Only single-line string benchmark directives are supported")
            removed.append(number)
        else:
            kept.append(line)
    generated = b"".join(kept)
    destination.write_bytes(generated)
    record = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": str(source), "output": str(destination),
        "source_sha256": hashlib.sha256(original).hexdigest(),
        "output_sha256": hashlib.sha256(generated).hexdigest(),
        "removed_benchmark_lines": removed,
        "transformation": "Remove only single-line benchmark directives; preserve all other bytes",
        "scope": "Execution workaround; no CPU/memory benchmark or new scientific result",
    }
    destination.with_suffix(destination.suffix + ".provenance.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("workflow/Snakefile"))
    parser.add_argument("--output", type=Path, default=Path("workflow/local_validation.smk"))
    args = parser.parse_args()
    record = prepare(args.source, args.output)
    print(f"Wrote {args.output}; omitted {len(record['removed_benchmark_lines'])} benchmarks")
