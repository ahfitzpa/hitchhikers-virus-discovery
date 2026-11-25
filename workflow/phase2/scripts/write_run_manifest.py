#!/usr/bin/env python3
"""Write an atomic, machine-readable run provenance manifest (run_manifest.json).

Records workflow version, git commit SHA (if inside a git checkout), Python and
OS details, UTC timestamp, CLI invocation, and SHA-256 checksums for input
manifests and configuration files so every HPC or local run is auditable.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import logging
import os
import platform
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

__version__ = "1.0.0"
LOGGER = logging.getLogger(__name__)


def sha256_file(path: str | Path) -> str:
    """Compute hex SHA-256 digest of a file using streaming 64 KiB chunks."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit_sha(repo_dir: str | Path) -> str | None:
    """Return HEAD git commit SHA if available, otherwise None."""
    try:
        res = subprocess.run(
            ["git", "-C", str(repo_dir), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
        return res.stdout.strip() if res.returncode == 0 and res.stdout.strip() else None
    except (OSError, subprocess.SubprocessError):
        return None


def build_manifest(
    workflow: str,
    input_files: list[str | Path],
    output_dir: str | Path,
    command: str | None = None,
    repo_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Assemble a JSON-serializable provenance dictionary."""
    checksums: dict[str, str] = {}
    for item in input_files:
        p = Path(item)
        if not p.is_file():
            raise FileNotFoundError(f"Input file not found for provenance checksum: {p}")
        checksums[p.name] = sha256_file(p)

    return {
        "schema_version": "1.0",
        "workflow": workflow,
        "workflow_version": __version__,
        "git_commit": git_commit_sha(repo_dir or Path.cwd()),
        "timestamp_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "command": command or " ".join(sys.argv),
        "output_dir": str(Path(output_dir)),
        "input_sha256": checksums,
    }


def write_manifest_atomic(manifest: dict[str, Any], dest: str | Path) -> None:
    """Atomically write provenance JSON to dest via sibling temp file + os.replace()."""
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{dest_path.name}.", suffix=".tmp", dir=dest_path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=2, sort_keys=True)
            fh.write("\n")
        os.replace(tmp_name, dest_path)
    except BaseException:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for writing a run provenance manifest."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workflow", required=True, help="Workflow name (e.g. phase1, phase2)")
    ap.add_argument("--input", action="append", default=[], help="Input manifest/config to SHA-256 checksum (repeatable)")
    ap.add_argument("--output-dir", required=True, help="Run output directory")
    ap.add_argument("--out", required=True, help="Destination JSON file (e.g. results/run_manifest.json)")
    ap.add_argument("--command", default=None, help="Optional command string to record")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args(argv)

    logging.basicConfig(level=getattr(logging, args.log_level), format="[%(levelname)s] %(message)s")
    manifest = build_manifest(
        workflow=args.workflow,
        input_files=args.input,
        output_dir=args.output_dir,
        command=args.command,
    )
    write_manifest_atomic(manifest, args.out)
    LOGGER.info("Wrote provenance manifest to %s", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
