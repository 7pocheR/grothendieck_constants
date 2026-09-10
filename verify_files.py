"""Check the SHA-256 hashes of the distributed files; no mathematical replay."""

import hashlib
from pathlib import Path
import sys


def main():
    root = Path(__file__).resolve().parent
    manifest = root / "SHA256SUMS"
    seen = set()
    failures = []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", 1)
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("Invalid digest in SHA256SUMS")
        target = (root / relative).resolve()
        if not target.is_relative_to(root) or relative in seen:
            raise ValueError("Invalid or duplicate path in SHA256SUMS")
        seen.add(relative)
        if not target.is_file():
            failures.append(f"Missing: {relative}")
            continue
        actual = hashlib.sha256()
        with target.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                actual.update(chunk)
        if actual.hexdigest() != digest:
            failures.append(f"Mismatch: {relative}")
    if not seen:
        raise ValueError("Empty SHA256SUMS")
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print(f"PASS: {len(seen)} file hashes match. Mathematical verification is separate.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
