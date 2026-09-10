"""Run the fixed independent vertex tasks with a bounded number of workers."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evidence", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--seconds-per-part", type=int, default=1200)
    args = ap.parse_args()
    assert 1 <= args.workers <= 16 and 1 <= args.seconds_per_part <= 3600
    checker = Path(__file__).with_name("check_fibers.py")
    common = ["--evidence", args.evidence, "--output", args.output, "--parts", "16"]

    def part(i):
        cmd = [sys.executable, str(checker), "part", *common, "--part", str(i),
               "--seconds", str(args.seconds_per_part)]
        subprocess.run(cmd, check=True, timeout=args.seconds_per_part+120)

    subprocess.run([sys.executable, str(checker), "inventory", *common], check=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(part, range(16)))
    subprocess.run([sys.executable, str(checker), "finish", *common], check=True)


if __name__ == "__main__":
    main()
