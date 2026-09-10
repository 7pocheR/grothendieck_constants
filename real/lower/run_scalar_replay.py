"""Bounded resident-allocation replay of the complete scalar partitions."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--seconds-per-part", type=int, default=1200)
    opt = ap.parse_args()
    assert 1 <= opt.workers <= 16 and 1 <= opt.seconds_per_part <= 3600
    checker = Path(__file__).with_name("check_scalars.py")
    common = ["--campaign", opt.campaign, "--output", opt.output, "--parts", "16"]
    subprocess.run([sys.executable, str(checker), "noise", *common], check=True, timeout=120)

    def part(i):
        subprocess.run([sys.executable, str(checker), "part", *common, "--part", str(i),
                        "--seconds", str(opt.seconds_per_part)], check=True,
                       timeout=opt.seconds_per_part+120)

    with ThreadPoolExecutor(max_workers=opt.workers) as pool:
        list(pool.map(part, range(16)))
    subprocess.run([sys.executable, str(checker), "finish", *common], check=True, timeout=120)


if __name__ == "__main__":
    main()
