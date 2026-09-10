"""Serialize the bounded exact rational tail sums with large denominators."""
from pathlib import Path
import runpy
import sys

sys.set_int_max_str_digits(40000)
runpy.run_path(str(Path(__file__).with_name("check_rational_tails.py")),run_name="__main__")
