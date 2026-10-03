\
"""Run all command-line experiments in order.

The QoolQit experiment may take longer than the pure NumPy/SciPy reference
experiments because it performs quantum emulation.

Usage:
    python run_all.py
"""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent

SCRIPTS = [
    "experiments/01_two_site_reference.py",
    "experiments/02_three_site_reference.py",
    "experiments/03_qoolqit_three_atom.py",
    "experiments/04_compare_three_site.py",
]


def main():
    for script in SCRIPTS:
        path = ROOT / script
        print("\n" + "=" * 72)
        print(f"RUNNING: {script}")
        print("=" * 72)

        subprocess.run(
            [sys.executable, str(path)],
            cwd=ROOT,
            check=True,
        )


if __name__ == "__main__":
    main()
