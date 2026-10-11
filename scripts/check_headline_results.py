"""Compare headline estimates in results/ with the values reported in the paper."""
import sys
from pathlib import Path

import pandas as pd

RESULTS = Path("results")

# (label, file, column, selector, expected, decimals)
CHECKS = [
    ("UK coefficient", "uk/primary/M0_primary_coefficient.csv", "estimate", 0.843, 3),
    ("UK standard error", "uk/primary/M0_primary_coefficient.csv", "se", 0.174, 3),
    ("UK average marginal effect", "uk/primary/M0_average_marginal_effect.csv", "estimate", 0.149, 3),
    ("Italy coefficient", "italy/models/italy_M0_primary_coefficient.csv", "estimate", 1.047, 3),
    ("Italy standard error", "italy/models/italy_M0_primary_coefficient.csv", "se", 0.074, 3),
    ("Italy average marginal effect", "italy/models/italy_M0_average_marginal_effect.csv", "estimate", 0.251, 3),
]


def main() -> int:
    failures = 0
    for label, rel, column, expected, decimals in CHECKS:
        path = RESULTS / rel
        if not path.exists():
            print(f"MISSING  {label}: {path}")
            failures += 1
            continue
        value = float(pd.read_csv(path)[column].iloc[0])
        ok = round(value, decimals) == expected
        failures += not ok
        print(f"{'OK' if ok else 'DIFFERS':8} {label}: {value:.{decimals + 3}f} (paper: {expected})")
    print("All headline values match." if failures == 0 else f"{failures} check(s) failed.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
