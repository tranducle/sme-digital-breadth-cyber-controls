from pathlib import Path
import hashlib
import json
import pandas as pd
import numpy as np

ROOT = Path(r".")
EXEC = Path("results/next_wave")

panel_path = ROOT / "results/panel/analysis_ready_panel.csv"
raw_path = ROOT / "data/raw/UKDA-9356-stata/stata/stata13/csls_waves_1-4.dta"

panel = pd.read_csv(panel_path)
raw = pd.read_stata(
    raw_path,
    columns=["incident12", "incident13", "sector", "bus_reg2"],
    convert_categoricals=False,
)
raw["source_row"] = np.arange(1, len(raw) + 1)
d = panel.merge(raw, on="source_row", how="left", validate="one_to_one")
d["I13"] = d["incident13"].map({1: 1, 2: 0})
d["I12"] = d["incident12"].map({1: 1, 2: 0})

adjacent = []
gap2 = []

for oid, g in d.sort_values(["id", "wave"]).groupby("id"):
    rec = {int(row["wave"]): row for row in g.to_dict("records")}

    for w in [1, 2, 3]:
        if w in rec and w + 1 in rec:
            a = rec[w]
            b = rec[w + 1]
            adjacent.append({
                "id": oid,
                "wave_t": w,
                "wave_next": w + 1,
                "source_row_t": int(a["source_row"]),
                "source_row_next": int(b["source_row"]),
                "E": float(a["E"]),
                "R": float(a["R"]),
                "R3": float(a["R3"]),
                "E_score": int(a["E_score"]),
                "R_score": int(a["R_score"]),
                "R3_score": int(a["R3_score"]),
                "I13_t": int(a["I13"]),
                "I13_next": int(b["I13"]),
                "I12_t": int(a["I12"]),
                "I12_next": int(b["I12"]),
                "web_mode_next": int(b["web_mode"]),
                "sector": int(a["sector"]),
                "region2": int(a["bus_reg2"]),
            })

    for w in [1, 2]:
        if w in rec and w + 2 in rec:
            a = rec[w]
            b = rec[w + 2]
            gap2.append({
                "id": oid,
                "wave_t": w,
                "wave_next": w + 2,
                "source_row_t": int(a["source_row"]),
                "source_row_next": int(b["source_row"]),
                "E": float(a["E"]),
                "R": float(a["R"]),
                "R3": float(a["R3"]),
                "E_score": int(a["E_score"]),
                "R_score": int(a["R_score"]),
                "R3_score": int(a["R3_score"]),
                "I13_t": int(a["I13"]),
                "I13_next": int(b["I13"]),
                "I12_t": int(a["I12"]),
                "I12_next": int(b["I12"]),
                "web_mode_next": int(b["web_mode"]),
                "sector": int(a["sector"]),
                "region2": int(a["bus_reg2"]),
            })

adj = pd.DataFrame(adjacent)
g2 = pd.DataFrame(gap2)

present = set(zip(d["id"], d["wave"]))
elig = d[d["wave"].isin([1, 2, 3])].copy()
elig["retain_adjacent"] = [
    int((oid, wave + 1) in present)
    for oid, wave in zip(elig["id"], elig["wave"])
]
elig = elig[["source_row", "id", "wave", "E", "R", "I13", "retain_adjacent"]].copy()

required_adj = [
    "E", "R", "R3", "I13_t", "I13_next", "I12_t", "I12_next",
    "web_mode_next", "sector", "region2",
]
required_elig = ["E", "R", "I13", "retain_adjacent"]

checks = {
    "adjacent_n": len(adj) == 447,
    "adjacent_orgs": adj["id"].nunique() == 308,
    "adjacent_i13_events": int(adj["I13_next"].sum()) == 268,
    "adjacent_i13_nonevents": int((1 - adj["I13_next"]).sum()) == 179,
    "s4_n": int((adj["wave_t"] == 3).sum()) == 89,
    "s4_events": int(adj.loc[adj["wave_t"] == 3, "I13_next"].sum()) == 55,
    "gap2_n": len(g2) == 154,
    "gap2_orgs": g2["id"].nunique() == 116,
    "gap2_events": int(g2["I13_next"].sum()) == 105,
    "elig_n": len(elig) == 686,
    "elig_retained": int(elig["retain_adjacent"].sum()) == 447,
    "adjacent_complete": int(adj[required_adj].isna().sum().sum()) == 0,
    "elig_complete": int(elig[required_elig].isna().sum().sum()) == 0,
    "sector_valid": bool(adj["sector"].between(1, 15).all()),
    "region_valid": bool(adj["region2"].between(1, 6).all()),
}

focal_cells = {}
for rscore in [3, 4]:
    for escore in [1, 3]:
        key = f"R{rscore}_E{escore}"
        focal_cells[key] = int(((adj["R_score"] == rscore) & (adj["E_score"] == escore)).sum())

expected_cells = {"R3_E1": 43, "R3_E3": 26, "R4_E1": 69, "R4_E3": 71}
checks["focal_cells_match"] = focal_cells == expected_cells

if not all(checks.values()):
    failed = [k for k, v in checks.items() if not v]
    raise RuntimeError("input check failed: " + ", ".join(failed))

adj.to_csv(EXEC / "inputs/configuration_adjacent_pairs.csv", index=False)
g2.to_csv(EXEC / "inputs/configuration_two_wave_gap_pairs.csv", index=False)
elig.to_csv(EXEC / "inputs/configuration_retention_eligibility.csv", index=False)

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

manifest = {
    "source_hashes": {
        "raw_csls": sha256(raw_path),
        "panel": sha256(panel_path),
    },
    "counts": {
        "adjacent_pairs": int(len(adj)),
        "adjacent_organizations": int(adj["id"].nunique()),
        "primary_events": int(adj["I13_next"].sum()),
        "primary_nonevents": int((1 - adj["I13_next"]).sum()),
        "strict_wave3_to4_pairs": int((adj["wave_t"] == 3).sum()),
        "two_wave_gap_pairs": int(len(g2)),
        "two_wave_gap_organizations": int(g2["id"].nunique()),
        "retention_eligible_rows": int(len(elig)),
        "retention_observed_adjacent": int(elig["retain_adjacent"].sum()),
    },
    "focal_cells": focal_cells,
    "checks": checks,
}
(EXEC / "inputs/configuration_input_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

print(json.dumps(manifest, indent=2))
