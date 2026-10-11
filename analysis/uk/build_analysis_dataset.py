from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

src = Path(r"data/raw/UKDA-8971-stata/stata/stata13/dcms_csbs_combined_2016-2022_banded_cost_data_only_v3_public.dta")
out_dir = Path(r"results/uk/build")
out_dir.mkdir(parents=True, exist_ok=True)

cols = [
    "uniqser","IMID2","year","weight","questtype","sizeb","region_comb","sector_comb2",
    "online4","online5","online9","mobile",
    "rules2","rules4","rules5","rules7",
    "manage1","manage2"
]

df = pd.read_stata(src, columns=cols, convert_categoricals=False)

counts = {"source_rows": int(len(df))}

# Population and years, using verified raw codes.
x = df[(df["questtype"] == 1) & (df["sizeb"].isin([2, 3])) & (df["year"].isin([2, 3, 4, 5, 6]))].copy()
counts["population_year_rows"] = int(len(x))
assert len(x) == 2915, f"population_year_rows expected 2915, got {len(x)}"

x = x[x["region_comb"].notna()].copy()
counts["region_complete_rows"] = int(len(x))
assert len(x) == 2890, f"region_complete_rows expected 2890, got {len(x)}"

x = x[x["weight"] > 0].copy()
counts["positive_weight_rows"] = int(len(x))
assert len(x) == 2889, f"positive_weight_rows expected 2889, got {len(x)}"

# Normalization occurs before item complete-case filtering.
x["w_norm"] = x.groupby("year")["weight"].transform(lambda s: s / s.mean())

# Verified year labels.
year_map = {2:2017, 3:2018, 4:2019, 5:2020, 6:2021}
x["year_calendar"] = x["year"].map(year_map)
x["size_medium"] = (x["sizeb"] == 3).astype(int)

# Mobile recode.
x["mobile_bin"] = np.where(x["mobile"] == 1, 1, np.where(x["mobile"] == 2, 0, np.nan))

foot_items = ["online4","online5","online9","mobile_bin"]
resp_items = ["rules2","rules4","rules5","rules7"]

def valid_binary(s):
    return s.isin([0, 1])

foot_valid = np.logical_and.reduce([valid_binary(x[c]).to_numpy() for c in foot_items])
resp_valid = np.logical_and.reduce([valid_binary(x[c]).to_numpy() for c in resp_items])

x = x.loc[foot_valid & resp_valid].copy()
counts["complete_construct_rows"] = int(len(x))
assert len(x) == 2859, f"complete_construct_rows expected 2859, got {len(x)}"

# Moderator validity is expected to be complete in the primary period.
assert x["manage1"].isin([0, 1]).all(), "manage1 contains invalid primary-period values"
assert x["manage2"].isin([0, 1]).all(), "manage2 contains invalid primary-period values"

x["footprint_score"] = x[foot_items].sum(axis=1).astype(int)
x["footprint_breadth"] = x["footprint_score"] / 4.0
x["response_score"] = x[resp_items].sum(axis=1).astype(int)
x["response_breadth"] = x["response_score"] / 4.0

# Structural assertions only, no association calculations.
assert x["footprint_score"].between(0, 4).all()
assert x["response_score"].between(0, 4).all()
assert x["footprint_breadth"].between(0, 1).all()
assert x["response_breadth"].between(0, 1).all()
assert np.isfinite(x["w_norm"]).all()
assert (x["w_norm"] > 0).all()

keep = [
    "uniqser","IMID2","year","year_calendar","weight","w_norm",
    "sizeb","size_medium","region_comb","sector_comb2",
    "online4","online5","online9","mobile","mobile_bin",
    "footprint_score","footprint_breadth",
    "rules2","rules4","rules5","rules7",
    "response_score","response_breadth",
    "manage1","manage2"
]

analysis = x[keep].copy()
csv_path = out_dir / "uk_analysis_input.csv"
analysis.to_csv(csv_path, index=False)

# Record normalization means at the stage where normalization was defined.
base = df[(df["questtype"] == 1) & (df["sizeb"].isin([2, 3])) & (df["year"].isin([2, 3, 4, 5, 6]))].copy()
base = base[base["region_comb"].notna() & (base["weight"] > 0)].copy()
base["w_norm"] = base.groupby("year")["weight"].transform(lambda s: s / s.mean())
weight_means = {
    str(year_map[int(y)]): {
        "n": int(len(g)),
        "raw_weight_mean": float(g["weight"].mean()),
        "normalized_weight_mean": float(g["w_norm"].mean())
    }
    for y, g in base.groupby("year")
}

manifest = {
    "source_file": str(src),
    "source_sha256": hashlib.sha256(src.read_bytes()).hexdigest(),
    "row_counts": counts,
    "weight_normalization_stage": "after population, year, region, and positive-weight filters; before item complete-case",
    "weight_means_by_wave": weight_means,
    "analysis_csv": str(csv_path),
    "analysis_csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
    "analysis_rows": int(len(analysis)),
    "analysis_columns": list(analysis.columns),
    "focal_models_executed": False,
    "cross_construct_association_computed": False
}

manifest_path = out_dir / "uk_analysis_dataset_manifest.json"
manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

print(json.dumps(manifest, indent=2))
