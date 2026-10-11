from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

source = Path(r"data/raw/UKDA-9356-stata/stata/stata13/csls_waves_1-4.dta")
out = Path(r"results/panel/csls_raw_selected_columns.csv")
manifest = Path(r"results/panel/io_bridge_manifest.json")

cols = [
    "wave","qtypdum","LongId","onlinea","onlineb","devices",
    "rulesb","rulesf","rulesh","rulesi","qsize","modetype"
]

df = pd.read_stata(source, columns=cols, convert_categoricals=False)
df.insert(0, "source_row", np.arange(1, len(df) + 1, dtype=np.int64))
df.to_csv(out, index=False)

back = pd.read_csv(out, dtype={"LongId":"string"}, keep_default_na=False)
assert len(back) == len(df)
assert list(back.columns) == list(df.columns)
assert np.array_equal(back["source_row"].to_numpy(dtype=np.int64), df["source_row"].to_numpy(dtype=np.int64))

for c in [x for x in cols if x != "LongId"]:
    a = pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float)
    b = pd.to_numeric(back[c], errors="coerce").to_numpy(dtype=float)
    assert np.allclose(a, b, equal_nan=True, rtol=0.0, atol=0.0), c

a = df["LongId"].astype("string").fillna("").to_numpy()
b = back["LongId"].astype("string").fillna("").to_numpy()
assert np.array_equal(a, b)

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()

obj = {
    "purpose": "I/O bridge only, no focal model and no filtering",
    "source": str(source),
    "source_sha256": sha256(source),
    "output": str(out),
    "output_sha256": sha256(out),
    "rows": int(len(df)),
    "columns": list(df.columns),
    "roundtrip_equal": True,
    "pandas_version": pd.__version__
}
manifest.write_text(json.dumps(obj, indent=2), encoding="utf-8")
print(json.dumps(obj, indent=2))
