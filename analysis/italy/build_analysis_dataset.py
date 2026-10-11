from pathlib import Path
import pandas as pd
import numpy as np
import hashlib
import json

src=Path(r"data/raw/istat_ict_2019_public_use/MICRODATI/ICT_Microdati_2019.txt")
out_dir=Path(r"results/italy/build")
out_dir.mkdir(parents=True,exist_ok=True)

d=pd.read_csv(src,sep="\t",dtype=str,keep_default_na=False)
for c in d.columns:
    d[c]=d[c].str.strip()

counts={"source_rows":int(len(d))}

x=d[d["clad3"].isin(["cl1","cl2"])].copy()
counts["sme_rows"]=int(len(x))
assert len(x)==16914

def num(col):
    return pd.to_numeric(x[col].replace({"":np.nan,".":np.nan}),errors="coerce")

for c in ["c5","c7","C8a","C10a","d1","D2a","D2b"]:
    x[c+"_num"]=num(c)

x["C8a_routed"]=x["C8a_num"]
struct=(x["c7_num"]==0) & x["C8a"].isin(["","."])
x.loc[struct,"C8a_routed"]=0

for suffix in list("abcdefghij"):
    c="e1"+suffix
    x[c+"_num"]=num(c)

for c in ["e3a","e3b"]:
    x[c+"_num"]=num(c)

tech4=["c5_num","c7_num","C8a_routed","C10a_num"]
tech6=["c5_num","c7_num","C10a_num","d1_num","D2a_num","D2b_num"]
resp10=["e1"+s+"_num" for s in list("abcdefghij")]
resp6=["e1a_num","e1b_num","e1d_num","e1e_num","e1f_num","e1h_num"]

tech_ok=np.logical_and.reduce([x[c].isin([0,1]).to_numpy() for c in tech4])
resp_ok=np.logical_and.reduce([x[c].isin([0,1]).to_numpy() for c in resp10])
counts["tech4_complete"]=int(tech_ok.sum())
counts["response10_complete"]=int(resp_ok.sum())

control_ok=(
    x["clad3"].isin(["cl1","cl2"]) &
    ~x["rip"].isin(["","."]) &
    ~x["mac"].isin(["","."]) &
    ~x["coeffin"].isin(["","."])
)
analysis=x.loc[tech_ok & resp_ok & control_ok].copy()
counts["primary_analysis_rows"]=int(len(analysis))
assert len(analysis)==16809

analysis["weight"]=pd.to_numeric(analysis["coeffin"],errors="coerce")
assert analysis["weight"].notna().all()
assert (analysis["weight"]>0).all()

analysis["size_medium"]=(analysis["clad3"]=="cl2").astype(int)
analysis["footprint_score"]=analysis[tech4].sum(axis=1).astype(int)
analysis["footprint_breadth"]=analysis["footprint_score"]/4.0
analysis["response_score"]=analysis[resp10].sum(axis=1).astype(int)
analysis["response_breadth"]=analysis["response_score"]/10.0

broad_ok=np.logical_and.reduce([analysis[c].isin([0,1]).to_numpy() for c in tech6])
analysis["digital_system_score_6"]=np.where(broad_ok,analysis[tech6].sum(axis=1),np.nan)
analysis["digital_system_breadth_6"]=analysis["digital_system_score_6"]/6.0

core_ok=np.logical_and.reduce([analysis[c].isin([0,1]).to_numpy() for c in resp6])
analysis["security_core_score_6"]=np.where(core_ok,analysis[resp6].sum(axis=1),np.nan)
analysis["security_core_breadth_6"]=analysis["security_core_score_6"]/6.0

p99=float(analysis["weight"].quantile(0.99))
analysis["high_weight_gt_p99"]=(analysis["weight"]>p99).astype(int)
counts["weight_gt_p99_rows"]=int(analysis["high_weight_gt_p99"].sum())

assert analysis["footprint_score"].between(0,4).all()
assert analysis["response_score"].between(0,10).all()
assert analysis["footprint_breadth"].between(0,1).all()
assert analysis["response_breadth"].between(0,1).all()

keep=[
    "codice_","clad3","size_medium","rip","mac","dom4","ricavi_cl",
    "coeffin","weight",
    "c5","c7","C8a","C10a","c5_num","c7_num","C8a_routed","C10a_num",
    "footprint_score","footprint_breadth",
    "d1_num","D2a_num","D2b_num","digital_system_score_6","digital_system_breadth_6",
    "e1a_num","e1b_num","e1c_num","e1d_num","e1e_num","e1f_num","e1g_num","e1h_num","e1i_num","e1j_num",
    "response_score","response_breadth","security_core_score_6","security_core_breadth_6",
    "e3a_num","e3b_num","high_weight_gt_p99"
]
out=analysis[keep].copy()
csv_path=out_dir/"italy_analysis_input.csv"
out.to_csv(csv_path,index=False)

w=analysis["weight"].astype(float)
manifest={
    "source_file":str(src),
    "source_sha256":hashlib.sha256(src.read_bytes()).hexdigest(),
    "row_counts":counts,
    "primary_controls":["size_medium","rip","mac"],
    "primary_predictor":"footprint_breadth",
    "primary_outcome":"response_breadth",
    "primary_weight":"weight",
    "factor_references":{
        "size_medium":0,
        "rip":sorted(analysis["rip"].unique())[0],
        "mac":sorted(analysis["mac"].unique())[0]
    },
    "weight":{
        "p99_threshold":p99,
        "gt_p99_rows":int((w>p99).sum()),
        "gt_p99_weight_share":float(w[w>p99].sum()/w.sum()),
        "kish_ess":float(w.sum()**2/(w.pow(2).sum())),
        "ess_fraction":float((w.sum()**2/(w.pow(2).sum()))/len(w))
    },
    "sensitivity_support":{
        "digital_system_breadth_6_complete":int(analysis["digital_system_breadth_6"].notna().sum()),
        "security_core_breadth_6_complete":int(analysis["security_core_breadth_6"].notna().sum()),
        "dom4_complete":int((~analysis["dom4"].isin(["","."])).sum())
    },
    "analysis_csv":str(csv_path),
    "analysis_csv_sha256":hashlib.sha256(csv_path.read_bytes()).hexdigest(),
    "focal_association_computed":False
}
(out_dir/"italy_analysis_input_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
print(json.dumps(manifest,indent=2))
