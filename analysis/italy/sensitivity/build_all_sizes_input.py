import pandas as pd, numpy as np
from pathlib import Path
src=Path(r"data/raw/istat_ict_2019_public_use/MICRODATI/ICT_Microdati_2019.txt")
out=Path(r"results/italy/models/italy_S8_all_sizes_input.csv")
d=pd.read_csv(src,sep="\t",dtype=str,keep_default_na=False)
for c in d.columns: d[c]=d[c].str.strip()
x=d[d.clad3.isin(["cl1","cl2","cl3"])].copy()
def num(c): return pd.to_numeric(x[c].replace({"":np.nan,".":np.nan}),errors="coerce")
for c in ["c5","c7","C8a","C10a"]: x[c+"_num"]=num(c)
x.loc[(x["c7_num"]==0)&x["C8a"].isin(["","."]),"C8a_num"]=0
for s in list("abcdefghij"): x["e1"+s+"_num"]=num("e1"+s)
tech=["c5_num","c7_num","C8a_num","C10a_num"]
resp=["e1"+s+"_num" for s in list("abcdefghij")]
ok=np.logical_and.reduce([x[c].isin([0,1]).to_numpy() for c in tech+resp])
ok=ok & ~x.rip.isin(["","."]) & ~x.mac.isin(["","."]) & ~x.coeffin.isin(["","."])
x=x.loc[ok].copy()
assert len(x)==19804
x["weight"]=pd.to_numeric(x.coeffin,errors="coerce")
x["footprint_score"]=x[tech].sum(axis=1).astype(int)
x["footprint_breadth"]=x.footprint_score/4
x["response_score"]=x[resp].sum(axis=1).astype(int)
x["response_breadth"]=x.response_score/10
x[["codice_","clad3","rip","mac","weight","footprint_score","footprint_breadth","response_score","response_breadth"]].to_csv(out,index=False)
print("IT_S8_INPUT_PASS",len(x),x.clad3.value_counts().to_dict())
