from pathlib import Path
import numpy as np
import pandas as pd

src=Path(r"data/raw/UKDA-8971-stata/stata/stata13/dcms_csbs_combined_2016-2022_banded_cost_data_only_v3_public.dta")
out=Path(r"results/uk/sensitivity/S1_unweighted_input.csv")

cols=["uniqser","IMID2","year","weight","questtype","sizeb","region_comb","sector_comb2",
      "online4","online5","online9","mobile","rules2","rules4","rules5","rules7","manage1","manage2"]

d=pd.read_stata(src,columns=cols,convert_categoricals=False)
x=d[(d.questtype==1)&d.sizeb.isin([2,3])&d.year.isin([2,3,4,5,6])&d.region_comb.notna()].copy()
x["mobile_bin"]=np.where(x.mobile==1,1,np.where(x.mobile==2,0,np.nan))
foot=["online4","online5","online9","mobile_bin"]
resp=["rules2","rules4","rules5","rules7"]
ok=np.logical_and.reduce([x[c].isin([0,1]).to_numpy() for c in foot+resp])
x=x.loc[ok].copy()
assert len(x)==2860, f"Expected 2860 rows, got {len(x)}"
assert x.manage1.isin([0,1]).all()
assert x.manage2.isin([0,1]).all()
x["year_calendar"]=x.year.map({2:2017,3:2018,4:2019,5:2020,6:2021})
x["size_medium"]=(x.sizeb==3).astype(int)
x["footprint_score"]=x[foot].sum(axis=1).astype(int)
x["footprint_breadth"]=x.footprint_score/4.0
x["response_score"]=x[resp].sum(axis=1).astype(int)
x["response_breadth"]=x.response_score/4.0
x["one_weight"]=1.0
x.to_csv(out,index=False)
print("S1_INPUT_PASS")
print("rows",len(x))
print("raw_zero_weight_rows",int((x.weight<=0).sum()))
