from pathlib import Path
import pandas as pd

d1 = pd.read_csv(r"results/uk/decomposition/D1_later_window_no_sector.csv")
s3 = pd.read_csv(r"results/uk/sensitivity/S3_sector_adjusted_key_results.csv")
primary = pd.read_csv(r"results/uk/primary/M1_interaction_coefficient.csv")

rows=[]
mapping=[("M0","D1_M0","S3_M0"),("M1","D1_M1","S3_M1"),("M2","D1_M2","S3_M2")]
for label,d1name,s3name in mapping:
    a=d1[d1.model==d1name].iloc[0]
    b=s3[s3.model==s3name].iloc[0]
    rows.append({
        "model":label,
        "later_window_no_sector_estimate":a.estimate,
        "later_window_no_sector_p":a.p_value,
        "later_window_with_sector_estimate":b.estimate,
        "later_window_with_sector_p":b.p_value,
        "sector_increment_estimate_change":b.estimate-a.estimate,
        "absolute_magnitude_ratio_with_sector_vs_no_sector":abs(b.estimate)/abs(a.estimate) if a.estimate!=0 else None
    })
out=pd.DataFrame(rows)

p_est=float(primary.estimate.iloc[0])
d1_m1=float(d1.loc[d1.model=="D1_M1","estimate"].iloc[0])
s3_m1=float(s3.loc[s3.model=="S3_M1","estimate"].iloc[0])

summary=pd.DataFrame([{
    "primary_m1_estimate":p_est,
    "later_window_no_sector_m1":d1_m1,
    "later_window_with_sector_m1":s3_m1,
    "primary_to_later_no_sector_magnitude_reduction_fraction":1-abs(d1_m1)/abs(p_est),
    "later_no_sector_to_sector_magnitude_reduction_fraction":1-abs(s3_m1)/abs(d1_m1),
    "sector_increment_as_fraction_of_primary_magnitude":abs(s3_m1-d1_m1)/abs(p_est)
}])

out.to_csv(r"results/uk/decomposition/D2_s3_decomposition_table.csv",index=False)
summary.to_csv(r"results/uk/decomposition/D2_m1_attenuation_summary.csv",index=False)
print(out.to_string(index=False))
print(summary.to_string(index=False))
