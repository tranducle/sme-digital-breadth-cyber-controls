from pathlib import Path
import pandas as pd
import numpy as np
import json
from itertools import combinations

UK=Path(r"data/raw/UKDA-8971-stata/stata/stata13/dcms_csbs_combined_2016-2022_banded_cost_data_only_v3_public.dta")
IT=Path(r"data/raw/istat_ict_2019_public_use/MICRODATI/ICT_Microdati_2019.txt")
OUT=Path(r"results/measurement")
OUT.mkdir(parents=True, exist_ok=True)

def wmean(x,w):
    x=pd.Series(x); w=pd.Series(w)
    m=x.notna() & w.notna() & (w>0)
    return float((x[m].astype(float)*w[m].astype(float)).sum()/w[m].sum())

def wcorr(x,y,w):
    x=pd.Series(x); y=pd.Series(y); w=pd.Series(w)
    m=x.notna() & y.notna() & w.notna() & (w>0)
    x=x[m].astype(float).to_numpy(); y=y[m].astype(float).to_numpy(); ww=w[m].astype(float).to_numpy()
    sw=ww.sum(); mx=(ww*x).sum()/sw; my=(ww*y).sum()/sw
    vx=(ww*(x-mx)**2).sum()/sw; vy=(ww*(y-my)**2).sum()/sw
    if vx<=0 or vy<=0: return np.nan
    return float((ww*(x-mx)*(y-my)).sum()/sw/np.sqrt(vx*vy))

# Build UK 2019 private SME eligible sample
u=pd.read_stata(UK, columns=[
    "year","questtype","sizeb","region_comb","sector_comb2","weight",
    "online2","online3","online4","rules1","rules4","rules13","rules14"
], convert_categoricals=False)
u=u[(u.year==4)&(u.questtype==1)&u.sizeb.isin([2,3])].copy()
u["website"]=u.online2
u["social"]=u.online3
u["online_tx"]=u.online4
u["updates"]=u.rules1
u["access"]=u.rules4
u["backup"]=np.where((u.rules13==1)|(u.rules14==1),1,np.where((u.rules13==0)&(u.rules14==0),0,np.nan))
u["weight_use"]=u.weight
uk_raw=len(u)

# Build Italy SME eligible sample
i=pd.read_csv(IT,sep="\t",dtype=str,keep_default_na=False)
for c in i.columns: i[c]=i[c].str.strip()
i=i[i.clad3.isin(["cl1","cl2"])].copy()
def num(df,c): return pd.to_numeric(df[c].replace({"":np.nan,".":np.nan}),errors="coerce")
i["website"]=num(i,"c7")
i["social"]=num(i,"C10a")
i["online_tx"]=num(i,"C8a")
i.loc[(i.website==0)&i.C8a.isin(["","."]),"online_tx"]=0
i["updates"]=num(i,"e1b")
i["access"]=num(i,"e1f")
i["backup"]=num(i,"e1e")
i["weight_use"]=num(i,"coeffin")
it_raw=len(i)

tech=["website","social","online_tx"]
resp=["updates","access","backup"]

def diagnostics(country,df,geo,sector):
    rows=[]
    for c in tech+resp:
        rows.append({
            "country":country,
            "item":c,
            "raw_n":int(len(df)),
            "valid_binary_n":int(df[c].isin([0,1]).sum()),
            "missing_or_invalid_n":int((~df[c].isin([0,1])).sum()),
            "unweighted_prevalence":float(df.loc[df[c].isin([0,1]),c].mean()) if df[c].isin([0,1]).any() else np.nan,
            "weighted_prevalence":wmean(df.loc[df[c].isin([0,1]),c],df.loc[df[c].isin([0,1]),"weight_use"]) if df[c].isin([0,1]).any() else np.nan
        })
    valid=np.logical_and.reduce([df[c].isin([0,1]).to_numpy() for c in tech+resp])
    controls=df[geo].notna() & df[sector].notna() & (df.weight_use>0)
    a=df.loc[valid & controls].copy()
    a["tech_score"]=a[tech].sum(axis=1).astype(int)
    a["resp_score"]=a[resp].sum(axis=1).astype(int)

    score=[]
    for construct,score_col,maxs in [("Technology","tech_score",3),("Response3","resp_score",3)]:
        total=a.weight_use.sum()
        for s in range(maxs+1):
            m=a[score_col]==s
            score.append({
                "country":country,"construct":construct,"score":s,
                "n":int(m.sum()),
                "unweighted_share":float(m.mean()),
                "weighted_share":float(a.loc[m,"weight_use"].sum()/total)
            })

    patterns=[]
    total=a.weight_use.sum()
    for vals,g in a.groupby(resp,dropna=False):
        if not isinstance(vals,tuple): vals=(vals,)
        mask=np.ones(len(a),dtype=bool)
        for c,v in zip(resp,vals): mask &= (a[c].to_numpy()==v)
        patterns.append({
            "country":country,
            "updates":int(vals[0]),"access":int(vals[1]),"backup":int(vals[2]),
            "pattern":f"{int(vals[0])}{int(vals[1])}{int(vals[2])}",
            "n":int(mask.sum()),
            "unweighted_share":float(mask.mean()),
            "weighted_share":float(a.loc[mask,"weight_use"].sum()/total)
        })

    corr=[]
    for construct,items in [("Technology",tech),("Response3",resp)]:
        for x,y in combinations(items,2):
            corr.append({
                "country":country,"construct":construct,
                "var1":x,"var2":y,
                "weighted_phi":wcorr(a[x],a[y],a.weight_use)
            })

    return a,pd.DataFrame(rows),pd.DataFrame(score),pd.DataFrame(patterns),pd.DataFrame(corr)

ua,ukitems,ukscores,ukpatterns,ukcorr=diagnostics("UK",u,"region_comb","sector_comb2")
ia,ititems,itscores,itpatterns,itcorr=diagnostics("Italy",i,"rip","mac")

item_audit=pd.DataFrame([
{"country":"UK","construct":"Technology","domain":"Website presence","source_item":"online2","exact_wording":"A website or blog","eligible_population":"All 2019 CSBS respondents; analysis restricted to private-sector SMEs with 10-249 employees","routing":"No item-specific skip","coding":"Yes=1, No=0","missing_rule":"Only valid binary values retained","weight":"CSBS released weight","reference_period":"Current use at interview","residual_mismatch":"Italy asks enterprise website, UK allows website or blog"},
{"country":"Italy","construct":"Technology","domain":"Website presence","source_item":"c7","exact_wording":"Indicare se l'impresa ha un sito Web, ovvero una o piu pagine su Internet","eligible_population":"Istat 2019 enterprises; analysis restricted to SMEs with 10-249 employees","routing":"If No, questionnaire skips C8 and proceeds to C10","coding":"Si=1, No=0","missing_rule":"Valid binary values retained","weight":"coeffin","reference_period":"Current status","residual_mismatch":"None material at domain level"},
{"country":"UK","construct":"Technology","domain":"Social-media presence/use","source_item":"online3","exact_wording":"Accounts or pages on social media sites (e.g. Facebook or Twitter)","eligible_population":"All 2019 CSBS respondents; private-sector SME analysis","routing":"No item-specific skip","coding":"Yes=1, No=0","missing_rule":"Only valid binary values retained","weight":"CSBS released weight","reference_period":"Current use at interview","residual_mismatch":"UK asks accounts/pages; Italy asks use of social networks"},
{"country":"Italy","construct":"Technology","domain":"Social-media presence/use","source_item":"C10a","exact_wording":"Social network (ad es. Facebook, LinkedIn, MySpace, Google+, Xing, Viadeo, Yammer)","eligible_population":"All Istat 2019 enterprises reaching C10; SME analysis","routing":"C10 follows website section and is asked regardless of website status","coding":"Si=1, No=0","missing_rule":"Valid binary values retained","weight":"coeffin","reference_period":"Current use","residual_mismatch":"Italy asks active use, not merely account/page presence"},
{"country":"UK","construct":"Technology","domain":"Online transaction capability","source_item":"online4","exact_wording":"The ability for customers to order, book or pay for products or services online","eligible_population":"All 2019 CSBS respondents; private-sector SME analysis","routing":"No website prerequisite in Q6 battery","coding":"Yes=1, No=0","missing_rule":"Only valid binary values retained","weight":"CSBS released weight","reference_period":"Current capability","residual_mismatch":"UK includes online payment; Italy wording is ordering/booking via website"},
{"country":"Italy","construct":"Technology","domain":"Online transaction capability","source_item":"C8a","exact_wording":"Possibilita di effettuare ordinazioni o prenotazioni online (es. carrello della spesa online)","eligible_population":"Enterprises with C7 website=Yes; SME analysis","routing":"C7=No skips C8; structural C8a missing is recoded 0 only when C7=0","coding":"Si=1, No=0; structural zero rule when no website","missing_rule":"Non-structural missing remains missing","weight":"coeffin","reference_period":"Current website capability","residual_mismatch":"Italy is website-conditional and excludes stand-alone external-platform capability"},
{"country":"UK","construct":"Response","domain":"Software updating","source_item":"rules1","exact_wording":"Applying software updates when they are available","eligible_population":"All 2019 CSBS respondents; private-sector SME analysis","routing":"No item-specific skip","coding":"Yes=1, No=0","missing_rule":"Only valid binary values retained","weight":"CSBS released weight","reference_period":"Controls currently in place","residual_mismatch":"Italy explicitly includes operating systems"},
{"country":"Italy","construct":"Response","domain":"Software updating","source_item":"e1b","exact_wording":"Mantenere aggiornato il software (inclusi i sistemi operativi)","eligible_population":"Istat 2019 enterprises; SME analysis","routing":"E1 battery, no item-specific skip identified","coding":"Si=1, No=0","missing_rule":"Valid binary values retained","weight":"coeffin","reference_period":"Measures currently used","residual_mismatch":"None material at domain level"},
{"country":"UK","construct":"Response","domain":"Access control","source_item":"rules4","exact_wording":"Restricting IT admin and access rights to specific users","eligible_population":"All 2019 CSBS respondents; private-sector SME analysis","routing":"No item-specific skip","coding":"Yes=1, No=0","missing_rule":"Only valid binary values retained","weight":"CSBS released weight","reference_period":"Controls currently in place","residual_mismatch":"Focuses privilege/user-right restriction"},
{"country":"Italy","construct":"Response","domain":"Access control","source_item":"e1f","exact_wording":"Controllo dell'accesso alla rete (gestione dell'accesso da parte di dispositivi e utenti alla rete aziendale)","eligible_population":"Istat 2019 enterprises; SME analysis","routing":"E1 battery, no item-specific skip identified","coding":"Si=1, No=0","missing_rule":"Valid binary values retained","weight":"coeffin","reference_period":"Measures currently used","residual_mismatch":"Includes device and user network access, broader than UK privilege restriction"},
{"country":"UK","construct":"Response","domain":"Secure backup","source_item":"rules13 OR rules14","exact_wording":"Backing up data securely via a cloud service OR backing up data securely via other means","eligible_population":"All 2019 CSBS respondents; private-sector SME analysis","routing":"Paired response options in Q31","coding":"Derived 1 if either item=1; 0 if both valid 0","missing_rule":"Missing unless presence established or both components valid zero","weight":"CSBS released weight","reference_period":"Controls currently in place","residual_mismatch":"Two implementation routes combined into one domain"},
{"country":"Italy","construct":"Response","domain":"Secure backup","source_item":"e1e","exact_wording":"Backup dei dati in una posizione separata dagli originali (incluso il backup nel cloud)","eligible_population":"Istat 2019 enterprises; SME analysis","routing":"E1 battery, no item-specific skip identified","coding":"Si=1, No=0","missing_rule":"Valid binary values retained","weight":"coeffin","reference_period":"Measures currently used","residual_mismatch":"Requires separate location; UK wording emphasizes secure backup route"}
])
item_audit.to_csv(OUT/"measurement_item_level_alignment_audit.csv",index=False)
pd.concat([ukitems,ititems],ignore_index=True).to_csv(OUT/"measurement_item_missingness_prevalence.csv",index=False)
pd.concat([ukscores,itscores],ignore_index=True).to_csv(OUT/"measurement_score_distributions.csv",index=False)
pd.concat([ukpatterns,itpatterns],ignore_index=True).sort_values(["country","pattern"]).to_csv(OUT/"measurement_response3_joint_patterns.csv",index=False)
pd.concat([ukcorr,itcorr],ignore_index=True).to_csv(OUT/"measurement_pairwise_phi.csv",index=False)

summary={
 "uk":{"raw_n":uk_raw,"aligned_analysis_n":int(len(ua)),"response3_max_weighted_score_share":float(ukscores[ukscores.construct=="Response3"].weighted_share.max()),"response3_top_score":3},
 "italy":{"raw_n":it_raw,"aligned_analysis_n":int(len(ia)),"response3_max_weighted_score_share":float(itscores[itscores.construct=="Response3"].weighted_share.max()),"response3_top_score":3},
 "response3_inferential_status":"BLOCKED_DESCRIPTIVE_ONLY",
 "common_technology_status":"PARTIAL_ITEM_CONTENT_CORRESPONDENCE_ONLY",
 "common_response_model_authorized":False,
 "focal_association_computed":False
}
(OUT/"measurement_alignment_diagnostics.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))
print("\nRESPONSE JOINT PATTERNS")
print(pd.concat([ukpatterns,itpatterns],ignore_index=True).sort_values(["country","weighted_share"],ascending=[True,False]).to_string(index=False))
print("\nPAIRWISE PHI")
print(pd.concat([ukcorr,itcorr],ignore_index=True).to_string(index=False))
