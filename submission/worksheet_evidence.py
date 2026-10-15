import sys, json, logging, os, warnings
sys.path.insert(0, "."); warnings.filterwarnings("ignore"); logging.disable(logging.INFO)
import numpy as np, pandas as pd
import pyarrow.feather as pf
import config
from src import fetch_fafb, fetch_malecns, eda_b
from src.populations import flag
from eda_shared.eda_helpers import load_neuropil_table, ALIASES, basic_features
OUT = sys.argv[1]; E = {}
T = load_neuropil_table()

fn, fe, fr = fetch_fafb.load()
mn, me = fetch_malecns.load_without_edgeroi()

# ---- Week 2 printouts
E["source_columns"] = {
  "fafb_annotations_v3.2.0": list(pd.read_csv("data/raw/fafb/Supplemental_file1_neuron_annotations_v3.2.0.tsv", sep="\t", nrows=1).columns),
  "malecns_body_annotations": pf.read_table("data/raw/malecns/body-annotations-male-cns-v1.0-minconf-0.5.feather").schema.names,
}
E["value_counts"] = {d: {c: n[c].value_counts(dropna=False).astype(int).to_dict() for c in ["side","super_class","nt_type"]}
                     for d, n in [("fafb", fn), ("malecns", mn)]}
E["null_counts"] = {d: n.isna().sum().astype(int).to_dict() for d, n in [("fafb", fn), ("malecns", mn)]}
E["unknown_counts"] = {d: {c: int((n[c].astype(str)=="unknown").sum()) for c in ["side","super_class","nt_type","region"]} for d, n in [("fafb", fn), ("malecns", mn)]}
E["edges"] = {d: {"tau1": len(e), "tau5": int((e.syn_count>=5).sum()), "syn_tau1": int(e.syn_count.sum())} for d, e in [("fafb", fe), ("malecns", me)]}
E["cache_files_MB"] = {p.name: round(p.stat().st_size/1e6,1) for p in sorted((config.resolve(config.CACHE_DIR) if hasattr(config,'resolve') else __import__('pathlib').Path('data/cache')).glob('*.parquet'))}
E["edgeroi_rows_fafb"] = len(fr)

# ---- D2
mx = int(fn.neuron_id.max()); E["D2"] = {"max_id": mx, "as_float64": int(float(mx)), "diff": int(float(mx)) - mx,
  "collisions_fafb": int(len(fn) - fn.neuron_id.astype("float64").nunique())}
# ---- D4
E["D4"] = {d: n.nt_confidence.astype(float).describe(percentiles=[.1,.25,.5,.75]).round(3).to_dict() for d,n in [("fafb",fn),("malecns",mn)]}

# ---- per-neuron per-neuropil synapses
def base(s):
    s = s.astype(str).str.replace(r"(_[LR]|\([LR]\))$", "", regex=True); return s.replace(ALIASES)
f_np = pd.concat([fr.groupby(["pre_id","neuropil"]).syn_count.sum().rename_axis(["nid","np"]),
                  fr.groupby(["post_id","neuropil"]).syn_count.sum().rename_axis(["nid","np"])]).groupby(level=[0,1]).sum().reset_index(name="syn")
roi = pd.read_parquet("data/raw/malecns/neuprint_traced_roiinfo_v1.0.parquet").dropna(subset=["roi"])
m_np = roi.assign(syn=roi.upstream+roi.downstream).rename(columns={"bodyId":"nid","roi":"np"})[["nid","np","syn"]]
m_np = m_np[m_np.syn>0]
def optic_of(np_):
    b = base(np_); o = b.map(T["is_optic"]).astype("boolean")
    o[b=="Optic-unspecified"] = True
    return o.fillna(False).astype(bool)
mf = flag(mn); mbrain = set(mf.loc[mf.in_brain | mf.region.isin(["CentralBrain-unspecified","Optic-unspecified(L)","Optic-unspecified(R)"]), "neuron_id"])
ff = flag(fn)
E["D13"] = {}; E["D15"] = {}
for d, P, ids in [("fafb", f_np, set(fn.neuron_id)), ("malecns_brain", m_np, mbrain)]:
    P = P[P.nid.isin(ids)].copy(); P["opt"] = optic_of(P.np)
    tot = P.groupby("nid").syn.sum(); osh = P[P.opt].groupby("nid").syn.sum().reindex(tot.index, fill_value=0)/tot
    E["D13"][d] = {"neurons_with_syn": len(tot), "optic_share_ge_0.5": int((osh>=.5).sum()),
                   "ambiguous_20_80": int(((osh>.2)&(osh<.8)).sum()), "ambiguous_10_90": int(((osh>.1)&(osh<.9)).sum())}
    s = P.sort_values(["nid","syn"], ascending=[True,False]); s["rk"] = s.groupby("nid").cumcount()
    top = s[s.rk==0].set_index("nid").syn; sec = s[s.rk==1].set_index("nid").syn.reindex(top.index, fill_value=0)
    E["D15"][d] = {"neurons": len(top), "exact_tie": int((sec==top).sum()), "second_ge_90pct": int((sec>=.9*top).sum()),
                   "second_ge_80pct": int((sec>=.8*top).sum()), "top_share_lt_50pct": int((top/P.groupby("nid").syn.sum().reindex(top.index)<.5).sum())}
E["D13"]["flag_share"] = {"fafb": float(ff.is_optic.mean()), "malecns_brain": float(mf[mf.in_brain].is_optic.mean())}
E["D13"]["fafb_superclass_optic_share"] = float(fn.super_class.isin(["optic","visual_projection","visual_centrifugal"]).mean())

# ---- D16 at tau=10
E["D16"] = {}
for d, n, e, pop in [("fafb", ff, fe, None), ("malecns_brain", mf[mf.in_brain], me, None)]:
    ids = set(n.neuron_id); ee = e[e.pre_id.isin(ids)&e.post_id.isin(ids)]
    bf = basic_features(n, ee, tau=10).merge(n[["neuron_id","is_optic"]], on="neuron_id")
    for p, sub in [("all", bf), ("central", bf[~bf.is_optic])]:
        E["D16"][f"{d}_{p}"] = {"neurons": len(sub), "zero_in": int((sub.in_degree==0).sum()), "zero_out": int((sub.out_degree==0).sum()),
                                "zero_both": int(((sub.in_degree==0)&(sub.out_degree==0)).sum()),
                                "zero_total_in_syn_tau1": int((sub.total_in_syn==0).sum()), "zero_total_out_syn_tau1": int((sub.total_out_syn==0).sum())}

# ---- D18: pair by cell_type, left vs right
E["D18"] = {}
for d, n in [("fafb", ff), ("fafb_central", ff[~ff.is_optic])]:
    typed = n[n.cell_type.notna()]
    lr = typed[typed.side.isin(["left","right"])].groupby(["cell_type","side"], observed=True).size().unstack(fill_value=0)
    lr = lr.reindex(columns=["left","right"], fill_value=0)
    E["D18"][d] = {"neurons": len(n), "untyped": int(n.cell_type.isna().sum()),
                   "typed_side_center_or_unknown": int((~typed.side.isin(["left","right"])).sum()),
                   "paired": int(2*np.minimum(lr.left, lr.right).sum()), "excess_unequal": int((lr.left-lr.right).abs().sum()),
                   "types": len(lr), "types_one_side_only": int(((lr.left==0)|(lr.right==0)).sum()), "types_unequal": int((lr.left!=lr.right).sum())}

# ---- D14: B7 on all vs central (fafb)
E["D14"] = {}
for p in ["central", "all"]:
    X = eda_b.b7_raw_features(ff, fe, fr, p, tau=10); r, _, info = eda_b.b7_pca(X)
    E["D14"][p] = {"info": info, "PC1": r.iloc[0].to_dict(), "max_abs_spearman_PC1_5": float(r.spearman_rho.abs().max())}

json.dump(E, open(OUT, "w"), indent=1, default=str); print("ok")
