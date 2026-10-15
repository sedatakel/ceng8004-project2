"""Workstream B computations (Week 3 EDA Plan): edge weights and thresholds.

Every function works on one dataset's EDA tables (src.populations.eda_tables).
Plots live in the notebook; these return numbers only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from src.populations import population_ids, restrict

TAUS = (1, 5, 10, 31)
FRACS = (0.009, 0.026)


def ccdf(values: pd.Series) -> pd.DataFrame:
    """Share of items with value >= x, for every distinct x."""
    counts = values.value_counts().sort_index()
    at_least = counts[::-1].cumsum()[::-1]
    return pd.DataFrame({"x": counts.index, "share": at_least.to_numpy() / len(values)})


# B1 -------------------------------------------------------------------------

def b1_counts(edges: pd.DataFrame, edgeroi: pd.DataFrame | None, self_logged: int) -> dict:
    out = {"edges": len(edges),
           "synapses": int(edges["syn_count"].sum()),
           "self_connections_in_table": int((edges["pre_id"] == edges["post_id"]).sum()),
           "self_connections_dropped_logged": self_logged}
    if edgeroi is not None:
        out["edgeroi_rows"] = len(edgeroi)
        out["edgeroi_rows_per_edge"] = len(edgeroi) / len(edges)
        out["edgeroi_synapses"] = int(edgeroi["syn_count"].sum())
    return out


# B2 -------------------------------------------------------------------------

def _describe(label: str, kept: pd.DataFrame, ids: pd.Index, base_edges: int, base_syn: int) -> dict:
    in_deg = kept.groupby("post_id").size().reindex(ids, fill_value=0)
    out_deg = kept.groupby("pre_id").size().reindex(ids, fill_value=0)
    connected = (in_deg + out_deg) > 0
    return {
        "criterion": label,
        "edges": len(kept),
        "edges_pct": 100 * len(kept) / base_edges,
        "synapses": int(kept["syn_count"].sum()),
        "synapses_pct": 100 * kept["syn_count"].sum() / base_syn,
        "neurons_connected": int(connected.sum()),
        "neurons_isolated": int((~connected).sum()),
        "isolated_pct": 100 * (~connected).mean(),
        "median_in_degree": float(in_deg.median()),
        "median_out_degree": float(out_deg.median()),
        "p99_out_degree": float(out_deg.quantile(0.99)),
    }


def b2_threshold_table(neurons: pd.DataFrame, edges: pd.DataFrame, population: str) -> pd.DataFrame:
    """Edges/synapses kept, isolated neurons and degrees for every threshold.

    Degrees count partners inside the population and include zeros (every
    neuron of the population is a row). `frac` uses the fixed denominator."""
    ids = population_ids(neurons, population)
    e = restrict(edges, ids)
    base_edges, base_syn = len(e), int(e["syn_count"].sum())
    rows = [_describe(f">= {t} syn", e[e["syn_count"] >= t], ids, base_edges, base_syn) for t in TAUS]
    rows += [_describe(f">= {f:.1%} of input", e[e["frac"] >= f], ids, base_edges, base_syn) for f in FRACS]
    out = pd.DataFrame(rows)
    out.insert(0, "population", population)
    out.insert(1, "neurons", len(ids))
    return out


# B3 -------------------------------------------------------------------------

def b3_crosstab(edges: pd.DataFrame, tau: int = 10, f: float = 0.009) -> pd.DataFrame:
    return pd.crosstab(edges["syn_count"] >= tau, edges["frac"] >= f,
                       rownames=[f"syn_count >= {tau}"], colnames=[f"frac >= {f:.1%}"])


def b3_disagreement(edges: pd.DataFrame, tau: int = 10, f: float = 0.009) -> dict:
    a, r = edges["syn_count"] >= tau, edges["frac"] >= f
    either = a | r
    return {"abs_only": int((a & ~r).sum()), "rel_only": int((~a & r).sum()),
            "both": int((a & r).sum()), "neither": int((~a & ~r).sum()),
            "disagree_pct_of_all": 100 * (a != r).mean(),
            "disagree_pct_of_passing_either": 100 * (a != r)[either].mean()}


# B4 -------------------------------------------------------------------------

def b4_lorenz(edges: pd.DataFrame, tau: int = 10) -> tuple[pd.DataFrame, dict]:
    w = np.sort(edges["syn_count"].to_numpy())[::-1]
    cum_edges = np.arange(1, len(w) + 1) / len(w)
    cum_syn = np.cumsum(w) / w.sum()
    k = int((w >= tau).sum())
    marker = {"tau": tau, "edges_share": k / len(w), "synapses_share": float(cum_syn[k - 1]) if k else 0.0}
    # thin the curve for plotting: 2,000 points are plenty
    idx = np.unique(np.linspace(0, len(w) - 1, 2000).astype(int))
    return pd.DataFrame({"cum_edges": cum_edges[idx], "cum_synapses": cum_syn[idx]}), marker


# B5, B6 ---------------------------------------------------------------------

def type_edges(neurons: pd.DataFrame, edges: pd.DataFrame, type_col: str = "cell_type",
               side: str | None = None) -> pd.DataFrame:
    """Type-level edges (EDA Plan Definitions): for types A, B on side s, the sum
    of syn_count over edges whose pre neuron is (A, s) and post neuron is (B, s).
    Neurons with side center/unknown or no type are left out.
    Returns columns side, pre_type, post_type, weight."""
    n = neurons.loc[neurons["side"].isin(["left", "right"]) & neurons[type_col].notna(),
                    ["neuron_id", type_col, "side"]]
    if side is not None:
        n = n[n["side"] == side]
    lab = n.set_index("neuron_id")
    e = edges[["pre_id", "post_id", "syn_count"]]
    e = e[e["pre_id"].isin(lab.index) & e["post_id"].isin(lab.index)]
    pre = lab.reindex(e["pre_id"].to_numpy())
    post = lab.reindex(e["post_id"].to_numpy())
    same = pre["side"].to_numpy() == post["side"].to_numpy()
    t = pd.DataFrame({"side": pre["side"].to_numpy()[same].astype(str),
                      "pre_type": pre[type_col].to_numpy()[same],
                      "post_type": post[type_col].to_numpy()[same],
                      "weight": e["syn_count"].to_numpy()[same]})
    return t.groupby(["side", "pre_type", "post_type"], observed=True)["weight"].sum().reset_index()


def _pair_agreement(a: pd.Series, b: pd.Series) -> tuple[pd.DataFrame, dict]:
    """a, b: weights indexed by (pre_type, post_type). Pairs present on at least one side."""
    both = pd.concat([a.rename("a"), b.rename("b")], axis=1).fillna(0)
    shared = both[(both["a"] > 0) & (both["b"] > 0)]
    la, lb = np.log10(shared["a"]), np.log10(shared["b"])
    ratio = np.maximum(shared["a"], shared["b"]) / np.minimum(shared["a"], shared["b"])
    stats = {"pairs_either": len(both), "pairs_both": len(shared),
             "only_a": int((both["b"] == 0).sum()), "only_b": int((both["a"] == 0).sum()),
             "pearson_r_log10": float(np.corrcoef(la, lb)[0, 1]),
             "share_diff_gt_30pct": float((ratio > 1.3).mean()),
             "median_log2_b_over_a": float(np.median(np.log2(shared["b"] / shared["a"])))}
    return both, stats


def b5_left_right(neurons: pd.DataFrame, edges: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    t = type_edges(neurons, edges).set_index(["pre_type", "post_type"])
    return _pair_agreement(t.loc[t["side"] == "left", "weight"], t.loc[t["side"] == "right", "weight"])


def matched_types(fafb: pd.DataFrame, malecns: pd.DataFrame) -> set:
    """fafb cell types that appear verbatim in malecns.flywire_type."""
    return set(fafb["cell_type"].dropna()) & set(malecns["flywire_type"].dropna())


def b6_across(fafb_n: pd.DataFrame, fafb_e: pd.DataFrame,
              mcns_n: pd.DataFrame, mcns_e: pd.DataFrame, side: str = "right") -> tuple[pd.DataFrame, dict]:
    """Matched types, one side: fafb weights (a) vs malecns weights (b)."""
    m = matched_types(fafb_n, mcns_n)
    f = fafb_n.assign(mtype=fafb_n["cell_type"].where(fafb_n["cell_type"].isin(m)))
    c = mcns_n.assign(mtype=mcns_n["flywire_type"].where(mcns_n["flywire_type"].isin(m)))
    a = type_edges(f, fafb_e, "mtype", side).set_index(["pre_type", "post_type"])["weight"]
    b = type_edges(c, mcns_e, "mtype", side).set_index(["pre_type", "post_type"])["weight"]
    both, stats = _pair_agreement(a, b)
    stats["matched_types"] = len(m)
    return both, stats


# B7 -------------------------------------------------------------------------

def b7_raw_features(neurons: pd.DataFrame, edges: pd.DataFrame, edgeroi: pd.DataFrame,
                    population: str = "central", tau: int = 10) -> pd.DataFrame:
    """Raw-count feature table: basic_features at tau plus input and output
    synapses per neuropil from edgeroi (all cached rows), one row per neuron."""
    from eda_shared.eda_helpers import basic_features
    ids = population_ids(neurons, population)
    feats = basic_features(neurons, edges, tau=tau).set_index("neuron_id").loc[ids]
    r = edgeroi[["pre_id", "post_id", "neuropil", "syn_count"]]
    out_np = (r[r["pre_id"].isin(ids)].groupby(["pre_id", "neuropil"], observed=True)["syn_count"].sum()
              .unstack(fill_value=0).add_prefix("out_"))
    in_np = (r[r["post_id"].isin(ids)].groupby(["post_id", "neuropil"], observed=True)["syn_count"].sum()
             .unstack(fill_value=0).add_prefix("in_"))
    X = feats.join(in_np).join(out_np).fillna(0)
    X = X.loc[:, X.sum() > 0]                         # neuropils nobody in the population uses
    return X


def b7_pca(X: pd.DataFrame, k: int = 5, normalise_rows: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Standardise, PCA, and Pearson r of PC1..k with log10(total synapses).

    normalise_rows=True is the Project 2 contrast (not part of the EDA item):
    each neuron's per-neuropil columns divided by its own total, degree and
    total columns dropped, so magnitude is removed before PCA."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    size = np.log10(X["total_in_syn"] + X["total_out_syn"])
    keep = X.index[size.notna() & np.isfinite(size)]
    M = X.loc[keep]
    if normalise_rows:
        np_cols = [c for c in M.columns if c.startswith(("in_", "out_"))
                   and c not in ("in_degree", "out_degree")]
        M = M[np_cols].div(M[np_cols].sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    Z = StandardScaler().fit_transform(M.to_numpy(float))
    pca = PCA(n_components=k, random_state=0).fit(Z)
    scores = pd.DataFrame(pca.transform(Z), index=keep, columns=[f"PC{i + 1}" for i in range(k)])
    s = size.loc[keep]
    # Pearson r is the EDA item; Spearman rho is reported beside it because a few
    # very large neurons dominate the standardised raw counts and the relation is not linear.
    r = pd.DataFrame({"component": scores.columns,
                      "r": [float(np.corrcoef(scores[c], s)[0, 1]) for c in scores.columns],
                      "spearman_rho": [float(scores[c].corr(s, method="spearman")) for c in scores.columns],
                      "explained_variance": pca.explained_variance_ratio_})
    info = {"neurons": len(keep), "columns": M.shape[1]}
    return r, scores.assign(log10_total_syn=s), info
