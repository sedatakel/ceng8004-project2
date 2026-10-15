"""MaleCNS (Janelia male CNS v1.0) tables to the Week 2 Data Conventions.

Sources (D6): Traced neurons only, synapse confidence >= 0.5.
  neurons  bulk body-annotations + body-neurotransmitters (flat-connectome, v1.0)
           + neuPrint male-cns:v1.0 for n_pre, n_post and roiInfo (needs the token)
  edges    bulk connectome-weights-...-minconf-0.5-traced-only.feather
  edgeroi  bulk syn-partners-...-minconf-0.5-traced-only.feather, per primary_post

Raw downloads live in data/raw/malecns/ (gitignored); the cache in data/cache/.
"""
from __future__ import annotations

import json
import logging
import os
import urllib.request

import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa

import config
from src.cache import cached, root
from src.schema import log_counts, recode, to_schema, top_region
from src.validate import validate_edgeroi, validate_edges, validate_neurons

log = logging.getLogger(__name__)

DATASET = "malecns"
RELEASE = config.MALECNS_DATASET.split(":")[1]                   # v1.0
SOURCE = config.MALECNS_STATUS.lower()                             # traced
CONF = config.MALECNS_SYN_CONFIDENCE

RAW_DIR = root() / "data" / "raw" / "malecns"
BUCKET = f"https://storage.googleapis.com/flyem-male-cns/{RELEASE}/connectome-data/flat-connectome"
F_ANN = f"body-annotations-male-cns-{RELEASE}-minconf-{CONF}.feather"
F_NT = f"body-neurotransmitters-male-cns-{RELEASE}.feather"
F_WEIGHTS = f"connectome-weights-male-cns-{RELEASE}-minconf-{CONF}-traced-only.feather"
F_PARTNERS = f"syn-partners-male-cns-{RELEASE}-minconf-{CONF}-traced-only.feather"
F_ROIINFO = f"neuprint_traced_roiinfo_{RELEASE}.parquet"            # our own extract

# D3: somaSide, else rootSide
SIDE_MAP = {"L": "left", "R": "right", "M": "center", "unknown": "unknown"}
# D4: predictedNt; "unclear" -> unknown
NT_MAP = {**{v: v for v in config.NT_VALUES}, "unclear": "unknown"}
# super_class: MaleCNS values mapped onto FAFB's. PROVISIONAL, reported to the
# group: MaleCNS-only levels are not yet in config.SUPER_CLASS_VALUES, so they
# map to unknown until the group adds them ("after EDA" in config.py).
SUPER_CLASS_MAP = {
    "ol_intrinsic": "optic", "cb_intrinsic": "central",
    "visual_projection": "visual_projection", "visual_centrifugal": "visual_centrifugal",
    "cb_sensory": "sensory", "ol_sensory": "sensory", "vnc_sensory": "sensory",
    "sensory_ascending": "sensory",
    "ascending_neuron": "ascending", "descending_neuron": "descending",
    "cb_motor": "motor", "vnc_motor": "motor",
    "cb_endocrine": "endocrine", "vnc_endocrine": "endocrine",
    # MaleCNS-only levels
    "vnc_intrinsic": "unknown", "vnc_efferent": "unknown", "cb_efferent": "unknown",
    "efferent_ascending": "unknown", "efferent_descending": "unknown",
    "sensory_descending": "unknown", "ENS": "unknown",
    # "_tbc" (to be confirmed) levels
    "vnc_tbc": "unknown", "vnc_sensory_tbc": "unknown", "cb_sensory_tbc": "unknown",
    "visual_projection_tbc": "unknown", "sensory_ascending_tbc": "unknown",
}


def _download(name: str) -> None:
    path = RAW_DIR / name
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    log.info("downloading %s", name)
    tmp = path.with_suffix(path.suffix + ".part")
    urllib.request.urlretrieve(f"{BUCKET}/{name}", tmp)
    tmp.rename(path)


def _feather_reader(name: str, columns: list[str]) -> pa.RecordBatchReader:
    """Stream a feather file batch by batch, so DuckDB never holds it whole."""
    f = pa.ipc.open_file(pa.memory_map(str(RAW_DIR / name)))
    schema = pa.schema([f.schema.field(c) for c in columns])
    batches = (f.get_batch(i).select(columns) for i in range(f.num_record_batches))
    return pa.RecordBatchReader.from_batches(schema, batches)


def _con() -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute("SET memory_limit='3GB'")
    con.execute(f"SET temp_directory='{root() / 'data' / 'raw' / 'duckdb_tmp'}'")
    return con


def fetch_edges() -> pd.DataFrame:
    """The bulk weights are already summed per pair (D8); drop self-connections."""
    _download(F_WEIGHTS)
    con = _con()
    w = _feather_reader(F_WEIGHTS, ["body_pre", "body_post", "weight"])
    df = con.sql("""SELECT body_pre AS pre_id, body_post AS post_id, SUM(weight) AS syn_count
                    FROM w GROUP BY 1, 2""").df()
    self_conn = df["pre_id"] == df["post_id"]
    log.info("malecns edges: dropped %d self-connections", int(self_conn.sum()))
    df = df[~self_conn & (df["syn_count"] >= config.FETCH_MIN_SYN)]
    df["dataset"], df["release"], df["min_syn_fetched"] = DATASET, RELEASE, config.FETCH_MIN_SYN
    return to_schema(df, config.EDGE_SCHEMA)


def fetch_edgeroi() -> pd.DataFrame:
    """One row per (pre, post, neuropil); neuropil is the synapse's primary_post ROI.
    Synapses outside every primary ROI (null or '<unspecified>') get neuropil 'unknown'."""
    _download(F_PARTNERS)
    con = _con()
    p = _feather_reader(F_PARTNERS, ["body_pre", "body_post", "primary_post"])
    df = con.sql("""SELECT body_pre AS pre_id, body_post AS post_id,
                           CASE WHEN primary_post IS NULL OR CAST(primary_post AS VARCHAR) = '<unspecified>'
                                THEN 'unknown' ELSE CAST(primary_post AS VARCHAR) END AS neuropil,
                           COUNT(*) AS syn_count
                    FROM p WHERE body_pre <> body_post GROUP BY 1, 2, 3""").df()
    df["dataset"], df["release"], df["min_syn_fetched"] = DATASET, RELEASE, config.FETCH_MIN_SYN
    return to_schema(df, config.EDGEROI_SCHEMA)


def _neuprint_client():
    from dotenv import load_dotenv
    from neuprint import Client
    load_dotenv(root() / ".env")
    return Client("neuprint.janelia.org", dataset=config.MALECNS_DATASET,
                  token=os.environ["NEUPRINT_APPLICATION_CREDENTIALS"])


def fetch_roiinfo(batch: int = 5000) -> pd.DataFrame:
    """n_pre, n_post and per-primary-ROI pre/post/upstream/downstream for every Traced neuron.
    Long format: bodyId, pre, post (neuron totals, repeated), roi, roi_pre, roi_post,
    upstream, downstream."""
    path = RAW_DIR / F_ROIINFO
    if path.exists():
        return pd.read_parquet(path)
    c = _neuprint_client()
    primary = set(c.primary_rois)
    ids = c.fetch_custom(f"MATCH (n:Neuron) WHERE n.status = '{config.MALECNS_STATUS}' "
                         "RETURN n.bodyId AS bodyId")["bodyId"].to_numpy()
    rows = []
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch].tolist()
        q = c.fetch_custom(f"""MATCH (n:Neuron) WHERE n.bodyId IN {chunk}
                               RETURN n.bodyId AS bodyId, n.pre AS pre, n.post AS post,
                                      n.roiInfo AS roiInfo""")
        for b, pre, post, info in q.itertuples(index=False):
            rois = {r: v for r, v in json.loads(info or "{}").items() if r in primary}
            if not rois:
                rows.append((b, pre, post, None, 0, 0, 0, 0))
            for r, v in rois.items():
                rows.append((b, pre, post, r, v.get("pre", 0), v.get("post", 0),
                             v.get("upstream", 0), v.get("downstream", 0)))
        log.info("neuPrint roiInfo: %d / %d neurons", min(i + batch, len(ids)), len(ids))
    df = pd.DataFrame(rows, columns=["bodyId", "pre", "post", "roi", "roi_pre", "roi_post",
                                     "upstream", "downstream"])
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    return df


def fetch_neurons() -> pd.DataFrame:
    _download(F_ANN)
    _download(F_NT)
    a = pd.read_feather(RAW_DIR / F_ANN)
    a = a[a["status"] == config.MALECNS_STATUS].copy()
    nt = pd.read_feather(RAW_DIR / F_NT, columns=["body", "predicted_nt", "predicted_nt_confidence"])
    a = a.merge(nt, how="left", left_on="bodyId", right_on="body")

    side_src = a["somaSide"].where(a["somaSide"].notna(), a["rootSide"])
    df = pd.DataFrame({
        "neuron_id":      a["bodyId"].astype("int64"),
        "dataset":        DATASET,
        "release":        RELEASE,
        "cell_type":      a["type"],
        "cell_class":     a["class"],
        "cell_sub_class": a["subclass"],
        "super_class":    recode(a["superclass"], SUPER_CLASS_MAP, "super_class"),
        "side":           recode(side_src, SIDE_MAP, "side"),
        "nt_type":        recode(a["predicted_nt"], NT_MAP, "nt_type"),
        "nt_confidence":  a["predicted_nt_confidence"],
        "flywire_type":   a["flywireType"],
        "manc_type":      a["mancType"],
        "hemibrain_type": a["hemibrainType"],
        "dimorphism":     a["dimorphism"],
        "status":         a["status"],
        "hemilineage":    a["itoleeHl"],
        "nt_known":       pd.NA,                       # FAFB-only column
    })
    df.loc[df["nt_type"] == "unknown", "nt_confidence"] = pd.NA

    roi = fetch_roiinfo()
    counts = roi.drop_duplicates("bodyId").set_index("bodyId")[["pre", "post"]]
    df["n_pre"] = df["neuron_id"].map(counts["pre"])
    df["n_post"] = df["neuron_id"].map(counts["post"])
    # D15 for MaleCNS: connections per ROI (upstream + downstream), the unit FAFB's
    # edgeroi counts (Guide §1: each sender-receiver contact is one synapse).
    per = roi.dropna(subset=["roi"]).assign(syn=lambda d: d["upstream"] + d["downstream"])
    per = per[per["syn"] > 0].rename(columns={"bodyId": "neuron_id", "roi": "neuropil"})
    df["region"] = df["neuron_id"].map(top_region(per)).fillna("unknown")

    log_counts("malecns superclass (source)", a["superclass"])
    log_counts("malecns side", df["side"])
    return to_schema(df, config.NEURON_SCHEMA)


def load(force_refresh: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """neurons, edges, edgeroi for MaleCNS, from the cache when present."""
    m = config.FETCH_MIN_SYN
    neurons = cached(config.cache_path(DATASET, "neurons", RELEASE, SOURCE),
                     fetch_neurons, validate_neurons, force_refresh)
    edges = cached(config.cache_path(DATASET, "edges", RELEASE, SOURCE, min_syn=m, conf=CONF),
                   fetch_edges, validate_edges, force_refresh)
    edgeroi = cached(config.cache_path(DATASET, "edgeroi", RELEASE, SOURCE, min_syn=m, conf=CONF),
                     fetch_edgeroi, validate_edgeroi, force_refresh)
    return neurons, edges, edgeroi


def load_without_edgeroi(force_refresh: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    m = config.FETCH_MIN_SYN
    neurons = cached(config.cache_path(DATASET, "neurons", RELEASE, SOURCE),
                     fetch_neurons, validate_neurons, force_refresh)
    edges = cached(config.cache_path(DATASET, "edges", RELEASE, SOURCE, min_syn=m, conf=CONF),
                   fetch_edges, validate_edges, force_refresh)
    return neurons, edges
