"""FAFB (FlyWire, materialisation 783) tables to the Week 2 Data Conventions.

Sources (D6), all public, no token needed:
  neurons  flywire_annotations v3.2.0, Supplemental_file1_neuron_annotations.tsv
  edges    Codex per-neuropil connection table, Buhmann synapses, no threshold
           (gs://flywire-data/codex/data/fafb/783/connections_buhmann_no_threshold.csv.gz)

Raw downloads live in data/raw/fafb/ (gitignored); the cache in data/cache/.
"""
from __future__ import annotations

import logging
import urllib.request

import duckdb
import pandas as pd

import config
from src.cache import cached, root
from src.schema import log_counts, recode, to_schema, top_region
from src.validate import validate_edgeroi, validate_edges, validate_neurons

log = logging.getLogger(__name__)

DATASET = "fafb"
RELEASE = str(config.FAFB_MATERIALIZATION)
ANN_SOURCE = "ann" + config.FAFB_ANNOTATIONS_TAG.lstrip("v")      # ann3.2.0
SYN_SOURCE = config.FAFB_SYNAPSE_TABLE                            # buhmann

RAW_DIR = root() / "data" / "raw" / "fafb"
ANN_URL = ("https://raw.githubusercontent.com/flyconnectome/flywire_annotations/"
           f"{config.FAFB_ANNOTATIONS_TAG}/supplemental_files/Supplemental_file1_neuron_annotations.tsv")
ANN_FILE = RAW_DIR / f"Supplemental_file1_neuron_annotations_{config.FAFB_ANNOTATIONS_TAG}.tsv"
CODEX_URL = "https://storage.googleapis.com/flywire-data/codex/data/fafb/{release}/{name}"
CONN_NAME = {"buhmann": "connections_buhmann_no_threshold.csv.gz",
             "princeton": "connections_princeton_no_threshold.csv.gz"}[SYN_SOURCE]
CONN_FILE = RAW_DIR / f"codex_{CONN_NAME}"

# D3: FAFB side values -> agreed vocabulary
SIDE_MAP = {"left": "left", "right": "right", "center": "center", "na": "unknown"}
# D4: FAFB top_nt is already spelled out
NT_MAP = {v: v for v in config.NT_VALUES}
# super_class: FAFB values are the reference list. sensory_ascending is not in
# config.SUPER_CLASS_VALUES; mapped to sensory (flow = afferent) PROVISIONALLY,
# reported to the group.
SUPER_CLASS_MAP = {v: v for v in config.SUPER_CLASS_VALUES}
SUPER_CLASS_MAP["sensory_ascending"] = "sensory"


def _download(url: str, path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    log.info("downloading %s", url)
    tmp = path.with_suffix(path.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.rename(path)


def _raw_files() -> None:
    _download(ANN_URL, ANN_FILE)
    _download(CODEX_URL.format(release=RELEASE, name=CONN_NAME), CONN_FILE)


def _con() -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute("SET memory_limit='3GB'")
    return con


def _roi_sql() -> str:
    """Per (pre, post, neuropil) synapse counts; the Codex rows already are that."""
    return f"""SELECT pre_root_id AS pre_id, post_root_id AS post_id, neuropil,
                      SUM(syn_count) AS syn_count
               FROM read_csv_auto('{CONN_FILE}')
               GROUP BY 1, 2, 3"""


def fetch_edgeroi() -> pd.DataFrame:
    _raw_files()
    con = _con()
    df = con.sql(f"SELECT * FROM ({_roi_sql()}) WHERE pre_id <> post_id").df()
    n_self = con.sql(f"SELECT COUNT(*) FROM ({_roi_sql()}) WHERE pre_id = post_id").fetchone()[0]
    log.info("fafb edgeroi: dropped %d self-connection rows", n_self)
    df["dataset"], df["release"], df["min_syn_fetched"] = DATASET, RELEASE, config.FETCH_MIN_SYN
    return to_schema(df, config.EDGEROI_SCHEMA)


def fetch_edges() -> pd.DataFrame:
    """D8: sum over neuropils first, then threshold; drop self-connections and log."""
    _raw_files()
    con = _con()
    agg = f"""SELECT pre_id, post_id, SUM(syn_count) AS syn_count
              FROM ({_roi_sql()}) GROUP BY 1, 2"""
    n_self = con.sql(f"SELECT COUNT(*) FROM ({agg}) WHERE pre_id = post_id").fetchone()[0]
    log.info("fafb edges: dropped %d self-connections", n_self)
    df = con.sql(f"""SELECT * FROM ({agg})
                     WHERE pre_id <> post_id AND syn_count >= {config.FETCH_MIN_SYN}""").df()
    df["dataset"], df["release"], df["min_syn_fetched"] = DATASET, RELEASE, config.FETCH_MIN_SYN
    return to_schema(df, config.EDGE_SCHEMA)


def _regions() -> pd.Series:
    con = _con()
    per = con.sql(f"""SELECT neuron_id, neuropil, SUM(syn) AS syn FROM (
                          SELECT pre_id AS neuron_id, neuropil, syn_count AS syn FROM ({_roi_sql()})
                          UNION ALL
                          SELECT post_id, neuropil, syn_count FROM ({_roi_sql()}))
                      GROUP BY 1, 2""").df()
    return top_region(per)


def fetch_neurons() -> pd.DataFrame:
    _raw_files()
    a = pd.read_csv(ANN_FILE, sep="\t", low_memory=False, dtype={"root_id": "int64"})

    # Ids in the 783 edge table but not in the annotation file get unknown
    # annotations (Step 0). Expected 7; found 14 (and 7 annotation ids that are
    # not 783 ids). Reported to the group.
    con = _con()
    edge_ids = con.sql(f"""SELECT DISTINCT id FROM (
                               SELECT pre_root_id AS id FROM read_csv_auto('{CONN_FILE}')
                               UNION SELECT post_root_id FROM read_csv_auto('{CONN_FILE}'))""").df()["id"]
    extra = pd.Index(edge_ids).difference(pd.Index(a["root_id"]))
    log.info("fafb neurons: %d edge ids missing from annotations, added as unknown", len(extra))
    a = pd.concat([a, pd.DataFrame({"root_id": extra.astype("int64")})], ignore_index=True)

    df = pd.DataFrame({
        "neuron_id":      a["root_id"].astype("int64"),
        "dataset":        DATASET,
        "release":        RELEASE,
        "cell_type":      a["cell_type"],
        "cell_class":     a["cell_class"],
        "cell_sub_class": a["cell_sub_class"],
        "super_class":    recode(a["super_class"], SUPER_CLASS_MAP, "super_class"),
        "side":           recode(a["side"], SIDE_MAP, "side"),
        "nt_type":        recode(a["top_nt"], NT_MAP, "nt_type"),
        "nt_confidence":  a["top_nt_conf"],
        "n_pre":          pd.NA,                       # not provided by FAFB (D1)
        "n_post":         pd.NA,
        "flywire_type":   pd.NA,                       # MaleCNS-only columns
        "manc_type":      pd.NA,
        "hemibrain_type": a["hemibrain_type"],
        "dimorphism":     a["dimorphism"],
        "status":         a["status"],
        "hemilineage":    a["ito_lee_hemilineage"],
        "nt_known":       a["known_nt"],
    })
    df.loc[df["nt_type"] == "unknown", "nt_confidence"] = pd.NA
    df["region"] = df["neuron_id"].map(_regions()).fillna("unknown")
    log_counts("fafb super_class (source)", a["super_class"])
    log_counts("fafb side", df["side"])
    return to_schema(df, config.NEURON_SCHEMA)


def load(force_refresh: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """neurons, edges, edgeroi for FAFB, from the cache when present."""
    m = config.FETCH_MIN_SYN
    neurons = cached(config.cache_path(DATASET, "neurons", RELEASE, ANN_SOURCE),
                     fetch_neurons, validate_neurons, force_refresh)
    edges = cached(config.cache_path(DATASET, "edges", RELEASE, SYN_SOURCE, min_syn=m),
                   fetch_edges, validate_edges, force_refresh)
    edgeroi = cached(config.cache_path(DATASET, "edgeroi", RELEASE, SYN_SOURCE, min_syn=m),
                     fetch_edgeroi, validate_edgeroi, force_refresh)
    return neurons, edges, edgeroi
