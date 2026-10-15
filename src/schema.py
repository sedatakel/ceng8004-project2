"""Recode source columns into the agreed vocabularies and dtypes (config.py)."""
from __future__ import annotations

import logging

import pandas as pd

log = logging.getLogger(__name__)


def recode(series: pd.Series, mapping: dict, column: str) -> pd.Series:
    """Data Conventions recode: an unmapped value raises, missing becomes unknown."""
    out = series.map(mapping)
    unmapped = series[out.isna() & series.notna()].unique()
    if len(unmapped):
        raise ValueError(f"{column}: unmapped values {list(unmapped)[:10]}")
    return out.fillna("unknown")


def to_schema(df: pd.DataFrame, schema: dict) -> pd.DataFrame:
    """Select the schema's columns in order and cast each to its dtype."""
    out = pd.DataFrame(index=df.index)
    for col, (spec, _) in schema.items():
        out[col] = df[col].astype(spec)
    return out.reset_index(drop=True)


def strip_side(neuropil: pd.Series) -> pd.Series:
    return neuropil.str.replace(r"(_[LR]|\([LR]\))$", "", regex=True)


def top_region(syn_by_neuropil: pd.DataFrame) -> pd.Series:
    """Region rule (D15, provisional): the neuropil holding most of a neuron's
    synapses, inputs plus outputs, over all cached edges (tau = 1).
    Ties go to the alphabetically first neuropil. The tie rule is a week 3 decision.

    syn_by_neuropil: columns neuron_id, neuropil, syn.
    """
    s = syn_by_neuropil.sort_values(["neuron_id", "syn", "neuropil"],
                                    ascending=[True, False, True])
    top = s.drop_duplicates("neuron_id")
    return top.set_index("neuron_id")["neuropil"]


def log_counts(name: str, series: pd.Series) -> None:
    log.info("%s:\n%s", name, series.value_counts(dropna=False).to_string())
