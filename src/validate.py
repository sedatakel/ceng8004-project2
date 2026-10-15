"""The eight checks of the Week 2 Data Conventions, read from config.py.

Each check raises ValueError naming the offending column.
"""
from __future__ import annotations

import pandas as pd

import config


def _dtype(spec) -> str:
    return str(pd.api.types.pandas_dtype(spec))


def check_columns(df: pd.DataFrame, schema: dict, name: str) -> None:
    missing = [c for c in schema if c not in df.columns]
    extra = [c for c in df.columns if c not in schema]
    if missing or extra:
        raise ValueError(f"{name}: missing columns {missing}, unexpected columns {extra}")
    if list(df.columns) != list(schema):
        raise ValueError(f"{name}: column order differs from config.py")


def check_dtypes(df: pd.DataFrame, schema: dict, name: str) -> None:
    for col, (spec, _) in schema.items():
        want = pd.api.types.pandas_dtype(spec)
        got = df[col].dtype
        if isinstance(want, pd.CategoricalDtype):
            if not isinstance(got, pd.CategoricalDtype) or list(got.categories) != list(want.categories):
                raise ValueError(f"{name}.{col}: dtype {got} differs from the fixed CategoricalDtype")
        elif str(got) != str(want):
            raise ValueError(f"{name}.{col}: dtype {got}, expected {want}")


def check_nulls(df: pd.DataFrame, schema: dict, name: str) -> None:
    for col, (_, nullable) in schema.items():
        if not nullable and df[col].isna().any():
            raise ValueError(f"{name}.{col}: {int(df[col].isna().sum())} nulls in a non-nullable column")


def check_keys(df: pd.DataFrame, keys: list[str], name: str) -> None:
    dup = df.duplicated(keys).sum()
    if dup:
        raise ValueError(f"{name}: {dup} repeated keys on {keys}")


def check_references(edges: pd.DataFrame, neurons: pd.DataFrame, name: str) -> None:
    ids = pd.Index(neurons["neuron_id"].unique())
    for col in ("pre_id", "post_id"):
        bad = ~edges[col].isin(ids)
        if bad.any():
            raise ValueError(f"{name}.{col}: {int(bad.sum())} ids have no neuron row in the same dataset")


def check_vocabulary(df: pd.DataFrame, schema: dict, name: str) -> None:
    # Categorical columns cannot hold values outside their categories; a value that
    # fell outside becomes NaN, which check_nulls catches. Region has no list yet.
    if config.REGION_VALUES is not None and "region" in df:
        bad = ~df["region"].isin(config.REGION_VALUES)
        if bad.any():
            raise ValueError(f"{name}.region: values outside REGION_VALUES {df.loc[bad, 'region'].unique()[:10]}")


def check_threshold(edges: pd.DataFrame, name: str) -> None:
    bad = edges["syn_count"] < edges["min_syn_fetched"]
    if bad.any():
        raise ValueError(f"{name}.syn_count: {int(bad.sum())} rows below min_syn_fetched")


def check_plausibility(edges: pd.DataFrame, neurons: pd.DataFrame, name: str) -> None:
    """Input synapses in the edge table must not exceed the source's n_post."""
    n_post = neurons.set_index("neuron_id")["n_post"].dropna()
    if n_post.empty:
        return
    total_in = edges.groupby("post_id")["syn_count"].sum()
    joined = total_in.to_frame("total_in").join(n_post, how="inner")
    bad = joined["total_in"] > joined["n_post"]
    if bad.any():
        raise ValueError(f"{name}.n_post: {int(bad.sum())} neurons receive more synapses in edges than n_post")


def validate_neurons(neurons: pd.DataFrame) -> None:
    s = config.NEURON_SCHEMA
    check_columns(neurons, s, "neurons")
    check_dtypes(neurons, s, "neurons")
    check_nulls(neurons, s, "neurons")
    check_keys(neurons, ["neuron_id", "dataset"], "neurons")
    check_vocabulary(neurons, s, "neurons")


def validate_edges(edges: pd.DataFrame, schema: dict | None = None, name: str = "edges") -> None:
    s = schema or config.EDGE_SCHEMA
    keys = ["pre_id", "post_id", "dataset"] + (["neuropil"] if "neuropil" in s else [])
    check_columns(edges, s, name)
    check_dtypes(edges, s, name)
    check_nulls(edges, s, name)
    check_keys(edges, keys, name)
    check_threshold(edges, name)
    if (edges["pre_id"] == edges["post_id"]).any():
        raise ValueError(f"{name}.pre_id: self-connections present")


def validate_edgeroi(edgeroi: pd.DataFrame) -> None:
    validate_edges(edgeroi, config.EDGEROI_SCHEMA, "edgeroi")


def validate_dataset(neurons: pd.DataFrame, edges: pd.DataFrame, edgeroi: pd.DataFrame) -> None:
    """All eight checks on one dataset's three tables."""
    validate_neurons(neurons)
    validate_edges(edges)
    validate_edgeroi(edgeroi)
    check_references(edges, neurons, "edges")
    check_references(edgeroi, neurons, "edgeroi")
    check_plausibility(edges, neurons, "edges")


def assert_interchangeable(a: pd.DataFrame, b: pd.DataFrame) -> None:
    assert list(a.columns) == list(b.columns), "column sets differ"
    assert a.dtypes.equals(b.dtypes), "dtypes differ"
