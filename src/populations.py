"""Datasets and populations as the Week 3 EDA Plan defines them.

fafb      all neurons
malecns   its brain part: in_brain == True (add_provisional_flags with the
          nerve-cord neuropils passed as nerve_cord) and the edges among them
Populations: all, central (is_optic == False), optic (is_optic == True).
An edge belongs to a population when both of its neurons do.
"""
from __future__ import annotations

import pandas as pd

from eda_shared.eda_helpers import add_provisional_flags, relative_weight

# MaleCNS nerve-cord neuropils and nerves (neuPrint primary ROIs outside the brain)
# that occur as a neuron's region. Passed to add_provisional_flags as nerve_cord.
MALECNS_NERVE_CORD = frozenset({
    "ANm", "IntTct", "LTct", "VNC-unspecified",
    "LegNp(T1)(L)", "LegNp(T1)(R)", "LegNp(T2)(L)", "LegNp(T2)(R)", "LegNp(T3)(L)", "LegNp(T3)(R)",
    "WTct(UTct-T2)(L)", "WTct(UTct-T2)(R)", "HTct(UTct-T3)(L)", "HTct(UTct-T3)(R)",
    "NTct(UTct-T1)(L)", "NTct(UTct-T1)(R)", "Ov(L)", "Ov(R)",
    "mVAC(T1)(L)", "mVAC(T1)(R)", "mVAC(T2)(L)", "mVAC(T2)(R)", "mVAC(T3)(L)", "mVAC(T3)(R)",
    "DMetaN(L)", "DMetaN(R)", "AbN1(L)", "AbN1(R)", "AbN2(L)", "AbN2(R)", "AbN3(L)", "AbN3(R)",
    "AbN4(L)", "AbN4(R)", "AbNT(L)", "AbNT(R)", "DProN(L)", "DProN(R)",
    "MesoLN(L)", "MesoLN(R)", "MetaLN(L)", "MetaLN(R)", "ProLN(L)", "ProLN(R)",
    "MesoAN(L)", "MesoAN(R)", "ProAN(L)", "ProAN(R)", "MetaAN(L)", "MetaAN(R)",
    "PrN(L)", "PrN(R)", "CV-unspecified",
})

# Brain regions add_provisional_flags rejects; reported to the group, NOT added to
# eda_shared. Until the group decides, these neurons are treated as region unknown
# for the EDA flags only (so is_optic = False and in_brain = False).
MALECNS_UNRESOLVED = frozenset({"CentralBrain-unspecified",
                                "Optic-unspecified(L)", "Optic-unspecified(R)"})


def flag(neurons: pd.DataFrame) -> pd.DataFrame:
    """add_provisional_flags for either dataset."""
    if (neurons["dataset"] == "malecns").any():
        n = neurons.copy()
        n["region"] = n["region"].mask(n["region"].isin(MALECNS_UNRESOLVED), "unknown")
        out = add_provisional_flags(n, nerve_cord=MALECNS_NERVE_CORD)
        out["region"] = neurons["region"]          # keep the source value in region
        return out
    return add_provisional_flags(neurons)


def eda_tables(neurons: pd.DataFrame, edges: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Flagged neurons and edges for one dataset as the EDA uses it.

    malecns is cut to its brain part. Edges gain `frac`, the relative weight with
    the fixed denominator total_in_syn(post) taken from all of the dataset's
    edges at tau = 1, computed before any filtering.
    """
    n = flag(neurons)
    e = edges
    if (n["dataset"] == "malecns").any():
        n = n[n["in_brain"]].reset_index(drop=True)
        keep = pd.Index(n["neuron_id"])
        e = e[e["pre_id"].isin(keep) & e["post_id"].isin(keep)]
    e = e.reset_index(drop=True)
    e["frac"] = relative_weight(e)
    return n, e


def population_ids(neurons: pd.DataFrame, population: str) -> pd.Index:
    if population == "all":
        m = pd.Series(True, index=neurons.index)
    elif population == "central":
        m = ~neurons["is_optic"]
    elif population == "optic":
        m = neurons["is_optic"]
    else:
        raise ValueError(f"population: {population}")
    return pd.Index(neurons.loc[m, "neuron_id"])


def restrict(edges: pd.DataFrame, ids: pd.Index) -> pd.DataFrame:
    """Edges whose two neurons are both in `ids`."""
    return edges[edges["pre_id"].isin(ids) & edges["post_id"].isin(ids)]
