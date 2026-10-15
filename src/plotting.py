"""Shared figure style and file naming for the EDA figures.

File names carry item, dataset, release, threshold and population
(EDA Plan rule), e.g. B2_fafb-malecns_783-v1.0_grid_all-central-optic.png.
"""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

from src.cache import root

COLORS = {"fafb": "#2a78d6", "malecns": "#eb6834", "third": "#1baf7a", "muted": "#8a8984"}
LINESTYLES = {"all": "-", "central": "--", "optic": ":"}
RELEASES = {"fafb": "783", "malecns": "v1.0"}
LABELS = {"fafb": "FAFB 783 (female)", "malecns": "MaleCNS v1.0, brain (male)"}

FIG_DIR = root() / "figures" / "eda"
TABLE_DIR = root() / "reports" / "eda_B_tables"


def style() -> None:
    mpl.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 160, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
        "xtick.color": "#52514e", "ytick.color": "#52514e",
        "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6,
        "lines.linewidth": 2, "legend.frameon": False,
    })


def name(item: str, datasets: list[str], threshold: str, populations: list[str]) -> str:
    rel = "-".join(RELEASES[d] for d in datasets)
    return f"{item}_{'-'.join(datasets)}_{rel}_{threshold}_{'-'.join(populations)}"


def save(fig: plt.Figure, item: str, datasets: list[str], threshold: str,
         populations: list[str], caption: str) -> str:
    """Add the caption (dataset, release, threshold, population) and save."""
    fig.text(0.0, -0.02, caption, ha="left", va="top", fontsize=7.5, color="#52514e", wrap=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fname = name(item, datasets, threshold, populations) + ".png"
    fig.savefig(FIG_DIR / fname)
    return fname


def save_table(df, fname: str) -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(TABLE_DIR / fname, index=False)
