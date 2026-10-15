"""Cache helper (Week 2 Data Conventions, D7-D9).

Every fetch goes through `cached`: read the parquet file if it exists,
otherwise fetch, validate, write. Only force_refresh=True refetches.
"""
from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Callable

import pandas as pd

import config

log = logging.getLogger(__name__)


def root() -> Path:
    """Repository root (the directory holding config.py)."""
    return Path(config.__file__).resolve().parent


def resolve(path: Path) -> Path:
    """CACHE_DIR in config.py is relative to the repository root."""
    return path if path.is_absolute() else root() / path


def cached(path: Path, fetch_fn: Callable[[], pd.DataFrame],
           validate_fn: Callable[[pd.DataFrame], None] | None = None,
           force_refresh: bool = False) -> pd.DataFrame:
    path = resolve(path)
    if path.exists() and not force_refresh:
        log.info("cache hit  %s", path.name)
        return pd.read_parquet(path)
    log.info("cache miss %s: fetching", path.name)
    t0 = time.time()
    df = fetch_fn()
    if validate_fn is not None:
        validate_fn(df)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    log.info("wrote %s (%d rows) in %.0f s", path.name, len(df), time.time() - t0)
    return df
