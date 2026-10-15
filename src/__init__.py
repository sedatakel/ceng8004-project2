import pandas as pd

# Arrow-backed strings: tens of millions of neuropil / release values fit in 8 GB RAM.
# str(dtype) is still "string", so the config.py schema checks are unaffected.
pd.set_option("mode.string_storage", "pyarrow")
