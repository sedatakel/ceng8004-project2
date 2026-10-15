# CENG8004 — Project 2: Modules without magnitude

How much of the reported modular structure of the *Drosophila* connectome survives when
connection magnitude is removed from the similarity measure?

EDA workstream owned by this repository: **B — Edge weights and thresholds** (Week 3 EDA Plan).

## Pinned releases (D6)

| Dataset | Release | Neurons | Edges | Filter |
|---|---|---|---|---|
| `fafb` | materialisation 783, annotations v3.2.0, Buhmann synapses | flywire_annotations v3.2.0 | Codex `connections_buhmann_no_threshold.csv.gz` (per neuropil) | none |
| `malecns` | `male-cns:v1.0` | bulk `body-annotations` + `body-neurotransmitters`; neuPrint for `n_pre`, `n_post`, `roiInfo` | bulk `connectome-weights-…-minconf-0.5-traced-only`; edgeroi from `syn-partners-…-traced-only` | Traced, synapse confidence ≥ 0.5 |

All values live in `config.py` (copied unchanged from the Week 2 Data Conventions).

Bulk download addresses:

- FAFB annotations: `https://raw.githubusercontent.com/flyconnectome/flywire_annotations/v3.2.0/supplemental_files/Supplemental_file1_neuron_annotations.tsv`
- FAFB connections: `https://storage.googleapis.com/flywire-data/codex/data/fafb/783/connections_buhmann_no_threshold.csv.gz`
  (the public bucket behind codex.flywire.ai downloads)
- MaleCNS: `https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome/<file>`

## Setup

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh     # uv; Homebrew builds Rust from source on Intel macOS
uv sync                                             # Python 3.12 + locked packages (.venv)
cp .env.example .env                                # then paste CAVE_TOKEN and NEUPRINT_APPLICATION_CREDENTIALS
uv run pytest eda_shared/tests                      # 15 tests must pass
```

In VS Code select `.venv` (Python 3.12) as the notebook kernel. Tokens live only in the gitignored `.env`.

`cryptography<49` is pinned in `pyproject.toml`: newer versions ship no Intel-macOS wheel.

Machine used: Intel i5 macOS, 8 GB RAM. Heavy aggregation runs in DuckDB (streamed from disk); strings are
Arrow-backed (`src/__init__.py`).

## Building the cache

```python
import src
from src import fetch_fafb, fetch_malecns
from src.validate import validate_dataset, assert_interchangeable
fn, fe, fr = fetch_fafb.load()          # downloads raw files on first run
mn, me, mr = fetch_malecns.load()
```

| Step | Cold | Warm (cache) |
|---|---|---|
| FAFB neurons + edges + edgeroi | 76 s (after download), peak RSS 2.6 GB | 6 s, no network |
| MaleCNS neurons (neuPrint roiInfo for 165,122 neurons) | 862 s | — |
| MaleCNS edges | 20 s, peak RSS 2.9 GB | — |
| MaleCNS edgeroi (124M synapse rows → 29.4M, DuckDB) | 228 s, peak RSS 2.6 GB | — |
| `notebooks/eda_B.ipynb`, all B items, network blocked | — | 159–211 s |

Cache files (D7), in `data/cache/` (gitignored):

```
fafb_neurons_783_ann3.2.0.parquet
fafb_edges_783_buhmann_min1.parquet
fafb_edgeroi_783_buhmann_min1.parquet
malecns_neurons_v1.0_traced.parquet
malecns_edges_v1.0_traced_min1_conf0.5.parquet
malecns_edgeroi_v1.0_traced_min1_conf0.5.parquet
```

## Step 0 counts

| Table | Expected | This repository |
|---|---|---|
| fafb neurons | 139,248 | 139,262 (see note 1) |
| fafb edges ≥ 5 | 2,700,513 | **2,700,513** |
| fafb edges, τ = 1 | — | 15,091,983 (54,492,922 synapses) |
| fafb edgeroi rows | — | 16,847,997 |
| malecns edgeroi rows | — | 29,436,294 (per-pair sums equal edges for all 25,563,096 pairs) |
| malecns neurons | 165,122 | **165,122** |
| malecns edges ≥ 5 | one agreed number | **6,235,649** (self-connections dropped; 6,235,682 with them) |
| malecns edges, τ = 1 | 25.56M | 25,563,096 (25,563,197 rows before dropping self-connections) |
| self-connections dropped | fafb 0, malecns 101 | fafb 0, malecns 101 |

## Points reported to the group

1. **FAFB ids.** The v3.2.0 annotation file has 139,248 rows, but 7 of its root ids are not 783 ids (later edits:
   ORN_DM4, 2× L2, Mi9, 2× HBeyelet, one unlabelled optic neuron); they appear in no 783 edge. Conversely
   14 ids of the 783 proofread set (139,255) appear in edges but not in the annotations. 139,248 − 7 + 14 = 139,255.
   Here the 14 are added with unknown annotations (needed for referential integrity) and the 7 are kept,
   giving 139,262 rows. The Step 0 note ("add the 7 missing ids") undercounts.
2. **`sensory_ascending`** is a FAFB super_class (612 neurons; 537 in MaleCNS) missing from
   `config.SUPER_CLASS_VALUES`. Provisionally mapped to `sensory` (flow = afferent).
3. **MaleCNS-only super_class levels** (`vnc_intrinsic`, `vnc_efferent`, `cb_efferent`, `efferent_ascending`,
   `efferent_descending`, `sensory_descending`, `ENS`, and all `*_tbc`) map to `unknown` until the group adds them.
4. **MaleCNS edges ≥ 5**: 6,235,649 from the bulk traced-only weights at confidence ≥ 0.5, self-connections dropped.
5. **MaleCNS region names `add_provisional_flags` rejects** (not added locally):
   `CentralBrain-unspecified` (365 neurons) and `Optic-unspecified(L)` / `Optic-unspecified(R)` (359 + 37) are brain regions missing from `neuropils_provisional.csv`. Group decision: add both as
   brain neuropils (`is_optic` True / False, `region_group = unknown`); until the shared table has them they count as unknown. The other 34 rejected names are nerve-cord neuropils
   and nerves, passed as `nerve_cord`.

## Source column mapping

| Target | FAFB | MaleCNS |
|---|---|---|
| neuron_id | root_id | bodyId |
| cell_type / cell_class / cell_sub_class | cell_type / cell_class / cell_sub_class | type / class / subclass |
| super_class | super_class | superclass (mapped, see above) |
| side | side (`na` → unknown) | somaSide, else rootSide (L/R/M → left/right/center) |
| region | computed from edgeroi | computed from neuPrint roiInfo (primary ROIs) |
| nt_type / nt_confidence | top_nt / top_nt_conf | predicted_nt / predicted_nt_confidence (`unclear` → unknown) |
| n_pre / n_post | missing | neuPrint pre / post |
| flywire_type / manc_type | missing | flywireType / mancType |
| hemibrain_type | hemibrain_type | hemibrainType |
| dimorphism, status | dimorphism, status | dimorphism, status |
| hemilineage | ito_lee_hemilineage | itoleeHl |
| nt_known | known_nt | missing |
| pre_id / post_id / syn_count | pre_root_id / post_root_id / syn_count | body_pre / body_post / weight |
| neuropil (edgeroi) | neuropil (e.g. `AVLP_R`) | primary_post (e.g. `AVLP(R)`); null → `unknown` |

**Region rule (D15, provisional).** The neuropil holding most of the neuron's synapses, inputs plus outputs,
counted per connection, the Guide's unit (§1: each sender–receiver contact is one synapse).
FAFB: edgeroi at τ = 1; MaleCNS: roiInfo `upstream + downstream` per primary ROI. Not `pre + post`, which counts
output release sites; FAFB has no such count (its per-neuropil input and output totals are equal).
Ties go to the alphabetically first neuropil. Region keeps the source spelling, side marker included;
`eda_shared.canonical_neuropil` maps both spellings to FlyWire base names.

## EDA workstream B

`notebooks/eda_B.ipynb` regenerates every figure in `figures/eda/B*` and every table in `reports/eda_B_tables/`
from the cache, with no network calls. Summary: `reports/eda_B.md`. Numbers: `reports/eda_B_numbers.json`.
