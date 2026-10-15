# EDA Workstream B — Edge weights and thresholds (Project 2)

fafb = FlyWire 783, annotations v3.2.0, Buhmann synapses. malecns = male CNS v1.0, Traced, conf ≥ 0.5, **brain part** (141,996 neurons).
Central = `is_optic == False` (fafb 44,494; malecns 43,623). Figures: `figures/eda/B*`; tables: `reports/eda_B_tables/`; notebook: `notebooks/eda_B.ipynb`.

| Item | Number | What it means |
|---|---|---|
| **B1** Edge strength | 9.1% (fafb) and 11.9% (malecns) of central edges reach 10 synapses; 1.8% / 2.8% reach 31. Aggregation check: edgeroi sums equal edge sums exactly (fafb 54.49M, malecns-brain 97.31M synapses); self-connections 0 in both tables (logged drops: fafb 0, malecns 101) | Most edges are weak (Guide §1). D8 is implemented correctly |
| **B2** Threshold table | Central, τ = 10: isolated neurons **11.0% fafb vs 1.6% malecns**; edges kept 9.1% / 11.9%; synapses kept 53.7% / 60.8%; median in/out degree 4/5 vs 9/13. τ = 31 isolates 47.5% / 16.9%. Relative 0.9% isolates 1.5% / 0.1% | τ = 10 leaves both central graphs usable, but **one absolute threshold is a different cut in the two datasets** (D10, D11) |
| **B3** Absolute vs relative | Central, ≥ 10 syn vs ≥ 0.9% of input: disagree on **74% (fafb) / 67% (malecns)** of edges passing either (15.6% / 11.6% of all edges). fafb: 206k pass only the count, 685k only the share, 315k both | The criteria are not interchangeable; one must be primary (D12) |
| **B4** Concentration | Edges ≥ 10 syn: 9% of central edges carry 54% of synapses (fafb); 12% carry 61% (malecns) | Lower than Schlegel et al.'s 16% / 79%, which is **type-level**; ours is neuron-level |
| **B5** Left–right (fafb) | Type-level, central, τ = 10: **R(log10) = 0.79** on 76,603 pairs present on both sides; **56%** of them differ by > 30%; 65,249 pairs on one side only. τ = 1: R = 0.77, 68% | Well below the published R = 0.97: within one brain, weights of the same type pair are noisy (Guide §6) |
| **B6** Across animals | 7,026 matched types, right side, τ = 10: R = 0.77; 72% of shared pairs differ by > 30%; **median log2(malecns / fafb) = +0.58 (×1.5)** | A systematic pipeline offset (malecns records more synapses per connection); extra spread over B5 is small |
| **B7** Magnitude in features | fafb central, 43,871 neurons × 162 raw-count columns: PC1 vs log10 total synapses **Pearson r = 0.43, Spearman ρ = 0.86**; PC3 r = −0.39 (ρ = −0.61). Row-normalised profiles: PC1 r = 0.00, ρ = −0.05 | **Neuron size orders the raw feature space**; Pearson understates it because a few huge neurons dominate PC1. Normalising removes it. Project 2's first measurement |

**What surprised me.** The two datasets differ by almost a factor of two in synapses (central: 46.3M synapses on 8.65M edges in malecns vs 24.2M on 5.73M in fafb, for about the same number of neurons, 43.6k vs 44.5k). So "10 synapses" isolates 11% of fafb central neurons but under 2% of malecns, and type-pair weights are 1.5× higher in malecns.

## Proposals for the worksheet (Project 2's view)

- **D10.** τ = 10 as default, weak confidence. At τ = 10, 11.0% (fafb) / 1.6% (malecns) of central neurons are isolated and 54% / 61% of synapses are kept. A connection below threshold is *unknown*, not absent: 56% of left–right type pairs already differ by > 30% (B5).
- **D11.** 1 synapse (permissive) and 2.6% of input (strict). They bracket the default from both sides, and both have measured persistence (42% and 99%, Schlegel et al. 2024).
- **D12.** Relative primary for cross-dataset claims, weak confidence. The absolute criterion is not comparable across pipelines (B2, B6: ×1.5 offset), and the two criteria disagree on 67–74% of passing edges (B3).
- **Bias disclosure.** Project 2's argument is that magnitude is unreliable. That makes me favour normalised or relative choices, so treat my D12 preference with that in mind.

## Choices made where the plan was silent

- **B5/B6 threshold.** τ = 10 is applied to neuron edges before summing to type level, following the plan's default. τ = 1 is reported alongside.
- **"Differ by > 30%".** Defined as larger / smaller > 1.3.
- **B2 degrees and isolated neurons.** Degrees count partners inside the population and include zeros. "Isolated" means no edge left inside the population.
- **B7 per-neuropil columns.** Use all cached edgeroi rows (τ = 1), the same basis as `total_in_syn`. 623 neurons with no synapses are dropped (log undefined). Spearman ρ is added next to the specified Pearson r.
- **MaleCNS brain part.**
  - Excludes 338 neurons with no region.
  - Region rule (D15): the neuropil with most synapses counted per connection, as in FAFB's edgeroi; for MaleCNS roiInfo `upstream + downstream` (group decision; `pre + post` would count output release sites, which FAFB does not).
  - Excludes, for now, 761 neurons whose region `eda_shared` rejects (`CentralBrain-unspecified` 365, `Optic-unspecified` 359 + 37). The group agreed to add both as brain neuropils (`region_group = unknown`); until the shared table has them, they stay out and these numbers will move slightly.
