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

Proposals for D10–D12, the bias disclosure and the choices made where the plan was silent: `reports/eda_B_notes.md`.
