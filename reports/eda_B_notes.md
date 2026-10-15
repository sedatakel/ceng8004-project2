# EDA Workstream B — notes behind the summary

Companion to `reports/eda_B.md`, kept separate so the summary stays on one page.

## Proposals for the worksheet (Project 2's view)

The completed worksheet is in `submission/`.

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
