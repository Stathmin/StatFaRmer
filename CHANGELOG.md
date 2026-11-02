# Changelog

All notable changes to StatFaRmer are documented here.

## [Unreleased]

### Added
- Multi-trait DR on **named growth-summary metrics** (`auc` / `peak` / `growth_rate` / `last`) instead of silent across-time raw means (`collapseManyTraitsBySelectedTimes`, `require_named_temporal_metrics`).
- Shared HortControl trait labels (`shiny/common/hortcontrol_labels.R`) — e.g. `GLI bin [0, 0.15] · auc`, `Digital Biomass · last`.
- Biology readout module (`shiny/common/biology_readout.R`) with HortControl §4.3 spectral glossary, GLI-bin polarity, size↔greenness archetypes; Clippy/pipeline coach wiring.
- Soy-325 (~325 G) demo suite scripts under `scripts/render_soy325_*.R` + `soy325_attach_cultivar_names.R`.
- Growth Summaries **multi-trait plant panel** UI (select many PlantEye traits → `trait__auc` / `__peak` / `__growth_rate` / `__last`) so Trait DR matches the offline soy-325 multi-trait PCA path.
- Puppeteer smoke harness `scripts/browser/smoke_gs_dr_blup.mjs` (+ `run_app_smoke.R`) for GS → Trait DR → BLUP.
- Regression test that bare `validate(need…)` is masked by `jsonlite::validate` (`test-shiny-validate-namespace.R`).
- Cursor rule forbidding across-time raw means for posthoc / multi-trait DR.
- Cursor rule: never commit/push `data/project_*` or experiment `outputs/` (program only).
- Tests: biology readout, temporal DR metrics, high-card plots, trait follow-up.

### Changed
- DR score biplots: biomass + NDVI arrows, dense-cloud-aware scaling, full score range on axes, short subtitle + caption.
- Loadings / biplot axes use HortControl-style parameter names.
- BLUP / Trait DR / follow-up: `shiny::validate` / `shiny::need` so tidyverse’s `jsonlite::validate` mask cannot crash Clippy.
- Clippy coach: swallow sidecar reactive errors instead of blanking the floating panel.
- Trait DR enforces `require_named_temporal_metrics`; BLUP prefers named `*__auc` / `*__last` over bare legacy metric names.

### Warnings
- **!!! THIS CODE SHOULD BE CHECKED BY A HUMAN PLZ !!!** — biology text, GLI-bin maps, and polarity heuristics need human review before citation.

### Known issues (tracked on GitHub)
- Main `distPlot` with ~300+ cultivars + several timestamps / timeseries: legend/plot overflow and cropping.
- Single time-cluster (no timeseries): main tab shows a small subset of cultivars as boxplots instead of violins.
