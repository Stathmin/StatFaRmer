# Changelog

All notable changes to StatFaRmer are documented here.

## [Unreleased]

### Added
- Multi-trait DR on **named growth-summary metrics** (`auc` / `peak` / `growth_rate` / `last`) instead of silent across-time raw means (`collapseManyTraitsBySelectedTimes`, `require_named_temporal_metrics`).
- Shared HortControl trait labels (`shiny/common/hortcontrol_labels.R`) — e.g. `GLI bin [0, 0.15] · auc`, `Digital Biomass · last`.
- Biology readout module (`shiny/common/biology_readout.R`) with HortControl §4.3 spectral glossary, GLI-bin polarity, size↔greenness archetypes; Clippy/pipeline coach wiring.
- Soy-325 (~325 G) demo suite scripts under `scripts/render_soy325_*.R` + `soy325_attach_cultivar_names.R`.
- Cursor rule forbidding across-time raw means for posthoc / multi-trait DR.
- Cursor rule: never commit/push `data/project_*` or experiment `outputs/` (program only).
- Tests: biology readout, temporal DR metrics, high-card plots, trait follow-up.

### Changed
- DR score biplots: biomass + NDVI arrows, dense-cloud-aware scaling, full score range on axes, short subtitle + caption.
- Loadings / biplot axes use HortControl-style parameter names.

### Warnings
- **!!! THIS CODE SHOULD BE CHECKED BY A HUMAN PLZ !!!** — biology text, GLI-bin maps, and polarity heuristics need human review before citation.
