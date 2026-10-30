---
cursor:
  subagentId: "bc-ac726eba-e2d0-569c-8d26-7ade819d48ce"
---

# StatFaRmer feature seams — DR + BLUPs

Companion to `docs/statfarmer-dr-and-blups-design.md`. Code-path map for implementers. Design English; no app code changed in this pass.

## Current analysis spine

```text
wizard (data/project_*/ → merged_table.rds, config.json, vector_of_groups)
    ↓
app ui.R (project, factors, cultivars, treatments, single trait, model_selection_method)
    ↓ Submit
reactive_outputs.R::analysisResults
    → formulate(out_var, anova_factors)
    → performANOVA(...)          # stats.R override
        → resolveModel(...)      # model_resolver.R
            → aov | lmer(time RE) | spline | high-card→aov
    → performTukey if max_card ≤ 20 else skip message
    ↓
tabs: Config | Raw | Descriptive | ANOVA | Tukey | Assumptions | Letters
      + effectSizesUI | growthSummariesUI
```

### `resolveModel` decision facts (do not break)

| Condition | Behavior |
|-----------|----------|
| `force_method == 'aov'` | Classic aov |
| `dbscan_cluster` ≤10, auto | Fast-path aov |
| High-card **time-like** only (`dbscan_cluster`, `timestamp_group`) | Candidates for `(1\|…)` |
| Auto + **any** factor max levels >20 | **Disable** lmer/spline → aov |
| Forced lmer/spline + card >20 | Try keep; Anova may timeout → aov fallback |
| Cultivar / genotype high card | Remains **fixed**; does **not** become RE today |

Blocking validation: `random_effects_validator.R` (`buildValidatedRandomEffects`). Unit RE permanently disabled.

### Existing additive module pattern (copy this)

| Piece | Effect sizes | Growth summaries |
|-------|--------------|------------------|
| UI fn | `effectSizesUI(id)` → `tabPanel` | `growthSummariesUI(id)` |
| Server | `effectSizesServer(id, filteredData, analysisResults, combined_inputs)` | `growthSummariesServer(id, filteredData, combined_inputs, projectData)` |
| source | `ui.R` + `server.R` | same |
| Depends on resolver? | Reads `analysisResults` optionally | **No** — own calc |

**Implication:** BLUP + DR should follow Growth Summaries more than ANOVA: own fit, own tab, no `force_method` required.

## Where Feature A (BLUPs) plugs in

```text
NEW shiny/common/blup_ranking.R
  fitBlupModels(data, trait, genotype_var, treatment_var, block_vars, ...)
  → list(genotype_table, treatment_table, varcomp, formula, meta)

NEW shiny/app/module_blup_ranking.R
  blupRankingUI / blupRankingServer
  inputs: filteredData, combined_inputs (, projectData)
  does NOT call resolveModel

WIRE
  ui.R tabsetPanel: blupRankingUI("blup_ranking")
  server.R: blupRankingServer("blup_ranking", ...)

SOFT HOOK (phase 2)
  reactive_outputs.R tukey_skip_reason text → mention BLUP tab
  optional help under Advanced Model Options (no new force_method yet)
```

### Factor mapping for soy-shaped data

| Concept | Likely columns | Role in Feature A |
|---------|----------------|-------------------|
| Genotype | `cultivar`, `genotype` | `(1\|genotype)` or fixed BLUE |
| Treatment | `treatment` (from handmade `Treatment`) | Fixed if ≥2 levels |
| Rep | `repeat` / handmade `Repeat` | `(1\|repeat)` |
| Run block T | `table_no` (T in `T:X:Y`) | `(1\|table_no)` — **not spatial** |
| Table position | X, Y from translation (may be absent in RDS) | Phase 4+ residual only |
| Time | `dbscan_cluster`, `timestamp_group` | Subset or separate from genotype RE path |

**Keep time-block lmer/spline path untouched** for growth questions; BLUP module is for genotype/treatment ranking on a chosen slice or summary.

## Where Feature B (DR) plugs in

```text
NEW shiny/common/trait_reduction.R
  buildGenotypeTraitMatrix(data|blups, traits, row_keys, value_source)
  runTraitPca(matrix, scale=TRUE)
  clusterTraits(cor_matrix)   # optional companion

NEW shiny/app/module_trait_reduction.R
  traitReductionUI / Server
  inputs: filteredData, combined_inputs, optional genotypeTraitMatrix reactive from BLUP batch

WIRE same tabset pattern
```

**Does not** go through `resolveModel`. May call Feature A helper to fill cells with BLUPs (phase 3), else means.

Trait discovery heuristic: numeric columns minus meta set  
`{unit, dbscan_cluster, timestamp, cultivar, genotype, treatment, repeat, table_no, timestamp_group, outlier, var, label, v_t_r, g_alias, sort_name, …}`.

## What not to touch (phase 1–3)

- `resolveModel` auto cardinality thresholds
- Spline / `fitSplineModel` / time RE validation defaults
- Launcher files (`install.R`, `just_shiny.R`, `launch_statfarmer.R`)
- Tukey letter algorithms (except skip messaging)
- Wizard DBSCAN (time clustering ≠ genotype/trait clustering)

## Optional phase 5 resolver awareness (gated)

If product wants auto-detect only:

```text
if (genotype_levels > THRESH && is.null(force_method))
  # do NOT rewrite formula
  # set flag analysisResults$suggest_blup_module <- TRUE
```

Still fit aov for ANOVA tab compatibility; UI banner points to BLUP module.

## Contract sketches

```r
# blup_table columns
# effect_type, level, estimate, std_error, rank, model_role, trait,
# formula, converged, singular, H2, reliability, n

# trait_pc_result
# list(scores, loadings, variance, source_matrix, traits_used, scaling)
```

## Test anchors (local only)

- Project: `data/project_2025-07-14_b51-91_soy/`
- Expect: 325 cultivars, 3 reps, 1 treatment, ~41 `table_no`, merged ~60 numeric traits, ~4 time clusters
- Assert: Feature A returns 325 genotype rows; treatment table empty-state; existing ANOVA submit still returns `model_type` consistent with today

## Git readiness note

Repo git-initialized on `master`; **no remote** observed this session. Working tree dirty with **data/** only — leave uncommitted. Feature work later: branch `cursor/<feature>-48ce` off clean `master` (or stash data noise); PR when implementing, not for design-only Agent Store docs.
