# BLUP + Trait DR Sidecars Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship additive **BLUP** and **Trait DR** Shiny sidecars (approach B): plant-curve collapse → optional per-table spatial residualize → `lme4` genotype BLUP on one trait or plant PC hyper-trait; soft Tukey hint; debug on soy-325 via RDS artifacts.

**Architecture:** Two-stage sidecar (not `resolveModel`). Growth Summaries honors selected times and yields one row/plant. Spatial stage (`mgcv` per-`table_no` X/Y smooths) residualizes when X/Y exist. BLUP/`DR` modules mirror Effect Sizes / Growth Summaries. Scripts and tests `saveRDS`/`readRDS` results — agents parse RDS, do not `tail -f` loops.

**Tech Stack:** R 4.x, Shiny, `lme4`, `emmeans`, `psych`, `mgcv` (already in lock transitive; ensure project use), `testthat` 3.x, tidyverse-style dplyr.

**Spec:** `docs/superpowers/specs/2026-10-30-blup-dr-sidecars-design.md` (approved)

## Global Constraints

- **No** edits to `install.R`, `just_shiny.R`, `launch_statfarmer.R` unless Alana names them.
- **No** `resolveModel` rewrite; Tukey path = soft copy only.
- BLUP only on **one** plant-level trait or **one** plant PC hyper-trait — never mass-BLUP-all-traits-then-PCA.
- DR = **plant × traits** only; factors are design columns, not PCA rows/vars.
- Collapse time only via Growth Summaries over **selected** clusters; block silent multi-time BLUP on raw traits.
- When X/Y present: always spatial stage with **own smooth per `table_no`** (T|X + T|Y, coarse).
- Tab title: **`BLUP`**. Checks: top-k if no check vector; no new wizard check UX.
- Git: no data dirt; feature commits only.
- **REQUIRED project skills** (read before coding in that area):
  - `.cursor/skills/r-package-development/SKILL.md` — `Rscript`, load helpers
  - `.cursor/skills/r-tidyverse-style/SKILL.md` — dplyr/style
  - `.cursor/skills/r-testthat/SKILL.md` + `references/fixtures.md` — fixtures, `readRDS` expected outputs
  - `.cursor/skills/shiny-r-optimize/SKILL.md` + `references/data-loading.md` — load once; precompute → `saveRDS`/`readRDS`

## Review Focus

- Multi-time raw trait + BLUP click → **block + message**, never silent mean.
- Soy Treatment 1-level → treatment panel **empty**, genotype table still ranks.
- Missing X/Y → BLUP still runs **without** spatial stage (no crash).
- DR never accepts cultivar as matrix rows; plant PC scores align 1:1 with plant IDs for BLUP.
- 41-table spatial on soy → finishes or times out with clear error object in RDS (no hung session without message).

---

## File map

| File | Responsibility |
|------|----------------|
| `tests/testthat/helper-statfarmer.R` | source common R without Shiny; paths via `here` |
| `tests/testthat/fixtures/plant_curve_mini.rds` | tiny synthetic plants × times × traits + T:X:Y |
| `tests/testthat/fixtures/*.expected.rds` | golden outputs per task |
| `shiny/common/txy_parse.R` | parse `T:X:Y` → `table_no`, `pos_x`, `pos_y` |
| `shiny/common/plant_collapse.R` | collapse filtered times → one row/plant metrics |
| `shiny/common/spatial_residualize.R` | per-table `mgcv` smooths; return residuals + meta |
| `shiny/common/blup_ranking.R` | `fitBlupModels(...)` → genotype/treatment tables |
| `shiny/common/trait_reduction.R` | plant×trait PCA + cluster; attach PC columns |
| `shiny/app/module_blup_ranking.R` | BLUP tab UI/server |
| `shiny/app/module_trait_reduction.R` | DR tab UI/server |
| `shiny/app/module_growth_summaries.R` | use selected times; expose plant summary reactive |
| `shiny/app/ui.R` / `server.R` | source + wire modules |
| `shiny/app/reactive_outputs.R` | Tukey skip copy → BLUP tab |
| `scripts/probe_blup_soy325.R` | offline soy probe; writes `tmp/probe_*.rds` (gitignored) |

---

### Task 1: testthat scaffold + mini fixture

**Files:**
- Create: `tests/testthat.R`, `tests/testthat/helper-statfarmer.R`, `tests/testthat/fixtures/plant_curve_mini.rds`
- Create: `scripts/build_plant_curve_mini_fixture.R` (regenerates fixture; commit the `.rds`)

**Interfaces:**
- Produces: fixture list/df with columns at least `unit`, `V.T.R`, `cultivar`, `treatment`, `repeat`, `table_no` or `T:X:Y`, `timestamp_group`/`dbscan_cluster`, ≥2 numeric traits, ≥2 times per plant, ≥2 tables with X/Y variation. Small: ~12 plants, 4 cultivars, 2 reps, 3 tables.

**Skills:** `r-testthat` fixtures; `r-package-development` for `Rscript` regenerate.

- [ ] **Step 1: Write failing test** `tests/testthat/test-fixture-mini.R`

```r
test_that("plant_curve_mini fixture has plant-time structure", {
  fix <- readRDS(testthat::test_path("fixtures", "plant_curve_mini.rds"))
  expect_true(all(c("unit", "cultivar") %in% names(fix)))
  expect_gte(dplyr::n_distinct(fix$unit), 8)
  expect_gte(dplyr::n_distinct(fix$dbscan_cluster), 2)
})
```

- [ ] **Step 2: Run** `Rscript -e "testthat::test_file('tests/testthat/test-fixture-mini.R')"` — expect FAIL (missing fixture)

- [ ] **Step 3: Build + save fixture** via `scripts/build_plant_curve_mini_fixture.R` → `saveRDS` under `tests/testthat/fixtures/`

- [ ] **Step 4: Re-run test** — PASS

- [ ] **Step 5: Commit** `test: add plant_curve_mini fixture and testthat scaffold`

---

### Task 2: `txy_parse` + `plant_collapse`

**Files:**
- Create: `shiny/common/txy_parse.R`, `shiny/common/plant_collapse.R`
- Test: `tests/testthat/test-plant-collapse.R`
- Fixture out: `tests/testthat/fixtures/plant_collapsed_mini.expected.rds`

**Interfaces:**
- Produces:
  - `ensureTxyColumns(data) -> data.frame` with `table_no`, `pos_x`, `pos_y` (from `T:X:Y` if needed)
  - `collapsePlantsBySelectedTimes(data, trait, time_col, time_values, metrics) -> data.frame` one row per `unit`/`V.T.R` with metric columns + design cols preserved
- Consumes: fixture from Task 1

**Skills:** `r-tidyverse-style`; `r-testthat` expected RDS.

- [ ] **Step 1: Failing tests** for parse + collapse (selected 2 clusters → 1 row/plant; full times vs subset differ)

- [ ] **Step 2: Run** — FAIL

- [ ] **Step 3: Implement** `ensureTxyColumns`, `collapsePlantsBySelectedTimes` (AUC/peak/mean as needed; reuse Growth Summaries formulas where practical)

- [ ] **Step 4: `saveRDS` actual collapse → compare to expected** (create expected on first green run, commit)

- [ ] **Step 5: Commit** `feat: parse T:X:Y and collapse plants over selected times`

---

### Task 3: Growth Summaries curve bridge (selected times)

**Files:**
- Modify: `shiny/app/module_growth_summaries.R` (stop using full unfiltered `merged_table` for metrics; filter to selected timestamp groups / `dbscan_cluster` from `combined_inputs` / `filteredData`)
- Modify: `shiny/app/server.R` if return value must expose summaries to siblings (prefer shared reactive or `session$userData` / returned reactive — pick one pattern and document in Interfaces)
- Test: `tests/testthat/test-growth-selected-times.R` (unit-test the collapse helper path Growth Summaries will call; optional shinytest2 out of scope)

**Interfaces:**
- Consumes: `collapsePlantsBySelectedTimes`
- Produces: plant-level summary table available to BLUP/DR (exact mechanism: `growthSummariesServer` returns `reactive` of summaries, or writes to a `plantLevelData` reactiveValues passed in — **builder picks one; prefer return reactive list**)

- [ ] **Step 1: Failing test** — collapse with selected clusters ≠ all clusters

- [ ] **Step 2: Implement filter** in Growth Summaries calc path

- [ ] **Step 3: Wire expose** plant summaries to server for later modules

- [ ] **Step 4: Commit** `fix: Growth Summaries use selected time clusters only`

---

### Task 4: Spatial residualize (per-table smooths)

**Files:**
- Create: `shiny/common/spatial_residualize.R`
- Test: `tests/testthat/test-spatial-residualize.R`
- Artifact: `tests/testthat/fixtures/spatial_meta_mini.expected.rds`

**Interfaces:**
- Produces: `residualizeSpatial(data, y_col, table_col="table_no", x_col="pos_x", y_pos_col="pos_y", ...) -> list(data, meta)` where `data[[paste0(y_col, "_spat_resid")]]` holds residuals; `meta` has formula, converged, per-table n. If X/Y absent → identity (residuals = y) + meta flag `spatial_applied=FALSE`.
- Uses `mgcv`; **one smooth per table** for X and Y (coarse `k`).

**Skills:** `shiny-r-optimize` data-loading (saveRDS fit meta); do not hang tests — use mini fixture only.

- [ ] **Step 1: Failing test** — with X/Y, residual variance &lt; raw; without X/Y, passthrough

- [ ] **Step 2: Implement** `residualizeSpatial`

- [ ] **Step 3: Run tests; saveRDS meta artifact for agent inspection**

- [ ] **Step 4: Commit** `feat: per-table mgcv spatial residualize helper`

---

### Task 5: `fitBlupModels` core

**Files:**
- Create: `shiny/common/blup_ranking.R`
- Test: `tests/testthat/test-blup-ranking.R`
- Artifact: `tests/testthat/fixtures/blup_geno_mini.expected.rds`

**Interfaces:**
- Produces: `fitBlupModels(data, trait, genotype_var="cultivar", treatment_var="treatment", repeat_var="repeat", table_var="table_no", apply_spatial=TRUE, checks=NULL, top_k=20) -> list(genotype_table, treatment_table, varcomp, formula, meta, h2)`
- Behavior: optional spatial first; then `lme4::lmer` with `(1|genotype)` (+ treatment FE if ≥2 levels) + `(1|repeat)` + `(1|table_no)` unless spatial already absorbed table (document choice in `meta`). Genotype table: level, blup, se (if avail), rank. Treatment empty df + message if &lt;2 levels. Top-k always; Dunnett only if `checks` non-NULL.
- Guard: if `meta$n_time_per_plant > 1` on raw uncollapsed input → `success=FALSE` with reason (tests pass this).

- [ ] **Step 1: Failing tests** — ranks length = n genotypes; 1-level treatment empty; multi-time guard

- [ ] **Step 2: Implement** `fitBlupModels`

- [ ] **Step 3: `saveRDS` genotype_table; compare expected

- [ ] **Step 4: Commit** `feat: fitBlupModels genotype BLUP helper`

---

### Task 6: BLUP Shiny module + wire

**Files:**
- Create: `shiny/app/module_blup_ranking.R`
- Modify: `shiny/app/ui.R` (source + `blupRankingUI("blup_ranking")` in tabset)
- Modify: `shiny/app/server.R` (source + `blupRankingServer(...)` with plant-level data / filtered + combined_inputs)
- Style: match Effect Sizes / Growth Summaries; title **BLUP**

**Interfaces:**
- Consumes: `fitBlupModels`, plant collapsed data (from Task 3), `combined_inputs$out_variables` or user pick among plant traits / PCs
- Produces: DT genotype + treatment panels, download CSV, formula/varcomp text

**Skills:** `shiny-r-optimize` — heavy fit behind actionButton, not on every reactive tick.

- [ ] **Step 1: Scaffold UI/server** with Calculate button calling `fitBlupModels`

- [ ] **Step 2: Wire ui/server**

- [ ] **Step 3: Manual smoke** optional; unit tests already cover fit

- [ ] **Step 4: Commit** `feat: BLUP tab module wired`

---

### Task 7: Trait DR (`trait_reduction`)

**Files:**
- Create: `shiny/common/trait_reduction.R`
- Create: `shiny/app/module_trait_reduction.R`
- Modify: `ui.R` / `server.R`
- Test: `tests/testthat/test-trait-reduction.R`

**Interfaces:**
- Produces: `reduceTraitsPCA(plant_df, trait_cols, n_comp=NULL, scale=TRUE) -> list(plant_df_with_pcs, loadings, variance, cluster, meta)`
- Rows = plants only; `trait_cols` user-selected phenotypes (default all eligible numeric minus IDs/design).
- Attach `PC1`… to same plant rows for BLUP response dropdown.
- Never take genotype-aggregated matrix.

- [ ] **Step 1: Failing tests** — nrow(scores)==n_plants; cultivar not in PCA vars; PC columns bind

- [ ] **Step 2: Implement** helper + module (actionButton; scree/loadings/scores DT; download)

- [ ] **Step 3: Wire**; BLUP trait chooser can see PC columns when present in shared plant table

- [ ] **Step 4: Commit** `feat: plant×trait PCA DR tab`

---

### Task 8: Soft Tukey hint

**Files:**
- Modify: `shiny/app/reactive_outputs.R` (~tukey_skip_reason and high-card model_info copy)

**Interfaces:**
- Exact copy fragment: mention **BLUP** tab for high-cardinality genotype ranking (keep growth summaries mention if useful).

- [ ] **Step 1: Update** skip strings to include BLUP tab guidance

- [ ] **Step 2: Commit** `fix: point high-card Tukey skip to BLUP tab`

---

### Task 9: Soy-325 offline probe (RDS artifacts)

**Files:**
- Create: `scripts/probe_blup_soy325.R`
- Ensure `tmp/` in `.gitignore` if not present

**Skills:** `r-package-development` Rscript; `shiny-r-optimize` precompute; agents **`readRDS`** probe outputs.

- [ ] **Step 1: Script** loads `data/project_2025-07-14_b51-91_soy` merged RDS; collapse selected/all clusters via helpers; optional spatial; `fitBlupModels` on one trait; `saveRDS` list to `tmp/probe_blup_soy325.rds` with genotype head, n, formula, timings, errors

- [ ] **Step 2: Run** `Rscript scripts/probe_blup_soy325.R` — expect RDS written

- [ ] **Step 3: Agent reads RDS** (`str`, top ranks) — no log tail loop

- [ ] **Step 4: Commit script only** `chore: soy325 BLUP probe script` (not tmp RDS, not data dirt)

---

### Task 10: Full static interface tests

**Files:**
- Create: `tests/testthat/test-interface-blup-dr.R`
- Create: `tests/testthat/fixtures/interface_contract.expected.rds` (optional summary of signatures / column contracts)

**Goal:** After helpers+modules exist, freeze the **public R interface** without launching Shiny: every exported helper signature, return-list names, genotype_table columns, DR plant-row invariant, multi-time guard, 1-level treatment empty panel shape, spatial passthrough when no X/Y.

**Skills:** `r-testthat`; parse via `readRDS` of mini + optional soy probe head.

- [ ] **Step 1: Write comprehensive interface tests** covering `ensureTxyColumns`, `collapsePlantsBySelectedTimes`, `residualizeSpatial`, `fitBlupModels`, `reduceTraitsPCA` contracts from Tasks 2–7

- [ ] **Step 2: Run** full `testthat::test_dir('tests/testthat')` — all green

- [ ] **Step 3: Commit** `test: static interface contract suite for BLUP/DR`

---

### Task 11: Shiny polish

**Files:**
- Modify: `shiny/app/module_blup_ranking.R`, `shiny/app/module_trait_reduction.R`, `shiny/app/module_growth_summaries.R` as needed
- Modify: UI copy / empty states / progress (`withProgress` or disable button while fitting)
- Optionally light help text under BLUP / DR matching Effect Sizes tone

**Goal:** After static tests green — UX polish only: empty states, loading, CSV buttons, factor labels, soft errors in-panel, no resolver changes.

**Skills:** `shiny-r-optimize` (actionButton gates; no heavy work on idle reactives).

- [ ] **Step 1: Polish** BLUP + DR tabs (and Growth Summaries bridge messaging if thin)

- [ ] **Step 2: Manual checklist note in report** (soy empty treatment; PC selectable after DR)

- [ ] **Step 3: Commit** `polish: BLUP/DR Shiny UX`

---

## Self-review (spec coverage)

| Spec item | Task |
|-----------|------|
| BLUP tab + empty treatment | 5–6 |
| Spatial per-table when X/Y | 4–5 |
| DR plant×traits → hyper-traits | 7 |
| Curve bridge selected times | 2–3 |
| Soft Tukey hint | 8 |
| No resolveModel / no launchers | Global + tasks avoid |
| Soy debug via artifacts | 9 |
| No mass-BLUP-then-PCA | 5, 7 tests |

---

## Execution

Recommended: **Subagent-driven** (tasks have hard interface contracts; spatial + BLUP easy to regress). Use vendored skills every task. After plan OK → `superpowers:subagent-driven-development`.
