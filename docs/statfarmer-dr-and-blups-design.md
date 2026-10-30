---
cursor:
  subagentId: "bc-ac726eba-e2d0-569c-8d26-7ade819d48ce"
---

# StatFaRmer design: trait dimensionality reduction + side-by-side BLUPs

**Status:** design only (no feature implementation).  
**Drivers:** 325-cultivar soybean screen (`project_2025-07-14_b51-91_soy`); many correlated PlantEye traits; `T:X:Y` translation where **T is a run/batch factor on the same physical table**, not field spatial autocorrelation.  
**Related:** `docs/soybean-325-analysis-brief.md`, `internal/soybean-325-analysis-proposal.md`, `internal/statfarmer-feature-seams.md`.

---

## 1. Problem

### 1.1 High-cardinality genotype ranking
StatFaRmer’s analysis path is built around classical ANOVA / time-block LMM / spline, then Tukey letters. For a genotype factor with hundreds of levels:

- Auto `resolveModel()` treats only **time-like** high-cardinality factors (`dbscan_cluster`, `timestamp_group`) as random. Cultivar stays **fixed**.
- Any factor with **>20 levels** forces the auto path to **classic ANOVA** (performance guardrail).
- Tukey is skipped when max cardinality >20; even if run, all-pairs among 325 genotypes (~52k comparisons) is unusable.

Breeding-style questions need **shrunk genotype effects (BLUPs)** or fixed BLUEs, optional heritability, top-k / Dunnett-vs-checks — not letter plots.

### 1.2 Many correlated traits
PlantEye exports dozens of traits (biomass, height, NDVI family, hue/lightness bins, etc.). Users currently pick **one** trait (`out_variables`, single select). Correlated traits inflate analysis surface and hide multivariate structure. Need an **intelligent reduction** path (PCA / sparse PCA / clustering of traits or of genotype profiles) that can consume genotype×trait summaries — preferably BLUPs, not raw means alone.

### 1.3 Design columns `T:X:Y`
Local soy project (counts only): 325 cultivars × 3 reps = 975 plants; **Treatment = 1 level**; translation `T:X:Y` with **T ≈ 41** distinct values (= `table_no` in handmade), **X ≈ 8**, **Y ≈ 3**.

**Modeling rule (Alana):** treat **T as blocking / run / batch** on the same physical table. Do **not** default to field-plot spatial autocorrelation on T. Still decide whether residual structure on **X:Y (table position)** helps vs ignore.

---

## 2. Constraints

| Constraint | Implication |
|------------|-------------|
| Do not break existing aov / lmer / spline paths | Prefer **additive modules** + optional soft UI hints; avoid changing auto `resolveModel` defaults in v1 |
| Launchers off-limits unless explicitly requested | No edits to `install.R`, `just_shiny.R`, `launch_statfarmer.R` |
| No case-registration plugin API | Mirror `module_effect_sizes.R` / `module_growth_summaries.R` pattern |
| Single-trait core today | DR needs multi-trait selection UI; BLUP module may still be per-trait first |
| Treatment may be degenerate | Soy screen has one Treatment; module must degrade gracefully and shine on multi-treatment projects |
| Performance | 325 G × ~40 traits × few time clusters; avoid `car::Anova` / emmeans explosion; guardrails similar to `MAX_FIT_SECONDS` |
| New R packages OK | Prefer packages already in `renv.lock` when sufficient (`lme4`, `emmeans`, `psych` present as deps); add only for clear gaps |

---

## 3. Proposed UX

### 3.1 Feature A — Side-by-side BLUPs (genotype + treatments)

**Placement:** new tab alongside Effect Sizes / Growth Summaries (does not replace ANOVA tab).

**Controls:**
- Genotype factor (default: `cultivar` / `genotype` if present)
- Treatment factor (default: `treatment`; disabled / message if <2 levels)
- Blocking factors: `repeat` / `Repeat`, **`table_no` / T (run block)** — checked by default when present
- Optional: include `dbscan_cluster` / time as additional block or fit **per time cluster** / on growth summary
- Mode: Genotype **random (BLUP)** vs **fixed (BLUE)**; Treatment **fixed** (default) vs random when levels justify
- Checks (optional multi-select) → Dunnett-style contrasts
- Top-k / bottom-k
- Trait: reuse current single `out_variables`, plus “run for all selected traits” batch later

**Outputs (side-by-side panels):**
1. **Genotype table:** level, BLUP/BLUE, SE (if available), rank, optional reliability / H² (when estimable)
2. **Treatment table:** same schema when ≥2 levels; otherwise empty state explaining single-level treatment
3. Model formula + variance components summary
4. Downloads (CSV)

**Messaging hooks in existing ANOVA/Tukey UI (non-breaking):** when cardinality of genotype-like factor > threshold, show link/text: “Use BLUP Ranking tab instead of all-pairs Tukey.”

### 3.2 Feature B — Intelligent trait dimensionality reduction

**Placement:** second new tab (or sub-tabs under one “Multivariate” tab). Depends logically on genotype×trait matrix.

**Controls:**
- Trait multi-select (defaults: numeric phenotyping columns minus meta IDs)
- Observation unit for matrix rows: **genotype** (default), or genotype×treatment
- Cell values: **BLUPs** (preferred, if Feature A available) / **BLUE** / **group means**
- Method: PCA (default) | correlation heatmap + hierarchical trait clustering | optional sparse PCA later
- Scaling: unit variance (default)
- Missing-data policy: pairwise complete / impute (simple median or `missMDA` if added)

**Outputs:**
- Scree / cumulative variance
- Loading table (traits × components) + biplot
- Genotype scores table (download)
- Suggested “reduced trait set” (high-loading unique traits per PC) for feeding back into single-trait ANOVA/BLUP

**Non-goals v1:** AMMI/GGE (needs multi-env), full GWAS export, deep autoencoder DR.

---

## 4. Statistical approach

### 4.1 Per-trait mixed model (Feature A)

Default single-environment form (illustrative):

```text
y ~ treatment + (1|genotype) + (1|repeat) + (1|table_no)
```

Notes:
- Drop `treatment` if <2 levels.
- **`table_no` / T = run block**, not spatial field trend.
- Genotype fixed alternative: `y ~ genotype + treatment + (1|repeat) + (1|table_no)` → BLUEs via `emmeans` / `fixef`.
- Time: either subset one `dbscan_cluster`, add `(1|dbscan_cluster)` only when scientifically a block, or analyze growth-summary endpoints from existing Growth Summaries module.
- **X:Y residual structure (phase 2+):** optional correlated residual or smooth on table coordinates; default **off** until residual diagnostics show position pattern after T-block. Do not invent field AR1/exp spatial on T.

**Heritability (when design supports):** broad-sense from variance components for genotype vs residual (+ blocks as appropriate). Report only when RE genotype model converged non-singularly; else mark “not estimable.”

**Inference:** overall genotype variance / LR or parametric bootstrap optional later; v1 prioritizes ranking table + optional Dunnett vs checks + FDR on that small contrast set. Never enable all-pairs Tukey for genotype > ~20 levels.

### 4.2 Multivariate reduction (Feature B)

1. Build matrix **G × P** (genotypes × traits) from BLUP/BLUE/means.
2. Standardize columns; PCA via `stats::prcomp` (baseline).
3. Trait clustering on absolute correlation distance as companion view for “which traits are redundant.”
4. Optional later: `psych::fa` / sparse PCA if loadings unstable or P large relative to G.

Link to Feature A: DR should **prefer BLUP matrix** so reduction reflects genetic ranking signal, not raw phenotyping noise.

---

## 5. Package candidates

| Need | Prefer first | Optional add |
|------|--------------|--------------|
| LMM BLUP/BLUE | **`lme4`** (already used) + `broom.mixed` if helpful | `sommer` only if multi-trait G-structure needed later |
| Means / Dunnett | **`emmeans`** + `multcomp` (already used) | — |
| PCA / FA | **`stats::prcomp`**, **`psych`** (already in lock as dependency) | `FactoMineR`, `missMDA` |
| Sparse PCA | — | `elasticnet` / `sparsepca` (defer) |
| Heritability helpers | hand-rolled from `VarCorr` | `heritability` package (defer unless reuse pays off) |

**v1 recommendation:** no new hard dependency for BLUPs; add at most one multivariate helper (`FactoMineR` or `missMDA`) if base+`psych` proves awkward in UI.

---

## 6. Seams (modules, contracts, UI hooks)

Full call-graph notes: `internal/statfarmer-feature-seams.md`.

### 6.1 Do not invade auto `resolveModel` in phase 1
Current chain: UI Submit → `analysisResults` → `performANOVA` → `resolveModel` → aov | time-block lmer | spline → optional Tukey.

**Feature A/B are sidecar modules** (same pattern as Effect Sizes / Growth Summaries):
- Files: `shiny/app/module_blup_ranking.R`, `shiny/app/module_trait_reduction.R`
- Stats helpers: `shiny/common/blup_ranking.R`, `shiny/common/trait_reduction.R` (sourced from `stats.R` or module)
- Wire: `ui.R` tabset + `server.R` `*Server(...)` with `filteredData`, `combined_inputs`, optionally `analysisResults` / shared BLUP cache

### 6.2 Soft, non-breaking UI hooks
| Location | Change |
|----------|--------|
| `ui.R` tabset | Add two module UIs |
| `server.R` | Register servers |
| `reactive_outputs.R` Tukey skip message | Point users to BLUP tab when genotype-like card high |
| Advanced Model Options | **Do not** add `force_method='blup'` until phase 3; keeps aov/lmer/spline semantics intact |

### 6.3 Data contracts (tibble schemas)

**`blup_table`**
- `effect_type` ∈ {`genotype`,`treatment`}
- `level`, `estimate`, `std_error`, `rank`, `model_role` ∈ {`BLUP`,`BLUE`}
- `trait`, `formula`, `converged`, `singular`
- optional: `reliability`, `H2`, `n`

**`trait_pc_result`**
- `scores`: genotype (+ optional treatment) × PC
- `loadings`: trait × PC
- `variance`: PC, prop, cumulative
- `source_matrix`: `BLUP` | `BLUE` | `mean`
- `traits_used`, `scaling`

Shared reactive (optional): `genotypeTraitMatrix()` produced by BLUP batch or means fallback, consumed by DR module.

### 6.4 Wizard / data layer
No wizard changes required for v1 if `table_no` / `treatment` / `cultivar` already land in merged tables (raw pipeline joins handmade + translation). **Open check:** ensure processed merged RDS retains `treatment` and exposes T (as `table_no`) consistently; X/Y may need parsing from translation if not already columns — additive parser only if Feature A phase 2 spatial diagnostics need them.

---

## 7. Phased delivery

| Phase | Scope | Risk to existing paths |
|-------|--------|-------------------------|
| **0 — Design** (this doc) | Seams + contracts | None |
| **1 — BLUP module scaffold** | `module_blup_ranking.R` + `blup_ranking.R`; genotype random + repeat + table_no; treatment panel empty-state; CSV download; fixture on soy project | None (additive) |
| **2 — Treatment BLUPs + Dunnett/top-k** | Side-by-side when treatment ≥2; checks UI; Tukey skip copy update | Copy-only in reactive_outputs |
| **3 — Trait DR module** | Multi-trait select + PCA on means; then switch source to BLUP matrix | Additive |
| **4 — Polish** | H², batch-all-traits, X:Y residual option behind advanced toggle, docs in `app_guide.md` / `common_guide.md` | Still avoid resolver default changes |
| **5 — Optional resolver awareness** | Detect genotype card >20 → suggest module; still no silent formula rewrite | Careful; gated feature flag |

**Tiny stub allowed before Alana OK:** empty Shiny module tabs that render “Design approved?” — **not required**; prefer docs-only until approval.

---

## 8. Risks

- Fitting 325-level RE with singular fits / timeout → need clear fallback to BLUE or means + warning.
- Users confuse **T run block** with field spatial models → UI help text must state Alana rule.
- DR on raw means overstates technical/table artifacts; BLUP-first matrix mitigates.
- Single-treatment soy project makes treatment panel look “broken” unless empty-state is explicit.
- Touching `resolveModel` auto thresholds could regress time-series ANOVA performance → keep out of critical path.
- Package creep (`sommer`, asreml-style) without MET data.

---

## 9. Open questions (block coding until Alana answers)

1. **Primary goal:** rank genotypes for selection, publish tests, or both?
2. **Check cultivars:** names/IDs for Dunnett?
3. **Confirm T:** always run/batch block (yes per latest guidance) — any experiment where T *is* field spatial?
4. **X:Y:** ignore for v1, or require residual diagnostics toggle in phase 4?
5. **Priority traits** for default multi-select in DR?
6. **Time handling default:** per `dbscan_cluster`, last cluster only, or growth-summary endpoint?
7. **Treatment** in soy: intentionally single-level forever, or metadata fix expected?
8. Product name / tab labels: “BLUP Ranking” vs `multi_genotype_ranking` case language from prior brief?

---

## 10. Recommended next implement slice (after Alana OK)

**Phase 1 only:** additive `shiny/common/blup_ranking.R` + `shiny/app/module_blup_ranking.R`, wired like Effect Sizes; fit  
`y ~ (1|cultivar) + (1|repeat) + (1|table_no)` on current trait; genotype BLUP table + variance components; treatment panel empty-state when <2 levels; no `resolveModel` changes; no new packages; test mentally / with local soy RDS — no cloud raw dump.

---

## 11. Evidence index (this machine)

- Checkout: `C:\Users\user\StatFaRmer_new\StatFaRmer`
- Resolver: `shiny\common\model_resolver.R`
- ANOVA override: `shiny\common\stats.R`
- Tukey skip: `shiny\app\reactive_outputs.R` (`max_card > 20`)
- Module pattern: `shiny\app\module_effect_sizes.R`, `module_growth_summaries.R`
- Wiring: `shiny\app\ui.R`, `server.R`
- Example project: `data\project_2025-07-14_b51-91_soy\` (325 G, 3 reps, Treatment=1, T≈41, X≈8, Y≈3)
