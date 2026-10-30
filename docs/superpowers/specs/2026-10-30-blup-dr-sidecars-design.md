# Design: BLUP + trait DR sidecars (approach B)

**Status:** approved (Alana, 2026-10-30)  
**Date:** 2026-10-30  
**Repo:** StatFaRmer (`C:\Users\user\StatFaRmer_new\StatFaRmer`)  
**Debug fixture:** `data/project_2025-07-14_b51-91_soy` (325 cultivars × 3 reps, Treatment = 1 level, T≈41 tables)  
**Related:** `docs/statfarmer-dr-blups-handoff.md`, `docs/statfarmer-dr-and-blups-design.md`, `internal/statfarmer-feature-seams.md`

**Doc precedence:** This spec + the 2026-10-30 Alana conversation win over older handoff/design notes when they conflict (notably: X/Y spatial is in scope under approach B; DR must not mass-BLUP-then-PCA; collapse via Growth Summaries over selected times).

---

## Intent and why

Classical ANOVA + all-pairs Tukey cannot rank ~325 genotypes. Users need breeding-style **BLUP rank + planned tests**, and a path to compress many correlated PlantEye traits into **hyper-traits** without treating repeated plant scans as independent rows.

Observations are **growth curves of specific plants** (`V.T.R`). Time dependence must be collapsed at the plant before genotype mixed models. Table geometry matters: **T = run/batch**; **X/Y = spatial cofactors** with a **per-table** smooth.

---

## Goals

1. **BLUP tab** — genotype ranking (BLUP) + planned tests (Dunnett/FDR when checks exist); treatment side panel when ≥2 levels; empty state when 1.
2. **Trait DR tab** — PCA (+ trait clustering) on **plant × traits**; user selects which traits; emit plant-level hyper-traits for BLUP.
3. **Curve bridge** — Growth Summaries over **user-selected** time points/clusters → one y per plant → feed BLUP/DR.
4. **Soft Tukey hint** when genotype cardinality is high → point at BLUP tab.
5. Debug and harden on soy-325 locally.

## Non-goals

- Rewriting `resolveModel` / main ANOVA auto path in this project.
- All-pairs Tukey for huge genotype.
- AMMI/GGE / multi-env.
- Field AR1/exp spatial treating T as a field plot.
- One-shot joint `gamm4` genotype+spatial (approach A) in v1 — deferred; packages may still be added for the spatial stage.

## Anti-goals (failure even if “works”)

- BLUP-fitting **all raw traits** then PCA on that matrix.
- Treating multi-time rows as i.i.d. in the genotype LMM.
- Silently auto-collapsing time without going through Growth Summaries over selected points.
- Dropping X/Y when present, or using one shared smooth instead of **per-table** smooths.
- Editing launchers (`install.R`, `just_shiny.R`, `launch_statfarmer.R`) unless Alana asks.
- Committing data dirt / new project folders.

---

## Constraints

| Constraint | Implication |
|------------|-------------|
| Additive modules only | Mirror `module_effect_sizes.R` / `module_growth_summaries.R`; wire `ui.R` / `server.R` only |
| User formula / factors | Follow existing `formulate(out_var, anova_factors)` + selected time clusters; do not invent a parallel time system |
| Packages | Prefer `lme4`, `emmeans`, `psych`; **spatial stage may add** `mgcv` (and helpers if needed) |
| Performance | 325 G × ~40 traits × ~4 clusters; guard slow paths; DR must not require BLUP-all-traits |
| Git | Feature commits without data junk; `.gitignore` already hides `data/project_*/**` untracked dirt |

---

## Approaches considered

| Approach | Idea | Tradeoff |
|----------|------|----------|
| **A — One-shot gamm4/mgcv** | Joint model: per-table spatial smooths + `(1\|cultivar)` | Statistically nicest; heaviest on soy (41 smooths × 325 G); harder Dunnett/`ranef` UX |
| **B — Two-stage sidecar (chosen)** | Collapse → spatial residualize per table → `lme4` BLUP on residuals; DR → hyper-traits → same BLUP | Debuggable; matches existing modules; joint uncertainty not fully shared |
| **C — mgcv-only RE smooths** | Genotype as `bs='re'` inside `gam` | Awkward breeding tables / emmeans vs `lme4` |

**Choice:** **B**, with Alana’s note that A is attractive later. Spatial libs are still allowed for stage 2 of B.

---

## Decided how

### Pipeline (strict order)

```text
User selects time points/clusters (existing UI)
        ↓
Growth Summaries over those selected points
  → one row per plant (V.T.R) with curve traits / endpoint
        ↓
Optional: Trait DR
  → user picks **traits** (default: all eligible numeric phenotypes)
  → matrix rows = **plants** (collapsed phenotype vectors) — never factor levels
  → PCA + trait correlation clustering on the trait space
  → hyper-trait columns (PC scores) on the **same plant rows**
        ↓
BLUP on exactly one of:
  (1) a single plant-level trait (raw endpoint or curve summary), or
  (2) a plant-level DRed hyper-trait (PC score)
Never: BLUP every raw PlantEye trait to feed DR.
Never: DR “on factors” (cultivar / treatment / … are design columns, not phenotypes).
```

### Time / dependence

- Biological unit = plant (`V.T.R`).
- Multiple `dbscan_cluster` / timestamps = dependent curve, not independent ANOVA rows for BLUP.
- Collapse **only** via Growth Summaries over **selected** points (AUC / peak / rate / etc. already in that module; extend as needed for curve-bridge phase).
- **Curve-bridge wiring:** today’s Growth Summaries uses full `projectData()$merged_table` and ignores the main time filter — that is wrong for this design. Bridge must summarize **only user-selected** timestamp groups / clusters (same selection as the main analysis filter), one row per plant/`unit`/`V.T.R`.
- If user requests BLUP on a raw trait while multiple time clusters remain uncollapsed: **block with guidance** to run Growth Summaries first (builder may refine copy; do not silent-mean).

### Spatial cofactors

- When X/Y present (parse from `T:X:Y` if not already columns): **always** include in the spatial stage.
- Structure: **table-specific** smooths — **its own smooth per `table_no` (T)** for X and for Y (notation T|X + T|Y), coarse only (tilted plane / saddle-scale), not fine nugget field models on T.
- **T / `table_no`** remains run/batch block in the BLUP stage (`(1|table_no)` and/or already absorbed via per-table spatial stage — builder picks non-double-counting form; mark as builder call below).
- Preferred implementation under B: `mgcv` smooths `s(X, by=table_no)` + `s(Y, by=table_no)` (or equivalent), residualize y, then `lme4` on residuals. Exact smooth basis = builder call within “coarse” constraint.

### BLUP stage (`lme4`)

Illustrative after residualization / when spatial absent:

```text
y ~ treatment + (1|genotype) + (1|repeat) + (1|table_no)   # drop treatment if <2 levels
```

- Genotype factor default: `cultivar` / `genotype`.
- Mode: random → BLUP (`ranef`); optional fixed → BLUE via `emmeans` later if cheap.
- Outputs: genotype table (effect, SE if available, rank), H² when estimable, formula + varcomps, CSV.
- Treatment panel: side-by-side when ≥2 levels; empty/disabled when 1 (soy).
- Tests: Dunnett vs checks when a check set is already present in project config/metadata; else **top-k only**. FDR on that small set. **No** all-pairs Tukey for high-card genotype.
- Check IDs: consume wizard-clean factors only — **no new check-ID UX and no wizard work in this project**. As of design time, soy-325 / app inputs expose no dedicated check list → ship top-k; wire Dunnett when/if a check vector already exists elsewhere.

### DR stage

- **What is reduced:** the **trait space** (correlated PlantEye / curve phenotypes). Factors are not DR inputs.
- Default traits: all eligible numeric phenotyping columns (minus IDs/meta); user can deselect.
- **Rows:** always **plants** (`unit` / `V.T.R`) after time collapse. Cultivar is not an observed phenotype row — do not build a cultivar×trait mean matrix as the DR default path, and do not offer factor levels as “row grain” for PCA.
- Cell values: plant-level collapsed trait values (means only if multiple rows still exist for one plant after collapse — not genotype BLUPs).
- Design columns (`cultivar`, `treatment`, `repeat`, `table_no`, X, Y) stay attached to plant rows for later BLUP/spatial; they are not PCA variables.
- Methods: PCA (`prcomp` / `psych` as needed) + trait correlation clustering.
- Outputs: scree, loadings, **plant** scores, suggested reduced trait set, download; hyper-traits (plant PC scores) selectable as BLUP response.

### Soft Tukey hint

- When Tukey path skips or genotype-like card high: soft copy pointing to **BLUP** tab (tab title **“BLUP”**, not “BLUP Ranking”).
- No resolver behavior change.

### Packages

- Core: existing `lme4`, `emmeans`, `psych`.
- Spatial stage: add `mgcv` (justified). Further spatial helpers only if `mgcv` insufficient for per-table coarse smooths.
- Do not add packages solely for approach A joint fit in this project.

### Seams / files (expected)

| Piece | Path |
|-------|------|
| Spatial + BLUP fit helpers | `shiny/common/blup_ranking.R` (and small spatial helper if cleaner) |
| BLUP module | `shiny/app/module_blup_ranking.R` — `blupRankingUI` / `blupRankingServer` |
| DR helpers | `shiny/common/trait_reduction.R` |
| DR module | `shiny/app/module_trait_reduction.R` |
| Curve bridge | extend Growth Summaries and/or thin adapter used by BLUP/DR |
| Wire | `shiny/app/ui.R`, `shiny/app/server.R` |
| Tukey hint | `shiny/app/reactive_outputs.R` (copy only) |
| X/Y parse | additive in common helper if merged RDS lacks X/Y columns |

Launchers off-limits.

### Phased delivery (all in this project)

1. BLUP tab on plant-collapsed single trait; blocks; treatment empty if 1 level; tests stub/Dunnett when checks exist.
2. Spatial residual stage (per-table X/Y smooths) always when X/Y present; treatment panel for multi-level projects.
3. DR tab (plant × selected traits; plant PC scores → BLUP).
4. Curve-trait bridge polished: Growth Summaries selections feed BLUP/DR as first-class traits.
5. **Out of v1:** joint approach-A model; fancy residual diagnostics UI beyond coarse smooth.

### Debug

- Primary: `project_2025-07-14_b51-91_soy`.
- Confirm: 325 ranks produce; single-level treatment empty panel; spatial runs with 41 tables without freezing the session (timeouts/progress messaging OK).

---

## Left to the builder (explicit)

- Exact `mgcv` smooth specification (basis, `k`, X+Y separate vs `s(X,Y,by=table)`) within coarse per-table requirement.
- Whether `(1|table_no)` stays in the BLUP stage after per-table residualization (avoid double-counting).
- How hyper-trait columns are stored in-session (reactive table vs cbind onto plant-level frame).
- H² formula details when design is unbalanced.
- UI layout details inside tabs (match Effect Sizes / Growth Summaries density).
- Whether BLUE mode ships in phase 1 or right after BLUP-only.

---

## Builder calls (settled here without further Alana questions)

1. **Block, don’t silent-mean** if BLUP requested on uncollapsed multi-time raw trait.
2. **Tab name:** `BLUP`.
3. **DR:** plant rows × selected traits only; never mass-BLUP-then-PCA; never PCA with cultivar/factor levels as rows.
4. **Spatial package default:** `mgcv` for per-table smooths under approach B.
5. **Approach A** remains a documented future option, not this build.
6. **Doc precedence** as stated in the header.
7. **No check list in repo today** → top-k first; Dunnett only if checks already exist in config.
8. **Growth Summaries must honor selected times** for the curve bridge (fix full-series behavior).

### Closed — Alana correction

DR does not run on factors. Cultivar is not an observed phenotype. PCA rows = plants; columns = traits; BLUP then uses plant PC scores with the usual design.

---

## Approval gate

No product code until Alana approves this document. After approval → `writing-plans` → implement phases 1–4 on branch; keep data dirt out of commits.
