---
status: handoff
audience: implementer
updated: 2026-10-30
related:
  - docs/statfarmer-dr-and-blups-design.md
  - docs/soybean-325-analysis-brief.md
  - internal/statfarmer-feature-seams.md
---

# Handoff: StatFaRmer BLUP + trait DR

Alana lock-ins + implementer brief. Prefer this over open questions in the design doc when they conflict.

## Decisions (locked)

| Topic | Decision |
|-------|----------|
| Goal | **Both** — rank *and* test (ranking table + planned tests / Dunnett-vs-checks / FDR; never all-pairs Tukey for huge genotype). |
| Check IDs | Wizard already validates project IDs upstream — **do not re-invent check-ID UX** in these modules; consume wizard-clean factors. |
| X:Y | **Ignore in v1.** Meaning: do **not** fit table-position residual / spatial smooth on X:Y yet. Still use **T / `table_no` as run/batch block**. (X:Y was only a later optional residual idea — skipped for first cut.) |
| DR default traits | **All** eligible numeric phenotyping traits (minus IDs/meta). User can deselect; default = full set. |
| Time | Follow **user-entered / resolved formula** and existing time path — not a single hard-coded “last cluster.” Prefer feeding DR/BLUP from **curve-derived traits** when available (see below). |
| Treatment | Soy may be **1 level** → empty/disabled treatment panel OK. Architecture must still support **many treatments** (side-by-side BLUPs when ≥2). |
| Tab naming | **“BLUP”** (not “BLUP Ranking”) — covers ranking + test-ish outputs. |

## Features to build

### 1. Tab **BLUP** (sidecar; do not break `resolveModel`)
- Per-trait (then batch) LMM: genotype random → BLUP; treatment fixed/random when levels allow → side panel.
- Blocks: `repeat` + **T/`table_no` as run block** (same physical table; not field spatial).
- Outputs: genotype table + treatment table; SE/rank; H² when estimable; Dunnett/planned tests + FDR; CSV.
- Soft hint from Tukey path when genotype cardinality high → “use BLUP tab.”

Default illustrative form (drop treatment if <2 levels):

```text
y ~ treatment + (1|genotype) + (1|repeat) + (1|table_no)
```

### 2. Trait intelligent DR (sidecar)
- Default: **all** numeric traits → G×P matrix (prefer **BLUP cells**, else BLUE/means).
- PCA (+ trait correlation clustering); scree, loadings, scores, suggested reduced set.
- No AMMI/GGE until multi-env.

## Curve → BLUP inputs (important follow-on)

Time should not force one snapshot. Prefer **generalized timecurve traits** as BLUP/DR inputs when growth modules / bio-object libs exist, e.g.:
- growth-curve parameters (asymp, rate, lag)
- AUC / final size
- time of max growth / **highest derivative**
- other PlantEye / bio-object curve summaries already in repo or easy package adds

**Implementer task:** inventory existing Growth Summaries / curve helpers and R packages already used or addable; expose curve-derived columns as first-class traits into BLUP + DR. Do not invent a second time system parallel to the formula path.

## Seams (do not invade)

- Add `module_blup*.R` + `module_trait_reduction*.R` like Effect Sizes / Growth Summaries.
- Wire `ui.R` / `server.R` only; **no** `resolveModel` rewrite in phase 1–2.
- Soft copy in Tukey skip messaging only.
- Packages: prefer `lme4`, `emmeans`, `psych`; new pkgs OK if justified; PR-ready git (repo already init; add remote when pushing).

## Phased slice

1. **BLUP tab** — genotype BLUPs + run-block T; treatment empty-state if 1 level; tests stub/Dunnett when checks exist via wizard factors.
2. **Treatment side panel** — multi-level projects.
3. **DR tab** — all traits default; means → then BLUP matrix.
4. **Curve-trait bridge** — auto growth/derivative features into BLUP/DR inputs.
5. Optional: X:Y residual diagnostics (explicitly **not** v1).

## Git / PR

- Checkout: `C:\Users\user\StatFaRmer_new\StatFaRmer`
- Git initialized; historically **no remote** — set remote + branch before PR; keep data dirt out of feature commits.
- Launchers (`install.R`, `just_shiny.R`, …) off-limits unless Alana asks.

## Caveman summary

BLUP tab + all-trait DR sidecars. T=run block. Skip X:Y v1. Rank+test. Curve traits feed BLUPs. Treatment panel scales 1→many.
