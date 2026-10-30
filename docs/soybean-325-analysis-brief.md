---
cursor:
  subagentId: "bc-f2d53306-9ad8-5cfc-a28e-bc627b8129ed"
---

# 325 soybean cultivars — analysis brief

**Verdict:** Do not rely on classical all-pairs ANOVA + Tukey for 325 cultivars. Rank with mixed-model BLUP/BLUE; use Dunnett/FDR only for planned contrasts.

## Local experiment shape (counts only)
Project `project_2025-07-14_b51-91_soy`: 325 cultivars × 3 reps (975 plants), one Treatment level, ~4 time clusters. Single-environment screen, not a multi-location MET.

## StatFaRmer today
Cases = ANOVA / LMM (time blocks) / spline + effect-sizes + growth modules. Auto path falls back to classic ANOVA when any factor has >20 levels. Cultivar is not treated as a random genetic effect. Pairwise Tukey is filtered above 1000 comparisons — useless as a 325-way letter plot.

## Recommended stack
1. QC + descriptives (wizard/app).
2. Per-trait LMM: genotype random (BLUP rank) or fixed (BLUE); include rep/block if present; report H² when estimable.
3. Checks: Dunnett; otherwise top-k. FDR on planned contrasts only.
4. Multivariate clustering / PCA on cultivar×trait matrix; growth/spline for time.
5. AMMI/GGE only after multi-env data exist.

## Need from Alana
Goal (test vs rank vs select); check cultivars; CRD/RCBD/lattice + spatial?; priority traits; multi-env plans; confirm single Treatment.

## Next StatFaRmer case (if approved)
`multi_genotype_ranking`: BLUP/BLUE table, H², top-k, Dunnett-vs-checks; auto-suppress all-pairs Tukey when genotype levels exceed ~20.

Full notes: `internal/soybean-325-analysis-proposal.md`.
