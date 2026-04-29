# NucleusGrid

**Multi-Skyrme microscopic predictions for the superheavy stability island.**

A reproducible HFBTHO-based survey of the Mc (Z=115), Nh (Z=113), and Ts (Z=117) isotope chains across the predicted N=184 island, with deformation-corrected ground states, Bayesian model averaging across four Skyrme functionals, and α-decay half-life predictions calibrated against AME2020.

| | |
|---|---|
| **Status** | Phase 3 compute complete (1,512 / 1,520 converged, 99.5%) |
| **Solver** | HFBTHO v3.00 — Navarro Perez et al., *Comput. Phys. Commun.* **220** (2017) 363 |
| **Functionals** | UNE1, SLY4, SKM\*, HFB9 |
| **Calibration** | AME2020 — Wang et al., *Chinese Phys. C* **45** (2021) 030003 |
| **Hardware** | Consumer-grade — AMD Ryzen 9 5950X + Ryzen 7 3800X |
| **Paper draft** | [`docs/NucleusGrid_Paper_Draft.docx`](docs/NucleusGrid_Paper_Draft.docx) |

---

## Coverage

| Element | Z | Mass range | β₂ sweep | Calculations |
|---|---:|---|---|---:|
| Mc | 115 | A = 285–311 | 8-point: −0.20, −0.10, 0.00, 0.10, 0.15, 0.20, 0.25, 0.30 | 528 |
| Nh | 113 | A = 281–307 | 8-point (same grid) | 528 |
| Ts | 117 | A = 291–311 | spherical only (β₂ = 0) — deformation extension planned | 312 |
| Rg / Mt / Bh / Db | 105–111 | chain daughters | β₂ = 0 | 144 |

Total: 1,512 self-consistent HFB calculations across 4 functionals.

## Methodology

For each (Z, A, functional) triple, the ground state is the minimum binding energy across the eight β₂ inputs (deformation correction). Per-functional binding-energy offsets are calibrated against the AME2020 Q_α evaluation; residuals are converted to Bayesian weights via a Gaussian likelihood with σ = 1 MeV, and α-decay half-lives are estimated via the Viola-Seaborg systematics from the weighted Q_α.

The deformation correction is the dominant systematic improvement: cohort RMSE drops from 0.91–1.53 MeV (spherical-basis) to 0.16–0.19 MeV (deformation-corrected), eliminating the +1.4 MeV Q_α bias that all four functionals carry under the spherical assumption.

## Headline findings (preliminary, paper draft)

- **Oblate ground state at N = 179** for Mc-294 and Nh-292: all four functionals agree on β₂ ≈ −0.07 ± 0.03, in contrast to spherical solutions ~1–2 MeV less bound.
- **N = 184 spherical shell closure** corroborated at Mc-299, Nh-297, and Ts-301.
- **Millennial-scale chain endpoints**: the Mc-299 → Db-279 α-cascade terminates at Bh→Db transitions whose deformation-corrected, weighted half-lives reach t½ ≈ 3 × 10³ y (Bh-283 → Db-279). Macroscopic timescales potentially observable through chemical-isolation methods, independent of the prompt α-correlation chain.
- **Functional weight rebalance**: the spherical-basis weighting is dominated by SKM\* (w = 0.46); the deformation-corrected weighting is near-parity (UNE1 = 0.30, SKM\* = 0.26, SLY4 = 0.22, HFB9 = 0.22).

See `docs/NucleusGrid_Paper_Draft.docx` for the full results, tables, and references.

## Repository layout

```
hfbtho/                          HFBTHO v3.00 source + Makefile + binary
hfbtho_wrapper.py                Python wrapper for single-shot runs
HFBTHO_SETUP.md                  Build instructions for the solver

scripts/
  machine1_run.sh                Mc/Nh PES sweep driver (5950X)
  machine2_run.sh                Mc/Nh extended + α-chain driver (3800X)
  tennessine_run.sh              Ts (Z=117) driver
  highprecision_run.sh           N_shell=20, acc=1e-5 refinement runs
  kewpie2_setup.sh               KEWPIE2 fusion-cross-section setup
  synthesis_setup.sh             TALYS synthesis-channel setup

analysis/
  analyze_results.py             Bayesian-weighted analysis (any element)

results/phase3/
  machine1_all.json              Mc/Nh PES sweep raw results (888)
  machine2_all.json              Mc/Nh chains + extended (312)
  tennessine_all.json            Ts raw results (312)
  combined_analysis.json         Full Mc-cohort analysis (spherical baseline)
  mc_analysis.json               Mc-cohort, deformation-corrected
  tennessine_analysis.json       Ts-cohort
  *_output.txt, *.jsonl          run logs and streaming JSONL companions

docs/
  NucleusGrid_Paper_Draft.docx   Phase 3 paper draft
  build_paper.py                 Reproducible builder (reads from results/phase3/)
```

## Reproducing

1. Build HFBTHO v3.00 — see `HFBTHO_SETUP.md`. The binary lands at `hfbtho/hfbtho_main`.
2. Run the drivers (long-running; resume-aware on the JSONL companion files):
   ```bash
   nohup scripts/machine1_run.sh   > results/phase3/machine1_output.txt 2>&1 &
   nohup scripts/machine2_run.sh   > results/phase3/machine2_output.txt 2>&1 &
   nohup scripts/tennessine_run.sh > results/phase3/tennessine_output.txt 2>&1 &
   ```
3. Run the analysis:
   ```bash
   python3 analysis/analyze_results.py \
       --single       results/phase3/machine2_all.json \
       --calibration  results/phase3/machine1_all.json \
       --out          results/phase3/combined_analysis.json
   ```
4. Rebuild the paper draft from the analysis JSONs:
   ```bash
   python3 docs/build_paper.py
   ```

## Open issues

- **KEWPIE2** (2015 release) hits a convergence issue on the ⁴⁸Ca + ²⁴³Am → Mc compound nucleus. Setup is in `scripts/kewpie2_setup.sh`; contact with the authors for an updated release is planned.
- The **Ts deformation sweep** is the principal outstanding extension to this work.

## Citation

Paper draft in preparation. In the meantime, the analysis JSONs in `results/phase3/` are the canonical record. Repository: <https://github.com/DesilatorRx/nucleusgrid>.
