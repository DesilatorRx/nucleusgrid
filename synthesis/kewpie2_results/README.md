# KEWPIE2 — ⁴⁸Ca + ²⁴³Am for Mc-291..287 synthesis

## Status

| Stage | Status |
|---|---|
| GSL prerequisite | ✅ installed (`libgsl-dev`) |
| Source clone | ✅ from `Sonal-Jaiswal/Research_phy` (final_kewpie2/) |
| Compilation | ✅ `g++ -O3 -ansi`, builds cleanly |
| Reference smoke test (¹⁶O + ²⁰⁸Pb) | ✅ runs in <1 s, produces expected output |
| **⁴⁸Ca + ²⁴³Am SHE calculation** | ❌ **hangs in this environment** |

## What works

The shipped sample input — ¹⁶O on ²⁰⁸Pb forming ²²⁴Th\* — runs in **under 1 second** and emits 21 energy points × 5 channels (0n..4n) of fusion + evaporation residue cross sections. The
output starts with the "Compound nucleus:" diagnostic block and matches the
shipped reference `Results/cross-sections.dat`.

So the build and the input-parsing path are correct; the binary is functional
for the calibration/sample case.

## What doesn't work — the SHE case

For the requested reaction `⁴⁸Ca + ²⁴³Am → ²⁹¹Mc*`, KEWPIE2 enters the
`InitTable` / `Element::Nucleus` / fusion-barrier setup path and **never
prints the "Compound nucleus" diagnostic**, regardless of how aggressively the
input is reduced. Tested configurations that all hang for >5 minutes producing
zero bytes of stdout (and still at 100%+ CPU when killed):

| L_max | E\* range | E\* step | Max neutrons | Fusion model | Result |
|---|---|---|---|---|---|
| 199 | 30–55 MeV | 0.1 | 4 | EBD | hung at ≥22 min |
| 80  | 30–55 MeV | 0.1 | 4 | EBD | hung at ≥22 min |
| 40  | 30–45 MeV | 0.5 | 4 | EBD | hung at ≥22 min |
| 30  | 35–38 MeV | 1.0 | 3 | WKB | hung at ≥3 min |
| 30  | 35–36 MeV | 1.0 | 1 | WKB | hung at ≥5.5 min |

Memory was not exhausted (peak RSS <10 MB), CPU was pinned at ~110% (single
core). The mass table (`Input/MassTab/MollerNix.dat`) does contain Mc-291
through Mc-272, Am-243, and Ca-48, so the table-lookup-then-cin.get() trap in
`Element.cpp` is unlikely to be the issue. Feeding the process unlimited
Enters via `yes "" |` did not unblock it, also consistent with the hang
being in actual computation rather than input.

The most plausible diagnosis (not verified) is a numerical pathology in the
SHE region — perhaps the partial-wave / WKB integration in `Transcoef.cpp` or
the deeply sub-barrier fusion (E\* ≈ 30 MeV corresponds to E_cm ≈ 32 MeV
*below* the Bass barrier of ~232 MeV for Ca+Am) takes an enormous number of
integration steps in the Schrödinger-equation solver.

## Files in this directory

- `Ca48_Am243_Mc291.inp` — the production input deck (full E\* sweep,
  L_max=80, max 4 neutrons, EBD fusion). Drop into
  `synthesis/kewpie2/Input/input.inp` and run.
- `README.md` — this file.

## Next steps (if pursuing further)

1. **Try a different KEWPIE2 distribution** — the GitHub mirror used here is
   the September 2015 CPC release. There may be later versions in the
   nuclear physics community with SHE-region fixes; ask the original authors
   (H. Lü, A. Boilley, Y. Abe).
2. **Add timing/progress prints to `Cascade.cpp` / `Transcoef.cpp`** —
   recompile and run to identify which routine is the bottleneck. The hang
   is BEFORE the "Compound nucleus" print in `main.cpp`, so it's in
   `Table::InitTable` or `Element::Nucleus` or one of their dependencies.
3. **Run on a server without WSL2** — even if KEWPIE2 needed 1–2 hours, the
   WSL2 host hibernated mid-run on this laptop, killing background processes
   silently. A non-hibernating Linux box would be more reliable.
4. **Cross-check with a different SHE code** — NRV (web-based,
   http://nrv.jinr.ru/) has a free fusion-evaporation calculator for SHE
   that can serve as a reference for what KEWPIE2 *should* produce here.
   Experimental: Oganessian 2003 reported ~3 pb for the 3n channel of this
   exact reaction.

## Citation when this works

Lü, A. Boilley, Y. Abe, et al., *KEWPIE2: A cascade code for the study of
dynamical decay of excited nuclei*, Comput. Phys. Commun. 200 (2016) 381.
