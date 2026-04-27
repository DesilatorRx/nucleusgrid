#!/usr/bin/env python3
"""NucleusGrid Phase 3 — Machine 2 analysis.

Loads HFBTHO results, finds ground states per (Z, A, functional), applies
per-functional calibration offsets, derives Bayesian model weights from RMSE
against experimental Q_alpha for Mc-286..290, and prints:
  - calibration RMSE & weights
  - full stability surface (BE per functional, weighted average, spread)
  - per-chain Q_alpha tables (Mc-291..299, 5 alpha steps each)
  - model agreement / disagreement summary
  - N=184 magic and N=179 minimum highlights
The full structured analysis is dumped to JSON.
"""
from __future__ import annotations
import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

E_ALPHA = 28.2957  # MeV

CALIBRATION_OFFSETS = {
    "UNE1": 57.146,
    "SLY4": 57.091,
    "SKM*": 55.935,
    "HFB9": 56.850,
}

# Experimental Q_alpha for the calibration cohort (MeV)
EXP_Q_ALPHA = {286: 10.33, 287: 10.59, 288: 10.01, 289: 9.95, 290: 9.78}

SYMBOL = {115: "Mc", 113: "Nh", 111: "Rg", 109: "Mt", 107: "Bh", 105: "Db"}

FUNCTIONALS = ["UNE1", "SLY4", "SKM*", "HFB9"]

N_MAGIC = 184
N_MIN = 179
AGREEMENT_THRESHOLD_MEV = 0.5


def load_runs(path: Path) -> list[dict]:
    with open(path) as f:
        return json.load(f)


def ground_states(runs: list[dict]) -> dict[tuple[int, int, str], dict]:
    grouped: dict[tuple[int, int, str], list[dict]] = defaultdict(list)
    for r in runs:
        if r.get("status") != "converged" or r.get("BE_mev") in (None, "null"):
            continue
        grouped[(r["Z"], r["A"], r["functional"])].append(r)
    return {k: min(v, key=lambda r: float(r["BE_mev"])) for k, v in grouped.items()}


def calibrate(machine1_gs: dict) -> dict[str, tuple[float, list]]:
    out = {}
    for f in FUNCTIONALS:
        residuals = []
        for A_mc, q_exp in EXP_Q_ALPHA.items():
            kp, kd = (115, A_mc, f), (113, A_mc - 4, f)
            if kp not in machine1_gs or kd not in machine1_gs:
                continue
            be_p = float(machine1_gs[kp]["BE_mev"])
            be_d = float(machine1_gs[kd]["BE_mev"])
            q_pred = be_p - be_d - E_ALPHA + CALIBRATION_OFFSETS[f]
            residuals.append((A_mc, q_pred, q_exp, q_pred - q_exp))
        if residuals:
            rmse = math.sqrt(sum(d * d for _, _, _, d in residuals) / len(residuals))
            out[f] = (rmse, residuals)
    return out


def bayesian_weights(rmse_data: dict) -> dict[str, float]:
    if not rmse_data:
        return {f: 1.0 / len(FUNCTIONALS) for f in FUNCTIONALS}
    raw = {f: 1.0 / (v[0] ** 2) for f, v in rmse_data.items()}
    s = sum(raw.values())
    weights = {f: raw[f] / s for f in raw}
    for f in FUNCTIONALS:
        weights.setdefault(f, 0.0)
    return weights


def calibrated_q(gs, Z_p, A_p, func) -> float | None:
    p = gs.get((Z_p, A_p, func))
    d = gs.get((Z_p - 2, A_p - 4, func))
    if p is None or d is None:
        return None
    return float(p["BE_mev"]) - float(d["BE_mev"]) - E_ALPHA + CALIBRATION_OFFSETS[func]


def viola_seaborg_log(Q: float | None, Z_parent: int, N_parent: int) -> float | None:
    """Viola-Seaborg 1966 with Sobiczewski hindrance for odd-Z / odd-N parents.
    log10(T_1/2 / s) = (aZ + b)/sqrt(Q) + (cZ + d) + h(parity)."""
    if Q is None or Q <= 0:
        return None
    a, b, c, d = 1.66175, -8.5166, -0.20228, -33.9069
    log_t = (a * Z_parent + b) / math.sqrt(Q) + (c * Z_parent + d)
    odd_z, odd_n = Z_parent % 2, N_parent % 2
    if odd_z and odd_n:
        log_t += 0.3522
    elif odd_z:
        log_t += 0.2392
    elif odd_n:
        log_t += 0.1130
    return log_t


def fmt_halflife(log_t: float | None) -> str:
    if log_t is None:
        return "n/a"
    if log_t < -12:
        return "<1 ps"
    if log_t < -9:
        return f"{10**(log_t+12):.1f} ps"
    if log_t < -6:
        return f"{10**(log_t+9):.1f} ns"
    if log_t < -3:
        return f"{10**(log_t+6):.1f} µs"
    if log_t < 0:
        return f"{10**(log_t+3):.1f} ms"
    if log_t < math.log10(60):
        return f"{10**log_t:.1f} s"
    if log_t < math.log10(3600):
        return f"{10**log_t/60:.1f} min"
    if log_t < math.log10(86400):
        return f"{10**log_t/3600:.1f} h"
    if log_t < math.log10(86400 * 365):
        return f"{10**log_t/86400:.1f} d"
    return f"{10**log_t/86400/365:.2g} y"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--single", required=True, help="Machine 2 results JSON")
    p.add_argument("--calibration", default=None,
                   help="Path to Machine 1 results JSON for calibration "
                        "(default: <repo>/nucleusgrid_results.json)")
    p.add_argument("--out", default=None,
                   help="Output JSON path (default: <input dir>/machine2_analysis.json)")
    args = p.parse_args()

    src = Path(args.single).resolve()
    runs = load_runs(src)
    gs = ground_states(runs)

    cal_path = Path(args.calibration).resolve() if args.calibration else (
        src.parent.parent.parent / "nucleusgrid_results.json").resolve()
    rmse_data = {}
    if cal_path.exists():
        cal_runs = load_runs(cal_path)
        rmse_data = calibrate(ground_states(cal_runs))
    weights = bayesian_weights(rmse_data)

    n_failed = sum(1 for r in runs if r.get("status") != "converged")

    print(f"\n{' NucleusGrid Phase 3 — Machine 2 Analysis ':=^82}")
    print(f"  Input:        {src}")
    print(f"  Calibration:  {cal_path if cal_path.exists() else '— (no source — uniform weights)'}")
    print(f"  Runs:         {len(runs)}  ({n_failed} failed)")
    print(f"  Functionals:  {', '.join(FUNCTIONALS)}")
    print(f"  Offsets:      " + ", ".join(f"{f}={CALIBRATION_OFFSETS[f]}" for f in FUNCTIONALS))

    if rmse_data:
        print(f"\n--- Calibration RMSE vs experimental Q_α (Mc-286..290) ---")
        print(f"  exp: " + ", ".join(f"Mc-{A}={q}" for A, q in EXP_Q_ALPHA.items()) + "  MeV")
        print(f"\n  {'functional':<10} {'RMSE (MeV)':>12} {'weight':>10}   residuals (Q_pred − Q_exp, MeV)")
        for f in FUNCTIONALS:
            if f not in rmse_data:
                print(f"  {f:<10} {'—':>12} {'—':>10}   (insufficient calibration data)")
                continue
            rmse, resid = rmse_data[f]
            res_str = " ".join(f"{d:+.3f}" for _, _, _, d in resid)
            print(f"  {f:<10} {rmse:>12.4f} {weights[f]:>10.4f}   {res_str}")
    else:
        print(f"\n  (no calibration source found — using uniform weights = 0.25)")

    # ---- Stability surface ----
    print(f"\n{' STABILITY SURFACE: BE per (Z, A) per functional ':=^82}")
    print(f"\n  {'sym':<3} {'Z':>3} {'A':>4} {'N':>4}  " +
          "  ".join(f"{f:>11}" for f in FUNCTIONALS) +
          f"  {'<BE>_w':>11}  {'spread':>7}  flags")
    print("  " + "-" * 100)
    surface = []
    all_ZA = sorted({(Z, A) for (Z, A, _) in gs}, key=lambda p: (-p[0], p[1]))
    for (Z, A) in all_ZA:
        N = A - Z
        sym = SYMBOL.get(Z, f"Z{Z}")
        be = {f: float(gs[(Z, A, f)]["BE_mev"]) if (Z, A, f) in gs else None
              for f in FUNCTIONALS}
        present = [(f, v) for f, v in be.items() if v is not None]
        if not present:
            continue
        wsum = sum(weights[f] for f, _ in present)
        be_avg = sum(v * weights[f] for f, v in present) / wsum
        spread = max(v for _, v in present) - min(v for _, v in present)
        flags = []
        if N == N_MAGIC:
            flags.append(f"N={N_MAGIC} MAGIC")
        if N == N_MIN:
            flags.append(f"N={N_MIN} min")
        be_strs = "  ".join(f"{be[f]:>11.3f}" if be[f] is not None else "          —"
                            for f in FUNCTIONALS)
        print(f"  {sym:<3} {Z:>3} {A:>4} {N:>4}  {be_strs}  "
              f"{be_avg:>11.3f}  {spread:>7.3f}  {' / '.join(flags)}")
        surface.append({
            "Z": Z, "A": A, "N": N, "symbol": sym,
            "BE_per_functional": be,
            "BE_weighted": round(be_avg, 4),
            "BE_spread": round(spread, 4),
            "is_magic_184": N == N_MAGIC,
            "is_minimum_179": N == N_MIN,
        })

    # ---- Decay chains ----
    print(f"\n{' DECAY CHAINS: Q_α (calibrated) per step, 5 α steps per chain ':=^82}")
    chains = {}
    for A0 in range(291, 300):
        chain_label = f"Mc-{A0}"
        print(f"\n  ─── Chain {chain_label} (Mc-{A0} → Db-{A0-20}) ───")
        print(f"  {'parent':<10} → {'daughter':<10} {'N_p':>4}  " +
              "  ".join(f"{f:>9}" for f in FUNCTIONALS) +
              f"  {'Q_w':>8}  {'spread':>7}  {'t½(w)':>10}  flags")
        steps = []
        for k in range(5):
            Z_p, A_p = 115 - 2 * k, A0 - 4 * k
            label_p = f"{SYMBOL[Z_p]}-{A_p}"
            label_d = f"{SYMBOL[Z_p-2]}-{A_p-4}"
            N_p = A_p - Z_p
            qs = {f: calibrated_q(gs, Z_p, A_p, f) for f in FUNCTIONALS}
            present = [(f, q) for f, q in qs.items() if q is not None]
            if present:
                wsum = sum(weights[f] for f, _ in present)
                q_w = sum(q * weights[f] for f, q in present) / wsum
                spread = max(q for _, q in present) - min(q for _, q in present)
                t_log = viola_seaborg_log(q_w, Z_p, N_p)
            else:
                q_w = spread = t_log = None
            flags = []
            if N_p == N_MAGIC:
                flags.append("N=184")
            if N_p == N_MIN:
                flags.append("N=179")
            if (A_p - 4) - (Z_p - 2) == N_MAGIC:
                flags.append("d:N=184")
            q_strs = "  ".join(f"{qs[f]:>9.4f}" if qs[f] is not None else "        —"
                               for f in FUNCTIONALS)
            qw_s = f"{q_w:>8.4f}" if q_w is not None else "       —"
            sp_s = f"{spread:>7.3f}" if spread is not None else "      —"
            t_s = f"{fmt_halflife(t_log):>10}"
            print(f"  {label_p:<10} → {label_d:<10} {N_p:>4}  {q_strs}  "
                  f"{qw_s}  {sp_s}  {t_s}  {' '.join(flags)}")
            steps.append({
                "step": k,
                "parent": label_p, "daughter": label_d,
                "Z_parent": Z_p, "A_parent": A_p, "N_parent": N_p,
                "Q_per_functional": {f: round(qs[f], 4) if qs[f] is not None else None
                                     for f in FUNCTIONALS},
                "Q_weighted": round(q_w, 4) if q_w is not None else None,
                "spread_mev": round(spread, 4) if spread is not None else None,
                "halflife_log10_s": round(t_log, 4) if t_log is not None else None,
                "halflife_human": fmt_halflife(t_log),
                "flags": flags,
            })
        chains[chain_label] = {
            "members": [s["parent"] for s in steps] + [steps[-1]["daughter"]],
            "steps": steps,
        }

    # ---- Model agreement ----
    print(f"\n{' MODEL AGREEMENT (Q_α spread across functionals) ':=^82}")
    agreeing, disagreeing = [], []
    for ch in chains.values():
        for step in ch["steps"]:
            sp = step["spread_mev"]
            if sp is None:
                continue
            entry = {"transition": f"{step['parent']} → {step['daughter']}",
                     "spread_mev": sp,
                     "Q_per_functional": step["Q_per_functional"]}
            (agreeing if sp <= AGREEMENT_THRESHOLD_MEV else disagreeing).append(entry)
    agreeing.sort(key=lambda e: e["spread_mev"])
    disagreeing.sort(key=lambda e: -e["spread_mev"])
    print(f"  agreement threshold: {AGREEMENT_THRESHOLD_MEV} MeV")
    print(f"  agreeing transitions:    {len(agreeing)}")
    print(f"  disagreeing transitions: {len(disagreeing)}")
    if disagreeing:
        print(f"\n  Top 5 DISAGREEMENTS:")
        for e in disagreeing[:5]:
            qs = ", ".join(f"{f}={v}" for f, v in e["Q_per_functional"].items() if v is not None)
            print(f"    {e['transition']:<27} spread={e['spread_mev']:.3f} MeV  [{qs}]")
    if agreeing:
        print(f"\n  Top 5 AGREEMENTS:")
        for e in agreeing[:5]:
            qs = ", ".join(f"{f}={v}" for f, v in e["Q_per_functional"].items() if v is not None)
            print(f"    {e['transition']:<27} spread={e['spread_mev']:.3f} MeV  [{qs}]")

    # ---- Highlights ----
    print(f"\n{' KEY HIGHLIGHTS ':=^82}")
    n184 = [s for s in surface if s["is_magic_184"]]
    n179 = [s for s in surface if s["is_minimum_179"]]
    if n184:
        print(f"\n  N=184 magic-shell isotopes ({len(n184)}):")
        for s in n184:
            print(f"    {s['symbol']}-{s['A']:<4} (Z={s['Z']:>3})   "
                  f"<BE>_w = {s['BE_weighted']:>11.3f} MeV   spread = {s['BE_spread']:.3f}")
    if n179:
        print(f"\n  N=179 minimum-stability isotopes ({len(n179)}):")
        for s in n179:
            print(f"    {s['symbol']}-{s['A']:<4} (Z={s['Z']:>3})   "
                  f"<BE>_w = {s['BE_weighted']:>11.3f} MeV   spread = {s['BE_spread']:.3f}")

    all_steps = [s for ch in chains.values() for s in ch["steps"]
                 if s["halflife_log10_s"] is not None]
    if all_steps:
        all_steps.sort(key=lambda s: -s["halflife_log10_s"])
        print(f"\n  Longest predicted α half-lives (top 5, weighted):")
        for s in all_steps[:5]:
            flag = f"  [{' '.join(s['flags'])}]" if s["flags"] else ""
            print(f"    {s['parent']:<8} → {s['daughter']:<8}  Q_w = {s['Q_weighted']:>7.4f} MeV  "
                  f"t½ ≈ {s['halflife_human']}{flag}")
        print(f"\n  Shortest predicted α half-lives (top 5, weighted):")
        for s in all_steps[-5:][::-1]:
            flag = f"  [{' '.join(s['flags'])}]" if s["flags"] else ""
            print(f"    {s['parent']:<8} → {s['daughter']:<8}  Q_w = {s['Q_weighted']:>7.4f} MeV  "
                  f"t½ ≈ {s['halflife_human']}{flag}")

    # ---- Save JSON ----
    out_path = Path(args.out).resolve() if args.out else (src.parent / "machine2_analysis.json")
    analysis = {
        "metadata": {
            "input_file": str(src),
            "calibration_source": str(cal_path) if cal_path.exists() else None,
            "calibration_offsets": CALIBRATION_OFFSETS,
            "experimental_q_alpha_mev": {f"Mc-{A}": q for A, q in EXP_Q_ALPHA.items()},
            "rmse_per_functional": {f: round(v[0], 4) for f, v in rmse_data.items()},
            "bayesian_weights": {f: round(weights[f], 4) for f in FUNCTIONALS},
            "n_runs": len(runs),
            "n_failed": n_failed,
            "agreement_threshold_mev": AGREEMENT_THRESHOLD_MEV,
        },
        "calibration_residuals": {
            f: [{"A": A, "Q_pred": round(qp, 4), "Q_exp": qe, "residual": round(d, 4)}
                for (A, qp, qe, d) in v[1]]
            for f, v in rmse_data.items()
        },
        "stability_surface": surface,
        "decay_chains": chains,
        "model_agreement": {
            "agreement_threshold_mev": AGREEMENT_THRESHOLD_MEV,
            "n_agreeing": len(agreeing),
            "n_disagreeing": len(disagreeing),
            "top_disagreements": disagreeing[:10],
            "top_agreements": agreeing[:10],
        },
    }
    with open(out_path, "w") as f:
        json.dump(analysis, f, indent=2)
    print(f"\n  Analysis JSON saved → {out_path}")


if __name__ == "__main__":
    main()
