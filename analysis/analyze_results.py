#!/usr/bin/env python3
"""NucleusGrid Phase 3 — calibrated SHE chain analysis.

Loads HFBTHO results, picks ground states per (Z, A, functional), derives
per-functional calibration offsets from a user-specified experimental cohort,
applies inverse-variance Bayesian model weighting, and emits the full
stability surface + per-chain Q_α tables for the requested parent element.

Defaults reproduce the Mc (Z=115) Machine-2 analysis. Override with:
  --parent-Z, --cohort-A, --cohort-Q, --chain-A, --cohort-source
to analyze a different element (e.g. Ts at Z=117).
"""
from __future__ import annotations
import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

E_ALPHA = 28.2957  # MeV

# Z → element symbol for the SHN region (Z=104..118).
SYMBOL_FULL = {
    104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs",
    109: "Mt", 110: "Ds", 111: "Rg", 112: "Cn", 113: "Nh",
    114: "Fl", 115: "Mc", 116: "Lv", 117: "Ts", 118: "Og",
}

FUNCTIONALS = ["UNE1", "SLY4", "SKM*", "HFB9"]

N_MAGIC = 184
N_MIN = 179
AGREEMENT_THRESHOLD_MEV = 0.5


def load_runs(path):
    with open(path) as f:
        return json.load(f)


def ground_states(runs):
    grouped = defaultdict(list)
    for r in runs:
        if r.get("status") != "converged" or r.get("BE_mev") in (None, "null"):
            continue
        grouped[(r["Z"], r["A"], r["functional"])].append(r)
    return {k: min(v, key=lambda r: float(r["BE_mev"])) for k, v in grouped.items()}


def derive_offsets(cohort_q, parent_Z, gs):
    """For each functional, offset = mean(Q_exp − (BE_p − BE_d − E_α))
    over cohort points where both parent and daughter are present."""
    out = {}
    for f in FUNCTIONALS:
        diffs = []
        for A_p, q_exp in cohort_q.items():
            kp = (parent_Z, A_p, f)
            kd = (parent_Z - 2, A_p - 4, f)
            if kp in gs and kd in gs:
                be_p = float(gs[kp]["BE_mev"])
                be_d = float(gs[kd]["BE_mev"])
                diffs.append(q_exp - (be_p - be_d - E_ALPHA))
        if diffs:
            out[f] = sum(diffs) / len(diffs)
    return out


def calibrate_residuals(cohort_q, parent_Z, gs, offsets):
    out = {}
    for f in FUNCTIONALS:
        residuals = []
        for A_p, q_exp in cohort_q.items():
            kp, kd = (parent_Z, A_p, f), (parent_Z - 2, A_p - 4, f)
            if kp not in gs or kd not in gs or f not in offsets:
                continue
            be_p = float(gs[kp]["BE_mev"])
            be_d = float(gs[kd]["BE_mev"])
            q_pred = be_p - be_d - E_ALPHA + offsets[f]
            residuals.append((A_p, q_pred, q_exp, q_pred - q_exp))
        if residuals:
            rmse = math.sqrt(sum(d * d for _, _, _, d in residuals) / len(residuals))
            out[f] = (rmse, residuals)
    return out


def bayesian_weights(rmse_data):
    if not rmse_data:
        return {f: 1.0 / len(FUNCTIONALS) for f in FUNCTIONALS}
    # Floor at 0.01 MeV so a single-cohort-point RMSE of 0 doesn't blow up.
    raw = {f: 1.0 / max(v[0], 0.01) ** 2 for f, v in rmse_data.items()}
    s = sum(raw.values())
    weights = {f: raw[f] / s for f in raw}
    for f in FUNCTIONALS:
        weights.setdefault(f, 0.0)
    return weights


def calibrated_q(gs, Z_p, A_p, func, offsets):
    p = gs.get((Z_p, A_p, func))
    d = gs.get((Z_p - 2, A_p - 4, func))
    if p is None or d is None or func not in offsets:
        return None
    return float(p["BE_mev"]) - float(d["BE_mev"]) - E_ALPHA + offsets[func]


def viola_seaborg_log(Q, Z_parent, N_parent):
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


def fmt_halflife(log_t):
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


def parse_int_list(s):
    return [int(x) for x in s.split(",")]


def parse_float_list(s):
    return [float(x) for x in s.split(",")]


def parse_a_range(s):
    if "-" in s and "," not in s:
        a, b = s.split("-")
        return list(range(int(a), int(b) + 1))
    return parse_int_list(s)


def find_cohort_source(args, gs):
    """If --cohort-source is given, use it. Else if cohort BE values are
    available in --single, use --single. Else fall back to the project-root
    nucleusgrid_results.json (preserves prior Mc behaviour)."""
    if args.cohort_source:
        return Path(args.cohort_source).resolve()
    cohort_A = parse_int_list(args.cohort_A)
    f0 = FUNCTIONALS[0]
    needed = [(args.parent_Z, A, f0) for A in cohort_A] + \
             [(args.parent_Z - 2, A - 4, f0) for A in cohort_A]
    if all(n in gs for n in needed):
        return Path(args.single).resolve()
    fallback = Path(args.single).resolve().parent.parent.parent / "nucleusgrid_results.json"
    if fallback.exists():
        return fallback
    return Path(args.single).resolve()


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--single", required=True,
                   help="Primary HFBTHO results JSON to analyse")
    p.add_argument("--parent-Z", type=int, default=115,
                   help="Top of decay chain (default 115=Mc; 117=Ts; 113=Nh; etc.)")
    p.add_argument("--cohort-A", default="286,287,288,289,290",
                   help="Comma-sep parent A values for the experimental cohort")
    p.add_argument("--cohort-Q", default="10.33,10.59,10.01,9.95,9.78",
                   help="Comma-sep experimental Q_α (MeV), one per --cohort-A")
    p.add_argument("--chain-A", default="291-299",
                   help="Range or comma-sep parent A values for decay chains")
    p.add_argument("--cohort-source", default=None,
                   help="JSON to read cohort BE values from "
                        "(default: --single if it has them, else "
                        "<repo>/nucleusgrid_results.json)")
    p.add_argument("--out", default=None,
                   help="Output JSON (default: <input dir>/<symbol>_analysis.json)")
    args = p.parse_args()

    src = Path(args.single).resolve()
    runs = load_runs(src)
    gs = ground_states(runs)

    cohort_A = parse_int_list(args.cohort_A)
    cohort_Q = parse_float_list(args.cohort_Q)
    if len(cohort_A) != len(cohort_Q):
        print(f"ERROR: --cohort-A has {len(cohort_A)} entries but --cohort-Q has {len(cohort_Q)}",
              file=sys.stderr)
        sys.exit(1)
    cohort = dict(zip(cohort_A, cohort_Q))
    chain_A = parse_a_range(args.chain_A)

    cal_path = find_cohort_source(args, gs)
    cohort_gs = gs if cal_path == src else ground_states(load_runs(cal_path))

    offsets = derive_offsets(cohort, args.parent_Z, cohort_gs)
    rmse_data = calibrate_residuals(cohort, args.parent_Z, cohort_gs, offsets)
    weights = bayesian_weights(rmse_data)

    z_set = sorted({z for (z, _, _) in gs}, reverse=True)
    symbols = {z: SYMBOL_FULL.get(z, f"Z{z}") for z in z_set}
    parent_sym = symbols.get(args.parent_Z, f"Z{args.parent_Z}")

    n_failed = sum(1 for r in runs if r.get("status") != "converged")

    print(f"\n{f' NucleusGrid Analysis — {parent_sym} (Z={args.parent_Z}) ':=^82}")
    print(f"  Input:        {src}")
    print(f"  Cohort src:   {cal_path}")
    print(f"  Runs:         {len(runs)}  ({n_failed} failed)")
    print(f"  Functionals:  {', '.join(FUNCTIONALS)}")
    print(f"  Z range:      {z_set[-1]}..{z_set[0]} ({', '.join(symbols[z] for z in z_set)})")
    print(f"  Cohort:       {parent_sym}-" +
          "/".join(str(a) for a in cohort_A) + " (Q_exp = " +
          ", ".join(f"{q}" for q in cohort_Q) + " MeV)")
    print(f"  Chain heads:  {parent_sym}-{chain_A[0]}..{parent_sym}-{chain_A[-1]} ({len(chain_A)} chains)")

    print(f"\n--- Derived per-functional calibration offsets ---")
    for f in FUNCTIONALS:
        if f in offsets:
            print(f"  {f:<10} offset = {offsets[f]:>9.4f} MeV")
        else:
            print(f"  {f:<10} (insufficient cohort data)")

    if rmse_data:
        print(f"\n--- Calibration RMSE (Q_pred − Q_exp after offset) ---")
        for f in FUNCTIONALS:
            if f not in rmse_data:
                continue
            rmse, resid = rmse_data[f]
            res_str = " ".join(f"{d:+.3f}" for _, _, _, d in resid)
            print(f"  {f:<10} RMSE={rmse:.4f}  weight={weights[f]:.4f}   residuals: {res_str}")

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
        sym = symbols.get(Z, f"Z{Z}")
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
    print(f"\n{f' DECAY CHAINS: Q_α (calibrated) per step ':=^82}")
    chains = {}
    for A0 in chain_A:
        chain_label = f"{parent_sym}-{A0}"
        # How deep does data let us go?
        max_steps = 0
        while True:
            Z_check = args.parent_Z - 2 * (max_steps + 1)
            A_check = A0 - 4 * (max_steps + 1)
            if Z_check < 1 or A_check < 1:
                break
            if not any((Z_check, A_check, f) in gs for f in FUNCTIONALS):
                break
            max_steps += 1
        if max_steps == 0:
            continue
        end_sym = symbols.get(args.parent_Z - 2 * max_steps,
                              f"Z{args.parent_Z - 2*max_steps}")
        end_A = A0 - 4 * max_steps
        print(f"\n  ─── Chain {chain_label} ({chain_label} → {end_sym}-{end_A}) ───")
        print(f"  {'parent':<10} → {'daughter':<10} {'N_p':>4}  " +
              "  ".join(f"{f:>9}" for f in FUNCTIONALS) +
              f"  {'Q_w':>8}  {'spread':>7}  {'t½(w)':>10}  flags")
        steps = []
        for k in range(max_steps):
            Z_p = args.parent_Z - 2 * k
            A_p = A0 - 4 * k
            label_p = f"{symbols.get(Z_p, f'Z{Z_p}')}-{A_p}"
            label_d = f"{symbols.get(Z_p-2, f'Z{Z_p-2}')}-{A_p-4}"
            N_p = A_p - Z_p
            qs = {f: calibrated_q(gs, Z_p, A_p, f, offsets) for f in FUNCTIONALS}
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

    # Q_α minimum across chains (per-Z, scanning N)
    print(f"\n  Q_α minima across chains (where each chain's Q_w is smallest):")
    chain_min = []
    for label, ch in chains.items():
        steps_with_q = [s for s in ch["steps"] if s["Q_weighted"] is not None]
        if not steps_with_q:
            continue
        m = min(steps_with_q, key=lambda s: s["Q_weighted"])
        chain_min.append((label, m))
    chain_min.sort(key=lambda x: x[1]["Q_weighted"])
    for label, s in chain_min[:5]:
        print(f"    {label:<8}   min Q_w = {s['Q_weighted']:>7.4f} MeV at "
              f"{s['parent']:<8} (N_p={s['N_parent']:>3})  "
              f"t½ ≈ {s['halflife_human']}")

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

    out_path = Path(args.out).resolve() if args.out else (
        src.parent / f"{parent_sym.lower()}_analysis.json")
    analysis = {
        "metadata": {
            "input_file": str(src),
            "cohort_source": str(cal_path),
            "parent_symbol": parent_sym,
            "parent_Z": args.parent_Z,
            "experimental_cohort": {f"{parent_sym}-{a}": q for a, q in cohort.items()},
            "calibration_offsets": {f: round(offsets[f], 4) for f in offsets},
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
