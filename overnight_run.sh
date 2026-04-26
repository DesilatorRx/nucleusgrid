#!/bin/bash
# ============================================================
# NucleusGrid Phase 2 — Overnight Batch Runner
# Runs all Mc-285 to Mc-303 + daughter nuclei
# All 4 functionals: UNE1, SLY4, SKM*, HFB9
# Results saved to nucleusgrid_results.log and results.json
# ============================================================
# Usage: bash overnight_run.sh
# Takes: 4-8 hours on 5950X
# ============================================================

HFBTHO="./hfbtho_main"
LOGFILE="../nucleusgrid_results.log"
JSONFILE="../nucleusgrid_results.json"
STARTTIME=$(date)

echo "============================================================" | tee -a $LOGFILE
echo "NucleusGrid Phase 2 Overnight Run" | tee -a $LOGFILE
echo "Started: $STARTTIME" | tee -a $LOGFILE
echo "============================================================" | tee -a $LOGFILE

# Initialize JSON
echo "[" > $JSONFILE

FIRST_ENTRY=1

# ── Function: run one HFBTHO calculation ──
run_hfbtho() {
    local Z=$1
    local N=$2
    local FUNCTIONAL=$3
    local BETA2=$4
    local LABEL=$5
    local A=$((Z + N))

    cat > hfbtho_NAMELIST.dat << NAMELIST
&HFBTHO_GENERAL
 number_of_shells = 20, oscillator_length = -1.0, basis_deformation = 0.0,
 proton_number = $Z, neutron_number = $N, type_of_calculation = 1 /
&HFBTHO_INITIAL
 beta2_deformation = $BETA2, beta4_deformation = 0.0 /
&HFBTHO_ITERATIONS
 number_iterations = 600, accuracy = 1.E-5, restart_file = 1 /
&HFBTHO_FUNCTIONAL
 functional = '$FUNCTIONAL', add_initial_pairing = F, type_of_coulomb = 2 /
&HFBTHO_PAIRING
 user_pairing = F, vpair_n = -250.0, vpair_p = -250.0,
 pairing_cutoff = 60.0, pairing_feature = 0.5 /
&HFBTHO_CONSTRAINTS
 lambda_values = 1, 2, 3, 4, 5, 6, 7, 8,
 lambda_active = 0, 0, 0, 0, 0, 0, 0, 0,
 expectation_values = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /
&HFBTHO_BLOCKING
 proton_blocking = 1, 0, 0, 0, 0, neutron_blocking = 0, 0, 0, 0, 0 /
&HFBTHO_PROJECTION
 switch_to_THO = 0, projection_is_on = 0,
 gauge_points = 1, delta_Z = 0, delta_N = 0 /
&HFBTHO_TEMPERATURE
 set_temperature = F, temperature = 0.0 /
&HFBTHO_FEATURES
 collective_inertia = F, fission_fragments = F, pairing_regularization = F,
 localization_functions = F /
&HFBTHO_NECK
 set_neck_constrain = F, neck_value = 0.00 /
&HFBTHO_DFRAG
 set_dfrag_constrain = F, dfrag_value = 0.00 /
&HFBTHO_CSI
 set_csi_constrain = F, csi_value = 0.00 /
&HFBTHO_DEBUG
 number_Gauss = 80, number_Laguerre = 80, number_Legendre = 80,
 compatibility_HFODD = F, number_states = 2000, force_parity = F,
 print_time = 0 /
NAMELIST

    echo -n "  Running $LABEL (Z=$Z N=$N A=$A $FUNCTIONAL β₂=$BETA2)... " | tee -a $LOGFILE

    # Run and capture output
    OUTPUT=$($HFBTHO 2>&1)
    EXIT=$?

    # Parse results
    BE=$(echo "$OUTPUT" | grep "tEnergy: ehfb" | awk '{print $NF}')
    BETA2_OUT=$(echo "$OUTPUT" | grep "deformation beta2" | awk '{print $NF}')
    CONVERGED=$(echo "$OUTPUT" | grep -c "iteration converged")
    ITERS=$(echo "$OUTPUT" | grep "iteration converged" | grep -oP 'after\s+\K[0-9]+')
    PAIR_N=$(echo "$OUTPUT" | grep "pairing energy" | awk '{print $3}')
    PAIR_P=$(echo "$OUTPUT" | grep "pairing energy" | awk '{print $4}')
    RMS=$(echo "$OUTPUT" | grep "rms-radius" | awk '{print $NF}')
    CPU=$(echo "$OUTPUT" | grep "Total CPU time" | awk '{print $5}')

    if [ "$CONVERGED" -eq 1 ]; then
        STATUS="converged"
        echo "OK  BE=$BE MeV  β₂=$BETA2_OUT  iters=$ITERS  CPU=${CPU}min" | tee -a $LOGFILE
    else
        STATUS="failed"
        BE="null"
        echo "FAILED (exit=$EXIT)" | tee -a $LOGFILE
    fi

    # Append to JSON
    if [ $FIRST_ENTRY -eq 0 ]; then
        echo "," >> $JSONFILE
    fi
    FIRST_ENTRY=0

    cat >> $JSONFILE << JSON
  {
    "label": "$LABEL",
    "Z": $Z, "N": $N, "A": $A,
    "functional": "$FUNCTIONAL",
    "beta2_input": $BETA2,
    "beta2_output": "${BETA2_OUT:-null}",
    "BE_mev": ${BE:-null},
    "pairing_n": "${PAIR_N:-null}",
    "pairing_p": "${PAIR_P:-null}",
    "rms_fm": "${RMS:-null}",
    "status": "$STATUS",
    "iterations": "${ITERS:-null}",
    "cpu_minutes": "${CPU:-null}"
  }
JSON
}

# ============================================================
# FUNCTIONALS TO RUN
# ============================================================
FUNCTIONALS=("UNE1" "SLY4" "SKM*" "HFB9")

# ============================================================
# ISOTOPES: Mc (Z=115) parent nuclei A=285 to A=303
# Plus daughters: Nh (Z=113) A=281 to A=299
# Plus alpha: He-4 (already known = -28.296 MeV)
# ============================================================

echo "" | tee -a $LOGFILE
echo "PHASE 1: Mc isotopes (Z=115, A=285-303)" | tee -a $LOGFILE
echo "------------------------------------------------------------" | tee -a $LOGFILE

for FUNC in "${FUNCTIONALS[@]}"; do
    echo "" | tee -a $LOGFILE
    echo "  Functional: $FUNC" | tee -a $LOGFILE
    for A in 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303; do
        N=$((A - 115))
        run_hfbtho 115 $N "$FUNC" 0.0 "Mc-${A}"
    done
done

echo "" | tee -a $LOGFILE
echo "PHASE 2: Nh daughter nuclei (Z=113, A=281-299)" | tee -a $LOGFILE
echo "------------------------------------------------------------" | tee -a $LOGFILE

for FUNC in "${FUNCTIONALS[@]}"; do
    echo "" | tee -a $LOGFILE
    echo "  Functional: $FUNC" | tee -a $LOGFILE
    for A in 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299; do
        N=$((A - 113))
        run_hfbtho 113 $N "$FUNC" 0.0 "Nh-${A}"
    done
done

# Close JSON
echo "]" >> $JSONFILE

# ============================================================
# COMPUTE Q_ALPHA VALUES
# ============================================================
echo "" | tee -a $LOGFILE
echo "============================================================" | tee -a $LOGFILE
echo "COMPUTING Q_ALPHA VALUES" | tee -a $LOGFILE
echo "============================================================" | tee -a $LOGFILE

python3 << 'PYEOF'
import json, math

E_ALPHA = 28.2957  # MeV

with open("../nucleusgrid_results.json") as f:
    data = json.load(f)

# Build BE lookup: (Z, A, functional) -> BE
be = {}
for r in data:
    if r["status"] == "converged" and r["BE_mev"] is not None:
        key = (r["Z"], r["A"], r["functional"])
        be[key] = float(r["BE_mev"])

print(f"\n{'Isotope':<10} {'Functional':<10} {'BE_parent':>12} {'BE_daughter':>13} {'Q_alpha':>10} {'Est t½'}")
print("-"*72)

results = []
for func in ["UNE1","SLY4","SKM*","HFB9"]:
    for A_mc in range(285, 304):
        A_nh = A_mc - 4
        key_mc = (115, A_mc, func)
        key_nh = (113, A_nh, func)

        if key_mc in be and key_nh in be:
            BE_mc = be[key_mc]
            BE_nh = be[key_nh]
            Q = BE_mc - BE_nh - E_ALPHA

            # WKB half-life estimate (Royer formula for odd-Z SHN)
            # log10(t_half/s) = a + b*Z/sqrt(Q) where a,b from Royer 2000
            if Q > 0.5:
                # Simplified Viola-Seaborg for odd-Z
                log_t = (1.64 * 115 / math.sqrt(Q)) - 8.5
                t_s = 10**log_t
                if t_s < 0.001:      hl = f"{t_s*1000:.2f} ms"
                elif t_s < 1:        hl = f"{t_s*1000:.0f} ms"
                elif t_s < 60:       hl = f"{t_s:.1f} s"
                elif t_s < 3600:     hl = f"{t_s/60:.1f} min"
                elif t_s < 86400:    hl = f"{t_s/3600:.1f} hr"
                else:                hl = f"{t_s/86400:.1f} days"
            else:
                hl = "SF/β-"

            print(f"Mc-{A_mc:<6} {func:<10} {BE_mc:>12.3f} {BE_nh:>13.3f} {Q:>10.4f}   {hl}")
            results.append({
                "A": A_mc, "functional": func,
                "BE_mc": BE_mc, "BE_nh": BE_nh,
                "Q_alpha": round(Q,4),
                "halflife_est": hl
            })

# Save Q_alpha results
with open("../q_alpha_results.json","w") as f:
    json.dump(results, f, indent=2)

print(f"\nQ_alpha results saved to q_alpha_results.json")

# Summary: best candidates
print("\n" + "="*72)
print("TOP STABILITY CANDIDATES (highest Q gives longest t½ near island)")
print("="*72)
by_A = {}
for r in results:
    A = r["A"]
    if A not in by_A: by_A[A] = []
    by_A[A].append(r)

summary = []
for A, rs in by_A.items():
    qs = [r["Q_alpha"] for r in rs if r["Q_alpha"] < 15]
    if qs:
        summary.append((A, min(qs), len(qs)))

summary.sort(key=lambda x: x[1])
print(f"\n{'Isotope':<10} {'Min Q_α (MeV)':>14} {'Models':>8}  Note")
print("-"*50)
for A, minQ, n_models in summary[:8]:
    N = A - 115
    note = " ← N=184 magic" if N==184 else (" ← N=178" if N==178 else "")
    print(f"Mc-{A:<6}  {minQ:>12.4f}  {n_models:>6}/4{note}")

PYEOF

echo "" | tee -a $LOGFILE
echo "============================================================" | tee -a $LOGFILE
echo "ALL DONE. Finished: $(date)" | tee -a $LOGFILE
echo "Results: nucleusgrid_results.json" | tee -a $LOGFILE
echo "Q values: q_alpha_results.json" | tee -a $LOGFILE
echo "Log: nucleusgrid_results.log" | tee -a $LOGFILE
echo "============================================================" | tee -a $LOGFILE
echo ""
echo "Paste q_alpha_results.json to Claude for analysis."
