#!/bin/bash
# ============================================================
# NucleusGrid Phase 3 — Machine 2
# Part 1: alpha decay chains for Mc-291..299 (6 nuclei deep, Z 115->105)
# Part 2: extended Mc-300..311 + Nh daughters Nh-296..307
# Functionals: UNE1, SLY4, SKM*, HFB9
# ============================================================

PROJECT_ROOT="/home/jordandesilets/nucleusgrid"
HFBTHO_DIR="$PROJECT_ROOT/hfbtho"
HFBTHO_BIN="$HFBTHO_DIR/hfbtho_main"
RESULTS_DIR="$PROJECT_ROOT/results/phase3"
JSONFILE="$RESULTS_DIR/machine2_all.json"

mkdir -p "$RESULTS_DIR"
cd "$HFBTHO_DIR" || { echo "FATAL: cannot cd to $HFBTHO_DIR"; exit 1; }
[ -x "$HFBTHO_BIN" ] || { echo "FATAL: $HFBTHO_BIN is not executable"; exit 1; }

declare -A SYMBOL=( [115]=Mc [113]=Nh [111]=Rg [109]=Mt [107]=Bh [105]=Db )

STARTTIME=$(date)
echo "============================================================"
echo "NucleusGrid Phase 3 — Machine 2"
echo "Started: $STARTTIME"
echo "Binary:  $HFBTHO_BIN"
echo "Output:  $JSONFILE"
echo "============================================================"

echo "[" > "$JSONFILE"
FIRST_ENTRY=1
RUN_COUNT=0

run_hfbtho() {
    local Z=$1 N=$2 FUNC=$3 BETA2=$4 LABEL=$5 PHASE=$6 CHAIN_PARENT=$7
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
 functional = '$FUNC', add_initial_pairing = F, type_of_coulomb = 2 /
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

    RUN_COUNT=$((RUN_COUNT + 1))
    echo -n "[$RUN_COUNT] $LABEL (Z=$Z N=$N A=$A $FUNC β₂=$BETA2)... "

    OUTPUT=$("$HFBTHO_BIN" 2>&1)
    EXIT=$?

    BE=$(echo "$OUTPUT"        | grep "tEnergy: ehfb"        | awk '{print $NF}')
    BETA2_OUT=$(echo "$OUTPUT" | grep "deformation beta2"    | awk '{print $NF}')
    CONVERGED=$(echo "$OUTPUT" | grep -c "iteration converged")
    ITERS=$(echo "$OUTPUT"     | grep "iteration converged"  | grep -oP 'after\s+\K[0-9]+')
    PAIR_N=$(echo "$OUTPUT"    | grep "pairing energy"       | awk '{print $3}')
    PAIR_P=$(echo "$OUTPUT"    | grep "pairing energy"       | awk '{print $4}')
    RMS=$(echo "$OUTPUT"       | grep "rms-radius"           | awk '{print $NF}')
    CPU=$(echo "$OUTPUT"       | grep "Total CPU time"       | awk '{print $5}')

    if [ "$CONVERGED" -eq 1 ]; then
        STATUS="converged"
        echo "OK BE=$BE β₂=$BETA2_OUT iters=$ITERS CPU=${CPU}min"
    else
        STATUS="failed"
        BE="null"
        echo "FAILED (exit=$EXIT)"
    fi

    if [ "$FIRST_ENTRY" -eq 0 ]; then
        echo "," >> "$JSONFILE"
    fi
    FIRST_ENTRY=0

    cat >> "$JSONFILE" << JSON
  {
    "label": "$LABEL",
    "Z": $Z, "N": $N, "A": $A,
    "functional": "$FUNC",
    "beta2_input": $BETA2,
    "beta2_output": "${BETA2_OUT:-null}",
    "BE_mev": ${BE:-null},
    "pairing_n": "${PAIR_N:-null}",
    "pairing_p": "${PAIR_P:-null}",
    "rms_fm": "${RMS:-null}",
    "status": "$STATUS",
    "iterations": "${ITERS:-null}",
    "cpu_minutes": "${CPU:-null}",
    "phase": "$PHASE",
    "chain_parent": "$CHAIN_PARENT"
  }
JSON
}

FUNCTIONALS=("UNE1" "SLY4" "SKM*" "HFB9")

# ============================================================
# PART 1: alpha decay chains Mc-291..299, 6 nuclei deep
# Z steps: 115 -> 113 -> 111 -> 109 -> 107 -> 105
# A steps: A0 -> A0-4 -> A0-8 -> A0-12 -> A0-16 -> A0-20
# ============================================================
echo ""
echo "------------------------------------------------------------"
echo "PART 1: alpha decay chains Mc-291..299 (6 nuclei deep)"
echo "------------------------------------------------------------"

for FUNC in "${FUNCTIONALS[@]}"; do
    echo ""
    echo "  Functional: $FUNC"
    for A0 in 291 292 293 294 295 296 297 298 299; do
        CHAIN="Mc-${A0}"
        for STEP in 0 1 2 3 4 5; do
            Z=$((115 - 2*STEP))
            A=$((A0 - 4*STEP))
            N=$((A - Z))
            SYM="${SYMBOL[$Z]}"
            run_hfbtho "$Z" "$N" "$FUNC" 0.0 "${SYM}-${A}" "part1_chain" "$CHAIN"
        done
    done
done

# ============================================================
# PART 2: extended Mc-300..311 + Nh daughters Nh-296..307
# ============================================================
echo ""
echo "------------------------------------------------------------"
echo "PART 2: extended Mc-300..311 + Nh-296..307"
echo "------------------------------------------------------------"

for FUNC in "${FUNCTIONALS[@]}"; do
    echo ""
    echo "  Functional: $FUNC (Mc-300..311)"
    for A in 300 301 302 303 304 305 306 307 308 309 310 311; do
        N=$((A - 115))
        run_hfbtho 115 "$N" "$FUNC" 0.0 "Mc-${A}" "part2_extended" ""
    done
    echo ""
    echo "  Functional: $FUNC (Nh-296..307)"
    for A in 296 297 298 299 300 301 302 303 304 305 306 307; do
        N=$((A - 113))
        run_hfbtho 113 "$N" "$FUNC" 0.0 "Nh-${A}" "part2_extended" ""
    done
done

echo "]" >> "$JSONFILE"

echo ""
echo "============================================================"
echo "Machine 2 finished: $(date)"
echo "Total runs: $RUN_COUNT"
echo "Results:    $JSONFILE"
echo "============================================================"
