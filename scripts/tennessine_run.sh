#!/bin/bash
# ============================================================
# NucleusGrid Phase 3 — Tennessine (Z=117)
# Part 1: alpha decay chains for Ts-291..299 (6 nuclei deep, Z 117 -> 107)
# Part 2: extended Ts-300..311 + Mc daughters Mc-296..307
# Functionals: UNE1, SLY4, SKM*, HFB9 ; β₂ = 0.0 throughout
#
# RESUME: previously-converged runs are detected from tennessine_all.jsonl
# and skipped. Final .json is rebuilt from .jsonl after the run.
# ============================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( dirname "$SCRIPT_DIR" )"
HFBTHO_DIR="$PROJECT_ROOT/hfbtho"
HFBTHO_BIN="$HFBTHO_DIR/hfbtho_main"
RESULTS_DIR="$PROJECT_ROOT/results/phase3"
JSONFILE="$RESULTS_DIR/tennessine_all.json"
JSONLFILE="$RESULTS_DIR/tennessine_all.jsonl"

mkdir -p "$RESULTS_DIR"
cd "$HFBTHO_DIR" || { echo "FATAL: cannot cd to $HFBTHO_DIR"; exit 1; }
[ -x "$HFBTHO_BIN" ] || { echo "FATAL: $HFBTHO_BIN is not executable"; exit 1; }

declare -A SYMBOL=( [117]=Ts [115]=Mc [113]=Nh [111]=Rg [109]=Mt [107]=Bh )

# 9 chains × 6 deep = 54  +  Part 2: 12 Ts + 12 Mc = 24  →  78 unique nuclei × 4 funcs
TOTAL=$(( (9*6 + 12 + 12) * 4 ))

# Resume: build DONE set from existing .jsonl
declare -A DONE
SKIPPED=0
if [ -f "$JSONLFILE" ] && [ -s "$JSONLFILE" ]; then
    DONE_TMP=$(mktemp)
    python3 - "$JSONLFILE" > "$DONE_TMP" <<'PYEOF'
import json, sys
with open(sys.argv[1]) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        if r.get("status") == "converged":
            print(f'{r["Z"]}|{r["A"]}|{r["functional"]}|{r["beta2_input"]}')
PYEOF
    while IFS= read -r key; do
        [ -z "$key" ] && continue
        DONE["$key"]=1
        SKIPPED=$((SKIPPED + 1))
    done < "$DONE_TMP"
    rm -f "$DONE_TMP"
fi

STARTTIME=$(date)
echo "============================================================"
echo "NucleusGrid Phase 3 — Tennessine (Z=117)"
echo "Started:        $STARTTIME"
echo "Binary:         $HFBTHO_BIN"
echo "Output (JSONL): $JSONLFILE"
echo "Output (JSON):  $JSONFILE  (rebuilt from JSONL at end)"
echo "Total planned:  $TOTAL runs"
echo "Resumed from:   $SKIPPED previously-converged runs"
echo "Remaining:      $((TOTAL - SKIPPED)) runs to execute"
echo "============================================================"

RUN_COUNT=$SKIPPED
NEW_THIS_SESSION=0

run_hfbtho() {
    local Z=$1 N=$2 FUNC=$3 BETA2=$4 LABEL=$5 PHASE=$6 CHAIN_PARENT=$7
    local A=$((Z + N))
    local KEY="$Z|$A|$FUNC|$BETA2"

    if [ -n "${DONE[$KEY]+x}" ]; then
        return 0
    fi

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
    NEW_THIS_SESSION=$((NEW_THIS_SESSION + 1))
    echo -n "[$RUN_COUNT/$TOTAL $((100*RUN_COUNT/TOTAL))%] $LABEL (Z=$Z N=$N A=$A $FUNC β₂=$BETA2)... "

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
        echo "OK BE=$BE β₂=$BETA2_OUT iters=$ITERS"
    else
        STATUS="failed"
        BE="null"
        echo "FAILED (exit=$EXIT)"
    fi

    {
        cat <<JSONL | tr -d '\n'
{"label":"$LABEL","Z":$Z,"N":$N,"A":$A,"functional":"$FUNC","beta2_input":$BETA2,"beta2_output":"${BETA2_OUT:-null}","BE_mev":${BE:-null},"pairing_n":"${PAIR_N:-null}","pairing_p":"${PAIR_P:-null}","rms_fm":"${RMS:-null}","status":"$STATUS","iterations":"${ITERS:-null}","cpu_minutes":"${CPU:-null}","phase":"$PHASE","chain_parent":"$CHAIN_PARENT"}
JSONL
        echo
    } >> "$JSONLFILE"

    if [ "$STATUS" = "converged" ]; then
        DONE["$KEY"]=1
    fi
}

FUNCTIONALS=("UNE1" "SLY4" "SKM*" "HFB9")

# ============================================================
# PART 1: alpha decay chains Ts-291..299, 6 nuclei deep
# Z steps: 117 -> 115 -> 113 -> 111 -> 109 -> 107
# ============================================================
echo ""
echo "------------------------------------------------------------"
echo "PART 1: alpha decay chains Ts-291..299 (6 nuclei deep)"
echo "------------------------------------------------------------"

for FUNC in "${FUNCTIONALS[@]}"; do
    echo ""
    echo "  Functional: $FUNC"
    for A0 in 291 292 293 294 295 296 297 298 299; do
        CHAIN="Ts-${A0}"
        for STEP in 0 1 2 3 4 5; do
            Z=$((117 - 2*STEP))
            A=$((A0 - 4*STEP))
            N=$((A - Z))
            SYM="${SYMBOL[$Z]}"
            run_hfbtho "$Z" "$N" "$FUNC" 0.0 "${SYM}-${A}" "part1_chain" "$CHAIN"
        done
    done
done

# ============================================================
# PART 2: extended Ts-300..311 + Mc daughters Mc-296..307
# ============================================================
echo ""
echo "------------------------------------------------------------"
echo "PART 2: extended Ts-300..311 + Mc-296..307"
echo "------------------------------------------------------------"

for FUNC in "${FUNCTIONALS[@]}"; do
    echo ""
    echo "  Functional: $FUNC (Ts-300..311)"
    for A in 300 301 302 303 304 305 306 307 308 309 310 311; do
        N=$((A - 117))
        run_hfbtho 117 "$N" "$FUNC" 0.0 "Ts-${A}" "part2_extended" ""
    done
    echo ""
    echo "  Functional: $FUNC (Mc-296..307)"
    for A in 296 297 298 299 300 301 302 303 304 305 306 307; do
        N=$((A - 115))
        run_hfbtho 115 "$N" "$FUNC" 0.0 "Mc-${A}" "part2_extended" ""
    done
done

# ============================================================
# Rebuild canonical .json from .jsonl
# ============================================================
python3 - "$JSONLFILE" "$JSONFILE" <<'PYEOF'
import json, sys
src, dst = sys.argv[1], sys.argv[2]
data = []
with open(src) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        try:
            data.append(json.loads(line))
        except json.JSONDecodeError:
            pass
with open(dst, "w") as f:
    json.dump(data, f, indent=2)
print(f"  Rebuilt {dst} ({len(data)} entries) from {src}")
PYEOF

echo ""
echo "============================================================"
echo "Tennessine run finished: $(date)"
echo "Total runs:        $RUN_COUNT/$TOTAL"
echo "  resumed:         $SKIPPED"
echo "  this session:    $NEW_THIS_SESSION"
echo "Streaming:         $JSONLFILE"
echo "Final results:     $JSONFILE"
echo "============================================================"
