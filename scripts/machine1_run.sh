#!/bin/bash
# ============================================================
# NucleusGrid Phase 3 — Machine 1
# Stability island core: PES deformation sweep
#   Mc-285..299 (Z=115) and Nh daughters Nh-281..295 (Z=113)
# beta2 sweep: -0.2, -0.1, 0.0, 0.1, 0.15, 0.2, 0.25, 0.3
# Functionals: UNE1, SLY4, SKM*, HFB9
#
# beta2=0.0 is SKIPPED for Mc-291..299 and Nh-287..295 because
# Machine 2's Part 1 alpha-chain run already covers those nuclei
# at beta2=0.0 (chain step 0 = Mc-291..299, step 1 = Nh-287..295).
# Mc-285..290 and Nh-281..286 are unique to Machine 1 and run
# the full 8-point sweep including beta2=0.0.
# ============================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( dirname "$SCRIPT_DIR" )"
HFBTHO_DIR="$PROJECT_ROOT/hfbtho"
HFBTHO_BIN="$HFBTHO_DIR/hfbtho_main"
RESULTS_DIR="$PROJECT_ROOT/results/phase3"
JSONFILE="$RESULTS_DIR/machine1_all.json"
JSONLFILE="$RESULTS_DIR/machine1_all.jsonl"
OUTPUTFILE="$RESULTS_DIR/machine1_output.txt"

# 4 functionals × (6 nuclei × 8 β₂ points + 9 nuclei × 7 β₂ points) × 2 sides (Mc, Nh)
TOTAL=$((  4 * (6*8 + 9*7) * 2 ))

# Push progress to GitHub every PUSH_INTERVAL completed runs.
# Set PUSH_INTERVAL=0 to disable progress pushes.
PUSH_INTERVAL=${PUSH_INTERVAL:-50}

mkdir -p "$RESULTS_DIR"
cd "$HFBTHO_DIR" || { echo "FATAL: cannot cd to $HFBTHO_DIR"; exit 1; }
[ -x "$HFBTHO_BIN" ] || { echo "FATAL: $HFBTHO_BIN is not executable"; exit 1; }

declare -A SYMBOL=( [115]=Mc [113]=Nh )

STARTTIME=$(date)
echo "============================================================"
echo "NucleusGrid Phase 3 — Machine 1"
echo "Started:        $STARTTIME"
echo "Binary:         $HFBTHO_BIN"
echo "Output (JSON):  $JSONFILE  (finalised at end of run)"
echo "Output (JSONL): $JSONLFILE (always-valid streaming)"
echo "Total planned:  $TOTAL runs"
echo "Progress push:  every $PUSH_INTERVAL runs (set PUSH_INTERVAL=0 to disable)"
echo "============================================================"

# ============================================================
# Resume support: if $JSONFILE already exists with converged
# entries, keep them, drop failed entries, and skip those keys
# during this run. Otherwise start a fresh JSON array. The
# streaming JSONL is truncated on fresh start and preserved on
# resume so it stays aligned with whatever JSONFILE entries we
# keep.
# ============================================================
SKIP_KEYS_FILE="$(mktemp)"
trap 'rm -f "$SKIP_KEYS_FILE"' EXIT

FIRST_ENTRY=1
RUN_COUNT=0
SKIP_COUNT=0

if [ -s "$JSONFILE" ]; then
    echo ""
    echo "Found existing $JSONFILE — scanning for converged entries..."
    if python3 - "$JSONFILE" "$SKIP_KEYS_FILE" "$JSONFILE.resume.tmp" << 'PYEOF'
import json, sys
from pathlib import Path

src = Path(sys.argv[1])
keys_out = Path(sys.argv[2])
json_out = Path(sys.argv[3])

raw = src.read_text().rstrip()
# The script writes entries as "  {...}\n,\n  {...}\n" with no trailing ].
# Recover by trimming any trailing comma/whitespace and appending ].
recovered = raw.rstrip(',').rstrip()
if not recovered.endswith(']'):
    recovered += '\n]'

try:
    entries = json.loads(recovered)
except json.JSONDecodeError:
    # Last entry may be partially written; truncate at last complete '}'
    last = recovered.rfind('}')
    if last < 0:
        sys.exit(2)
    truncated = recovered[:last+1].rstrip().rstrip(',') + '\n]'
    entries = json.loads(truncated)

converged = [e for e in entries if e.get('status') == 'converged']
dropped   = len(entries) - len(converged)

with keys_out.open('w') as f:
    for e in converged:
        b2 = e['beta2_input']
        # Match the bash literal: "0.0", "-0.2", "0.15", etc.
        b2_str = repr(b2) if isinstance(b2, float) else str(b2)
        f.write(f"{e['Z']}|{e['N']}|{e['functional']}|{b2_str}\n")

# Re-emit JSON with only converged entries, no trailing ] (so bash can append).
with json_out.open('w') as f:
    f.write('[')
    for i, e in enumerate(converged):
        if i > 0:
            f.write('\n,')
        be = e["BE_mev"]
        be_str = "null" if be is None else repr(be)
        f.write(f'''
  {{
    "label": "{e['label']}",
    "Z": {e['Z']}, "N": {e['N']}, "A": {e['A']},
    "functional": "{e['functional']}",
    "beta2_input": {repr(e['beta2_input'])},
    "beta2_output": "{e['beta2_output']}",
    "BE_mev": {be_str},
    "pairing_n": "{e['pairing_n']}",
    "pairing_p": "{e['pairing_p']}",
    "rms_fm": "{e['rms_fm']}",
    "status": "{e['status']}",
    "iterations": "{e['iterations']}",
    "cpu_minutes": "{e['cpu_minutes']}",
    "phase": "{e['phase']}",
    "chain_parent": "{e['chain_parent']}"
  }}''')

print(f"  resume scan: {len(converged)} converged kept, {dropped} failed dropped (will retry)")
PYEOF
    then
        mv "$JSONFILE.resume.tmp" "$JSONFILE"
        SKIP_COUNT=$(wc -l < "$SKIP_KEYS_FILE")
        if [ "$SKIP_COUNT" -gt 0 ]; then
            FIRST_ENTRY=0
            echo "  will skip $SKIP_COUNT already-converged calculations on this run"
            # JSONL preserved as-is (streaming continues from prior session)
        else
            FIRST_ENTRY=1
            echo "[" > "$JSONFILE"
            > "$JSONLFILE"
            echo "  no converged entries found — starting fresh"
        fi
    else
        echo "  WARNING: resume scan failed — starting fresh"
        echo "[" > "$JSONFILE"
        > "$JSONLFILE"
        FIRST_ENTRY=1
    fi
else
    echo "[" > "$JSONFILE"
    > "$JSONLFILE"
fi

push_progress() {
    [ "$PUSH_INTERVAL" -eq 0 ] && return 0
    local PR="$PROJECT_ROOT" msg
    git -C "$PR" add results/phase3/machine1_all.jsonl 2>/dev/null || true
    [ -f "$OUTPUTFILE" ] && git -C "$PR" add results/phase3/machine1_output.txt 2>/dev/null || true
    if git -C "$PR" diff --cached --quiet 2>/dev/null; then
        return 0
    fi
    msg="Machine 1 progress: $RUN_COUNT/$TOTAL ($((100*RUN_COUNT/TOTAL))%)"
    git -C "$PR" commit -m "$msg" --quiet 2>/dev/null || { echo "  [progress: commit failed]"; return 1; }
    if ! git -C "$PR" pull --rebase --autostash origin main >/dev/null 2>&1; then
        git -C "$PR" rebase --abort >/dev/null 2>&1 || true
        echo "  [progress: rebase failed, skipping push this cycle]"
        return 1
    fi
    if ! git -C "$PR" push origin main >/dev/null 2>&1; then
        echo "  [progress: push failed (auth?), will retry next cycle]"
        return 1
    fi
    echo "  [progress: pushed $RUN_COUNT/$TOTAL to GitHub]"
}

run_hfbtho() {
    local Z=$1 N=$2 FUNC=$3 BETA2=$4 LABEL=$5 PHASE=$6 CHAIN_PARENT=$7
    local A=$((Z + N))

    # Resume: skip if this exact (Z, N, functional, beta2) already converged.
    local KEY="$Z|$N|$FUNC|$BETA2"
    if [ "$SKIP_COUNT" -gt 0 ] && grep -qxF "$KEY" "$SKIP_KEYS_FILE"; then
        RUN_COUNT=$((RUN_COUNT + 1))
        echo "[$RUN_COUNT] $LABEL (Z=$Z N=$N A=$A $FUNC β₂=$BETA2)... SKIP (already converged)"
        return 0
    fi

    cat > hfbtho_NAMELIST.dat << NAMELIST
&HFBTHO_GENERAL
 number_of_shells = 16, oscillator_length = -1.0, basis_deformation = 0.0,
 proton_number = $Z, neutron_number = $N, type_of_calculation = 1 /
&HFBTHO_INITIAL
 beta2_deformation = $BETA2, beta4_deformation = 0.0 /
&HFBTHO_ITERATIONS
 number_iterations = 600, accuracy = 1.E-4, restart_file = 1 /
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

    # Always-valid streaming companion: one compact JSON object per line.
    {
        cat << JSONL | tr -d '\n'
{"label":"$LABEL","Z":$Z,"N":$N,"A":$A,"functional":"$FUNC","beta2_input":$BETA2,"beta2_output":"${BETA2_OUT:-null}","BE_mev":${BE:-null},"pairing_n":"${PAIR_N:-null}","pairing_p":"${PAIR_P:-null}","rms_fm":"${RMS:-null}","status":"$STATUS","iterations":"${ITERS:-null}","cpu_minutes":"${CPU:-null}","phase":"$PHASE","chain_parent":"$CHAIN_PARENT"}
JSONL
        echo
    } >> "$JSONLFILE"

    if [ "$PUSH_INTERVAL" -gt 0 ] && [ $((RUN_COUNT % PUSH_INTERVAL)) -eq 0 ]; then
        push_progress
    fi
}

FUNCTIONALS=("UNE1" "SLY4" "SKM*" "HFB9")
BETA2_FULL=(-0.2 -0.1 0.0 0.1 0.15 0.2 0.25 0.3)
BETA2_NOZERO=(-0.2 -0.1 0.1 0.15 0.2 0.25 0.3)

# ============================================================
# Mc isotopes (Z=115): PES deformation sweep
# ============================================================
echo ""
echo "------------------------------------------------------------"
echo "Mc PES sweep: Mc-285..290 (full 8-pt) + Mc-291..299 (7-pt, no β₂=0)"
echo "------------------------------------------------------------"

for FUNC in "${FUNCTIONALS[@]}"; do
    echo ""
    echo "  Functional: $FUNC (Mc-285..290 full sweep)"
    for A in 285 286 287 288 289 290; do
        N=$((A - 115))
        for B2 in "${BETA2_FULL[@]}"; do
            run_hfbtho 115 "$N" "$FUNC" "$B2" "Mc-${A}" "pes_sweep" ""
        done
    done

    echo ""
    echo "  Functional: $FUNC (Mc-291..299 sweep, β₂=0 skipped — covered by Machine 2)"
    for A in 291 292 293 294 295 296 297 298 299; do
        N=$((A - 115))
        for B2 in "${BETA2_NOZERO[@]}"; do
            run_hfbtho 115 "$N" "$FUNC" "$B2" "Mc-${A}" "pes_sweep" ""
        done
    done
done

# ============================================================
# Nh isotopes (Z=113): PES deformation sweep
# ============================================================
echo ""
echo "------------------------------------------------------------"
echo "Nh PES sweep: Nh-281..286 (full 8-pt) + Nh-287..295 (7-pt, no β₂=0)"
echo "------------------------------------------------------------"

for FUNC in "${FUNCTIONALS[@]}"; do
    echo ""
    echo "  Functional: $FUNC (Nh-281..286 full sweep)"
    for A in 281 282 283 284 285 286; do
        N=$((A - 113))
        for B2 in "${BETA2_FULL[@]}"; do
            run_hfbtho 113 "$N" "$FUNC" "$B2" "Nh-${A}" "pes_sweep" ""
        done
    done

    echo ""
    echo "  Functional: $FUNC (Nh-287..295 sweep, β₂=0 skipped — covered by Machine 2)"
    for A in 287 288 289 290 291 292 293 294 295; do
        N=$((A - 113))
        for B2 in "${BETA2_NOZERO[@]}"; do
            run_hfbtho 113 "$N" "$FUNC" "$B2" "Nh-${A}" "pes_sweep" ""
        done
    done
done

echo "]" >> "$JSONFILE"

echo ""
echo "============================================================"
echo "Machine 1 finished: $(date)"
echo "Total runs: $RUN_COUNT/$TOTAL"
echo "Results:    $JSONFILE  (final)"
echo "Streaming:  $JSONLFILE"
echo "============================================================"

# Final push including the finalised .json array
if [ "$PUSH_INTERVAL" -gt 0 ]; then
    git -C "$PROJECT_ROOT" add results/phase3/machine1_all.json results/phase3/machine1_all.jsonl 2>/dev/null || true
    [ -f "$OUTPUTFILE" ] && git -C "$PROJECT_ROOT" add results/phase3/machine1_output.txt 2>/dev/null || true
    if ! git -C "$PROJECT_ROOT" diff --cached --quiet 2>/dev/null; then
        git -C "$PROJECT_ROOT" commit -m "Machine 1 finished: $RUN_COUNT/$TOTAL runs complete" --quiet || true
        git -C "$PROJECT_ROOT" pull --rebase --autostash origin main >/dev/null 2>&1 \
            || git -C "$PROJECT_ROOT" rebase --abort >/dev/null 2>&1 || true
        git -C "$PROJECT_ROOT" push origin main >/dev/null 2>&1 \
            && echo "  [final push: results delivered to GitHub]" \
            || echo "  [final push: failed — push manually]"
    fi
fi
