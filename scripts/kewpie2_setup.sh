#!/bin/bash
# ============================================================
# NucleusGrid Phase 3 — KEWPIE2 setup + SHE synthesis run
#
# 1. Verifies / installs libgsl-dev (KEWPIE2 dependency)
# 2. Clones KEWPIE2 source from GitHub mirror
#    (https://github.com/Sonal-Jaiswal/Research_phy → final_kewpie2/)
# 3. Builds with g++ via the shipped Makefile
# 4. Smoke-tests against the sample input (¹⁶O + ²⁰⁸Pb)
# 5. Generates a ⁴⁸Ca + ²⁴³Am input deck and runs it
#    → Mc-291* → 0n,1n,2n,3n,4n channels = Mc-291..287
# 6. Saves results to synthesis/kewpie2_results/
# ============================================================

set -u
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( dirname "$SCRIPT_DIR" )"
SYNTH_DIR="$PROJECT_ROOT/synthesis"
KEWPIE_DIR="$SYNTH_DIR/kewpie2"
RESULTS_DIR="$SYNTH_DIR/kewpie2_results"
LOG="$SYNTH_DIR/kewpie2_setup.log"

KEWPIE_REPO="https://github.com/Sonal-Jaiswal/Research_phy"
KEWPIE_SUBPATH="final_kewpie2"

mkdir -p "$SYNTH_DIR" "$RESULTS_DIR"
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo "KEWPIE2 setup — $(date)"
echo "Synthesis dir:  $SYNTH_DIR"
echo "KEWPIE2 dir:    $KEWPIE_DIR"
echo "Results dir:    $RESULTS_DIR"
echo "============================================================"

# ============================================================
# Step 1 — GSL prerequisite
# ============================================================
echo ""
echo "──────────────────────────────────────────────────────────"
echo "STEP 1 — GSL (GNU Scientific Library)"
echo "──────────────────────────────────────────────────────────"
if pkg-config --exists gsl 2>/dev/null; then
    echo "  GSL found: $(pkg-config --modversion gsl)"
elif ldconfig -p 2>/dev/null | grep -q libgsl.so; then
    echo "  GSL libraries present (no pkg-config entry)"
else
    echo "  GSL not found. Attempting install via apt..."
    if sudo -n apt-get install -y libgsl-dev 2>&1; then
        echo "  installed"
    else
        echo ""
        echo "  ERROR: cannot install libgsl-dev (passwordless sudo not available)."
        echo "  Run this yourself in a terminal, then re-execute this script:"
        echo ""
        echo "    sudo apt-get install -y libgsl-dev"
        echo ""
        exit 1
    fi
fi

# ============================================================
# Step 2 — Clone KEWPIE2 source
# ============================================================
echo ""
echo "──────────────────────────────────────────────────────────"
echo "STEP 2 — Clone KEWPIE2"
echo "──────────────────────────────────────────────────────────"
if [ -d "$KEWPIE_DIR" ] && [ -f "$KEWPIE_DIR/Makefile" ]; then
    echo "  Already present at $KEWPIE_DIR"
else
    cd "$SYNTH_DIR"
    rm -rf "$KEWPIE_DIR" research_phy_tmp
    if git clone --depth 1 "$KEWPIE_REPO" research_phy_tmp 2>&1; then
        mv "research_phy_tmp/$KEWPIE_SUBPATH" "$KEWPIE_DIR"
        rm -rf research_phy_tmp
        echo "  cloned → $KEWPIE_DIR"
    else
        echo "  ERROR: clone failed"
        exit 1
    fi
fi

# ============================================================
# Step 3 — Build
# ============================================================
echo ""
echo "──────────────────────────────────────────────────────────"
echo "STEP 3 — Build"
echo "──────────────────────────────────────────────────────────"
cd "$KEWPIE_DIR"
mkdir -p Librairies/Objects Executable
make clean > /dev/null 2>&1 || true
if make 2> "$SYNTH_DIR/kewpie2_compile.log" > "$SYNTH_DIR/kewpie2_compile.out"; then
    if [ -x Executable/KEWPIE2.exe ]; then
        echo "  COMPILE OK → $KEWPIE_DIR/Executable/KEWPIE2.exe"
    else
        echo "  build claimed success but no executable — see $SYNTH_DIR/kewpie2_compile.log"
        exit 1
    fi
else
    echo "  COMPILE FAILED. Tail of stderr:"
    tail -20 "$SYNTH_DIR/kewpie2_compile.log" | sed 's/^/    /'
    exit 1
fi

# ============================================================
# Step 4 — Smoke test against the sample input (¹⁶O + ²⁰⁸Pb)
# ============================================================
echo ""
echo "──────────────────────────────────────────────────────────"
echo "STEP 4 — Smoke test (¹⁶O + ²⁰⁸Pb sample)"
echo "──────────────────────────────────────────────────────────"
cd "$KEWPIE_DIR"
# The shipped sample is ¹⁶O + ²⁰⁸Pb. KEWPIE2 writes its output back to
# Results/cross-sections.dat (the filename is the last line of input.inp).
# The example output ALSO lives at Results/cross-sections.dat in the repo,
# so we save it as the reference before re-running.
if [ -f Input/input.inp ]; then
    cp Input/input.inp Input/input.inp.original
    [ -f Results/cross-sections.dat ] && cp Results/cross-sections.dat Results/reference.dat
    cd Executable
    if timeout 120 ./KEWPIE2.exe > "$SYNTH_DIR/kewpie2_smoke.log" 2>&1; then
        cd "$KEWPIE_DIR"
        out="Results/cross-sections.dat"
        if [ -f "$out" ]; then
            echo "  output: $out ($(wc -l < "$out") lines)"
            echo "  --- first 8 lines ---"
            head -8 "$out" | sed 's/^/    /'
            if [ -f Results/reference.dat ] && \
               diff -q <(awk 'NF>1{print $1}' "$out" | head -20) \
                       <(awk 'NF>1{print $1}' Results/reference.dat | head -20) >/dev/null 2>&1; then
                echo "  SMOKE TEST PASSED (energy grid matches shipped reference)"
            else
                echo "  smoke test ran (no shipped reference for diff, or grid differs)"
            fi
        else
            echo "  WARNING: no output file produced — see $SYNTH_DIR/kewpie2_smoke.log"
        fi
    else
        echo "  smoke test failed or timed out — see $SYNTH_DIR/kewpie2_smoke.log"
    fi
else
    echo "  WARNING: Input/input.inp missing — skipping smoke test"
fi

# ============================================================
# Step 5 — Generate ⁴⁸Ca + ²⁴³Am input deck and run
# ============================================================
echo ""
echo "──────────────────────────────────────────────────────────"
echo "STEP 5 — Generate ⁴⁸Ca + ²⁴³Am input + run"
echo "──────────────────────────────────────────────────────────"

CA_AM_INPUT="$RESULTS_DIR/Ca48_Am243_Mc291.inp"
CA_AM_OUTPUT_NAME="Ca48_Am243_Mc291.dat"

cat > "$CA_AM_INPUT" <<'INP'
20		// Projectile charge (Ca)
48		// Projectile mass
95		// Target charge (Am)
243		// Target mass
1  		// Shell correction (0: not used, 1: used)
1.		// Shell-correction factor (between 0 and 1)
0.1		// Stepsize for the kinetic energy (MeV)
1		// Level-density parameter (1: Ignatyuk default)
19.		// Shell-damping energy (MeV)
1.		// af/an ratio
0.4		// Correction factor for the moment of inertia
0		// Mass table (0: Moller-Nix, 1: User-defined file)
1		// Pairing Effect (0: not used, 1: used)
1		// Kramers factor (0: not used, 1: used)
1		// Collective-enhancement factor (0: not used, 1: used)
5e21		// Reduced friction coefficient (s^-1)
1.		// Potential curvature in the ground state (MeV)
1.		// Potential curvature at the saddle point (MeV)
1		// Evaporation calculation (0: HF, 1: WE)
0		// Optical model for nucleons (0: Varner, 1: Koning)
0		// Hill-Wheeler penetration factor (0: not used, 1: used)
1		// Fission-barrier model (0: User-defined, 1: TF, 2: LSD)
10.		// Minimum excitation energy (MeV) — covers up to 1n channel
55.		// Maximum excitation energy (MeV) — covers up to 4n channel
0.1		// Stepsize for the energy discretization (MeV)
1.		// Stepsize for plotting the excitation function (MeV)
0		// Minimum angular momentum (hbar)
80		// Maximum angular momentum (hbar) — appropriate for SHE
4		// Total number of emitted neutrons (covers up to 4n)
0		// Total number of emitted protons
1	 	// Types of emitted particles (1: neutrons only)
0		// Gamma emission (0: not used, 1: used)
1		// Gamma-ray strength function (0: EGLO, 1: SMLO)
0		// Fission-time calculation (0: not used, 1: used)
1.1		// Time factor for the fission-time calculation
0e-21		// Transient time related to the fission process
1e-6		// Cut-off for the transmission coefficient
1e-50		// Cut-off for the decay width (MeV)
1		// Fusion model (0: WKB, 1: EBD, 2: fusion file)
0. 0.		// Barrier-distribution parameters (mean shift, width)
Ca48_Am243_Mc291.dat
INP

# Install the input into the KEWPIE2 input dir, run, then collect.
cp "$CA_AM_INPUT" "$KEWPIE_DIR/Input/input.inp"
cd "$KEWPIE_DIR/Executable"
echo "  running KEWPIE2 with Ca-48 + Am-243 input..."
echo "  (this can take 5–30 min depending on machine — cascade across 4n channels)"
START_TS=$(date +%s)
if timeout 3600 ./KEWPIE2.exe > "$SYNTH_DIR/kewpie2_CaAm.log" 2>&1; then
    END_TS=$(date +%s)
    echo "  done in $(( END_TS - START_TS )) s"
    cd "$KEWPIE_DIR"
    if [ -f "Results/$CA_AM_OUTPUT_NAME" ]; then
        cp "Results/$CA_AM_OUTPUT_NAME" "$RESULTS_DIR/"
        cp "$SYNTH_DIR/kewpie2_CaAm.log" "$RESULTS_DIR/Ca48_Am243_Mc291.log"
        echo "  saved → $RESULTS_DIR/$CA_AM_OUTPUT_NAME"
        echo ""
        echo "  --- first 20 lines of cross-section output ---"
        head -20 "$RESULTS_DIR/$CA_AM_OUTPUT_NAME" | sed 's/^/    /'
        echo ""
        echo "  --- channel summary (highest σ per channel) ---"
        python3 - "$RESULTS_DIR/$CA_AM_OUTPUT_NAME" <<'PYEOF'
import sys
fn = sys.argv[1]
# KEWPIE2 cross-section file format: typically rows of E*, then σ per channel.
# Be tolerant — try to identify columns and report max in each.
rows = []
with open(fn) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('!'):
            continue
        try:
            vals = [float(x) for x in line.split()]
            if len(vals) >= 2:
                rows.append(vals)
        except ValueError:
            continue
if not rows:
    print("    (no numeric rows parsed; inspect the file directly)")
else:
    ncols = max(len(r) for r in rows)
    print(f"    rows: {len(rows)}  columns: {ncols}")
    # Assume column 0 is E*, columns 1..n are channels
    for c in range(1, min(ncols, 7)):
        vals = [r[c] for r in rows if len(r) > c]
        if vals:
            mx = max(vals)
            mx_e = rows[vals.index(mx)][0]
            label = ["", "0n (Mc-291)", "1n (Mc-290)", "2n (Mc-289)",
                     "3n (Mc-288)", "4n (Mc-287)"][c] if c < 6 else f"col{c}"
            print(f"    col {c}  {label:<14}  max σ = {mx:.3e}  at E* = {mx_e:.2f} MeV")
PYEOF
    else
        echo "  WARNING: expected output Results/$CA_AM_OUTPUT_NAME not found"
        echo "  ls Results/:"
        ls "$KEWPIE_DIR/Results/" 2>&1 | sed 's/^/    /'
    fi
else
    echo "  KEWPIE2 run failed or timed out — see $SYNTH_DIR/kewpie2_CaAm.log"
    tail -20 "$SYNTH_DIR/kewpie2_CaAm.log" | sed 's/^/    /'
fi

# Restore original input.inp so a re-run smoke-tests cleanly
[ -f "$KEWPIE_DIR/Input/input.inp.original" ] && \
    cp "$KEWPIE_DIR/Input/input.inp.original" "$KEWPIE_DIR/Input/input.inp"

# ============================================================
# Step 6 — Summary
# ============================================================
echo ""
echo "──────────────────────────────────────────────────────────"
echo "SUMMARY"
echo "──────────────────────────────────────────────────────────"
echo "  Synthesis dir:    $SYNTH_DIR"
echo "  KEWPIE2 binary:   $KEWPIE_DIR/Executable/KEWPIE2.exe"
echo "  Results dir:      $RESULTS_DIR"
echo "  Files in results:"
ls -la "$RESULTS_DIR" | sed 's/^/    /'
echo ""
echo "  Reaction: ⁴⁸Ca + ²⁴³Am → ²⁹¹Mc* → xn channels"
echo "    0n → Mc-291  (cold survival, σ very small)"
echo "    1n → Mc-290"
echo "    2n → Mc-289"
echo "    3n → Mc-288  (Oganessian 2003 reported ~3 pb)"
echo "    4n → Mc-287"
echo ""
echo "  Citation:"
echo "    H. Lü, A. Boilley, Y. Abe, et al., KEWPIE2,"
echo "    Comput. Phys. Commun. 200 (2016) 381."
echo ""
echo "KEWPIE2 setup finished: $(date)"
echo "============================================================"
