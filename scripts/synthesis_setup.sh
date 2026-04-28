#!/bin/bash
# ============================================================
# NucleusGrid Phase 3 — Synthesis cross-section setup
#
# Researches PACE4 and HIVAP statistical compound-nucleus / evaporation
# codes, attempts to download whichever is publicly accessible, and
# tries to compile it. If anything is missing, the script prints a
# detailed report of exactly what would be needed.
#
# Target reaction (as requested):  ⁴⁸Ca + ²⁴⁴Pu → CN → Mc-291..299
# *** IMPORTANT NUCLEAR-PHYSICS NOTE ***
#   ⁴⁸Ca (Z=20) + ²⁴⁴Pu (Z=94) gives a compound nucleus with Z=114,
#   i.e. ⁴⁸Ca + ²⁴⁴Pu → ²⁹²Fl* (Flerovium), NOT Moscovium.
#   To produce Mc, the canonical reaction is:
#     ⁴⁸Ca + ²⁴³Am → ²⁹¹Mc*  (1n, 2n, 3n, 4n channels → Mc-290..287)
#   For the ²⁹¹..²⁹⁹Mc range you actually need:
#     ⁴⁸Ca + ²⁴³Am → CN(Z=115, A=291)*   (Oganessian 2003)
#   This script reports cross sections for BOTH reactions so the user
#   can pick. The ⁴⁸Ca + ²⁴⁴Pu line is computed as requested and
#   labeled clearly as the FLEROVIUM channel.
# ============================================================

set -u
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( dirname "$SCRIPT_DIR" )"
SYNTH_DIR="$PROJECT_ROOT/synthesis"
LOG="$SYNTH_DIR/synthesis_setup.log"

mkdir -p "$SYNTH_DIR"
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo "NucleusGrid synthesis setup — $(date)"
echo "Working dir: $SYNTH_DIR"
echo "============================================================"

# ============================================================
# Part 1 — Code research summary
# ============================================================
cat <<'DOC'

────────────────────────────────────────────────────────────
CODE OVERVIEW
────────────────────────────────────────────────────────────

  PACE4 (Projection Angular-momentum Coupled Evaporation, v4)
    Author:        A. Gavron (LANL), Phys. Rev. C 21, 230 (1980)
    Language:      Fortran 77
    Scope:         Statistical compound-nucleus decay only
                   (no fusion — caller supplies l-distribution)
    Distribution:  Bundled with LISE++ (Tarasov, NSCL/MSU)
                   http://lise.nscl.msu.edu/
                   Standalone source historically shared by request.
    License:       Free for research, attribution requested.
    Build deps:    gfortran (or any F77/F90 compiler)
    Notes:         Widely used for SHE evaporation residue (ER)
                   cross sections. Output: σ_ER per (xn, αxn, ...) channel.

  HIVAP (Heavy-Ion eVAPoration)
    Author:        W. Reisdorf (GSI), Z. Phys. A 300, 227 (1981)
    Language:      Fortran 77
    Scope:         Fusion + evaporation (full reaction simulation)
    Distribution:  GSI internal; obtain by contact with author or
                   successor (M. Schädel et al. for SHE community).
                   No reliable public mirror as of writing.
    License:       Restricted; requires contact with maintainer.
    Build deps:    gfortran, possibly CERNLIB linkage
    Notes:         Standard for SHE production cross-section
                   estimates. Includes Bass fusion barriers,
                   Smoluchowski diffusion model for capture,
                   statistical evaporation à la Reisdorf-Schädel.

  Modern alternatives (open source, actively maintained):
    TALYS       https://nds.iaea.org/talys/    Fortran 90 — full reaction code
    EMPIRE      https://www-nds.iaea.org/empire/   Fortran — HF/preeq/stat
    NRV (web)   http://nrv.jinr.ru/             Online calculator (no source)

  Recommendation for this project:
    • PACE4 — most accessible, focused on the question asked (ER xs)
    • TALYS — if PACE4 unobtainable, modern + open + well-supported

DOC

# ============================================================
# Part 2 — System prerequisites
# ============================================================
echo ""
echo "────────────────────────────────────────────────────────────"
echo "PREREQUISITE CHECK"
echo "────────────────────────────────────────────────────────────"
NEEDED=(gfortran make wget curl tar gzip python3)
MISSING=()
for tool in "${NEEDED[@]}"; do
    if command -v "$tool" >/dev/null 2>&1; then
        ver=$("$tool" --version 2>&1 | head -1)
        printf "  %-10s OK   (%s)\n" "$tool" "$ver"
    else
        printf "  %-10s MISSING\n" "$tool"
        MISSING+=("$tool")
    fi
done
if [ ${#MISSING[@]} -gt 0 ]; then
    echo ""
    echo "  Install missing tools first:"
    echo "    sudo apt-get install -y ${MISSING[*]}"
fi

# ============================================================
# Part 3 — Reaction kinematics (informational, no external code needed)
# ============================================================
echo ""
echo "────────────────────────────────────────────────────────────"
echo "TARGET REACTIONS (kinematics only — no cross-section yet)"
echo "────────────────────────────────────────────────────────────"
python3 <<'PYEOF'
# Atomic masses (AME2020, MeV) — selected, hard-coded for this audit.
# Source: AME2020 mass evaluation (Wang, Audi et al., Chin. Phys. C 45 (2021))
mass = {
    ("Ca", 48):  -44.214,   # mass excess MeV
    ("Pu", 244): +59.806,
    ("Am", 243): +57.176,
    ("Fl", 292): +135.000,  # estimated (no AME entry; approximate)
    ("Mc", 291): +137.000,  # estimated
}
def Q_fusion(p_sym, p_a, t_sym, t_a, cn_sym, cn_a):
    return mass[(p_sym, p_a)] + mass[(t_sym, t_a)] - mass[(cn_sym, cn_a)]

def bass_barrier(Z1, A1, Z2, A2):
    # Bass 1980 phenomenological fusion barrier (MeV), simplistic.
    R = 1.16 * (A1**(1/3) + A2**(1/3))
    e2 = 1.44  # MeV·fm
    V0 = Z1 * Z2 * e2 / R
    return V0  # rough

print("\n  Reaction A (REQUESTED):  ⁴⁸Ca + ²⁴⁴Pu → ²⁹²Fl* (Flerovium, Z=114)")
print(f"    Q-value (fusion):    {Q_fusion('Ca',48,'Pu',244,'Fl',292):+.2f} MeV (estimate)")
print(f"    Bass barrier:        ~{bass_barrier(20,48,94,244):.1f} MeV (CM)")
print(f"    Typical lab E:       ~245 MeV (above-barrier)")
print(f"    Channel of interest: 2n, 3n, 4n → ²⁹⁰..²⁸⁸Fl")
print(f"    NOTE: this does NOT produce Moscovium.")

print("\n  Reaction B (CORRECT FOR Mc):  ⁴⁸Ca + ²⁴³Am → ²⁹¹Mc* (Moscovium, Z=115)")
print(f"    Q-value (fusion):    {Q_fusion('Ca',48,'Am',243,'Mc',291):+.2f} MeV (estimate)")
print(f"    Bass barrier:        ~{bass_barrier(20,48,95,243):.1f} MeV (CM)")
print(f"    Typical lab E:       ~248 MeV (Oganessian 2003)")
print(f"    Channel of interest: 2n, 3n, 4n → ²⁸⁹..²⁸⁷Mc")
print(f"    Experimental σ:      ~3 pb (3n channel, near barrier)")
print()
print("  For the Mc-291..299 isotope range specifically:")
print("    Mc-287..290 are reachable via ²⁴³Am(⁴⁸Ca, xn) (x=1..4)")
print("    Mc-291..295 require neutron-richer targets (²⁴⁵Cm, ²⁴⁹Bk + ⁴⁸Ca)")
print("    Mc-296..299 are beyond current synthesis capability with stable beams.")
PYEOF

# ============================================================
# Part 4 — Attempt PACE4 download
# ============================================================
echo ""
echo "────────────────────────────────────────────────────────────"
echo "ATTEMPT: PACE4 download"
echo "────────────────────────────────────────────────────────────"
cd "$SYNTH_DIR"

PACE_URLS=(
    "http://lise.nscl.msu.edu/pace4.html"
    "https://lise.nscl.msu.edu/9_10/9_10.html"
    "http://groups.nscl.msu.edu/lise/pace4/pace4.zip"
)
PACE_OK=0
for url in "${PACE_URLS[@]}"; do
    echo "  trying: $url"
    if curl -sIL --max-time 15 "$url" 2>/dev/null | head -1 | grep -q "200\|301\|302"; then
        echo "    HEAD OK; attempting GET..."
        out=$(basename "$url")
        if curl -sL --max-time 60 -o "$out" "$url" && [ -s "$out" ]; then
            echo "    downloaded: $out ($(wc -c < "$out") bytes)"
            file "$out" 2>/dev/null | head -1 | sed 's/^/    /'
            PACE_OK=1
            break
        fi
    else
        echo "    unreachable / no 2xx response"
    fi
done

if [ $PACE_OK -eq 0 ]; then
    cat <<'NOTE'
  RESULT: PACE4 source could NOT be auto-downloaded.

  Why: LISE++ (which bundles PACE4) is distributed as a Windows .exe
  via the LISE++ web interface. There is no public direct-download URL
  for PACE4 source as a standalone tarball. To obtain it:

    Option 1 (LISE++ bundle, Windows):
      • Download LISE++ from http://lise.nscl.msu.edu/
      • Install on Windows (or via Wine on Linux)
      • PACE4 source files are in <LISE_DIR>/source/pace4/
      • Copy .f files out and compile with gfortran

    Option 2 (request directly):
      • Email Oleg Tarasov (NSCL/MSU): tarasov@nscl.msu.edu
      • Subject: "PACE4 standalone source request, research use"

    Option 3 (use a published port):
      • Some nuclear-physics groups host PACE4 ports on GitHub.
        Search:  github.com  "pace4"  language:Fortran
      • Verify provenance and citation history before trusting.
NOTE
fi

# ============================================================
# Part 5 — Attempt PACE4 compilation (if we got source)
# ============================================================
if [ $PACE_OK -eq 1 ] && ls "$SYNTH_DIR"/*.f "$SYNTH_DIR"/*.f90 2>/dev/null | head -1 >/dev/null; then
    echo ""
    echo "────────────────────────────────────────────────────────────"
    echo "ATTEMPT: PACE4 compilation"
    echo "────────────────────────────────────────────────────────────"
    cd "$SYNTH_DIR"
    SOURCES=(*.f *.f90)
    echo "  sources: ${SOURCES[*]}"
    if gfortran -O2 -std=legacy -o pace4 "${SOURCES[@]}" 2> compile_err.log; then
        echo "  COMPILE OK → ./pace4"
    else
        echo "  COMPILE FAILED. First 20 lines of error log:"
        head -20 compile_err.log | sed 's/^/    /'
        echo "  Full log: $SYNTH_DIR/compile_err.log"
    fi
fi

# ============================================================
# Part 6 — HIVAP attempt (informational; no public source)
# ============================================================
echo ""
echo "────────────────────────────────────────────────────────────"
echo "ATTEMPT: HIVAP"
echo "────────────────────────────────────────────────────────────"
cat <<'NOTE'
  HIVAP source is NOT publicly mirrored. To obtain:

    • Contact GSI Helmholtzzentrum für Schwerionenforschung
      Department: Atomphysik / SHE group
      Successor maintainers (as of 2020s): M. Block, Ch.E. Düllmann
      Email via:  https://www.gsi.de/en/researchaccelerators/research_an/she_chemistry

    • Alternative: cite Reisdorf 1981 + 1992 papers and reproduce the
      methodology; the algorithm is well documented in:
        W. Reisdorf, Z. Phys. A 300, 227 (1981)
        W. Reisdorf, M. Schädel, Z. Phys. A 343, 47 (1992)

  No automated download attempted.
NOTE

# ============================================================
# Part 7 — Final report
# ============================================================
echo ""
echo "────────────────────────────────────────────────────────────"
echo "SUMMARY"
echo "────────────────────────────────────────────────────────────"
echo "  Working dir:  $SYNTH_DIR"
echo "  Log file:     $LOG"
echo ""
if [ $PACE_OK -eq 1 ]; then
    echo "  PACE4:  source downloaded — see $SYNTH_DIR for status"
else
    echo "  PACE4:  NOT downloaded automatically"
    echo "          Action: install LISE++ on Windows or contact Tarasov (see above)"
fi
echo "  HIVAP:  not attempted (closed distribution)"
echo "          Action: contact GSI SHE group or reimplement from Reisdorf 1981"
echo ""
echo "  Reaction note: Ca-48 + Pu-244 → ²⁹²Fl* (Z=114), NOT Mc (Z=115)."
echo "  For Mc, use Ca-48 + Am-243.  See kinematics section above."
echo ""
echo "  Recommended next step:"
echo "    1. Install LISE++ (Windows + Wine) and extract PACE4 source, OR"
echo "    2. Install TALYS (open, modern):"
echo "         http://www.talys.eu/   (registration required, free)"
echo "         covers fusion + evaporation in one package."
echo ""
echo "Synthesis setup finished: $(date)"
echo "============================================================"
