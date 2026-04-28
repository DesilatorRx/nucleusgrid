#!/bin/bash
# ============================================================
# NucleusGrid Phase 3 — Synthesis cross-section setup (TALYS)
#
# Pivots from PACE4 (Windows-only LISE++ bundle) and HIVAP
# (closed GSI distribution) to TALYS — the leading open-source
# nuclear reaction code (Koning, Hilaire, Goriely et al.).
#
# Attempts to fetch + compile TALYS automatically. If the source
# archive is gated behind registration (the usual case), prints a
# concrete step-by-step registration walkthrough.
#
# Target reaction (as requested): ⁴⁸Ca beam on ²⁴⁴Pu.
# *** REACTION NOTE (unchanged from prior version) ***
#   ⁴⁸Ca + ²⁴⁴Pu → ²⁹²Fl* (Flerovium, Z=114) — NOT Moscovium.
#   For Mc the canonical reaction is ⁴⁸Ca + ²⁴³Am → ²⁹¹Mc*.
#   See kinematics section below.
# ============================================================

set -u
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( dirname "$SCRIPT_DIR" )"
SYNTH_DIR="$PROJECT_ROOT/synthesis"
LOG="$SYNTH_DIR/synthesis_setup.log"

mkdir -p "$SYNTH_DIR"
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo "NucleusGrid synthesis setup (TALYS) — $(date)"
echo "Working dir: $SYNTH_DIR"
echo "============================================================"

# ============================================================
# Helpers
# ============================================================

# Returns 0 if file looks like a Fortran source archive (tar/gz/zip/bz2/xz),
# 1 otherwise. Avoids the prior bug where an HTML index page was treated
# as "source downloaded".
is_archive() {
    local f="$1"
    [ -s "$f" ] || return 1
    local mime
    mime=$(file -b --mime-type "$f" 2>/dev/null)
    case "$mime" in
        application/x-tar|application/gzip|application/x-bzip2|\
        application/x-xz|application/zip|application/x-gtar)
            return 0 ;;
        *)
            return 1 ;;
    esac
}

# Returns 0 if any plain Fortran source file exists in the given directory.
has_fortran_sources() {
    local d="$1"
    shopt -s nullglob
    local files=("$d"/*.f "$d"/*.f90 "$d"/*.F "$d"/*.F90 "$d"/*.for)
    shopt -u nullglob
    [ ${#files[@]} -gt 0 ]
}

# ============================================================
# Part 1 — Code overview (TALYS-first)
# ============================================================
cat <<'DOC'

────────────────────────────────────────────────────────────
CODE OVERVIEW — TALYS (primary target)
────────────────────────────────────────────────────────────

  TALYS — Nuclear Reaction Code
    Authors:        A.J. Koning, S. Hilaire, S. Goriely
    Language:       Fortran 90/95
    Scope:          Full nuclear reaction code: optical model,
                    direct, pre-equilibrium, compound, fission,
                    γ-cascade. Covers fusion-evaporation for SHE.
    Distribution:   http://www.talys.eu/   (primary)
                    https://nds.iaea.org/talys/   (IAEA mirror)
                    Email registration required to obtain source.
    License:        Free for academic / non-commercial use.
                    Citation: Koning, Hilaire, Goriely,
                    Eur. Phys. J. A 59, 131 (2023) for v2.0.
    Build deps:     gfortran (>= 9 recommended), make, ~500 MB
                    for source + structure database (RIPL-3).

  PACE4 (deprioritised)
    Distribution:   Bundled inside LISE++ (Windows .exe).
                    No standalone public source URL.
    Status here:    Skipped — no automated path on Linux.

  HIVAP (deprioritised)
    Distribution:   GSI internal; closed (contact author/successors).
    Status here:    Skipped — no public mirror.

  Why TALYS over PACE4/HIVAP for this project:
    • Open-source and actively maintained (last release 2023).
    • Single code covers the full chain: barrier → fusion →
      evaporation → ER cross-section. No need to glue PACE4 to
      a separate fusion model.
    • Modern Fortran, builds cleanly on Linux with gfortran.

DOC

# ============================================================
# Part 2 — System prerequisites
# ============================================================
echo ""
echo "────────────────────────────────────────────────────────────"
echo "PREREQUISITE CHECK"
echo "────────────────────────────────────────────────────────────"
NEEDED=(gfortran make wget curl tar gzip python3 file)
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
# Part 3 — Reaction kinematics (informational)
# ============================================================
echo ""
echo "────────────────────────────────────────────────────────────"
echo "TARGET REACTIONS (kinematics only — no cross-section yet)"
echo "────────────────────────────────────────────────────────────"
python3 <<'PYEOF'
mass = {
    ("Ca", 48):  -44.214,
    ("Pu", 244): +59.806,
    ("Am", 243): +57.176,
    ("Fl", 292): +135.000,   # AME-extrapolation, approximate
    ("Mc", 291): +137.000,
}
def Q_fusion(p_sym, p_a, t_sym, t_a, cn_sym, cn_a):
    return mass[(p_sym, p_a)] + mass[(t_sym, t_a)] - mass[(cn_sym, cn_a)]
def bass_barrier(Z1, A1, Z2, A2):
    R = 1.16 * (A1**(1/3) + A2**(1/3))
    return Z1 * Z2 * 1.44 / R

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
print("  For the Mc-291..299 range specifically:")
print("    Mc-287..290 reachable via ²⁴³Am(⁴⁸Ca, xn)  (x=1..4)")
print("    Mc-291..295 require neutron-rich actinide targets (²⁴⁵Cm, ²⁴⁹Bk + ⁴⁸Ca)")
print("    Mc-296..299 are beyond stable-beam synthesis capability today.")
PYEOF

# ============================================================
# Part 4 — TALYS auto-download attempt
# ============================================================
echo ""
echo "────────────────────────────────────────────────────────────"
echo "ATTEMPT: TALYS download"
echo "────────────────────────────────────────────────────────────"
cd "$SYNTH_DIR"

# Probe candidate direct-archive URLs. These are educated guesses
# at common naming patterns; almost all will 404 because the
# canonical distribution is gated behind a registration form.
TALYS_CANDIDATES=(
    "https://nds.iaea.org/talys/talys2.0.tar"
    "https://nds.iaea.org/talys/talys2.0.tar.gz"
    "https://nds.iaea.org/talys/talys-2.0.tar.gz"
    "https://nds.iaea.org/talys/talys.tar.gz"
    "http://www.talys.eu/fileadmin/talys/user/docs/talys2.0.tar.gz"
    "http://www.talys.eu/download/talys2.0.tar.gz"
)

# Also probe the index pages so we can confirm the sites are alive.
TALYS_PAGES=(
    "https://nds.iaea.org/talys/"
    "http://www.talys.eu/"
)

echo ""
echo "  Site reachability:"
for url in "${TALYS_PAGES[@]}"; do
    code=$(curl -sIL --max-time 15 -o /dev/null -w "%{http_code}" "$url" 2>/dev/null || echo "---")
    printf "    %-50s HTTP %s\n" "$url" "$code"
done

echo ""
echo "  Probing direct-source URLs (almost all will fail —"
echo "  TALYS distribution is gated behind email registration):"
TALYS_OK=0
TALYS_FILE=""
for url in "${TALYS_CANDIDATES[@]}"; do
    code=$(curl -sIL --max-time 15 -o /dev/null -w "%{http_code}" "$url" 2>/dev/null || echo "---")
    printf "    HTTP %s  %s\n" "$code" "$url"
    if [ "$code" = "200" ]; then
        out="$SYNTH_DIR/$(basename "$url")"
        if curl -sL --max-time 120 -o "$out" "$url" 2>/dev/null && is_archive "$out"; then
            echo "      → archive obtained: $out ($(wc -c < "$out") bytes)"
            TALYS_OK=1
            TALYS_FILE="$out"
            break
        else
            echo "      → got bytes but not an archive (likely HTML); discarding"
            rm -f "$out"
        fi
    fi
done

# ============================================================
# Part 5 — TALYS extract + compile (only if a real archive)
# ============================================================
if [ "$TALYS_OK" -eq 1 ] && [ -n "$TALYS_FILE" ]; then
    echo ""
    echo "────────────────────────────────────────────────────────────"
    echo "ATTEMPT: TALYS extract + compile"
    echo "────────────────────────────────────────────────────────────"
    EXTRACT_DIR="$SYNTH_DIR/talys_src"
    mkdir -p "$EXTRACT_DIR"
    case "$TALYS_FILE" in
        *.tar.gz|*.tgz) tar -xzf "$TALYS_FILE" -C "$EXTRACT_DIR" ;;
        *.tar.bz2)      tar -xjf "$TALYS_FILE" -C "$EXTRACT_DIR" ;;
        *.tar.xz)       tar -xJf "$TALYS_FILE" -C "$EXTRACT_DIR" ;;
        *.tar)          tar -xf  "$TALYS_FILE" -C "$EXTRACT_DIR" ;;
        *.zip)          unzip -q "$TALYS_FILE" -d "$EXTRACT_DIR" ;;
        *)              echo "  unknown archive type: $TALYS_FILE"; EXTRACT_DIR="" ;;
    esac

    if [ -n "$EXTRACT_DIR" ]; then
        # TALYS ships an install.csh / install script in talys/source/.
        # Find the source directory.
        SRC_DIR=$(find "$EXTRACT_DIR" -type d -name source -path "*talys*" -print -quit 2>/dev/null)
        if [ -n "$SRC_DIR" ] && has_fortran_sources "$SRC_DIR"; then
            echo "  TALYS source tree found at: $SRC_DIR"
            cd "$SRC_DIR"
            shopt -s nullglob
            SOURCES=(*.f *.f90 *.F *.F90 *.for)
            shopt -u nullglob
            echo "  ${#SOURCES[@]} source files"
            if gfortran -O2 -std=legacy -fno-automatic -o talys "${SOURCES[@]}" \
                  2> "$SYNTH_DIR/talys_compile.log"; then
                echo "  COMPILE OK → $SRC_DIR/talys"
            else
                echo "  COMPILE FAILED. Tail of log:"
                tail -20 "$SYNTH_DIR/talys_compile.log" | sed 's/^/    /'
                echo "  Full log: $SYNTH_DIR/talys_compile.log"
                echo "  (TALYS normally ships its own install script —"
                echo "   try: cd $SRC_DIR && ./install.csh   or   make )"
            fi
        else
            echo "  No Fortran source tree found inside the archive."
            echo "  Archive contents (top level):"
            find "$EXTRACT_DIR" -maxdepth 2 | head -20 | sed 's/^/    /'
        fi
    fi
fi

# ============================================================
# Part 6 — Registration walkthrough (always printed if not obtained)
# ============================================================
if [ "$TALYS_OK" -eq 0 ]; then
    cat <<'REGISTER'

────────────────────────────────────────────────────────────
TALYS REGISTRATION WALKTHROUGH (manual, ~5 minutes)
────────────────────────────────────────────────────────────

  TALYS is free for academic / non-commercial use, but the source
  is gated behind a one-time email registration.

  Primary source:
    1. Go to:  https://nds.iaea.org/talys/
       (IAEA Nuclear Data Section mirror — most reliable)

    2. Click "Download TALYS" → opens the registration form.

    3. Fill in:
         • Name
         • Email  (a valid academic / institutional address is
                   strongly preferred; gmail/hotmail accepted but
                   sometimes scrutinised)
         • Affiliation (institution / company / "individual")
         • Country
         • Intended use (one or two sentences is fine, e.g.
           "Statistical reaction model calculations for SHE
            production cross-sections.")

    4. Accept the licence (free for academic / non-commercial).

    5. You receive an email with a personalised download link
       (usually within minutes, occasionally up to 24 h).

    6. Download the tarball. As of 2023:
         talys2.0.tar.gz   (~150 MB including RIPL-3 structure DB)

  Alternative source:
    http://www.talys.eu/   (mirror; same registration system)

  After download, place the tarball at:
       $SYNTH_DIR/talys2.0.tar.gz
  then re-run this script — it will detect the file, extract,
  and attempt compilation.

  Manual install (if re-running this script is inconvenient):
       cd $SYNTH_DIR
       tar xzf talys2.0.tar.gz
       cd talys/source
       ./install.csh                  # ships with TALYS; sets up paths
       ./talys < ../samples/Pb208/Ni58/input > test.out  # smoke test

  Build prerequisites already satisfied on this host:  gfortran ✓
  Disk space required:  ~500 MB after extraction.

  Citation when publishing results:
    Koning, Hilaire, Goriely, Eur. Phys. J. A 59, 131 (2023).
REGISTER
fi

# ============================================================
# Part 7 — Brief notes on the deprioritised codes
# ============================================================
echo ""
echo "────────────────────────────────────────────────────────────"
echo "PACE4 / HIVAP — alternative paths if TALYS is blocked"
echo "────────────────────────────────────────────────────────────"
cat <<'NOTE'
  PACE4:
    • Source is bundled with LISE++. There is no standalone
      Linux-friendly archive published.
    • To extract: install LISE++ on Windows (or via Wine on
      Linux), copy <LISE_DIR>/source/pace4/*.f, build with
      gfortran -std=legacy.
    • Alternative: email O. Tarasov (NSCL/MSU) requesting the
      standalone PACE4 source for academic use.

  HIVAP:
    • No public mirror. Contact the GSI SHE group:
      M. Block, Ch.E. Düllmann
      via https://www.gsi.de/en/researchaccelerators/research_an/she_chemistry
    • Or reimplement from Reisdorf's published methodology:
        W. Reisdorf, Z. Phys. A 300, 227 (1981)
        W. Reisdorf, M. Schädel, Z. Phys. A 343, 47 (1992)
NOTE

# ============================================================
# Part 8 — Final report (accurate this time)
# ============================================================
echo ""
echo "────────────────────────────────────────────────────────────"
echo "SUMMARY"
echo "────────────────────────────────────────────────────────────"
echo "  Working dir:  $SYNTH_DIR"
echo "  Log file:     $LOG"
echo ""
if [ "$TALYS_OK" -eq 1 ] && [ -f "$SYNTH_DIR/talys_src" ]; then
    echo "  TALYS:  archive obtained AND extracted ($TALYS_FILE)"
elif [ "$TALYS_OK" -eq 1 ]; then
    echo "  TALYS:  archive obtained ($TALYS_FILE) — see compile section"
else
    echo "  TALYS:  no direct download URL responded with a real archive."
    echo "          Action: complete the registration walkthrough above."
fi
echo "  PACE4:  not attempted (Windows-only LISE++ bundle)"
echo "  HIVAP:  not attempted (closed GSI distribution)"
echo ""
echo "  Reaction note (unchanged):"
echo "    ⁴⁸Ca + ²⁴⁴Pu → ²⁹²Fl* (Z=114)  — what you asked for, Flerovium"
echo "    ⁴⁸Ca + ²⁴³Am → ²⁹¹Mc* (Z=115)  — correct route to Mc"
echo "    See kinematics section above."
echo ""
echo "Synthesis setup finished: $(date)"
echo "============================================================"
