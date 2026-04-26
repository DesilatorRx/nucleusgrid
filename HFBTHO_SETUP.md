# NucleusGrid Phase 2 — HFBTHO Setup Guide
## From Zero to Running Real HFB Calculations on Your 5950X

---

## What You Need

- Linux (Ubuntu 22/24, Debian, or WSL2 on Windows)
- gfortran compiler + BLAS/LAPACK
- ~15 minutes of setup
- ~5-8 hours of compute time for the full sweep

---

## Step 1: Install Dependencies (5 minutes)

```bash
sudo apt update
sudo apt install -y gfortran liblapack-dev libblas-dev \
                    libopenblas-dev git python3-pip

pip install numpy scipy matplotlib --break-system-packages
```

---

## Step 2: Get HFBTHO Source (5 minutes)

HFBTHO is GPL-licensed and published in Computer Physics Communications.
The official source is available via:

**Option A — CPC Program Library (official):**
```
https://www.sciencedirect.com/science/article/pii/S0010465522000868
# Download supplementary material → hfbtho_v4.zip
```

**Option B — GitHub fork (community maintained):**
```bash
git clone https://github.com/samarcanda93/hfbtho-newConstraints
# Note: this is v2/v3. For v4.0 features (SeaLL1), use Option A.
```

**Option C — OSTI (DOE open access):**
```
https://www.osti.gov/biblio/1310177
```

---

## Step 3: Compile HFBTHO (5 minutes)

```bash
cd hfbtho/   # wherever you extracted it

# Edit Makefile — set your compiler
# Look for COMPILER variable and set:
#   COMPILER = GFORTRAN

# Compile with OpenMP (uses all 32 threads on your 5950X)
make COMPILER=GFORTRAN USE_OPENMP=1

# Verify it compiled
ls -la hfbtho_main
./hfbtho_main --help 2>/dev/null || echo "Binary exists, ready to test"
```

**If you see LAPACK errors:**
```bash
# Try linking explicitly
make COMPILER=GFORTRAN USE_OPENMP=1 \
     LAPACK_LIBS="-L/usr/lib/x86_64-linux-gnu -llapack -lblas"
```

**Test the binary works:**
```bash
cd hfbtho/runs/
ln -sf hfbtho_NAMELIST_Cr50_N12.dat hfbtho_NAMELIST.dat
../hfbtho_main > test.out
grep "converged" test.out   # should say "HFBTHO converged"
grep "Total HFB" test.out   # should show a binding energy
```

---

## Step 4: Run the Wrapper (the fun part)

### Quick test — single isotope, single functional:
```bash
cd ~/nucleusgrid/phase2/

python3 hfbtho_wrapper.py \
    --binary /path/to/hfbtho_main \
    --isotope 299 \
    --functional UNEDF1 \
    --threads 32 \
    --verbose
```

### Validate against experiment first (recommended):
```bash
# Runs Mc-286 to Mc-290 vs AME2020 data
# Produces model_weights.json for Bayesian ensemble
python3 hfbtho_wrapper.py \
    --binary /path/to/hfbtho_main \
    --validate \
    --threads 32
```

Expected output:
```
  FRDM-2012   : RMSE=0.4xx MeV  weight=0.xxxx
  DRHBc       : RMSE=0.3xx MeV  weight=0.xxxx
  ...
```

### Full sweep — all isotopes, all functionals:
```bash
# ~5-8 hours on 5950X with 32 threads
# Saves incrementally every 10 WUs — safe to interrupt
python3 hfbtho_wrapper.py \
    --binary /path/to/hfbtho_main \
    --all \
    --threads 32 \
    --output hfbtho_results.json
```

### Convert results → Phase 2 Bayesian format:
```bash
python3 hfbtho_wrapper.py \
    --integrate \
    --output hfbtho_results.json
# Produces: phase2_hfbtho_integrated.json
```

---

## Step 5: Run Phase 2 Analysis on Real Results

```bash
# Once you have phase2_hfbtho_integrated.json:
python3 runner_phase2.py --use-hfbtho phase2_hfbtho_integrated.json
```

---

## Runtime Estimates on Your 5950X (32 threads)

| Scope | WUs | Est. Time |
|-------|-----|-----------|
| Validation only (Mc-286–290) | 50 | ~1 hr |
| Focused sweep (Mc-291–299, 3 β₂) | 135 | ~3 hrs |
| Full sweep (Mc-285–303, all β₂) | 1,425 | ~8 hrs |

One HFBTHO calculation for a superheavy nucleus (Z=115, N~184)
takes 60-180 seconds converged, with 20 HO shells and 600 max iterations.
Your 32-thread 5950X runs one at a time per calculation (OpenMP within),
so parallelism helps convergence speed, not WU throughput.

**For multi-machine speedup:** Run different isotope ranges on each machine,
combine JSON files before integration step.

---

## What the Output Tells You

After integration, `phase2_hfbtho_integrated.json` contains per isotope:

- `BE_mean_mev` — Bayesian-averaged binding energy (MeV)
- `Q_alpha_mean` — predicted alpha decay Q-value (MeV)
- `Q_alpha_std` — **inter-model uncertainty** (the publishable number)
- `shell_mean` — shell correction energy (MeV)
- `beta2_mean` — ground state deformation
- `by_functional` — per-model breakdown with weights

The `Q_alpha_std` field is the key result. Where it's large (>0.5 MeV),
models fundamentally disagree — those are the isotopes where
experimental data is most needed and where our paper makes its claim.

---

## Troubleshooting

**"HFBTHO not converged"** — Normal for highly deformed SHN. Try:
- Increase `max_iter` to 1000 in `generate_namelist()`
- Use `beta2=0.0` as starting point (spherical)
- Reduce `accuracy` to 1e-5

**Binding energy looks wrong** — Check:
- `n_shells=20` (must be large for SHN, default 14 is insufficient)
- Pairing cutoff 60 MeV (standard for SHN)

**Segfault on large nuclei** — Stack size issue:
```bash
ulimit -s unlimited
./hfbtho_main > output.out
```

---

## Citing the Tools

If this produces publishable results, cite:

1. HFBTHO v4.0: Perez et al., Comput. Phys. Commun. 276 (2022) 108367
2. UNEDF1: Kortelainen et al., Phys. Rev. C 85 (2012) 024304
3. SLy4: Chabanat et al., Nucl. Phys. A 635 (1998) 231
4. AME2020: Wang et al., Chin. Phys. C 45 (2021) 030003

And the methodology paper we're building is the output of this project.
