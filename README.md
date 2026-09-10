# quantum-linear-rank-cone

Machine-certified data and tools for the **five-party quantum linear rank cone**
(QLR₅): an exact H-representation, certified extreme rays, and a self-test
suite of verification gates. Self-contained: no `psitip` dependency, all input
data bundled, a static Normaliz binary included, resumable long runs.

> **Status: work in progress.** Everything in `data/` is machine-certified by
> the gates in `selftest.py`, but several transcription and literature claims
> are still awaiting independent human verification (see §7 and
> `CHANGELOG.md`, which also records all errata). Treat the repository as a
> research notebook with reproducible artefacts, not as a refereed publication.

The repository grew out of an internal project (code name *EPR-1*, "complete
characterisation of the quantum entropy cone", plan A: reduce the n = 5 case to
decisions about three cones). Internal shorthand that survives in file names
and the changelog is explained in the glossary at the end.

---

## 0. Status at a glance (read this first)

**Decided / completed (do not redo):**

- [x] **psitip's linear-rank list is incomplete (round R3, reproduced on two
  machines).** The set of five-variable linear rank inequalities shipped with
  `psitip` (`linear_bound`) is missing one class, **Ingleton(39)** (the
  "overlapping" substitution form of Ingleton, arXiv:0910.0284 eq. (39)). The
  verdict `PSITIP-INCOMPLETE` is re-derived by `scripts/s12_dfz_reconcile.py`.
  The fix is bundled: `CLR_H_fixed` (1,905 rows) is the classical cone's
  irredundant H-representation and `templates28` / **pure28** (23,355 rows) the
  quantum-side construction in current use.
- [x] **Headline 18/19.** Under pure28, 18 of the 19 orbits of extreme rays of
  the five-party holographic entropy cone (HEC₅) remain extreme; the only
  non-extreme one is **ray #19** (exact tight-rank 29), in agreement with the
  count stated in the literature. **R11 turned this into a canonical
  certificate:** the minimal face containing #19 is two-dimensional, both edges
  were computed exactly, and **r₁₉ = 1·q₁ + 2·q₂** with q₁ = HEC ray #4 and
  both edges of rank 30 (`data/cert19_exact.npz`, guarded by gate G14).
- [x] **31 ↔ 31 facet matching, machine-certified (R12).** The "31 QLR
  inequalities" of Bao–Cheng–Hernández-Cuenca–Su (BCHS) are the
  Dougherty–Freiling–Zeger (DFZ) list with classical monotonicity replaced by
  weak monotonicity (pinned to the paper's footnote). The naive family
  (DFZ 28 + subadditivity elementals + weak monotonicity) has **34** S₆-classes;
  exactly 3 weak-monotonicity classes are valid but redundant; the remaining
  **31 coincide, as a set, with the 31 facet classes of pure28**. No new
  transcription was needed.
- [~] **A1 campaign, day one (R15, 2026-09-10).** `scripts/s4c_qlr.py` runs a
  symmetry-aware adjacency decomposition on the exact irredundant
  H-representation of QLR₅ (`qlr5_H_facets10860`, 10,860 rows = the facet
  instance count found in R11), seeded with the 59 known orbits plus q₂. In
  170 s, expanding the 12 lightest vertex figures, it reached **1,924 certified
  extreme-ray orbits** — far from saturation (91 % of the last figure's 1,067
  neighbours were new orbits; coordinates up to 162). **Conclusion on scale:
  QLR₅ has orders of magnitude more extreme rays than CLR₅, and full
  enumeration is out of reach for this engine** — which is why BCHS reported
  the enumeration as computationally infeasible in 2021. Statable results: a
  **lower bound of ≥ 1,924 orbits (59 previously known)**, the certified partial
  catalogue `qlr5_orbits_partial` (gate G17), and the S7 priority test set
  `qlr5_new_small200` (the 200 new orbits with the smallest coordinates,
  maximum coordinate 3..13). **Interpretation guard:** none of the 1,864 new
  orbits is witnessed by any known realisable pool — this only says that the
  pools' entropy scale (≤ 12) cannot reach rays with coordinates up to 162;
  **it is not evidence about the S7 conjecture.**
- [x] **The 59 = 40 + 17 + 2 ledger reproduced end to end (R13).** Under the
  S₆ (purification) symmetry the overlap between HEC and CLR orbits is exactly
  HEC #1; among the 26 S₆-orbits of six-qubit graph-state vectors exactly 2 are
  new QLR-extreme orbits (G11 and G15, the latter identified by the minimum
  −1 of inequality (4.4) over its full orbit). The 59 representatives are
  bundled (`qlr59_reps`, gate G16). **Theorem note:** equality of the 31 facet
  classes implies **cone(pure28) = QLR₅**; q₂ is a genuine extreme ray of QLR₅
  that evades every realisable pool we have, i.e. a **concrete test case for
  the S7 conjecture** (stabilizer₅ = QLR₅) and the primary object for A2.
- [x] **The "18th ray" question settled (R12).** Table 3 of Hernández-Cuenca's
  PRD letter was machine-parsed and matched to our rays as an S₆-orbit-level
  bijection: **the paper's ray 18 is our #19** (gate G15). The headline, the
  exact certificate r₁₉ = HEC#4 + 2·q₂ and the literature now interlock; the
  earlier worry came from the HEC database and Table 3 swapping rows 18/19 and
  using different orbit representatives.
- [x] **CLR₅ enumeration complete (R7): all 162 orbits of extreme rays found;
  orbit sizes sum to 7,943 rays**, matching the literature on both counts;
  every representative certified extreme (tight rank 30). Engine history:
  Normaliz primal mode blew up (R4) → a full lrs reverse search was infeasible
  (12.7 h on the server, 97 M bases, 12 % of the rays; the estimator was off by
  a factor 68, erratum C3, R6) → **symmetry-aware adjacency decomposition
  (s4c) finished in about thirteen minutes**. Data bundled
  (`clr5_orbit_reps162`, `clr5_rays7943`, gate G12).
- [x] **s6 complete (R7): all 162 representatives lie in the quantum cone and
  exactly 40 stay QLR-extreme** — the ledger's headcount, certified with exact
  Gram-matrix ranks (non-extreme rank distribution 25×56, 26×25, 27×18, 28×23).
- [x] **Completeness cross-certification closed at the logical level (R9).**
  DFZ §4 states the counts verbatim (7,943 rays / 162 orbits); the set of
  extreme rays of a cone is unique; our 162 orbits each carry an exact
  certificate; hence the two sets coincide and a vector-by-vector diff against
  DFZ's `rays5` file is logically redundant. The three example rays printed in
  the paper all hit our 7,943-ray set coordinate by coordinate. The `rays5`
  file is an optional reinforcement only (see §7). 122 of the 162 vertex
  figures are fully closed and the total 7,943 agrees, which independently
  corroborates the published count; a from-scratch `mplrs` run is a luxury for
  publication only.
- [x] The 27-class mutual irredundancy certificates and the pure(27) 16/19
  result are historical baselines kept as regression gates.

**Open items, in order (details in §7):** continue the A1 campaign on a
server (`s4c_qlr.py`); human spot check of the 11 single-source lines in
`data/dfz_ref28.csv`; stabilizer realisability of the small-coordinate new
rays (S7 front line); optional `rays5` cross-check.

---

## 1. Requirements

- Python ≥ 3.8, `pip install -r requirements.txt` (numpy, scipy).
- `gcc` for the C helpers used by s7 (OpenMP optional but recommended; without
  it the tools run single-threaded and print a warning).
- `lrs` only for re-enumeration / seeding (optional): `apt-get install -y lrslib`
  (or conda-forge `lrslib`). For cross-validation a **static Normaliz binary
  is bundled at `bin/normaliz` (v3.11.1, GPL-3)**; scripts resolve
  `$NORMALIZ` → `PATH` → bundled binary. To pin the bundled one (avoid another
  version on `PATH` winning; the R4 server actually used conda's 3.11.0):
  `NORMALIZ=$PWD/bin/normaliz bash scripts/s4_run_normaliz.sh ...`
- Linux x86_64 only (the parallel layer relies on `fork`; the bundled binary is
  statically linked for Linux).

## 2. Quick start (reference timings on 25 vCPU / 90 GB)

```
# ---- A. Verify the completed results (optional, ~3 minutes in total) ----
python3 selftest.py --full            # 20 gates, ~35 s; must be all green before anything else
python3 scripts/s1_build_qlr.py --variant pure28 --workers 25     # seconds
python3 scripts/s2_judge.py QLR_H_pure28.npy                      # expect: violated 0
python3 scripts/s3_rank19.py QLR_H_pure28.npy                     # expect: 18/19, only #19 has rank 29
python3 scripts/s12_dfz_reconcile.py                              # expect: PSITIP-INCOMPLETE / ing39

# ---- B. Server session, one command ----
sh scripts/run_b_batch.sh
#   = s7 (two layers, 0 violations) -> s8 (#19 certificate; unconditional pool first,
#     automatic escalation to the conditional layer, labelled) -> s10 (facet classes N +
#     redundant M = 84) -> s11 packaging.  Commit or push the s11 tarball afterwards.

# ---- C. A1 campaign on the QLR5 cone (hours; resumable) ----
python3 scripts/s4c_qlr.py --init data/qlr_seeds60.npy --state qlr_adj.npy
python3 scripts/s4c_qlr.py --state qlr_adj.npy --budget 14400 --rep-timeout 1800 --max-tight 400

# ---- D. Optional reinforcement: vector-by-vector diff against DFZ's rays5 ----
python3 scripts/s5b_diff_rays.py rays5     # expect the last line to read MATCH
```

Expected numbers: s2 = 0 violations; s3 = 18/19 with only #19 at rank 29
(pure28 sha16 `827079728c051bfe`); s7 = 0 violations in both layers; s10 prints
the facet-class count for the BCHS-31 reconciliation; s5b ends with `MATCH`.
**The CLR₅ enumeration is finished:** 162 orbits / 7,943 rays ship with the
repository (gate G12); there is no need to rerun it. Optional re-derivation
path: `s4_clr_prepare` → a short `s4b` (lrs) run for seeds → `s4c_adjacency`
batches until closure.

## 3. Pipeline stages

| Stage | Command essentials | Purpose / expectation |
| --- | --- | --- |
| selftest | `selftest.py [--full] [--workers N]` | 20 regression gates; run first on any machine |
| s1 | `--variant pure28` (current) / `pure`, `v3`, `pure2` (historical) | build the H-representation; seconds |
| s2 | `s2_judge.py <H.npy>` | 760 graph-state judge; `pure*` variants must give 0 violations |
| s3 | `s3_rank19.py <H.npy> [--extra-rows X]` | tight-rank test of the 19 HEC rays; pure28 → 18/19 |
| s4 | **`s4c_adjacency.py`** (`--init` seeds → `--state` batches) | CLR₅ extreme-ray enumeration engine (symmetry-aware adjacency decomposition; 162 orbits in R7); `s4b` (lrs) / `s4` (Normaliz) for seeds and cross-validation |
| s4c_qlr | **`s4c_qlr.py`** (`--H`, `--sym s6`, `--max-tight`) | the same engine on the QLR₅ cone (A1 campaign) |
| s5 | `s5_orbits.py clr5.out [--raw] [--expect 162]` | rays → S₅ orbits; reconcile against 162 |
| s5b | `s5b_diff_rays.py rays5` | one-command cross-check against DFZ's published ray list |
| s6 | `s6_sweep_clr_orbits.py reps.npy --H pure28` | CLR rays through the quantum cone; ledger count 40 |
| s7 | `--H pure28 [--f3-limit N] [--gf2-samples N]` | grow the realisable pools (full F₃ enumeration + 12-vertex GF(2) sampling) |
| s8 | `--H pure28 --pool <s7 outputs...>` | face-restricted decomposition certificate for **#19** (default target `targets1`) |
| s9 | `s9_sixvar.py <hand-transcribed CSV>` | six-variable inequalities (fallback; transcribe from the literature only) |
| s10 | `--H pure28 [--limit N]` | facet classes → BCHS-31 reconciliation |
| s11 | `s11_collect.py <globs> --out tar.gz` | package results with canonical sha manifest |
| s12 | `s12_dfz_reconcile.py` (bundled CSV by default) | re-derive the R3 verdict (PSITIP-INCOMPLETE / ing39) |

## 4. Monitoring progress

- selftest: one PASS/FAIL line per gate.
- s1 and other parallel builds: a `k/n chunks + ETA` line at most every 5 s
  (ETA shows `--` for the first 2 s).
- s4c / s4c_qlr: one line per expanded representative
  (`tight=… neighbours / new orbits / total / expanded`), state saved after
  every expansion (`--state`), interruptible and resumable; giant figures that
  exceed `--rep-timeout` are skipped and reported explicitly.
- lrs / Normaliz (seeding or cross-validation): `wc -l *.lrs.out`,
  `tail -f clr5.log`; stop with `pkill -x lrs` / `pkill -x normaliz`
  (**never** `pkill -f`, which once killed the wrong process in R2).
- s7: both C tools report counts on stderr (F₃: one line per 2²⁰
  configurations with a percentage; GF(2): one line per 1 M samples);
  byte check: the full F₃ enumeration writes 14,348,907 × 31 ≈ 445 MB.
- s10: one line per 10 classes; a resumed shard reports done/remaining counts.
- Memory guard: `EPR1_MAX_RSS_GB` applies to this package's workers only,
  **not to Normaliz**.

## 5. Pitfalls (real incidents first)

1. **Wrong variant.** `--variant pure` is the historical 27-template
   construction (16/19). Always use `pure28`. This happened once (R3 server
   logs); if s3 reports 16/19 with #17/#18 at rank 29, check this first.
2. **No `normaliz` for s4.** The script fails with an explicit message and two
   remedies (bundled binary / conda). If the bundled binary reports a
   permission error: `chmod +x bin/normaliz`.
3. **Enumeration engine history (three verdicts; do not go back).** Normaliz
   primal mode exploded in the intermediate stages (91.5 M hyperplanes at
   generator 62/1905, R4); a full lrs reverse search has an infeasible tail
   (12.7 h / 97 M bases / 12 % of the rays; the depth-2 estimator was off by
   68×, erratum C3, R6); **symmetry-aware adjacency decomposition finished in
   about thirteen minutes (R7)**. The results ship with the repository; normal
   use needs no enumeration at all.
4. **s5 MISMATCH.** First confirm that the input came from `CLR_H_fixed`
   (1,905 rows; `clr5.in` header `inequalities 1905`). The old `CLR_H`
   (1,875 rows) falls short.
5. **DFZ's raw `rays5` data.** `code.ucsd.edu/zeger/linrank/` disallows
   automated access; download it manually in a browser (plain HTTP, port 80;
   some networks block it). Then `scripts/s5b_diff_rays.py rays5`.
6. **First ETA line.** The extrapolation after the first finished shard is
   unreliable; it converges within seconds.
7. **s7 without gcc / without OpenMP.** Explicit error, respectively a
   downgrade warning; see §1.
8. **s9 iron rule.** Six-variable inequalities may only be transcribed from the
   literature (template `data/sixvar_template.csv`); rows that fail the
   per-variable balance check are rejected — that is almost always a copying
   error.
9. **Estimators lie.** The lrs depth-2 tree-size estimate systematically
   underestimates on highly unbalanced trees (68× measured in R6). Base
   long-run decisions on checkpoint base counts and ray-output rates, never on
   the estimate.
10. **Giant vertex figures.** For hyper-symmetric rays with ≳ 1,000 tight
    rows the vertex figure is as hard as the original cone; s4c skips them on
    `--rep-timeout` and reports "closure incomplete" honestly.
11. **Hour-scale jobs belong on a server.** Interactive container/notebook
    sessions reap background processes between turns (observed twice in R5,
    both truncated mid-line). lrs is stateless — just restart; rays already
    streamed are genuine (the s5 parser drops a truncated last line).
12. **Batch scripts and `tee`.** `cmd | tee log` hides the exit status from
    `set -e` (R10 ran through two crashed steps and still printed DONE).
    `run_b_batch.sh` now checks dependencies up front and verifies each step's
    output file.
13. **Patching with `str.replace` must assert a hit.** A silent no-op patch
    went unnoticed for three rounds (erratum C4).

## 6. Data inventory (`data/`)

| File | Rows | Meaning |
| --- | :-: | --- |
| `M_LR.npy` | 1,790 | LR instance family extracted from psitip (27 classes; historical) |
| `CLR_H.npy` / `CLR_H_fixed.npy` | 1,875 / **1,905** | classical cone H-representation (fixed = + the 30 instances of the Ingleton(39) orbit) |
| `REF28.npy` + `dfz_ref28.csv` | 28 | reference family (24 DFZ inequalities + 4 Ingleton forms); the CSV carries source tags and the list of single-source lines awaiting human check |
| `templates27.npy` / `templates28.npy` | 27 / 28 | substitution templates (28 = current) |
| `ING39_orbit.npy` | 30 | S₅ instances of Ingleton(39) |
| `graphstate_vecs.npy` | 760 | entropy vectors of the six-vertex GF(2) graph states (judge) |
| `hec5_rays_maskorder.npy` | 19 | HEC₅ extreme-ray orbit representatives (mask coordinate order) |
| `targets1.npy` / `targets3.npy` / `targets7.npy` | 1 / 3 / 7 | s8 targets: #19 (current) / #17,18,19 / historical R2 target set |
| `sep_certs_27.npz` | 27 | R3 mutual-irredundancy exact integer certificates (historical) |
| `clr5_orbit_reps162.npy` | **162** | CLR₅ extreme-ray orbit representatives (R7; each certified rank 30) |
| `clr5_rays7943.npy` | **7,943** | all CLR₅ extreme rays (orbit union; both literature counts match; gate G12) |
| `cert19_exact.npz` | — | exact certificate r₁₉ = α·q₁ + β·q₂ (α = 1, β = 2; q₁ = HEC #4; gate G14) |
| `pure28_facet31_reps.npy` | 31 | the 31 facet classes of pure28 (= QLR₅), S₆ representatives |
| `qlr59_reps.npy` | 59 | the 40 + 17 + 2 ledger of known QLR₅ extreme-ray orbits (gate G16) |
| `a3_hyp_ineq44.npy` | 1 | BCHS inequality (4.4), the known hypergraph/stabilizer separator (unbalanced) |
| `shc_table3_rays.npy` + `shc_table3_map.json` | 19 | Table 3 of the PRD letter, machine-parsed, with the orbit-level bijection to our numbering (gate G15) |
| `qlr5_H_facets10860.npy` | 10,860 | exact irredundant H-representation of QLR₅ (all facet instances of the 31 classes) |
| `qlr_seeds60.npy` | 60 | A1 seeds: the 59 known orbits + q₂ |
| `qlr5_orbits_partial.npy` | **1,924** | certified partial catalogue of QLR₅ extreme-ray orbits (A1 day one; gate G17; wide int16 sha) |
| `qlr5_new_small200.npy` | 200 | S7 priority test set: new orbits with the smallest coordinates |
| `sixvar_template.csv` | — | transcription template for s9 |
| `manifest_shas.json` | — | canonical sha256 of every row family (`sha_rows`; `sha_rows_wide` for wide coordinates) |

Build products (`QLR_H_*.npy`) are not committed; s1 regenerates them in
seconds and checks the sha.

## 7. What next, and how

**⓪ A1 campaign on the QLR₅ cone (server, hours, resumable).** Commands in
§2-C. The goal is not closure (out of reach) but a tighter lower bound, more
small-coordinate test cases and a growth curve (record "total orbits vs.
expanded" per batch); commit the `qlr_adj.npy` state file. Per the project
plan's stop-loss clause, A1 is restated as **A1′: certified partial catalogue +
lower bound + structural statistics.**

**① Completeness cross-certification of the 162 (closed at the logical level,
R9; file optional).** Proposition: the extreme-ray set of a cone is an
intrinsic invariant. Our 162 orbit representatives are pairwise inequivalent
and each carries an exact rank-30 certificate (so they lie in ExtRays(C)); DFZ
§4 states that C has exactly 7,943 extreme rays in 162 orbits; hence the two
sets coincide. Spot reinforcement: the three explicit rays printed in the paper
(the F³ five-subspace example, U₂,₄, U₂,₅) all hit our 7,943-ray set.
Honest note: the count premise is DFZ's own cddlib computation (2009, 2–3 days
per 31-dimensional run); our 122/162 closed vertex figures and the matching
total corroborate it independently. Optional reinforcement: obtain `rays5`
from `http://code.ucsd.edu/zeger/linrank/rays5` (browser / `curl -O` from a
network that allows it, or by e-mail to the authors) and run
`scripts/s5b_diff_rays.py rays5`.

**② Server session (one command, ~25 minutes).** `sh scripts/run_b_batch.sh`.
Checkpoints: s7 prints `pure violations 0` for both layers and writes
`s7_out/realizable_gf2_12v.npy` (unconditional layer) and
`s7_out/realizable_f3_conditional.npy` (**conditional layer**; its validity
premise is documented, do not mix the two); s8 uses the unconditional pool
first — `certificate saved` writes `s8_certs.npz` (the λ weights are the
certificate); on `no certificate` the script retries with the conditional
layer merged in and stores the result **separately** as
`s8_certs_conditional.npz` with a premise note (the two certificates have
different rigour and must be reported separately); s10 ends with
`facet classes N  redundant M` (N + M = 84); s11 packages everything with the
three logs.

**③ Human checks pending.** (a) The 11 single-source lines of
`data/dfz_ref28.csv` (`dfz10, 11, 13, 14, 15, 16, 19, 20, 21, 22, 23`) against
DFZ eqs. (10)–(23) — the CSV comment lines print the original I-notation; the
ar5iv rendering of arXiv:0910.0284 is a convenient second source. (b) The
stabilizer realisability of the small-coordinate new orbits in
`qlr5_new_small200` (the practical front line of the S7 conjecture; 31 new
orbits have maximum coordinate ≤ 5). (c) Machine verification of the 64-row
contraction table proving inequality (4.4).

**④ Optional luxury: from-scratch completeness proof.** A parallel `mplrs`
run (mplrs + OpenMPI, ~25 processes, days, checkpointable) or a separate
attack on the ~40 giant vertex figures (tight ≳ 1,000). Only if a publication
requires it.

## 8. Versions

See `CHANGELOG.md` and `VERSION`. The changelog also records the errata
(C3: lrs estimator; C4: silent patch failure; C5: report dates of rounds
R4–R13).

## 9. Repository workflow

This repository is public at
`https://github.com/RuifengCao/quantum-linear-rank-cone`.

- **CI:** `.github/workflows/selftest.yml` installs `lrslib` and the Python
  requirements and runs `python3 selftest.py --full` on every push and pull
  request (all 20 gates must pass).
- **`.gitignore`** excludes rebuildable products (`QLR_H_pure28.npy`), the s7
  pools, logs and campaign state files. Large result files (> 50 MB pools)
  should go to Git LFS or to release assets; small results can be committed to
  a `results/` branch.
- **Server workflow:** `git clone` → `pip install -r requirements.txt` →
  `sh scripts/run_b_batch.sh` or the A1 campaign commands → commit the state
  and log files.
- **Citing:** `CITATION.cff` holds the metadata; a Zenodo DOI will be attached
  to the first tagged release.
- **Licences:** code MIT, data CC-BY-4.0, `bin/normaliz` is an unmodified
  GPL-3 binary of Normaliz (provenance in `bin/`).

## 10. Glossary of internal shorthand

- **Rn** — round n of the project iteration (see `CHANGELOG.md`).
- **Cn** — erratum n (C3, C4, C5 are recorded in the changelog).
- **Gn** — gate n of `selftest.py`.
- **A1 / A2 / A3** — the three sub-problems of the project plan: A1 enumerate
  the extreme rays of QLR₅ (now A1′); A2 decide whether the five-party
  stabilizer cone equals QLR₅ (the **S7 conjecture** of BCHS); A3 find balanced
  facets separating the five-party hypergraph cone from the stabilizer cone.
- **CLR₅ / QLR₅** — the classical and quantum linear rank cones on five
  parties (31 coordinates); QLR₅ is cut out by the DFZ inequalities with
  classical monotonicity replaced by weak monotonicity.
- **pure28** — the 23,355-row construction built by s1 from the 28 templates;
  its facets coincide with the 31 QLR₅ classes, so cone(pure28) = QLR₅.
- **HEC₅** — the five-party holographic entropy cone (19 extreme-ray orbits).
- **DFZ** — Dougherty, Freiling, Zeger, *Linear rank inequalities on five or
  more variables*, arXiv:0910.0284.
- **BCHS** — Bao, Cheng, Hernández-Cuenca, Su, *A gap between the hypergraph
  and stabilizer entropy cones*, arXiv:2006.16292.
- **SHC** — Hernández-Cuenca, *Holographic entropy cone for five regions*,
  Phys. Rev. D 100, 026004 (2019), arXiv:1903.09148; HEC database:
  `github.com/SergioHC95/Holographic-Entropy-Cone`.
- **psitip** — Cheuk Ting Li's Python symbolic information-theoretic
  inequality prover, source of the `M_LR` family.
- **ing39** — the Ingleton form of DFZ eq. (39), the class missing from
  psitip's list.
- **mask order** — coordinate index = bitmask of the subset minus one, with
  parties A..E = bits 0..4; the literature's lexicographic-by-size order is
  converted by the permutation used in `scripts/` and `shc_table3_map.json`.
