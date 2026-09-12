# CHANGELOG

## R26 (2026-09-12) -- second server session analysed; A1' closed
- Session 2 (207 workers, coordinate-focused queue, 228 min): 1,136,210 →
  **1,917,706 certified orbits (≥ 1,368,353,902 rays)**, all valid,
  1,500-sample extremality 100 %; stopped at the 1.9 M ceiling (39.9 MB state).
- Focus-mode lesson: 1.1 % of the 781k new orbits have max coordinate ≤ 12
  (38 with ≤ 5) — no better than the unfocused session per orbit and ~30×
  costlier per figure; the sandbox smoke's 46 % was a small-sample artefact.
  A1' is declared closed; further server time goes to A2 tools.
- Witness scan of the 8,538 new small-coordinate orbits: +2 unconditional
  (GF(2)); merged asset `qlr5_new_witnessed` (37; 35 unconditional), small
  catalogue `qlr5_small_orbits` 23,721, refreshed `qlr5_new_small200`;
  gate G18 now reads its expectations from the manifest.

## R25 (2026-09-11) -- pre-session rehearsal on the real state
- Rehearsal on the 1,136,210-orbit state: canonicalisation throughput ~3,100
  orbits/s (load ~6 min, peak ~1.3 GB); the main process was found to be the
  bottleneck in focus mode (heavy figures, hundreds of candidates each, all
  canonicalised in one process while 207 workers wait).
- `s4c_qlr.py`: workers now canonicalise and pre-certify their candidates
  against the fork-time snapshot of the catalogue and return only unknown
  ones with keys and rank verdicts; the main process does dictionary work
  only.  Sequential and parallel runs verified to produce identical orbit
  sets (2,635 orbits over two controlled rounds, `--max-rounds`).
- `--max-rounds N` for controlled experiments.

## R24 (2026-09-11) -- hand-off fix
- `run_job.py`: pass-through arguments are shell-quoted (`shlex.quote`), so
  values with spaces such as `--engine-args="--order coord --queue-max-coord 20"`
  reach the session intact (previously the quotes were dropped and argparse
  rejected the fragments). Both `--engine-args=...` and `--engine-args "..."`
  now work; the manual uses the `=` form.
- `run_session.py --dry-run`: prints the step commands the session would run
  and exits — use it before spending server time.
- Diagnosis of the failed second-session attempt: the server was on a pre-R21
  tree (no `--engine-args`, session still staged a `session-runner` folder);
  the manual now makes `cat VERSION` a hard gate.

## R23 (2026-09-11) -- second-session tuning
- `s4c_qlr.py`: `--order coord` and `--queue-max-coord C` (expand small-coordinate
  representatives first / only: in a sandbox smoke 46 % of the new orbits had
  max coordinate <= 12 versus 1.3 % in the unfocused catalogue); queue
  construction vectorised (chunked BLAS tight counts instead of a Python loop
  over 10^6 representatives, which cost ~6 min per queue round); state saves
  throttled (`--save-every`, default 120 s; final save always forced).
- `run_session.py --engine-args "..."` forwards engine options to the campaign.
- Recommended second session: `run session -- --hours 6 --skip s7-pools
  --engine-args "--order coord --queue-max-coord 20"` (README §2b).

## R22 (2026-09-11) -- chordality filter (A2 node 1)
- `epr1kit/chordal.py` + `scripts/s15_chordal.py`: correlation hypergraph,
  line-graph chordality (MCS + perfect elimination), irreducibility,
  Hubeny-Rota Algorithm 1 (clique tree -> simple tree -> weights) and an
  independent min-cut re-computation of the entropy vector as certificate.
- Sanity (gate G19, 22 gates): HEC5 rays -> 9 chordal, 6 irreducible-chordal,
  6/6 trees verified.  Catalogue: 15,183 small-coordinate orbits -> only the 9
  HEC seeds are chordal; 0 of the 15,123 new orbits (summary in
  `results/2026-09-11_s15-chordal/`).  Implication: the new QLR5 extreme rays
  are not simple-forest-holographic; A2 needs a non-tree stabilizer search.

## R21 (2026-09-11) -- first server session analysed
- **Session results** (208-core server, 207 workers): catalogue 32,106 →
  **1,136,210 certified orbits** (≥ 810,354,812 rays), all valid, 2,000-sample
  extremality 100 %, max coordinate 583; pools: 50 M GF(2) samples →
  1,130,413 realisable vectors (0 violations), F₃ 1,130.
- **S7 evidence.** New `scripts/s14_witness.py` (hash scan + exact
  re-verification): 35 newly found small-coordinate extreme-ray orbits are
  stabilizer-realisable (33 unconditional, 2 conditional); assets
  `qlr5_new_witnessed35`, `qlr5_small_orbits` (15,183), refreshed
  `qlr5_new_small200`; gate G18 (21 gates).
- **Design flaw fixed.** The campaign queue was built once at start-up, so the
  1.1 M orbits found during the session were never expanded and the session
  idled for 4.4 of its 6.5 h. The queue is now rebuilt whenever it runs dry
  (`queue round n` lines in the log). Not the operator's fault.
- `run_job` honours `"stage": false` (the `session` wrapper no longer stages a
  duplicate folder); `--max-orbits` default 1.9 M (~40 MB compressed state,
  keeps the file under the 50 MB commit limit).

## R20 (2026-09-10) -- one-shot server session
- `scripts/run_session.py` + `session` job: self-test → A1 campaign (budget =
  hours − 1.5 h, workers = cpus − 1, `--max-orbits`) → `s7-pools` (50 M GF(2)
  samples + full F₃) → `SESSION_SUMMARY.md`; every step staged on its own;
  sandbox rehearsal `run session --smoke` (~2 min) passed end to end.
- `s4c_qlr.py`: compressed `.npz` state (×2.9 smaller; ~22 MB per 10⁶ orbits),
  vectorised canonicalisation in the merge step (the main process no longer
  throttles 24 workers on light figures), `--max-orbits` safety valve;
  legacy `.npy` states still load. `s13` summarises `.npz` states.
- New job `s7-pools` (unconditional GF(2) pool at 50 M samples + F₃ layer,
  downcast to int8 on staging) to enlarge the realisability witness pools for
  the S7 test set.
- Sandbox batch 3 (2 workers, 110 s, new merge): 18,985 → **32,106 orbits**;
  staged as `results/2026-09-10_a1-batch3/qlr_adj.npz` = the server's resume
  point. The `data/` catalogue stays at 18,985 until the session returns.

## R19 (2026-09-10) -- pre-server code review
- `s4c_qlr.py`: explicit `fork` multiprocessing context (Python >= 3.14 would
  default to forkserver and lose the inherited globals); compact state
  format 2 (int16/int32 representatives + boolean done/skipped flags, ~half
  the size, scales to 10^5+ orbits) with automatic conversion of legacy states.
- `s13_results_commit.py`: oversized `.npy` outputs are downcast losslessly
  before staging (107 MB int64 pool -> 13 MB int8) and anything still above
  the limit is skipped with a manifest note instead of aborting the staging;
  summarises format-2 states.
- `run_job.py`: `-- EXTRA` pass-through of engine options, `--no-resume`,
  Ctrl-C handling (`-interrupted` staging), partial outputs staged on failure;
  fixed an argument-parsing bug (REMAINDER swallowed `--smoke`).
- Sandbox tests: legacy state -> format 2 round trip (18,985 -> 21,215 orbits
  over two short runs), s13 downcast, smoke and pass-through runs; 20 gates.

## R18 (2026-09-10)
- **Hand-off protocol.** `jobs.json` registry (tier S sandbox / tier X server,
  each X job with a smoke variant) and `scripts/run_job.py` (preflight, resume
  from the latest `results/` folder, direct exit-code capture, `JOB_STATUS.json`,
  automatic staging into `results/<date>_<tag>/`, `-failed` staging on error).
  README §2b documents the split; manual CI workflow `smoke.yml` runs the smokes.
- **Parallel campaign engine.** `s4c_qlr.py --workers N` solves N vertex
  figures concurrently (fork + one lrs each; merge, canonicalisation and
  certification stay in the main process; budget checked between chunks of 4N).
  Sandbox smokes passed: a1-campaign (2 workers) and server-batch (s7 small
  samples 0 violations; s10 first classes).

## R17 (2026-09-10)
- **Workflow.** `.gitignore` now ignores scratch files at the repository root
  only and always commits `results/`; new `results/README.md` convention
  (`results/<date>_<tag>/` + `manifest.json`) and `scripts/s13_results_commit.py`
  to stage result files with sha256 / row-family shas and a suggested commit
  message. Round trip with GitHub Desktop documented in README §9.
- **Engine.** `core.class_reps` gains a `width` parameter (`u8` historical,
  `u16` big-endian for coordinates up to 65,535, `None` = auto); the uint8
  guard correctly stopped the campaign at a coordinate of 256.
  `s4c_qlr.py` keys its state with `u16` (legacy states are re-keyed on load)
  and appends a growth curve to `<state>.growth.csv`.
- **A1 campaign, batch 2.** 165 s, 1,615 light vertex figures →
  **18,985 orbits, all certified rank 30, S₆-distinct; orbit sizes sum to
  13,147,989 extreme rays**; max coordinate 256; 89 new orbits with maximum
  coordinate ≤ 5. Catalogue stored as int16; manifest records the full
  checks; gate G17 samples S₆-distinctness (400 rows) to stay fast.
  State and growth log committed under `results/2026-09-10_a1-batch2/`.

Round-by-round history of the kit. Errata are numbered C3–C5 (C1/C2 live in
the round reports of the internal project). Rounds are dated by their actual
production date; see erratum C5 for the rounds that were originally mislabelled.

## R16 (2026-09-10)
- Repository made English-only: README and CHANGELOG rewritten in English,
  Chinese comment lines in `data/dfz_ref28.csv` and messages in
  `scripts/run_b_batch.sh` translated; `CITATION.cff` now points at the public
  repository. No code or data semantics changed; all 20 gates unchanged.

## R15 (2026-09-10)
- **A1 campaign, day one.** New engine `scripts/s4c_qlr.py` (parameterised
  H-representation, S₆ canonicalisation, giant figures skipped, `--max-tight`
  triage); exact irredundant H-representation of QLR₅
  `data/qlr5_H_facets10860` (= the facet instance count found in R11);
  60 seeds (the 59 known orbits + q₂); 170 s / 12 vertex figures →
  **1,924 certified extreme-ray orbits** (far from saturation; coordinates up
  to 162). Conclusion on scale: full enumeration is out of reach for this
  engine (the cause of the BCHS 2021 failure); A1 restated as A1′ (certified
  partial catalogue + lower bound + structural statistics) per the plan's
  stop-loss clause.
- Assets: `qlr5_orbits_partial` (1,924; gate G17; canonical sha via the new
  int16 `sha_rows_wide`), `qlr5_new_small200` (S7 priority test set).
  `core.sha_rows_wide` added — the int8 guard of `sha_rows` correctly rejected
  coordinates of 162.
- Interpretation guard written into the README: zero realisable-pool witnesses
  among the new orbits is **not** evidence about S7 (the pools' entropy scale
  cannot reach them).
- Repository scaffolding: `.gitignore`, `.github/workflows/selftest.yml`,
  `LICENSE`, `CITATION.cff`; README §9 repository workflow.
- **Erratum C5.** The reports and packages of rounds R4–R13 were labelled
  "2026-08-13" by inertia; the actual production date was 2026-08-19 (transcript
  timestamps). The two R13 files were renamed to 08-19; the other historical
  files keep their names as a record, and this entry is authoritative.

## R13 (2026-08-19)
- **The 59-ray ledger reproduced end to end.** Under the S₆ symmetry the
  overlap between the HEC orbits and the CLR orbits is exactly #1 (→ 17 net new
  orbits); among the 26 S₆-orbits of graph-state vectors exactly 2 are new
  QLR-extreme orbits (→ +2, namely G11 and G15); 40 + 17 + 2 = 59 assembled as
  `data/qlr59_reps.npy` (S₆-distinct, all certified rank 30).
- **A3 asset.** The BCHS separating inequality (4.4) transcribed
  (`a3_hyp_ineq44`) and validated behaviourally three ways: minimum −1 over the
  full orbit of G15 (verbatim agreement with the paper), +1 for G11, and
  minimum exactly 0 over the 7,943 CLR rays; unbalanced components C/D/E each
  +2 (the paper's footnote cites E as the witness).
- **Theorem notes.** Equality of the 31 facet classes ⇒ **cone(pure28) =
  QLR₅** (earned in R12, stated here); q₂ (violates all five classical
  monotonicity inequalities, evades 431,950 realisable vectors) confirmed as a
  genuine extreme ray of QLR₅ = **a concrete test case for the S7 conjecture**
  (stabilizer = QLR).
- Gate G16 (19 gates, all green).

## R12 (2026-08-19)
- **Literature decoding.** BCHS arXiv:2006.16292 pinned verbatim: the
  59 = 40 + 17 + 2 ledger formula fully explained (40 CLR∩QLR + HEC orbits
  except "the 18th of Table 3 of [2]", one of which overlaps the 40, leaving
  17 + graph states 11 and 15); the "31 QLR inequalities" = the DFZ list with
  weak monotonicity; graph-state Table 2 corroborated.
- **31 ↔ 31 machine certification.** The S₆-classes of DFZ28 + SUB + WM number
  34 (the "34" of the R2 mystery); 3 WM classes are valid but redundant; the
  remaining 31 coincide as a set with the 31 facet classes of pure28 —
  mystery ⑥ fully resolved.
- **The "18th ray" question settled.** Machine parse of SHC Table 3 +
  S₆-orbit-level bijection: **the paper's ray 18 = our #19** — interlocking
  with the R11 exact certificate; the HEC database and Table 3 swap rows 18/19,
  which was the root of the historical worry. Data `shc_table3_rays` /
  `shc_table3_map` added, gate G15 (18 gates).

## R11 (2026-08-19)
- Post-mortem of the R10 server batch: s7 perfect (GF(2) 431,950 / F₃ 1,130
  vectors, 0 violations in both layers); s8/s10 crashed because the environment
  lacked scipy, and `tee` swallowed the exit codes so `set -e` did not stop the
  batch — the batch script now checks dependencies (scipy/gcc) up front and
  verifies each step's output.
- s8 rerun in the container and upgraded: the face containing #19 is
  two-dimensional; both edges q₁, q₂ computed exactly (rank 30 each);
  **r₁₉ = 1·q₁ + 2·q₂ exactly** = an unconditional, theorem-grade non-extremality
  certificate; q₁ == HEC ray #4; q₂ is in no realisable pool (an A2/A3-side
  observation). Asset `data/cert19_exact.npz` + gate G14 (17 gates).
- s10 rerun in the container: **31 facet classes, exactly the BCHS count**
  (53 redundant classes; 10,860 facet instances, 54 % redundancy) — the 32/34
  mystery of the R2 honesty list resolved at the count level. Asset
  `data/pure28_facet31_reps.npy`.
- **Erratum C4.** The R7 report claimed that `rank_exact_gram` had entered
  `core`; the patch had silently missed (a `str.replace` anchor that did not
  match raised no error) and the function was absent from the package. No
  existing result is affected (the gates take the exact Bareiss path). Added
  for real in this round; exposed and fixed by G14.

## R10 (2026-08-19)
- s8 accepts several pools (`--pool` takes multiple `.npy` files, merged and
  deduplicated).
- New `scripts/run_b_batch.sh`: one-command server session
  (s7 → s8 → s10 → s11 with logs and automatic packaging); s8 certificate
  policy: unconditional GF(2) pool first, the F₃ conditional layer is merged in
  only on failure, with the certificate stored separately as
  `s8_certs_conditional.npz` and the premise noted.
- Fix: the F₃ layer file is `realizable_f3_conditional.npy` (the documentation
  previously said `realizable_f3.npy`, which would crash the batch).
- README: §0 cleaned of pre-R9 leftovers (completeness closure ticked, rays5
  downgraded), §2-B replaced by the one-command batch.

## R9 (2026-08-19)
- The rays5 download was blocked (user's network cannot reach UCSD on port 80).
  Resolution by moving up a level: DFZ §4 states the counts verbatim
  (7,943 rays / 162 orbits); "uniqueness of the extreme-ray set + our exact
  certificates + the published count" closes the completeness
  cross-certification logically; the vector-by-vector diff is redundant and the
  rays5 file is downgraded to an optional reinforcement.
- Spot check: the three explicit rays printed in the paper (F³ five-subspace
  example / U₂,₄ / U₂,₅) hit our 7,943-ray set 3/3. README §7① rewritten.

## R8 (2026-08-19)
- Documentation switched to the "after the 162" era: §2 quick start rewritten
  as verification / current work / rays5 cross-check (the lrs full-run
  guidance condemned by C3 removed); pitfall 3 rewritten as the engine-history
  verdict; §7 rewritten as "what and how" (including the three-level
  escalation when s8 finds no certificate and the human spot-check list).
- New `scripts/s5b_diff_rays.py`: one-command rays5 cross-check (robust
  parser + set-level comparison + orbit counts; rehearsal prints MATCH);
  selftest gate G13 (parser), 16 gates in total.
- §1 lrs downgraded to optional; §6 data inventory gains the 162 / 7,943 rows.

## R7 (2026-08-19)
- **Result.** All 162 CLR₅ extreme-ray orbits found; orbit sizes sum to 7,943,
  matching the literature on both counts; s6 ledger: 40 / 162 QLR-extreme,
  exactly certified (non-extreme rank distribution 25×56 / 26×25 / 27×18 /
  28×23).
- **Erratum C3.** The R5 lrs tree-size estimate underestimated by a factor 68
  on a highly unbalanced tree (server: 12.7 h / 97,026,210 bases / 935 rays
  before interruption; interrupting was right).
- New engine `scripts/s4c_adjacency.py` (symmetry-aware adjacency
  decomposition, lrs as the vertex-figure sub-solver, exact integer ratio-test
  lift, state saved after every representative); data assets
  `clr5_orbit_reps162` / `clr5_rays7943`, gate G12 (15 gates).
- (Claimed here: `core.rank_exact_gram`; see erratum C4 — it did not actually
  enter the package until R11.)
- README: status, roadmap and pitfalls rewritten (pitfall 9 estimator, 10
  giant figures).

## R6 (2026-08-12)
- Documentation overhaul: status section as decided / to-do lists; new §7
  roadmap (nine items with owners and criteria); pitfall 9 (hour-scale jobs
  only on a server: interactive containers reap background processes between
  turns — observed twice in R5); gate count corrected 13 → 14.
- s4b completion markers documented (`end` + `*Totals` line); partial-stream
  pre-check recorded (125/125 genuine extreme rays covering 51/162 orbits).

## R5 (2026-08-12)
- Diagnosis of the R4 server s4 incident: Normaliz primal mode exploded in the
  intermediate hyperplane count (91.5 M at generator 62/1905, ~1.31× per
  generator); inputs verified by sha, the 1,905 rows are exactly the
  literature's facet list.
- New `s4b_run_lrs.sh` as the s4 engine (lrs reverse search; `estimate` mode
  first; measured estimate ~1.4 M bases / ~1 h / MB-scale memory — later
  found wrong, see C3); `s4_run_normaliz.sh` gains a `dual` mode and is
  downgraded to cross-validation.
- s5 auto-detects lrs V-representation output (skips `*` comments and the
  origin vertex), gcd-primitive normalisation on all paths.
- selftest 14 gates: G11 (lrs parser unit test, no lrs binary needed).
- README: quick start and pitfalls updated after the real incident; Normaliz
  version pinning documented.

## R4 (2026-08-12)
- Bundled static Normaliz 3.11.1 (`bin/`, with the GPL-3 licence and a
  provenance note); `s4_run_normaliz.sh` resolves `$NORMALIZ` → PATH → bundled
  binary and fails explicitly when none is found.
- New `data/targets1.npy` (ray #19); s8 default target changed from targets3
  to targets1.
- Assertion hardening: int8 range checks in `sha_rows` and the substitution
  builder, uint8 domain check in `class_reps`, divisibility check in exact
  Bareiss (prevents silent float truncation errors).
- selftest 8 → 13 gates: G0 (shard checkpoint/resume smoke test), G2b
  (`CLR_H_fixed` sha), G9 extended (REF28 sha, ing39 orbit consistency,
  targets1 consistency).
- s7: explicit error without gcc, downgrade warning without OpenMP,
  `--f3-limit` smoke parameter. s10: `--limit` smoke parameter, dead code
  removed. s12: docstrings rewritten as "decided / re-verify", dead function
  removed.
- Every script prints its module docstring under `--help`; s1 docstring
  escape fix; s5 value-range validation; s6 help wording for pure28; the
  parallel layer shows `--` as ETA for the first two seconds.
- README rewritten (status / quick start / stage table / pitfalls / data
  inventory); `requirements.txt`, `VERSION`, `CHANGELOG.md` added.

## R3 (2026-08-12)
- s12 verdict: PSITIP-INCOMPLETE, missing Ingleton(39); `CLR_H_fixed`,
  `templates28`, pure28 (18/19) bundled; `dfz_ref28.csv` transcribed
  (two sources cross-checked, four-fold validation).
- Progress output (pmap / shards / C tools); s12 verdict wording corrected.

## R2 / R1
- See the corresponding round reports of the internal project.
