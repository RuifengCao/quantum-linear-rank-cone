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
- [x] **A1′ campaign closed (R15–R26).** `scripts/s4c_qlr.py` (symmetry-aware
  adjacency decomposition on the exact irredundant H-representation of QLR₅,
  seeded with the 59 known orbits plus q₂) ran in two server sessions
  (207 workers): 2026-09-11 → 1,136,210 orbits; 2026-09-12 (coordinate-focused
  queue) → **1,917,706 certified extreme-ray orbits = at least 1,368,353,902
  extreme rays** (orbit sizes summed; 98.2 % free S₆-orbits), all valid,
  extremality sampled 100 %, coordinates up to 583; the second session stopped
  at the `--max-orbits 1.9 M` ceiling (state file 39.9 MB, the largest we
  commit). **Statable results: ≥ 1.92 M orbits and ≥ 1.37 B rays (59 orbits
  previously known); full catalogue in `results/2026-09-12_a1-campaign/qlr_adj.npz`;
  the small-coordinate part `qlr5_small_orbits` (23,721 orbits with maximum
  coordinate ≤ 12; 389 with ≤ 5) and the S7 priority set `qlr5_new_small200`
  are bundled.** Growth structure: the generic region saturates quickly while
  degenerate figures (43–200 tight rows) yield ~100 new orbits each.
  **Lesson from the focused session:** expanding small-coordinate
  representatives first does not harvest small-coordinate orbits efficiently —
  only 1.1 % of the 781k new orbits had maximum coordinate ≤ 12 (38 with ≤ 5),
  the same rate as the unfocused session, at ~30× the cost per figure (the
  46 % seen in a sandbox smoke was a small-sample artefact). Further growth of
  the catalogue has diminishing scientific value; the A1′ campaign is closed
  and server time goes to A2 tools from here on. **Interpretation guard:** an
  unwitnessed orbit is not an unrealisable one.
- [x] **R32 (2026-10-02): stabilizer entropy cones are linear rank cones in disguise; they depend on the
  local dimension** (`docs/R32_reduction_and_local_dimension.md`). *Theorem 1 (proved):* for every n and
  every prime p, Stab⁽ᵖ⁾_n = {S : h_norm(S) ∈ LR⁽ᵖ⁾_{n+1}} = π(LR⁽ᵖ⁾_{n+1}) = the CSS cone (π = matroid
  connectivity function); every stabilizer entropy vector, doubled, is a CSS vector. Consequences: **S7′
  holds** (Stab₅ is cut out by the six-variable linear rank inequalities, or by the balanced ones applied to
  S); the 1,155 undecided orbits are classical representability questions; for n = 4 the DFZ five-variable
  inequalities are redundant on the self-dual slice and all 46 extreme rays of h⁻¹(Shannon + Ingleton) are
  qubit-realised (`data/s26_stab4_certs.npz`, gates G29/G29b), reproducing the four-party theorem of LMRW.
  *Theorem 2 (proved modulo DFZ arXiv:1311.4601, Thms 8.6–8.9, checked against the PDF):* from seven parties
  on, the qubit cone and the cone of any odd prime are **incomparable**. Witnesses: the "cube states" (one qudit
  per vertex of {0,1}³, code of affine functions; over GF(2) the extended Hamming code), whose entropy vectors
  differ in a single 4|4 cut (2 versus 4) and violate the pulled-back DFZ inequalities (65) resp. (91) by −2
  (`scripts/s27_cube_states.py`, `data/s27_cube_witnesses.npz`, `data/dfz_char_dependent.csv`, gate G30).
  *Proposition 3:* doubled codes turn any characteristic-dependent n-variable linear rank inequality into an
  n-party separation, so the dimension dependence of Stab_n is sandwiched between those of LR_n and LR_{n+1}.
  Five and six parties remain open. Novelty: AI pre-check only (§7 ⓪).
- [x] **A2, R31 (2026-10-02): no qubit/qutrit difference found.** The 16 orbits that R29 realised
  only with qutrits are realised by qubits at λ = 2: GF(4) graph states at λ = 1, field-reduced to qubit
  graph states with 24–32 qubits (`scripts/s24_galois_search.py`, `cc/s24_anneal_gfq.c`,
  `data/s24_qubit_l2_certs.npz`, gate G28). Conversely the 80 qubit-realised extreme rays of Stab₅ with
  the fewest qubits are realised by qutrits (73 at λ = 1, 7 at λ = 2, `data/s24_qutrit_certs.npz`); for
  four of them (q₂ among them) an exhaustive enumeration shows that neither qubits nor qutrits work at
  λ = 1 (`cc/s24_exhaustive6.c`). The R29 qubit stalls were a small-λ effect (like U₂,₄, binary only after
  doubling the ranks), not a characteristic-2 obstruction; no extreme ray of Stab₅ is known that
  distinguishes p = 2 from p = 3. Two simultaneous common informations (`scripts/s25_two_ci.py`) and
  GF(4) searches decide none of the remaining small orbits (`results/2026-10-02_r31-sandbox/`).
- [x] **A2, R30 (2026-10-01): CI + Ingleton on every undecided orbit; 641 extreme rays of Stab₅.**
  The fourth server session (`a2-ingleton`, 208 cores, 40 min,
  `results/2026-10-02_a2-ingleton/`) ran the CI+Ingleton LP at every pair on the
  2,070 undecided orbits: **904 are excluded**, each by an exact certificate
  (`data/s22_ingleton_exclusions_a2i.json`, re-derived by
  `stabcert.verify_exclusion_ext`, gate G26); their inequalities form 628 new
  classes up to S₆ (`data/s22_ineq_classes_a2i.npy`), none violated by any known
  realisable vector (all certificates, the graph states, the 1.13 M GF(2) pool
  reduced to 7,929 S₆ orbits, 46,000 random qubit/qutrit graph states with up to
  5 qubits per party) and none cutting an orbit the server found feasible.
  Qubit annealing with a new seed realised **11 more orbits** (22–28 qubits);
  qutrit annealing realised none. **Catalogue: 1,915,910 excluded (1,914,541
  single-CI + 1,369 CI+Ingleton) / 641 realised / 1,155 undecided**
  (`results/2026-10-02_a2-status/`); **641 known extreme-ray orbits of Stab₅**.
  The exclusion rate grows with the number N of qubits a λ = 1 realisation
  would need: 10 % at N = 21–30, 88 % at 51–64, all 101 orbits with N > 64.
  The 10 orbits with N ≤ 20 (frontier ray 294 among them) survive, also
  against CI + Ingleton + all 128,040 DFZ five-variable instances involving Z
  (`scripts/s23_ci_dfz.py`, floating-point verdicts; instances checked by
  gate G27) and against longer qubit/qudit searches in the sandbox; a random
  sample of 16 larger undecided orbits is not excluded by the DFZ rows either
  (`results/2026-10-02_r30-sandbox/`).
- [x] **A2, R29 (2026-10-01): S7″ refuted, qudit realisations, 630 extreme rays of Stab₅.**
  (i) **Single common information does not suffice.** 18 extreme-ray orbits of
  QLR₅ pass *every* single-CI test exactly (rational witnesses for all 24,057
  non-trivial pairs, `data/s22_single_ci_witnesses.npz`) and are nevertheless
  outside the stabilizer cone: exact certificates combine Shannon inequalities
  with conditional Ingleton instances on the CI extension A..F,Z
  (`scripts/s17_ci_test.py --extension ingleton`, `scripts/s22_certify_ext.py`,
  `data/s22_ingleton_exclusions.json`, 18 new classes `data/s22_ineq_classes.npy`;
  independent checker `epr1kit.stabcert.verify_exclusion_ext`, gate G25). So
  Stab₅ ≠ QLR₅ ∩ (single-CI cone): question S7″ is answered in the negative.
  The 18 classes cut 447 further undecided orbits; none of them is violated by
  any known realisable vector (778 certificates, 760 graph states, the 1.13 M
  GF(2) pool). (ii) **Qutrits realise what qubits missed.** 16 of the 45
  undecided orbits with at most 20 qubits — including 12 of the 17 R27 frontier
  rays (4 more are among the 18 CI+Ingleton exclusions; frontier ray 294 is the
  only one left) — are entropy vectors of explicit qutrit graph states at λ = 1
  (`scripts/s21_qudit_search.py`, `cc/s21_anneal_gfp.c`, rank formula checked
  against exact state vectors, certificates re-verified by an independent GF(p)
  elimination, gate G24). Qubit annealing had stalled at one unit of residual on
  all of them; whether some qubit multiple λ realises them is open. (iii) The
  third server session (a2-survivors: complete single-CI test on 2,566 orbits,
  none excluded; 14 new qubit realisations) and one sandbox qubit realisation
  bring the catalogue to **1,915,006 excluded (1,914,541 single-CI + 465
  CI+Ingleton) / 630 realised / 2,070 undecided**
  (`results/2026-10-01_a2-status/`), and **630 known extreme-ray orbits of
  Stab₅**. Facets: still 119 certified; the 18 new classes are valid but not
  certified facets (tight ranks 14–27). Novelty of the CI+Ingleton inequalities
  relative to DFZ's six-variable lists is unchecked.
- [x] **A2 at catalogue scale (R28, 2026-10-01).** A cutting-plane pass over the
  whole A1′ catalogue (`scripts/s18_cutting_plane.py` logic; single common
  information, disjoint pairs) classifies all 1,917,706 extreme-ray orbits of
  QLR₅: **1,914,541 (99.84 %) are outside the stabilizer cone** (each by an exact
  certificate; 747 inequality classes up to S₆, `data/s17_ineq_classes.npy`;
  the 4,665 new certificates are in `data/s18_exclusions.npz`), **599 are
  realised** by explicit graph states or are BCHS ledger rays — every realised
  extreme ray of QLR₅ is an extreme ray of Stab₅, so **599 extreme-ray orbits of
  the five-party stabilizer cone are now known (59 before)**,
  `data/stab5_extreme_reps.npy` — and **2,566 are undecided** (feasible for every
  disjoint-pair single-CI test, not yet realised;
  `results/2026-10-01_s18-cutting-plane/`). **119 of the 747 classes are
  certified facets** of the qubit stabilizer cone (`scripts/s19_facet_check.py`,
  `data/s17_ineq_tight_rank.npy`). Validation: no violation of any class on the
  1.13 M pool, the 760 six-qubit and 1,130 qutrit graph states, or any of the 747
  realisation certificates. Literature anchor: DFZ's printed six-variable
  inequalities (44)–(48), with the purifier as sixth variable, already exclude
  45,292 catalogue orbits. Write-up with full proofs:
  `docs/A2_stabilizer_vs_QLR.md` (human review pending).
- [x] **A2 resolved in the negative: the five-party stabilizer cone is strictly
  smaller than the QLR₅ cone (R27, 2026-09-30).** QLR₅ here is the cone of
  Bao–Cheng–Hernández-Cuenca–Su (SciPost Phys. 9 (2020); arXiv:2002.05317 §2.3):
  the DFZ five-variable linear rank inequalities with monotonicity replaced by
  weak monotonicity and orbits completed under S₆ — exactly our 31 facet
  classes. BCHS (arXiv:2006.16292) wrote that there is "significant evidence
  for and none against" the coincidence of the stabilizer and QLR cones.
  **Counterexample:** extreme ray #14 of the frontier set (in
  `data/s17_exclusions.json`; valid for all 10,860 facet instances, exact
  tight rank 30) violates
  `S(ABD)+2S(ACD)+S(BCD)+2S(BCE)+S(ABCE)+2S(BDE) ≥ S(AD)+2S(CD)+S(ABCD)+2S(BE)+S(ACE)+S(BCDE)+S(ABCDE)`,
  which holds for every stabilizer state of prime local dimension. Proof chain
  (details in the R27 report): (i) normal form — each party's marginal of a
  stabilizer state is maximally mixed on a stabilizer subspace, so a party-local
  Clifford discards surplus qubits and party p may be assumed to own exactly
  S_p qubits; (ii) with stabilizer ↔ Lagrangian subspace L ⊂ F_p^{2N},
  h(X) = dim π_X(L) = S(X) + |X|, a GF(p)-linear polymatroid on the six parties
  (the observation behind Gross–Walter, arXiv:1302.6990); (iii) linear
  arrangements admit common information Z = U_X ∩ U_Y; (iv) an exact rational
  Farkas certificate: a nonnegative combination of Shannon inequalities on the
  seven elements A..F,Z whose Z-terms cancel under the three CI equalities
  (`scripts/s17_certify.py`, re-checked by `epr1kit/stabcert.py`, gate G20).
  The inequality is a genuinely six-variable linear rank inequality (the
  purifier is one of the variables), which the five-variable QLR construction
  does not contain. **Scale of the gap:** single-common-information LP tests
  (`s17_ci_test.py`) exclude 139 of the 300 frontier rays (exact certificates
  for all 139; 60 inequality classes up to S₆, of which 16 are certified facets
  of the qubit stabilizer cone by 30 linearly independent tight realisable
  vectors); screening the full A1′ catalogue with these 60 classes excludes
  **1,794,496 of the 1,917,706 orbits (93.6 %)**. Validation: none of the
  38,880 inequality images is violated by any of 1,130,413 + 760 + 1,130
  (qutrit) pool vectors, 200,000 random graph states with arbitrary party sizes,
  or any realisation certificate.
- [x] **A2 positive side: targeted stabilizer search (R27).** `cc/s16_anneal.c`
  + `scripts/s16_stab_search.py` anneal over qubit graph states in the normal
  form above (for fixed λ this space loses no generality), every hit re-verified
  independently; `cc/s16_exhaust.c` decides tiny cases exhaustively. Frontier
  (all 300 orbits with maximum coordinate ≤ 5 that were neither seeds nor
  pool-witnessed, plus q₂): **realised** by explicit graph states for λ = 1
  (126), λ = 2 (17) or λ = 3 (1) — 144 in total (`data/s16_realisations.npz`,
  `data/a2_frontier_status.json`) — **excluded** by exact certificates (139),
  **undecided** 17. Realised and excluded sets are disjoint, as they must be
  (two independent methods). **q₂ is stabilizer-realisable**:
  2·q₂ is the entropy vector of a 12-qubit graph state with two qubits per
  party (also checked by state-vector SVD), while q₂ itself is not (λ = 1 is
  excluded exhaustively: with one qubit per party the realisable set is exactly
  the 760 six-vertex graph-state vectors, recomputed independently). The two
  F₃-conditional witnesses of R21 are now realised by qubit graph states at
  λ = 2, so all 37 pool witnesses are unconditional.
- [x] **A2 pipeline, node 1: the Hubeny–Rota chordality filter (R22).**
  `epr1kit/chordal.py` implements the correlation hypergraph, the chordality
  test of its line graph (Theorem 1 of arXiv:2512.24490: chordal ⟺
  realisable by a holographic simple forest ⟹ inside the stabilizer cone) and
  Algorithm 1 (simple-tree construction) with an independent min-cut
  re-computation as certificate. Sanity: of the 19 HEC₅ rays, 9 are chordal,
  6 irreducible-chordal, and all 6 constructed trees reproduce the ray exactly
  (gate G19). **Result on the catalogue: of the 15,183 small-coordinate
  orbits only the 9 HEC seeds are chordal — none of the 15,123 new orbits
  is.** So the new extreme rays of QLR₅ (including the 35 realised by
  12-qubit graph states) are not simple-forest-holographic; if they are
  stabilizer-realisable at all, non-tree models or non-holographic stabilizer
  states are needed. The shortcut therefore adds no new S7-positive verdicts;
  the next A2 tool is a targeted (non-tree) stabilizer search.
- [x] **Pool witnesses of realisable new orbits (R21, R26; superseded by R27).** With the
  50 M-sample GF(2) pool (1,130,413 realisable vectors) the witness scan
  (`s14_witness.py`) finds **37 newly discovered small-coordinate extreme-ray
  orbits realised by stabilizer states** — 35 by 12-qubit graph states
  (unconditional), 2 via the conditional F₃ layer — every hit re-verified
  coordinate by coordinate (`qlr5_new_witnessed`, gate G18). All have maximum
  coordinate ≤ 5. The ledger of known realisable extreme-ray orbits of QLR₅
  grows from 59 to **96**. (Erratum C6: the number of unwitnessed orbits with
  maximum coordinate ≤ 5 is 299, not 292 — only 53 of the 60 seeds have maximum
  coordinate ≤ 5; see CHANGELOG.)
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

**Open items, in order (details in §7):** human novelty check and contact
with the BCHS authors about the R27 result; the undecided frontier rays;
human check of R32 (proofs of Theorems 1 and 2; DFZ arXiv:1311.4601
eqs. (65)/(91) against the paper; novelty); human spot check of the 11
single-source lines in `data/dfz_ref28.csv`; optional `rays5` cross-check.
(S7′ was settled in R32.)

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
python3 selftest.py --full            # 34 gate lines + ALL PASS, ~5 min; must be all green first
python3 scripts/s1_build_qlr.py --variant pure28 --workers 25     # seconds
python3 scripts/s2_judge.py QLR_H_pure28.npy                      # expect: violated 0
python3 scripts/s3_rank19.py QLR_H_pure28.npy                     # expect: 18/19, only #19 has rank 29
python3 scripts/s12_dfz_reconcile.py                              # expect: PSITIP-INCOMPLETE / ing39

# ---- B. Server session, one command ----
sh scripts/run_b_batch.sh
#   = s7 (two layers, 0 violations) -> s8 (#19 certificate; unconditional pool first,
#     automatic escalation to the conditional layer, labelled) -> s10 (facet classes N +
#     redundant M = 84) -> s11 packaging.  Commit or push the s11 tarball afterwards.

# ---- C. Server: one session runs everything pending (see §2b) ----
python3 scripts/run_job.py run session -- --hours 8

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

## 2b. Job tiers and the hand-off protocol (sandbox vs. server)

Work is split by cost. **Tier S** jobs run in the analysis sandbox (one CPU,
minutes, no MPI): the self-test, the core verification, small analyses and
short campaign batches. **Tier X** jobs (parallel, hours, MPI or large
memory) are handed to a server — and, because server time is rented, they are
bundled into **one session** that runs everything pending end to end:

```
python3 scripts/run_job.py --list                    # registry with tiers and runtimes
python3 scripts/run_job.py run session --smoke       # sandbox rehearsal of the whole session (~2 min)
python3 scripts/run_job.py run session -- --hours 8  # the real thing, once, on the server
# second session (R23 recommendation): focus the ceiling-limited growth on the S7 front line
python3 scripts/run_job.py run session -- --hours 6 --skip s7-pools --engine-args "--order coord --queue-max-coord 20"
```

`run_session.py` runs, in order: the full self-test (aborts if not green),
the A1 campaign (budget = hours − 1.5 h, `--workers` = cpus − 1, safety valve
`--max-orbits 1.5 M` ≈ 33 MB of state), the pool job `s7-pools` (50 M GF(2)
samples + the full F₃ layer, ~1 h on 25 cores), and finally writes
`SESSION_SUMMARY.md`. Every step goes through `run_job.py`, so each one is
staged into its own `results/<date>_<tag>/` folder (also on failure or
Ctrl-C, with partial outputs), and the session continues past a failed step.
Commit `results/` once at the end.

`run_job.py` does the preflight (Python deps, gcc / lrs / mplrs as declared),
copies the latest resume file from `results/` when the job declares one,
captures the log and the exit code directly (no `tee` pitfall), writes
`JOB_STATUS.json` (host, cores, timings, exit code) and stages the declared
outputs + log + status through `s13_results_commit.py`. Oversized `.npy`
outputs are downcast losslessly (int64 → int8/int16) before staging; anything
still above 50 MB is skipped with a manifest note (attach it to a GitHub
Release instead). Engine options can be overridden after `--`.

| Job | Tier | Runtime | Smoke variant (sandbox) |
| --- | :-: | --- | --- |
| `selftest` | S | ~40 s | — |
| `verify-core` | S | ~1 min | — |
| **`session`** | X | one rented session, default 8 h | self-test + both smokes below (~2 min) |
| `a1-campaign` | X | hours, resumable, compressed `.npz` state | 2 workers, 30 s budget from the seeds |
| `s7-pools` | X | ~1 h on 25 cores | small samples (~1 min) |
| `server-batch` | X | ~25 min (only if s8/s10 must be redone) | s7 small samples + s10 first 2 classes |
| `mplrs-completeness` | X | days, MPI (optional luxury) | prepare the `.ine` only |
| `s17-ci` | S | ~7 s per ray (about 1,300 LPs) | — |
| `a2-frontier` | X | ~1 h on 25 cores (17 frontier rays, λ = 2..4, 10 min each) | 2 rays, λ = 2, 5 s |
| `a2-survivors` | X | done 2026-10-01 (28 min on 208 cores): complete single-CI test on 2,566 orbits + annealing at λ = 1, 2 | 4 orbits, 3 s annealing (~40 s) |
| `a2-ingleton` | X | done 2026-10-02 (40 min on 208 cores): complete CI+Ingleton test on 2,070 orbits (904 excluded), qutrit and second-seed qubit annealing at λ = 1, 60 s (0 + 11 realised) | 2 orbits, 3 s annealing (~3 min) |

**Rule:** every tier-X job is handed over only after its smoke variant has
passed in the sandbox (and, optionally, in the manual CI workflow
`.github/workflows/smoke.yml`). The analysis side then pulls `main`, reads
`results/`, and refreshes the catalogue assets in `data/` at milestones.

## 3. Pipeline stages

| Stage | Command essentials | Purpose / expectation |
| --- | --- | --- |
| selftest | `selftest.py [--full] [--workers N]` | 34 regression gate lines (G0–G30, G29b) + ALL PASS; run first on any machine |
| s1 | `--variant pure28` (current) / `pure`, `v3`, `pure2` (historical) | build the H-representation; seconds |
| s2 | `s2_judge.py <H.npy>` | 760 graph-state judge; `pure*` variants must give 0 violations |
| s3 | `s3_rank19.py <H.npy> [--extra-rows X]` | tight-rank test of the 19 HEC rays; pure28 → 18/19 |
| s4 | **`s4c_adjacency.py`** (`--init` seeds → `--state` batches) | CLR₅ extreme-ray enumeration engine (symmetry-aware adjacency decomposition; 162 orbits in R7); `s4b` (lrs) / `s4` (Normaliz) for seeds and cross-validation |
| s4c_qlr | **`s4c_qlr.py`** (`--H`, `--sym s6`, `--max-tight`, `--workers N`, `--order coord`, `--queue-max-coord C`) | the same engine on the QLR₅ cone (A1 campaign); parallel vertex-figure solves; optional focus on small-coordinate representatives |
| s14 | `s14_witness.py --state <npz> --gf2 <pool> --f3 <pool>` | S7 witness scan of the small-coordinate orbits against the realisable pools; exact re-verification of hits |
| s15 | `s15_chordal.py <rays.npy> [--skip-sa] [--out models.npz]` | Hubeny–Rota chordality filter + Algorithm 1 simple-tree construction with min-cut certificate |
| s16 | `s16_stab_search.py <rays.npy> --lams 1 2 --seconds 20 --out certs.npz` | targeted stabilizer search: annealing over normal-form qubit graph states (`cc/s16_anneal.c`), independent re-verification of every hit; `cc/s16_exhaust.c` for exhaustive small cases |
| s17 | `s17_ci_test.py <rays.npy> [--first-only] [--pairs disjoint] [--workers N]` → `s17_certify.py --ray-file … --index i --X x --Y y` | common-information LP test of stabilizer realisability (six-variable linear rank; `--pairs disjoint` is the fast screen, `all` the complete single-CI test); exact rational Farkas certificate for every exclusion |
| s18 | `s18_cutting_plane.py --catalogue … --classes … --realised … --workers N` | catalogue-wide cutting-plane passes: test unresolved orbits, certify, add classes, re-screen |
| s19 | `s19_facet_check.py classes.npy stab5_extreme_reps.npy [--pool …]` | exact rank of known realisable vectors on each inequality's hyperplane (30 = certified facet of the stabilizer cone) |
| s17 + Ingleton | `s17_ci_test.py <rays.npy> --extension ingleton [--workers N]` → `s22_certify_ext.py --rays … --ci … --out …` | R29: the CI extension must also satisfy every conditional Ingleton instance involving Z; exact certificates with typed terms, checked by `stabcert.verify_exclusion_ext` |
| s23 | `s23_ci_dfz.py <rays.npy> [--index …] [--out …]` | R30: CI + Ingleton + the 128,040 DFZ five-variable instances involving Z, as cutting planes, at every pair (floating-point verdicts; no certificate path yet) |
| s24 | `s24_galois_search.py <rays: .npy or .npz> --q 4 --lams 1 --seconds 60 [--out …]` ; `--self-test` ; `--exhaustive --q 3 --index …` | R31: graph states over GF(q), q = 2, 3, 4, 5, 7, 8, 9 (`cc/s24_anneal_gfq.c`); a hit over GF(p^k) is field-reduced to a GF(p) graph state at k·λ and verified by `stabcert`; `--exhaustive`: all one-qudit-per-party graph states (`cc/s24_exhaustive6.c`) |
| s25 | `s25_two_ci.py <rays.npy> [--index …] [--maxsize 2]` | R31: two simultaneous common informations on small pairs, Shannon on the 8-element extension (floating-point verdicts; exploratory) |
| s26 | `s26_reduction.py --self-test` ; `--n4 [--lrs] [--realise --out …]` | R32: Theorem 1 checks -- CSS entropy formula and the doubling 2S against state vectors; n = 4: DFZ rows redundant on the self-dual slice (LP), exact extreme rays (lrs), qubit realisations re-verified |
| s27 | `s27_cube_states.py [--state-vector] [--random N] [--out …]` ; `--probe6 SECONDS` | R32: Theorem 2 -- cube states over GF(2) and GF(p), self-duality, contraction to the DFZ Fano / non-Fano configurations, pulled-back DFZ (65)/(91) with value −2; `--probe6`: six-party checks (exploratory) |
| s20 | `s20_residuals.py <rays.npy> --lam 1 --seconds 20 --runs 3` | where qubit annealing gets stuck: residual vector of the best state of several independent runs (diagnostic; also keeps any hit as a certificate) |
| s21 | `s21_qudit_search.py <rays.npy> --p 3 --lams 1 --seconds 60` ; `--self-test` | qudit (p = 3, 5, 7) graph-state search (`cc/s21_anneal_gfp.c`), independent GF(p) re-verification; the self-test checks S(X) = rank_GF(p) W[X, X^c] against exact state vectors |
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
14. **HiGHS threads on many-core machines (R28.1).** scipy's LP solver starts
    about nproc/2 threads in every process that solves an LP, so one LP
    worker per core means ~nproc²/2 threads (about 21,000 on 208 cores). Every
    LP in the kit passes `options={'threads': 1}`; keep it that way in new
    code, because within one process an LP with a different thread setting
    after the first one fails (scipy 1.17: status 4). scipy < 1.11 ignores the
    option; `s17_ci_test.py` then caps its workers and prints a NOTE. A 2-core
    sandbox cannot show any of this: emulate the server's CPU count with an
    `LD_PRELOAD` shim that overrides `get_nprocs`/`sysconf` (as done in R28.1).
15. **A preset `OMP_NUM_THREADS`.** Some server images export it (the third
    server session found `OMP_NUM_THREADS=25` on 208 cores); an OpenMP kernel
    then silently uses only that many threads (with 1: single-threaded). `s16_stab_search.py`
    sets the thread count explicitly and logs `kernel: ... OpenMP threads N`;
    check that N is the core count.
16. **Qubit-only and single-seed searches under-report realisability (R29).**
    16 orbits that stalled at one unit of residual for qubits were realised by
    qutrits at once, and one was realised by a qubit run with a different seed.
    An unrealised orbit is not evidence against stabilizer realisability.
    R31: the qutrit-only orbits are qubit-realisable at λ = 2, and four small
    extreme rays need λ = 2 for qubits and qutrits alike. Before reading a
    λ = 1 failure as an obstruction, search GF(4) (or GF(9)) graph states at
    λ = 1 with `s24_galois_search.py`: their search space is far smaller than
    that of qubits (qutrits) at λ = 2.
17. **A hung CI run is not a red gate (R30).** Run #17 (R29) hung for over an
    hour in `apt-get update`: the runner's Ubuntu mirror stopped answering and
    apt waited without an error, so `selftest.py` never started. The workflows
    now cap every step (each `apt-get update` attempt 150 s, three attempts;
    the whole job 40 min). A run that fails or hangs before the selftest step
    says nothing about the code: cancel it and use "Re-run all jobs".
18. **Staging twice on one day (R30).** `s13_results_commit.py` (and hence
    `run_job.py`) used to write into an existing `results/<date>_<tag>/` when
    called again with the same tag on the same day, overwriting the files
    staged there (it happened once in the R30 sandbox and was undone from git
    before anything was committed). It now stages into `<date>_<tag>-2`, `-3`,
    … and prints a NOTE.
19. **Web-tool transcriptions of inequalities (R31, R32).** A web tool's
    reading of the eight-variable inequality of DFZ arXiv:1401.2507 (Thm 3.1)
    was violated by random GF(2) arrangements, so it cannot be the theorem;
    it was not used. The inequalities of R32 (DFZ arXiv:1311.4601 eqs. (65)
    and (91)) were checked term by term against the rendered PDF pages and
    tested on random arrangements before use; `data/dfz_char_dependent.csv`
    records the pages. Note that Theorem 8.6 there lists "A, B, C, D, W, X,
    Y, Z" although D does not occur in (65).

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
| `qlr5_orbits_partial.npy` | **18,985** | certified partial catalogue of QLR₅ extreme-ray orbits (int16; all rank 30, S₆-distinct, ≥ 13.1 M rays; gate G17; wide int16 sha) |
| `qlr5_new_small200.npy` | 200 | S7 priority test set: new orbits with the smallest coordinates (max coordinate 2..5) |
| `qlr5_small_orbits.npy` | **23,721** | all catalogue orbits with maximum coordinate ≤ 12 (int8; S₆-distinct; from the 2026-09-12 session) |
| `qlr5_new_witnessed.npz` | 37 | new extreme-ray orbits witnessed realisable (rep, witness vector, multiple, permutation, source pool; 35 unconditional; gate G18). `qlr5_new_witnessed35.npz` is the R21 snapshot |
| `a2_frontier300.npy` + `a2_frontier_status.json` | 300 | the A2 frontier (max coordinate ≤ 5, not seed, not pool-witnessed; q₂ = index 299) and its status: realised 156 / excluded 143 / undecided 1 (R29; R27: 144 / 139 / 17) |
| `s16_realisations.npz` | 762 | graph-state realisation certificates (rep, λ, party sizes, adjacency rows as uint64; label frontier / q2 / witness / catalogue), gate G20 |
| `s17_exclusions.json` | 139 | exclusion certificates: ray, CI pair (X, Y), exact rational multipliers y of Shannon elemental inequalities on A..F,Z, S-form inequality F′, gate G20 |
| `s17_ineq_classes.npy` + `s17_ineq_tight_rank.npy` | **747** | the distinct (up to S₆) S-form inequalities from all exclusion certificates and the rank of their tight known realisable vectors (30 = certified facet of the qubit stabilizer cone; 119 classes) |
| `s18_exclusions.npz` | 4,665 | catalogue-scale exclusion certificates (R28), compact: ray, CI pair, S-form, integer multipliers (`epr1kit.stabcert.unpack_exclusion`), gate G21 |
| `stab5_extreme_reps.npy` + `stab5_extreme_src.npy` | **641** | known extreme-ray orbits of the five-party stabilizer cone (realised extreme rays of QLR₅ = the realised catalogue orbits; source: ledger59 / cert_r27 / cert_r28 / cert_r29 (server, qubit) / cert_r29s (sandbox, qubit) / cert_r29q3 (qutrit) / cert_r30 (server, qubit; the last 11 rows)) |
| `s21_qudit_realisations.npz` | 16 | qudit certificates (catalogue index, rep, λ, party sizes, weight matrix W over GF(p) padded to 64 × 64, p = 3), gate G24 |
| `s22_ingleton_exclusions.json` | 18 | CI+Ingleton exclusion certificates (catalogue index, ray, pair (X, Y), S-form, typed terms: ["E", k] Shannon elemental, ["I", a, b, c, d, K] conditional Ingleton), gate G25 |
| `s22_ineq_classes.npy` | 18 | their distinct S-form inequalities up to S₆ (none violated by any known realisable vector; tight ranks 14–27) |
| `s22_single_ci_witnesses.npz` | 24,057 | exact rational CI extensions (denominators ≤ 3) proving that the 18 excluded orbits pass every single-CI test (`stabcert.verify_single_ci_feasible`), gate G25 |
| `s22_ingleton_exclusions_a2i.json` | 904 | R30: CI+Ingleton exclusion certificates for the orbits the server found infeasible (row of `results/2026-10-01_a2-status/undecided_reps.npy`, catalogue index, ray, pair, S-form, class, typed terms), gate G26 |
| `s22_ineq_classes_a2i.npy` | 628 | their distinct S-forms up to S₆ (lexicographically smallest image; none violated by any known realisable vector) |
| `s24_qubit_l2_certs.npz` | 16 | R31: qubit graph states at λ = 2 (G = adjacency, sizes = qubits per party) for the 16 orbits of `s21_qudit_realisations.npz`, with the GF(4) matrices they come from (Wq), gate G28 |
| `s24_qutrit_certs.npz` | 80 | R31: qutrit graph states (λ = 1 for 73, 2 for 7) for the 80 qubit-realised rows of `stab5_extreme_reps.npy` with the fewest qubits (row in `stab5_row`), gate G28 |
| `dfz_char_dependent.csv` | 2 | R32: the characteristic-dependent seven-variable linear rank inequalities (65) [odd characteristic] and (91) [characteristic 2] of DFZ arXiv:1311.4601, transcribed and checked against the PDF pages (recorded in the header) |
| `s26_stab4_certs.npz` | 46 | R32: the extreme rays of h⁻¹(Shannon + Ingleton) on five elements (n = 4; 15 coordinates) with qubit graph-state certificates (λ, sizes, adjacency W), gates G29/G29b |
| `s27_cube_witnesses.npz` | 4 | R32: entropy vectors of the qubit and odd-p cube states (127 coordinates, parties A,B,C,W,X,Y,Z = vertices 100,…,111, purifier 000) and the pulled-back seven-party inequalities F65 (odd p) and F91 (p = 2), gate G30 |
| `sixvar_template.csv` | — | transcription template for s9 |
| `manifest_shas.json` | — | canonical sha256 of every row family (`sha_rows`; `sha_rows_wide` for wide coordinates) |

Build products (`QLR_H_*.npy`) are not committed; s1 regenerates them in
seconds and checks the sha.

## 7. What next, and how

**⓪ A2 after R30, R32.** (a) Human novelty check (≥ 3 communities: quantum
information / entropy cones; information theory / linear rank inequalities;
matroid theory / representability) before any claim leaves the repository —
now including the 646 CI+Ingleton classes against DFZ's six-variable lists
and the two R32 theorems (the converse direction of Theorem 1 and the
local-dimension dependence of Theorem 2 were not found in LMRW, Gross–Walter,
BCHS or Majenz's thesis by an AI pre-check); human reading of the R32 proofs
and of DFZ arXiv:1311.4601 eqs. (65)/(91); then write to the BCHS authors. (b) The 1,155 undecided orbits
(`results/2026-10-02_a2-status/`): most need 21–40 qubits at λ = 1, where 60 s
of annealing mostly ends far from a realisation, so stronger searches
(longer, λ = 2, restarts from the best states) are the obvious server job.
(c) The 10 orbits with at most 20 qubits survive everything tried (qubits at
λ = 1, 2, qutrits, p = 5, 7, GF(4), CI+Ingleton, CI+Ingleton+DFZ, two common
informations on small pairs), among them frontier ray 294. Untried: second
common informations of pairs that involve the first, larger pairs, DFZ
instances with unions of elements in one slot, and DFZ's own
six-variable list (arXiv:0910.0284 §6: 3,490 inequalities, 2,395,095 with
permuted forms; "hundreds" need two common informations) — if that list can
be obtained (the authors' site `code.ucsd.edu/zeger/linrank`, which also hosts
`rays5`), testing h_norm of the 10 orbits against it is cheap. (d) S7′ (all
six-variable linear rank inequalities) holds (R32, Theorem 1).
The A1′ catalogue campaign is closed (R26); `s4c_qlr.py` remains available.

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
orbits have maximum coordinate ≤ 5) — superseded by R27. (c) Inequality
(4.4): validity was machine-verified in R1 with a self-produced contraction
certificate (exhaustive check of 119,877,472 tuples) and its −1 violation by
graph states reproduced; only the optional transcription of the printed
64-row table remains (erratum C7).

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
  request (all 34 gate lines must pass). Every step has a time cap (job 40 min;
  each `apt-get update` attempt 150 s, three attempts), so a stalled package
  mirror fails the run within minutes instead of hanging for GitHub's 6-hour
  default; such a failure says nothing about the code -- re-run the job.
- **`.gitignore`** excludes rebuildable products (`QLR_H_pure28.npy`), the s7
  pools, logs and campaign state files. Large result files (> 50 MB pools)
  should go to Git LFS or to release assets; small results can be committed to
  a `results/` branch.
- **Server workflow:** `git clone` → `pip install -r requirements.txt` →
  `sh scripts/run_b_batch.sh` or the A1 campaign commands → stage the
  keepers with `python3 scripts/s13_results_commit.py --tag <tag> <files>`
  (creates `results/<date>_<tag>/` with a manifest) → commit and push.
- **Round trip with GitHub Desktop:** the analysis side delivers a zip of
  changed files with relative paths; unzip over the local clone, Desktop
  lists the changes, commit with the suggested message, push. Results flow
  back through `results/` (never through chat uploads). `.gitignore` ignores
  scratch files at the repository root only; everything under `results/` is
  committed. Files above 50 MB go to GitHub Releases.
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
