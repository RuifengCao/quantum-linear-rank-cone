# CHANGELOG

## R35 (2026-10-02) -- how far down the dimension dependence goes
- **Theorem 7 (proved).** Stab^(Z_d)_4 = h^-1(Shannon + Ingleton) for every d >= 2, so for n <= 4 the stabilizer
  cone does not depend on the local dimension.
  - Each of the 46 extreme rays (lrs, R32) has an integer graph state, the qubit certificate with signs on its
    edges, whose visible cut blocks have Smith invariants 0 and 1 only. It therefore realises the ray over every
    GF(p).
  - Signs were needed for rays 33, 36, 39, 41 and 45.
  - The lifts were re-verified with `stabcert` for p <= 13.
  - The inclusion in the other direction holds because Ingleton holds for subgroup rank functions.
- **Proposition 6 (proved).** The one-cut observation of R33 holds in general. Over a field of characteristic l
  that does not divide t - 1, the generalised cube matroid differs from the rational one exactly on the balanced
  sets X with:
  - c in X, 0 not in X, l | t - a, a != t; or
  - 0 in X, c not in X, l | a - 1, a != 1.
  - The drop is exactly 1. For l = t prime this is the single cut {A, C} | {B, 0}.
  - Short proof via the dependency equations; the exact rank functions agree for t = 2..7 and l <= 13.
- **Proposition 8 (proved, exact certificates).** Consider the six seven-variable characteristic-dependent
  inequalities DFZ (65), (91), (LA)_2, (LB)_2, Pena-Sarria (a), (b) with t = 2. For each, f + f o D is a
  nonnegative rational combination of Shannon elemental inequalities, where D is Kaced's duality
  (arXiv:1611.04109).
  - So, in every relabelling, they follow from strong subadditivity on the self-dual slice, hold for every pure state
    on 6 + 1 parties, and cannot separate local dimensions at n = 6.
  - A six-party separation needs a characteristic-dependent g with g o D = g.
  - Not covered: BKL's seven-variable inequalities (transcription unchecked) and substitution instances.
  - The exhaustive point-configuration probe (f(r + r*) >= 0 on all {0,1}-point configurations of GF(p)^3) is now an
    encoding check.
- **Five parties (pilot, not certified).** 40 of the first 60 known extreme rays of Stab_5 lift to all fields;
  the remaining 20 need a better search. No lifts are stored and no gate covers it (`--n5-pilot`).
- **Referee.** An independent referee re-derived Theorem 7, Proposition 6 and Proposition 8 with its own code:
  - its own double description and lrs gave the same 46 rays;
  - the one-cut prediction matched in all 30 cases, and the hypothesis that l does not divide t - 1 is shown to be
    needed;
  - the certificates were re-found;
  - 1,750 pure states x 5,040 relabellings x 6 forms gave minimum 0.

  Applied fixes: scope ("considered here", relabellings, BKL and substitutions not covered), pairwise statement of
  the dependence, the face argument for n < 4, the pilot marked as not certified, stale "open" remarks removed, and
  Kaced cited for D.
- **Code and data.**
  - `scripts/s30_dimension_profile.py` (`--n4`, `--n4-search`, `--one-cut`, `--six`, `--solve`, `--probe6`,
    `--n5-pilot`).
  - `data/s30_stab4_unimodular.npz` and `data/s30_selfdual_shannon_certs.npz`, both in the manifest.
- **Gates.**
  - G33: four-party lifts 46/46, six-party certificates 6/6, one cut for t = 2..5, manifest.
  - G33b (`--full`): one cut for t = 6, 7.
  - 40 gate lines.
- **Docs.** R32 note: section 4.7, with a header note, the one-cut paragraph of 4.4 and the open questions updated.
  README: R35 status, pipeline row s30, data inventory, gate counts. Raw output: `results/2026-10-02_r35-sandbox/`.

## R34 (2026-10-02) -- Z_d qudits: the stabilizer entropy cone remembers the local dimension
- **Theorem 5 (proved).** For d, d' >= 2 and n >= 2 max(d, d') + 3, Stab^(Z_d)_n is contained in Stab^(Z_d')_n
  exactly when d divides d' (entropies in bits). For d | d' the inclusion holds for every n.
  - Details in `docs/R32_reduction_and_local_dimension.md`, section 4.6.
  - Examples: qubits are strictly inside Z_4 from n = 7, and Z_4 strictly inside Z_8 from n = 11.
  - Galois qudits GF(4) give exactly the qubit cone, so the two standard theories of local dimension 4 differ.
- **Theorem 1(d) (proved).** Stab^(Z_d)_n = h^-1(Gamma^(d)_{n+1}) = pi(Gamma^(d)_{n+1}) = closed CSS cone over Z_d.
  Gamma^(d) is the cone of subgroup rank functions log|sum U_P| in abelian groups of exponent dividing d.
  - Entropy formula: Gross-Walter, Thm 2, valid for every d.
  - h_S is a subgroup rank function via the quotient Gamma of the Lagrangian L and Pontryagin duality.
- **Lemma 4(G) (proved).** Lemma 4 holds for subgroups of finite abelian groups:
  - (LB) when multiplication by t is injective;
  - (LA) when t V = 0.
  - (LB) gets a counting proof (no complements), which is also shorter for vector spaces.
- **Witnesses.** The generalised cube states over Z_N, N = p^{j+1}: free code, unit weights, C^perp = C diag(d),
  h_S = 2r, contraction = J - I over Z_N, and H(B) from the Smith form diag(1, ..., 1, t). Certified values:
  - Z_4, t = 2 (8 qudits): (LA)_2 = (LB)_2 = -2 bits, so it is outside every prime cone. Its entropy vector in
    bits equals the qubit cube vector plus the qutrit cube vector (inside their convex hull).
  - Z_8, t = 4: (LA)_4 = -2 bits.
  - Z_9, t = 3: (LA)_3 = -2 trits.
  - Z_4, t = 6: (LA)_6 = -2 bits.
  - The Z_4 cube state also matches its state vector (4^8 amplitudes, deviation 9e-16).
- **Remarks.**
  - For n = 4, Stab^(Z_d)_4 lies in h^-1(Shannon + Ingleton), with equality for even d, since subgroups have common
    information.
  - The A2 exclusion certificates use only Shannon, Ingleton and one common information, so they hold for every
    local dimension. In particular Stab_5 is strictly inside QLR_5 also for the cone over all d.
- **Open question added.** Can Z_d stabilizer states leave the convex hull of the prime cones? Up to one party this
  is the open inclusion Gamma^MixL subset Gamma^Abl of Khazaei, arXiv:2608.09543.
- **Referee.** An independent referee re-derived every item and stress-tested Lemma 4(G) (about 1.7 M evaluations
  on groups up to Z_12^6) with no violation of a claimed form. Applied fixes:
  - cite Gross-Walter for the entropy formula;
  - state the doubling explicitly;
  - the witnesses need t - 1 to be a unit;
  - the prime-power chain needs j >= 1;
  - GF(4) = qubits in both directions.
- **Code and data.** `scripts/s29_zd_stabilizer.py`: exact Z_N arithmetic through Smith forms, witnesses, state
  vector, Lemma 4(G) point tests. `data/s29_zd_witnesses.npz` with per-case int16 row sha in the manifest.
- **Gates.**
  - G32: Z_4 (t = 2), Z_8 (t = 4) and Z_9 (t = 3) witnesses, the partner value (LA)_4 = 0 on Z_4, the state vector,
    the cube-sum identity, stored data.
  - G32b (`--full`): Z_4 with t = 6, (LA)_6 >= 0 on GF(2) and GF(3), and Lemma 4(G) point tests. The claimed forms
    stay >= 0; (LA)_2 and (LB)_2 over Z_4 and (LA)_3 / (LA)_2 on J - I over Z_4 / Z_9 go negative.
  - 38 gate lines.
- **Docs.** R32 note: section 4.6, open questions, literature status (AI pre-check: no prior statement found that
  the cone depends on d), reproduce. README: R34 status, pipeline row s29, data inventory, gate counts, glossary.
  Raw output: `results/2026-10-02_r34-sandbox/`.

## R33 (2026-10-02) -- any two distinct primes give incomparable stabilizer cones
- **Theorem 2' (proved, self-contained).** For primes p < q and every n >= 2p + 3, Stab^(p)_n and Stab^(q)_n are
  incomparable. Details in `docs/R32_reduction_and_local_dimension.md`, section 4.4.
  - Witnesses: generalised cube states on the 2t+4 points e_i, c - e_i, c, 0 of F^{t+1} (t = p; the cube for
    t = 2), with the code of affine functions over GF(p) and over GF(q).
  - Self-duality by explicit weights (-1 on e_i, +1 on c - e_i, -(t-1) on c, t-1 on 0; G diag(d) G^T = 0 over Z).
    This holds whenever char does not divide t - 1, so h_norm = 2r.
  - Contracting the origin gives the configuration A_i = <e_i>, B_i = <c - e_i>, C = <c> (the matrix J - I).
- **Lemma 4 (proved in the note).** Two characteristic-dependent linear rank inequalities for that configuration,
  with explicit error terms:
  - (LB) valid when char does not divide t: M H(C) <= H(B) + M[gamma + M alpha + sum(delta_i + eps_i + zeta_i)].
  - (LA) valid when char divides t: H(B) <= (M-1) H(C) + M alpha + M gamma + sum(zeta_i + eta_i).
  - On the configuration each fails by 1 in the wrong characteristic; pulled back to the witnesses, the values are
    -2.
  - An independent referee re-derived every step and stress-tested both forms (billions of configurations, several
    families exhaustive up to GL, adversarial searches) with no violation; small textual fixes were applied.
  - Pena-Sarria arXiv:1905.00003v3 Example 6 (typed from the rendered page 4) gives the same values and is kept as
    a cross-check only.
- For t = 2, 3, 5 the two entropy vectors differ in exactly one coordinate: the cut {A_1..A_M, C} | {B_1..B_M, 0},
  with value t over GF(p) and t + 2 over GF(q).
- **Misprint found (pitfall 20).** DFZ arXiv:1401.2507v1 Theorem 3.1 (eight variables, characteristic != 3) is false
  as printed: C = X = a line and the rest zero give -4, over every field. An exhaustive check over GF(2)^2 / GF(3)^2
  finds 125 / 254 violations (`s28_generalized_cubes.py --dfz1401`). Not used.
- **Code and data.** `scripts/s28_generalized_cubes.py` (`lemma_forms`, `ps_forms`, witnesses, `--exhaustive01`,
  `--dfz1401`); `data/s28_gencube_t3.npz` (GF 3/5/7) and `data/s28_gencube_t5.npz` (GF 5/7) with manifest entries
  (int16 row sha).
- **Gates.**
  - G31: t = 3 and t = 5 witnesses (weights, self-duality, contraction, Lemma 4 values -2, Pena-Sarria agreement,
    single differing coordinate, stored data).
  - G31b (`--full`): (LA) over GF(2) and (LB) over GF(3) on all 2.1 M {0,1}-point configurations for t = 2, with
    the wrong-characteristic forms violated there (sensitivity).
  - 36 gate lines.
- **Docs.** R32 note: sections 4.4 (Lemma 4 with proof, Theorem 2') and 4.5, open questions, literature status,
  reproduce. README: R33 status, pipeline row s28, data inventory, gate counts, pitfall 20. Raw output:
  `results/2026-10-02_r33-sandbox/`.

## R32 (2026-10-02) -- reduction to linear rank cones; local-dimension dependence from seven parties
- **Theorem 1 (proved, `docs/R32_reduction_and_local_dimension.md`).** For every n and every prime p,
  Stab^(p)_n = {S : h_norm(S) in LR^(p)_{n+1}} = pi(LR^(p)_{n+1}) = CSS cone, where pi(r)(X) = r(X) + r(E\X) - r(E)
  is the connectivity function. The proof combines Lemma 2 of the A2 note with the CSS entropy formula
  S = pi(r) and the identity pi(h_norm(S)) = 2S. Every stabilizer entropy vector, doubled, is a CSS vector
  (the CSS state of the Lagrangian code). By balancing, Stab_n is cut out by the balanced (n+1)-variable linear
  rank inequalities applied to S. **S7' holds**; the 1,155 undecided orbits are classical representability
  questions.
- **n = 4 check.** The 24 DFZ five-variable inequalities pulled back through h_norm give 820 rows, all implied
  by the pulled-back Shannon + Ingleton rows (LP). The cone h^-1(Shannon + Ingleton) on five elements has
  exactly 46 extreme rays (lrs), all realised by qubit graph states (41 at lambda = 1, 5 at lambda = 2;
  `data/s26_stab4_certs.npz`, `stabcert.verify_realisation_n`). So Stab4 = h^-1(Shannon + Ingleton), consistent
  with LMRW.
- **Theorem 2 (proved, using DFZ arXiv:1311.4601 Thms 8.6-8.9).** For every n >= 7 and odd prime p, the qubit
  and p-qudit stabilizer cones are incomparable.
  - Witnesses: the cube states. One qudit per vertex of {0,1}^3 carries the code of affine functions (over
    GF(2) the extended Hamming code); the seven nonzero vertices are A,B,C,W,X,Y,Z in DFZ's order and the
    origin is the purifier.
  - Exact facts: the code's matroid is identically self-dual, so h_norm = 2r. Contracting the origin gives
    DFZ's Fano (GF(2)) resp. non-Fano (GF(p)) configuration. All 4x4 minors lie in {0, +-1, +-2}.
  - The pulled-back eq. (65) (valid for odd characteristic) is -2 on the qubit state; the pulled-back eq. (91)
    (valid for characteristic 2) is -2 on the odd-p state.
  - The two entropy vectors differ only on the cut ABCZ | WXY+purifier (2 versus 4).
  - (65) and (91) were transcribed into `data/dfz_char_dependent.csv` and checked against the PDF pages.
- **Proposition 3 (doubled codes, proved).** For any GF(p) configuration G, the CSS state of [G | G] (copies in
  the purifier) has h_norm = r_G + |X| on the visible parties. Hence (sandwich): LR^(p)_n not in LR^(q)_n implies
  Stab^(p)_n not in Stab^(q)_n, which implies LR^(p)_{n+1} not in LR^(q)_{n+1}. At n = 7 this gives a second proof
  of Theorem 2 (doubled Fano / non-Fano: balanced (65)/(91) = -1). Checked for the doubled Fano, non-Fano and T8
  codes in s27 and G30. Possible routes to "all pairs of primes" (DFZ arXiv:1401.2507, Pena-Sarria
  arXiv:1905.00003) are recorded as open; neither paper has been checked.
- **Code and data.**
  - `epr1kit.stabcert.verify_realisation_n`: n-party GF(p) graph-state certificates, p = 2 allowed.
  - `scripts/s26_reduction.py`: self-test against state vectors, `--n4 [--lrs] [--realise]`.
  - `scripts/s27_cube_states.py`: witnesses, `--state-vector`, `--random N`, `--probe6`.
  - `data/s27_cube_witnesses.npz` and manifest entries.
- **Gates.**
  - G29: CSS formula and doubling vs state vectors; n = 4 LP; the 46 certificates.
  - G30: the cube-state facts, both -2 values, state vectors, stored data, random-arrangement sanity.
  - `--full` adds G29b (exact lrs enumeration = the stored 46 rays).
  - 34 gate lines; quick run ~85 s, `--full` ~4.5 min in the sandbox.
- **Six-party probe (negative, exploratory).** The symmetrised Fano / non-Fano polymatroids (h_norm of the
  seven-qudit CSS states) never violate (65) / (91) under any of the 5,040 relabellings (minima 19 / 16).
  240 s of annealing over seven-party qutrit graph states also found no violation (best 4).
- **Docs.**
  - New `docs/R32_reduction_and_local_dimension.md`.
  - A2 note: header pointer; Sections 6, 8 and 9 updated (S7' answered; seven-variable characteristic-dependent
    inequalities exist; qubit/qudit status).
  - README status, open items, pipeline rows s26/s27, data inventory, gate counts, pitfall 19.
  - Raw output: `results/2026-10-02_r32-sandbox/`.

## R31 (2026-10-02) -- qubits and qutrits: no difference found; Galois-field search; two common informations
- **Qutrit-only orbits are qubit-realisable.** `scripts/s24_galois_search.py` + `cc/s24_anneal_gfq.c`
  search weighted graph states over GF(q), q = 2, 3, 4, 5, 7, 8, 9 (field tables), and turn a hit over
  GF(p^k) into a GF(p) certificate at k * lambda by standard field reduction (x-coordinates in the basis
  1, a, ..., a^(k-1), z-coordinates in the trace-dual basis), local Fourier transforms to graph form, and
  verification by `epr1kit.stabcert`; `--self-test` compares reduced ranks with k times the GF(q) ranks on
  random states (q = 4, 8, 9). All 16 orbits realised only by qutrits in R29 are GF(4) graph states at
  lambda = 1, hence qubit graph states at lambda = 2 with 24-32 qubits (`data/s24_qubit_l2_certs.npz`).
- **Reverse direction.** The 80 qubit-realised extreme rays of Stab5 with the fewest qubits at
  lambda = 1 (2-12) are realised by qutrits: 73 at lambda = 1, 7 at lambda = 2
  (`data/s24_qutrit_certs.npz`). For four of the seven (stab5 rows 43, 44, 195 = q2, 225; one qudit per
  party at lambda = 1) `cc/s24_exhaustive6.c` (`s24_galois_search.py --exhaustive`) enumerates all 2^15
  qubit and all 3^15 qutrit graph states: none works at lambda = 1, while both fields work at lambda = 2.
  Positive controls (rows 3 and 45) give hits. So the R29 qubit stalls were a small-lambda effect (as for
  U_{2,4}, binary only after doubling the ranks), not a characteristic-2 obstruction; no extreme ray of
  Stab5 is known that separates p = 2 from p = 3. The sharper route to "the five-party stabilizer cone
  depends on the local dimension" considered after R30 has no candidate left.
- **Gate G28**: the 16 qubit certificates (lambda = 2, `stabcert.verify_realisation`, three of them
  re-derived from the stored GF(4) matrices), the 80 qutrit certificates (`verify_realisation_gfp`, the
  selection rule recomputed), the field-reduction self-test. `--full`: 31 gate lines.
- **Two common informations.** `scripts/s25_two_ci.py`: Z1 and Z2 for every unordered pair of disjoint
  pairs with |X|, |Y| <= 2, Shannon on the 8-element extension (covers the structure of DFZ's
  two-common-information inequality (61), arXiv:0910.0284 Section 6). The 10 smallest undecided orbits
  and 6 larger ones are feasible for all of them (about 90,000 LPs). GF(4) at lambda = 1 realises none of
  the 10 (290 s each) nor 16 larger undecided orbits (55 s; residual 16-74). Raw output:
  `results/2026-10-02_r31-sandbox/`.
- Docs: `docs/A2_stabilizer_vs_QLR.md` new section 8 (R31); the R29 sentence "a qubit realisation at
  larger lambda is not excluded" now points to it; open questions updated. README status, pipeline rows
  s24/s25, data inventory, pitfall 16.

## R30 (2026-10-01) -- CI + Ingleton on every undecided orbit; 641 extreme rays of Stab5
- **Fourth server session** (`a2-ingleton`, 208 cores, 40 min, `results/2026-10-02_a2-ingleton/`,
  server date 2026-10-02): the CI+Ingleton LP at every pair on the 2,070 undecided orbits finds
  904 infeasible, 1,166 feasible, 0 solver errors (~17 min); qutrit annealing (p = 3, lambda = 1,
  60 s) realises none of 1,969 orbits; qubit annealing with a new seed realises 11 (22-28 qubits).
  `OMP_NUM_THREADS=25` was preset again and replaced by 208.
- **904 exact certificates** (`data/s22_ingleton_exclusions_a2i.json`, 94 s with
  `scripts/s22_certify_ext.py` in the sandbox), re-derived by `stabcert.verify_exclusion_ext`; every
  one uses conditional Ingleton rows (3-29 of 17-52 terms). 628 new classes up to S6
  (`data/s22_ineq_classes_a2i.npy`), disjoint from the 18 of R29. One-off validation: none of their
  444,600 images is violated by any certificate, graph state, the F3 vectors, the 1,130,413 GF(2)
  pool vectors (7,929 S6 orbits) or 46,000 random weighted graph states (qubits with 1-5 per party,
  qutrits with 1-3); none cuts an orbit the server found feasible.
- **Catalogue** (`results/2026-10-02_a2-status/`): 1,915,910 excluded (1,914,541 single-CI + 1,369
  CI+Ingleton) / 641 realised / 1,155 undecided. `data/stab5_extreme_reps.npy` grows to 641 rows
  (the last 11 are the server realisations, source `cert_r30`); its rows are exactly the realised
  catalogue orbits. Exclusion rate by qubit count N at lambda = 1: 0 % (N <= 20, 10 orbits),
  10 % (21-30), 33 % (31-40), 70 % (41-50), 88 % (51-64), 100 % (> 64, 101 orbits).
- **`scripts/s23_ci_dfz.py`**: CI + Ingleton + the DFZ five-variable linear rank inequalities
  (`data/dfz_ref28.csv`: eqs. (1)-(24) and Ingleton forms (36)-(39) of arXiv:0910.0284) on the
  extension, five distinct elements contracted by any subset of the other two: 149,520 distinct
  instances, 128,040 involving Z, used as cutting planes. All 10 orbits with N <= 20 stay feasible at
  every pair (at most 166 DFZ rows activated per orbit). A random sample of 16 larger undecided
  orbits is feasible as well (0 of 26 excluded): the DFZ rows add nothing measurable to CI + Ingleton
  here. Longer sandbox searches on the 10 small orbits (qubits lambda = 2, 150 s; qutrits lambda = 2,
  120 s; p = 7 lambda = 1 and p = 5 lambda = 2, 60 s) realise none. Raw output:
  `results/2026-10-02_r30-sandbox/`.
- **Gates**: G26 (904 certificates re-derived, classes, no realised vector or server-feasible ray cut,
  11 realisations re-verified, status and accounting, 641 extreme rays = realised orbits); G27 (the 28
  DFZ forms on random GF(2) arrangements, instance counts, the 21,480 instances on A..F hold for all
  realised vectors). G25 now checks the first 630 rows of `stab5_extreme_reps.npy` against the R29
  sha (`stab5_extreme_reps_r29` in `data/manifest_shas.json`). `--full`: 30 gate lines.
- **CI hardening.** Run #17 (R29) hung in `apt-get update` for over an hour: the runner's Ubuntu
  mirror (`azure.archive.ubuntu.com`) stopped answering and apt waited without an error, so
  `selftest.py` never started; the re-run passed in 2 min. Both workflows now cap the job (40 min for
  `selftest`, 60 min for the manual smokes) and every step; each `apt-get update` attempt is capped
  at 150 s and retried up to three times, with apt's own retry and timeout options set (checked
  against a fake apt-get that hangs on its first two calls: both attempts are killed, the third
  succeeds, no process is left over). `actions/checkout@v5` and `actions/setup-python@v6` (Node 24;
  run #18 warned that Node 20 actions are deprecated).
- **`s13_results_commit.py` never overwrites a staged result**: a second staging with the same tag
  on the same day goes to `<date>_<tag>-2`, `-3`, ... (README pitfall 18; it overwrote the R29
  `results/2026-10-01_a2-status/` in the sandbox once and was undone from git before any commit).

## R29 (2026-10-01) -- S7'' refuted; qutrit realisations; 630 extreme rays of Stab5
- **Third server session** (a2-survivors, 208 cores, 28 min, `results/2026-10-01_a2-survivors/`):
  the complete single-CI test excludes none of the 2,566 undecided orbits (all pairs feasible,
  no solver error); annealing realises 14 orbits at lambda = 1 (re-verified, gate G23). The
  server had `OMP_NUM_THREADS=25` preset; the R28.1 override used all 208 threads. The CI
  phase took ~11 min instead of the estimated 4 (~50 core-seconds per orbit on that machine).
- **Qudits.** `scripts/s21_qudit_search.py` + `cc/s21_anneal_gfp.c`: annealing over weighted
  graph states over GF(p), p = 3, 5, 7, in the same normal form as s16; `--self-test` checks
  S(X) = rank_GF(p) W[X, X^c] against exact state vectors; every hit is re-verified by an
  independent GF(p) elimination (`epr1kit.stabcert.verify_realisation_gfp`). Of the 45
  undecided orbits with N <= 20 -- all stalled at residual 1 for qubits -- 16 are realised by
  qutrit graph states at lambda = 1 (`data/s21_qudit_realisations.npz`, gate G24), including
  12 of the 17 R27 frontier rays; p = 5 realises none of the remaining 28 (30 s each).
- `scripts/s20_residuals.py` (+ `dump_best` option of `cc/s16_anneal.c`): residual vectors
  of the best qubit states of independent runs; a different seed realised one more orbit
  (N = 20, label catalogue in `data/s16_realisations.npz`); for most orbits the cut that stays
  off by one changes from run to run.
- **CI + Ingleton.** `s17_ci_test.py --extension ingleton` adds the 1,320 conditional Ingleton
  instances that involve Z to the 7-element extension; `s17_certify.certify(...,
  extension='ingleton')`; `scripts/s22_certify_ext.py`; independent checker
  `epr1kit.stabcert.ingleton_row` / `verify_exclusion_ext` (typed terms, rows rebuilt from
  their descriptors). 18 of the 28 remaining small orbits are excluded with exact certificates
  (`data/s22_ingleton_exclusions.json`), leaving 10 orbits with N <= 20; the R27 frontier is now
  156 realised / 143 excluded / 1 undecided (ray 294; `data/a2_frontier_status.json`). The 18
  classes (`data/s22_ineq_classes.npy`) cut 447
  more undecided orbits. Validation: no violation by the 778 realisation certificates, the 760
  graph states, the 1.13 M GF(2) pool or the F3 pool.
- **S7'' refuted.** The 18 orbits pass every single-CI test exactly: rational CI extensions
  (denominators <= 3) for all 24,057 non-trivial pairs, `data/s22_single_ci_witnesses.npz`,
  checked in integer arithmetic by `stabcert.verify_single_ci_feasible` (gate G25). Hence
  Stab5 is strictly smaller than QLR5 cut by all single-CI inequalities. Novelty of the 18
  classes relative to DFZ's six-variable inequalities is unchecked.
- **Catalogue** (`results/2026-10-01_a2-status/`): 1,915,006 excluded (1,914,541 single-CI,
  465 CI+Ingleton = status code 3) / 630 realised / 2,070 undecided; 630 known extreme-ray
  orbits of Stab5. Facets: still 119 (the 31 new extreme rays raise no tight rank to 30); the 18
  new classes are valid, tight ranks 14-27, not certified facets.
- Gates G23 (server merge), G24 (qudits), G25 (CI+Ingleton, S7'' witnesses, accounting):
  `--full` prints 28 gate lines. New tier-X job `a2-ingleton` for the fourth session.
- `s17_ci_test.py`: `--chunk` (rays per pool task) and an extension-dependent default
  `--result-timeout` (Ingleton LPs are ~4x slower: 1 ray per task, 30 min instead of 4 rays, 15 min).

## R28.1 (2026-10-01) -- pre-server hardening of the a2-survivors job
- HiGHS (scipy's LP solver) starts about nproc/2 threads in every process that
  solves an LP. The single-CI test runs one LP worker per core, i.e. about
  21,000 threads on the 208-core server (measured with the CPU count faked to
  208 by an LD_PRELOAD shim; invisible on the 2-core sandbox). Every LP in the
  kit now passes `options={'threads': 1}` (`s17_ci_test`, `s17_certify`,
  `core.decompose_certificate`, `s10`, `s12`): 3 threads per worker instead of
  106, same speed. scipy 1.11.4, 1.14.1 and 1.17.1 honour the option; 1.9.3
  and 1.10.1 ignore it, and `s17_ci_test` then caps the worker count (NOTE line).
- Mixed settings are a trap of their own: with scipy 1.17 an LP pinned to one
  thread after an unpinned LP in the same process fails with status 4, which
  the old `ci_feasible` would have read as "infeasible". Pinning only the CI
  test would have created exactly this situation in `s18_cutting_plane` on a
  many-core machine (certificate LPs in the parent, CI LPs in forked workers
  that inherit the parent's HiGHS state); hence the pin is applied to every LP,
  and the s18 pattern was re-run with the CPU count faked to 208 (statuses 0/2 only).
- `s17_ci_test`: HiGHS statuses other than 0/2 are reported as solver errors
  (ray UNDETERMINED), never as infeasibility; abort when 8 of the first 16 rays
  error (systemic); abort instead of waiting forever when no result arrives for
  15 min (a dead pool worker); final `summary:` line; `lp_errors` in the JSON.
- `s16_stab_search`: sets `OMP_NUM_THREADS` for the annealing kernel explicitly
  (`--threads`, else every CPU in the affinity mask) and prints a NOTE when it
  replaces a value from the environment -- some server images export
  `OMP_NUM_THREADS=1`, which would have serialised the annealing (up to ~53 h instead
  of ~17 min); prints the gcc flags and thread count, and a WARNING when the
  kernel had to be compiled without OpenMP (previously a silent fallback).
- `jobs.json` a2-survivors: deletes stale outputs first (files from the smoke
  run could otherwise be staged as partial outputs of a failed full run);
  declares gcc for the preflight; runtime estimate from the measured qubit
  counts (2,323 orbits with N <= 64 at lambda=1, 860 with 2N <= 64 at lambda=2).
- Gate G22 (single-CI LP statuses with the pinned solver): `--full` now prints 25 gate
  lines plus the ALL PASS line. Correction: the "25 gate lines" quoted for R28 counted
  the ALL PASS line; R28 had 24 gates.
- CI runners pinned to ubuntu-24.04 (ubuntu-latest moves to 26.04 on 2026-10-19).

## R28 (2026-10-01) -- A2 at catalogue scale
- Cutting plane over the whole A1' catalogue (1,917,706 orbits): 1,914,541
  excluded from the stabilizer cone by exact certificates (747 inequality
  classes up to S6; 4,665 new certificates in `data/s18_exclusions.npz`, all
  re-derived exactly by `epr1kit.stabcert`), 599 realised (747 graph-state
  certificates in `data/s16_realisations.npz`; every realised extreme ray of QLR5
  is an extreme ray of Stab5 -> `data/stab5_extreme_reps.npy`), 2,566 undecided
  (`results/2026-10-01_s18-cutting-plane/`).
- 119 of the 747 classes certified as facets of the qubit stabilizer cone
  (`scripts/s19_facet_check.py`).
- `s17_ci_test.py`: `--pairs disjoint` fast screen (every infeasible pair found
  so far is disjoint) and `--workers`; new `scripts/s18_cutting_plane.py`,
  `scripts/s19_facet_check.py`; compact certificate storage in `epr1kit.stabcert`.
- Literature anchor: DFZ's printed six-variable inequalities (44)-(48) (purifier
  as sixth variable) exclude 45,292 catalogue orbits.
- Headline counterexample changed to the smallest one: frontier ray 275 (10 qubit
  units, maximum coordinate 4), cut by a certified facet (class 29) with a
  26-term certificate; full write-up `docs/A2_stabilizer_vs_QLR.md`.
- Gate G21 (25 gate lines in `--full`); job `a2-survivors` (X).
- Incident: `pkill -f` matched its own shell and killed the controlling command
  (pitfall 3 of R2 repeated); use `kill <pid>` from `ps` instead.

## R27 (2026-09-30) -- A2: the five-party stabilizer cone is strictly smaller than QLR5
- **Main result.** Single-common-information LP tests (`scripts/s17_ci_test.py`)
  on the 300 frontier rays find 139 rays of QLR5 that are not in the stabilizer
  cone; each exclusion has an exact rational certificate (`scripts/s17_certify.py`,
  `data/s17_exclusions.json`) re-derived by an independent checker
  (`epr1kit/stabcert.py`, gate G20). The certificates are six-variable linear rank
  inequalities with the purifier as a variable (60 classes up to S6; 16 certified
  facets of the qubit stabilizer cone). Screening the 1.9M-orbit A1' catalogue with
  them excludes 1,794,496 orbits (93.6 %; `results/2026-09-30_s17-screen/`). This
  contradicts the coincidence of the stabilizer and QLR cones at five parties
  suggested by BCHS (arXiv:2006.16292). Validation: zero violations on 1,130,413
  + 760 + 1,130 pool vectors, 200,000 random graph states with arbitrary party
  sizes, and all realisation certificates.
- **Targeted stabilizer search.** `cc/s16_anneal.c` + `scripts/s16_stab_search.py`
  (normal-form qubit graph states, annealing, independent re-verification) and
  `cc/s16_exhaust.c` (exhaustive, small N). Frontier: 144 realised (lambda 1/2/3:
  126/17/1), 139 excluded, 17 undecided (`data/a2_frontier_status.json`). Controls:
  37/37 pool witnesses re-realised; the two F3-conditional ones are now realised by
  qubit graph states at lambda = 2 (all 37 unconditional).
- **q2 is stabilizer-realisable**: 2*q2 is the entropy vector of a 12-qubit graph
  state (checked also by state-vector SVD); q2 itself (lambda = 1) is not
  (exhaustive: the 760 six-vertex graph-state vectors, recomputed independently).
- New data: `s16_realisations.npz`, `s17_exclusions.json`, `s17_ineq_classes.npy`,
  `s17_ineq_tight_rank.npy`, `a2_frontier300.npy`, `a2_frontier_status.json`;
  gate G20 (24 gate lines in `--full`); jobs `a2-frontier` (X) and `s17-ci` (S).
- **Erratum C6.** R26 stated that 292 unwitnessed orbits with maximum coordinate
  <= 5 remain; the correct number is 299 (only 53 of the 60 seeds have maximum
  coordinate <= 5; the R26 figure subtracted all 60).
- **Erratum C7.** Rounds R14-R26 listed the machine verification of inequality
  (4.4) ("64-row contraction table") as overdue. It was completed in R1 with a
  self-produced contraction certificate (SAT-CEGAR, exhaustive check of
  119,877,472 tuples) plus the exhaustive -1 violation by graph states; only the
  optional transcription of the printed table remains.

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
