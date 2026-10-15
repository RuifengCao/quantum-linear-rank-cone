# CHANGELOG

## R43 (2026-10-15) -- from six parties on, the qubit cone differs from that of every other local dimension; erratum to R42
- **Result (proved; checked by an independent referee, an agent; like Theorem 12 not yet checked by hand; not in the
  paper).** `docs/R32_reduction_and_local_dimension.md` §4.11, Corollary 13.
  - The *seven-point state* over Z_d, the CSS state (one Z_d qudit per system) of the code spanned by the seven
    nonzero vectors of {0,1}³, has f12 = −log₂(d / gcd(d, 2)) bits. It is the Steane state for d = 2 (f12 = 0) and
    the non-Fano state for odd primes d (f12 = −1 in units of log d).
  - Hence Stab^(Z_d)_n ⊄ Stab⁽²⁾_n for every d ≥ 3 and n ≥ 6. Theorem 12 covered d with an odd prime factor; new
    are d = 4, 8, 16, …. With Theorem 5, Stab⁽²⁾_n ⊊ Stab^(Z_d)_n for even d ≥ 4 from n = 6 on (for d = 4 the Z_4 cube
    state gave n ≥ 7); with Theorem 9, incomparability for odd d.
  - Proof: Theorem 12, Theorem 1(d) step 3 (CSS entropies over Z_d), and the value of f12, which follows from the
    Chinese remainder theorem and the Smith invariants of the column sets (all 1 except at W, X, Y: 1, 1, 2).
    Equivalently S^(d) = log₂d·S_NF for odd d and S_F + (log₂d − 1)·S_NF for even d.
  - The statement is about Z_d qudits: Galois qudits GF(2^k) give exactly the qubit cone (§4.6).
- **Erratum (R42).** §5 said "n = 6 for Z_{2^k}, k ≥ 2: open in both directions against qubits and against odd d".
  Two of these four directions were known: Stab⁽²⁾ ⊆ Stab^(Z_{2^k}) (Theorem 5), and Stab^(Z_{2^k})_6 ⊄ Stab^(Z_d)_6
  for odd d (Theorem 9 with that inclusion). The item is struck through with a note; the open one left is whether
  Stab^(Z_d)_6 ⊆ Stab^(Z_{2^k})_6 for odd d. The R42 entry below is left as it was.
- **Referee.** A separate agent, with its own code and without reading the new functions, confirmed the closed form
  exactly for d = 2..40, the entropies against state vectors (d = 4, 6 dense; d = 2, 3, 5, 8, 9, 10, 12 sparse), the
  Smith invariants by two methods, and the proof. Its corrections are in the text: the statement restricted to Z_d
  qudits (not GF(2^k)), the step n > 6, why the number of invariants is r_NF(X), the remark on why Theorem 12 stops at
  characteristic 2 (Lemmas T1–T5 appear to carry over to Z_4 with log|·|; only (LA)₂ fails), S^(4) = S_F + S_NF lies
  in Stab⁽²⁾ + Stab⁽³⁾, and the range 6 ≤ n ≤ 2p + 2 of the remaining open question.
- **Checks.** Gate G39: 2^f12 = gcd(d, 2)/d exactly for d = 2..16 (brute force over the d³ messages), and the Smith
  invariants of all 127 column sets (sympy). Gate G39b: the seven-point states over Z_4, Z_6, Z_8 from their state
  vectors (largest deviation of the 126 cut entropies about 2·10⁻¹⁵). Sandbox: exact values for d = 2..30.
- **Code.** `scripts/s34_six_party_survey.py`: `seven_point_sizes`, `seven_point_f12`, `seven_point_smith`,
  `seven_point_state_vector`, called by `transfer_check`; gate texts of G39/G39b in `selftest.py`.
- **Results.** `results/2026-10-15_r43-sandbox/` (sandbox, no server time).
- **README.** Status, §0 (R43), §7 (⓪″). `VERSION`: R43; `CITATION.cff` unchanged. No GitHub Release.

## R42 (2026-10-14) -- six parties, the other direction for the closed cones: qubits and every odd local dimension are incomparable from n = 6
- **Result (proved; checked by an independent referee, an agent; not yet by hand; not in the paper).**
  `docs/R32_reduction_and_local_dimension.md` §4.11.
  - *Theorem 12.* Stab^(Z_d)_6 ⊄ Stab⁽²⁾_6 for every d ≥ 3 that is not a power of 2, hence for every n ≥ 6; for
    n ≥ 7 this was known (Theorem 2 with Theorem 5), so the new case is n = 6. With Theorem 9, Stab⁽²⁾_n and
    Stab^(Z_d)_n are incomparable for odd d from n = 6 on (before: n = 7), and Stab⁽²⁾_n ⊊ Stab^(Z_d)_n for even d
    that is not a power of 2 (before: n = 7). This answers the question of the paper's Discussion that R41 had left
    open for the closed cones.
  - *Witness and certificate.* The non-Fano state S_NF, and a new inequality f12 ≥ 0 (31 terms on the six-party cuts)
    valid for every qubit stabilizer state, with f12(S_NF) = −1 and f12(S_F) = 0.
  - *Proof ("symplectic transfer").* For a state reduced to S(P) qubits per system, let Z₀ be the sum of the six line
    pieces of the stabilizer group and D_P its restriction to P. Isotropy bounds dim D_P by linear forms in the
    entropies: at a corner, three pairwise orthogonal restrictions (Lemma T1); at a diagonal point, the rank of the
    form is transported along a line from its two corners (Lemmas T2, T3); and z(ABCZ) is bounded below (Lemma T4).
    Lemma 4 (LA)₂ applied to the subspace arrangement dual to Z₀ (Lemma T5) then gives f12 ≥ 0. Only (LA)₂ needs
    characteristic 2. This is the quantitative version of Proposition 10 that R41 proposed.
  - *Side results.* g12 = f12∘π is a self-dual characteristic-2 linear rank inequality on seven variables with
    g12(h_NF) = −2: the inequality Proposition 8 asks for in this direction, not a combination of Shannon and the
    eight forms of Proposition 8. Theorem 12 does not use the Gross–Walter entropy formula (for composite d the
    witness is a Z_d CSS state).
- **Checks.**
  - Gate G39: f12 built from its definition and pinned in the manifest (`s34_transfer_form`: sha256 of the 128
    coefficients and the 31 terms); f12(S_NF) = −1, f12(S_F) = 0; f12 ≥ 0 on the GF(2) list of G38 (minimum 0, 38
    zeros); on the GF(3) list f12 < 0 only at S_NF.
  - Gate G39b: every step of the proof (`transfer_chain`) on 40 states: random compressed qubit graph states,
    compressed CSS states of GF(2) families pushed towards the non-Fano dependencies, Steane copies after random
    local symplectic maps; the dual arrangement built explicitly; the compression (symplectic reduction). A sandbox
    run of 400 states found no failure.
  - Adversarial search [empirical]: `cc/s37_minform7.c` (`--transfer-anneal`) minimises f12 over qubit graph states;
    nine sandbox runs (1–3 qubits per system, two pulled towards 2·S_NF and 3·S_NF, about 6·10⁸ moves) never went
    below 0, and every reported minimiser passed every step of the proof.
- **Referee.** A separate agent, reading only the note and the data, recomputed the 31 terms (no difference) and
  the values at S_NF and S_F, checked Lemmas T1–T5, every term of (LA) and the direction of every bound, and with its
  own code checked every step on 26,549 qubit stabilizer states, all 2^21 seven-vertex qubit graph states, about 6 M
  random states and 78 M annealing proposals: no negative value; weakened forms do go negative. Its corrections are
  in the text: the new case is n = 6 (n ≥ 7 already followed from Theorems 2 and 5), no dependence on
  Gross–Walter, g12 over infinite fields of characteristic 2, notation, and the passages of the note that still
  called the question open (header, §4.7, §4.8, §4.10, §5).
- **Code.** `scripts/s34_six_party_survey.py`: `transfer_form`, `transfer_terms`, `transfer_folded`, `transfer_sha`,
  `transfer_chain`, `compress` (symplectic reduction to S(P) qubits per system), `transfer_check`, `transfer_ok`,
  `transfer_anneal`; CLI `--transfer [--full --trials N]`, `--transfer-anneal`, `--pull-nf`, `--seed`. New
  `cc/s37_minform7.c`. `selftest.py`: gates G39, G39b. `data/manifest_shas.json`: entry `s34_transfer_form`.
- **Results.** `results/2026-10-14_r42-sandbox/` (sandbox, no server time).
- **README.** Status, §0 (R42; R41, R36 and R35 entries annotated), the s34 row of the scripts table, §7 (⓪″).
  `VERSION`: R42; `CITATION.cff` unchanged. No GitHub Release.

## R41 (2026-10-14) -- six parties, the other direction: no qubit realisation of the non-Fano state at any multiplicity; erratum to R40
- **Result (proved; checked by an independent referee, an agent; not yet by hand; not in the paper).**
  `docs/R32_reduction_and_local_dimension.md` §4.10.
  - *Proposition 10, strengthened.* For every field F = GF(q) and every m ≥ 1, an F-linear stabilizer state of the
    seven systems whose entropy vector is m·S_NF has, at every point, one Lagrangian subspace as the restriction of
    all its line pieces. The new step 1′: every non-Fano line contains two of A, B, C, Z, where the restrictions are
    already Lagrangian (R40's step 1), so isotropy forces the restriction at the third point W, X or Y to be
    isotropic, hence Lagrangian, hence equal to that of the other line. R40's step 2 then gives an m-fold
    representation of the non-Fano matroid over F, so q is odd. Hence **no qubit stabilizer state has an entropy
    vector that is a positive multiple of S_NF**, with any number of qubits per system.
  - *Corollary 11 (uses the certified face ranks of §4.8, G35).* S_NF is not a finite nonnegative combination of qubit
    stabilizer entropy vectors, and h_NF is not a finite nonnegative combination of rank functions of GF(2) subspace
    families: a summand would lie on the face cone(r_NF, r_NF*), and the doubled CSS state of such a family would have
    entropy vector (a + b)·S_NF. So for every odd p the cone generated by six-party GF(p) stabilizer entropy vectors
    is not contained in the one generated by qubit ones. For the closed cones (the paper's open question) this
    remains open: S_NF could be a limit of such combinations.
- **Erratum (R40).** R40's Proposition 10 had (c) for m = 1 only and a part (d): "for q = 2, in every qubit state
  with entropy vector m·S_NF and m qubits per system, the two line pieces at one of W, X, Y restrict to different,
  non-Lagrangian subspaces". By step 1′ there is no such state, so (d) was true only vacuously, and the text built on
  it described a structure that cannot occur: the §4.8 item, "What this shows" in §4.10, §5, the README entry and
  ⓪″ (which proposed a targeted search at m = 2), and the R40 entry below. Neither R40 nor its referee noticed
  step 1′. The note and the README are corrected (the README keeps R40's sentence struck through with a note); the
  R40 entry below is left as it was.
- **Wording.** §4.8 "Faces": "the face of Γ₇ containing h_F" now reads "the smallest face"; §4.8 reverse-direction
  remark rewritten (exact points: none; only limits remain); §5 n = 6 item rewritten.
- **Check.** `line_pieces_check` in `scripts/s34_six_party_survey.py` (part of its self-test and shown in gate G38b):
  on the CSS states of the non-Fano code over GF(3) and of the Fano code over GF(2), with one and with two qudits
  per system, after random local symplectic maps, every line piece has dimension m and the restrictions at every
  point coincide and are Lagrangian. This checks the conventions on states that exist, not the proposition.
- **Referee.** A separate agent checked Proposition 10 (with step 1′), Corollary 11 and the erratum line by line,
  recomputed the face ranks, and confirmed steps 1, 1′ and 2 numerically (up to three qutrits per system for the
  non-Fano state). Its wording points are taken into the text.
- **Results.** `results/2026-10-14_r41-sandbox/` (sandbox, no server time).
- **README.** Status, §0 (R41; the R40 entry annotated), §7 (⓪″). `VERSION`: R41; `CITATION.cff` unchanged. No
  GitHub Release.

## R40 (2026-10-14) -- five parties: all 641 known extreme rays over every local dimension (roadmap item S5(a), finished); six parties, the other direction (roadmap item M3, steps 1-2)
- **Result 1 (certified; checked by an independent referee, an agent with its own code; not yet by hand; not in the
  paper).** All 641 known extreme-ray orbits of Stab₅ lie in Stab^(Z_d)_5 for every d ≥ 2, hence are extreme rays of
  it (`docs/R32_reduction_and_local_dimension.md` §4.9): 566 by unimodular realisations (unchanged) and 75 by covering
  certificates. R40 covered the five rows R39 left open (587, 589, 612, 630, 635).
  - **What closed them.** Rational searches of 1200 s from the qubit certificate and from R39's best near miss gave
    the first rational realisations of rows 587 and 589. Rounds of 600 s from the best near miss so far, at start
    temperature 0.3 instead of 1.0, gave those of rows 630, 635 and (in a later push) 612. Searches that are rational
    and right mod the uncovered primes, started from these realisations, then covered all five (60 s; row 612: 300 s,
    third round). R39's near misses had missed one to three cuts: the searches had not finished.
  - New pieces: 9 rational realisations (κ = 1) with bad primes, minimal covering subsets as in R39, and the 5 stored
    qubit certificates as GF(2) pieces. Every row's bad primes meet in {2}; row 630 needs only one rational piece.
  - Every one of the 641 rows now has an explicit qutrit realisation in the data (a unimodular certificate mod 3 or a
    covering piece right at 3).
  - Consequence (§5 of the R32 note): no known extreme ray of Stab₅ can witness a dimension dependence at n = 5. A
    witness would need an extreme ray not yet known or, for d divisible by a square, a vector of Stab^(Z_d)_5 outside
    Stab₅.
- **Result 2 (§4.10; status as marked; independent referee; not in the paper).** Six parties, the other direction: is
  Stab⁽ᵖ⁾_6 ⊄ Stab⁽²⁾_6 for some odd p? The only candidate so far was the non-Fano state h_NF (§4.8).
  - **Survey [certified: exhaustive].** Every stabilizer state of seven systems with at most one qudit per system is
    local-Clifford equivalent to a weighted graph state on seven vertices; `cc/s35_enum7.c` enumerates them up to
    scaling: 10,773 distinct entropy vectors over GF(2) (all 2^21 graphs) and 27,334 over GF(3) (255,602,200 normal
    forms). The 30 qubit-only vectors are the relabellings of the Steane state. The 16,591 qutrit-only vectors form 53
    orbits of S₇; for 52 of them `cc/s34_anneal7.c` found a 14-qubit graph state realising 2·S (exact ranks; none is a
    product of two states with one qubit per system). The remaining orbit is that of h_NF: among the states over Z_3
    with at most one qudit per system, it is the only possible witness.
  - **Proposition 10 [proved].** If a GF(q)-linear state of the seven systems with m qudits each has entropy vector
    m·S_NF (units of log q), then the line pieces at the four points on three lines restrict to one Lagrangian subspace
    each; if this also holds at the three diagonal points, the non-Fano matroid has an m-fold representation over
    GF(q), impossible in characteristic 2. For m = 1 it always holds. So S_NF is a GF(q)-linear entropy vector iff q is
    odd, no GF(2^k)-linear state has entropy vector k·S_NF, and every qubit realisation of m·S_NF (m qubits per
    system) has two different, non-Lagrangian line pieces at a diagonal point. Checked by exhaustive search:
    `cc/s36_target7.c` over GF(2), GF(3), GF(4), GF(5), GF(7), GF(8), GF(9) (S_NF only for odd q, S_F only for even
    q, random targets everywhere) and, for S_NF, GF(16) (none; 75 min); a separate Python implementation reproduces
    the counts for q ≤ 4 (gate G38b).
  - **Annealing [empirical].** No qubit state with entropy vector 2·S_NF or 3·S_NF was found (seven runs at m = 2 and
    three at m = 3, up to 600 s each). Every run stopped at residual m, the residual of m copies of the Steane state
    (S_F is the only GF(2) vector at L1 distance 1 from S_NF, and no sum of two GF(2) vectors is closer than 2 to
    2·S_NF; certified, G38b).
- **Code.** `cc/s35_enum7.c` (the enumeration, OpenMP), `cc/s34_anneal7.c` (annealing of qubit graph states with m
  qubits per system towards a target vector), `cc/s36_target7.c` (exhaustive search for graphs over GF(q),
  q = 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, with a given cut-rank vector), `scripts/s34_six_party_survey.py` (`--check
  [--full]`, `--self-test`, `--enumerate`, `--survey`, `--nonfano`, `--fields`, `--fields-python`).
- **Data and gates.** New `data/s34_seven_qudit_vectors.npz` (the two lists, the 53 orbit representatives with their
  GF(3) graphs, the 52 qubit graph states). `data/s32_stab5_cover.npz`: 75 rows, 199 pieces (125 rational
  realisations, 3 of them at κ = 2; 74 GF(2) / GF(3) realisations). Manifest: new entry `s34_seven_qudit_vectors`
  (with the counts of `cc/s36_target7.c`), `s32_stab5_cover` updated. New gates G38 (data checks) and G38b (`--full`:
  the GF(2) list recomputed in numpy, 20,000 random GF(3) graphs, the tensor-product and distance checks, the Python
  target search against the stored counts, the s34 self-test). `selftest.py --full` now has 49 gate lines.
- **Referees.** Two separate agents, each with its own code and without reading the repository's code for the result:
  - S5(a): re-verified the 14 new pieces (structure, ranks on all 31 cuts, bad primes by its own Smith form, by minors
    and by sympy) and the coverage of the five rows; the 185 earlier pieces and the unimodular file are unchanged, and
    every row 0-640 is certified one way or the other. No discrepancy.
  - M3, two rounds: re-enumerated both lists (its GF(3) normal form differs; same graph count) and found them equal to
    the stored ones; confirmed the orbits, sizes, representatives, the 52 qubit certificates, the Steane and non-Fano
    orbits, the tensor-product and distance statements, and with two exhaustive programs every count of
    `cc/s36_target7.c` from GF(2) to GF(9); checked Proposition 10 line by line (step 2 now uses its shorter argument).
    Its first round found the text wrong or overclaimed in five places, all corrected: §4.8 equated the exact points
    a = b = m with qubit states with m·S_NF (they give them; conversely an exact point gives 2m·S_NF); "all states"
    needed "over Z_3"; the residual was undefined and its baseline (m Steane copies) unstated; controls were missing
    or omitted; the search class and its size were not stated. Its second round tightened four sentences (the scope
    of the field-linear statement, the controls, the linear and non-linear conditions of the search, the n = 5
    sentence).
- **Results.** `results/2026-10-14_r40-sandbox/` (2-core sandbox, no server time): the five-party campaign log after
  the R39 cutoff with its raw kernel output, summary and cover check; the six-party enumeration, survey, non-Fano
  runs and controls, field runs, checks and self-tests; `selftest.py --full`.
- **Note and README.** R32 note: header (R40 paragraph), §4.7 (pilot), §4.8 (reverse direction, R40 items), §4.9
  (counts, round R40, the rows open after R39), §4.10 (new), §5 (n = 5, n = 6, five parties over every field), §7.
  README: status, §0 (R40; the R39 entry's open rows), §2, §2b (job table), §3 (selftest, s32, s34), §6 (data),
  §7 (⓪′ done, new ⓪″ for M3), §9 (CI gate count and time caps).
- `.github/workflows/selftest.yml`: time caps raised to 45 min for the job and 25 min for `selftest.py --full`
  (G38 and G38b add about three minutes: 8 min in the sandbox, where the lrs part of G29b is skipped; CI runs it).
- `jobs.json`: the server job `s32-unimodular` has nothing left to do (kept for rows added later). `VERSION`: R40;
  `CITATION.cff` unchanged. No GitHub Release for R40 (it would mint a new Zenodo DOI).

## R39 (2026-10-14) -- five parties: 636 of the 641 known extreme rays over every local dimension (roadmap item S5(a), continued)
- **Result (certified; checked by an independent referee, an agent with its own code; not yet by hand; not in the paper).** 636 of the 641 known extreme-ray orbits of Stab₅
  lie in Stab^(Z_d)_5 for every d ≥ 2 (577 in R38), hence are extreme rays of it (`docs/R32_reduction_and_local_dimension.md`
  §4.9): 566 by unimodular realisations (κ = 1: 515, 2: 50, 3: 1) and 70 by covering certificates. R39
  closed 59 of the 64 rows R38 left open: 2 by unimodular realisations, 57 by covering certificates.
  - **What closed them.** The covering search now starts each search "rational and right mod q" from a rational
    realisation of the row (the one with the fewest uncovered primes) and asks for all uncovered primes at once; R38
    started these searches from the qubit certificate, one prime at a time. This covered all 21 rows whose rational
    realisations from R38's runs left primes uncovered, within a few rounds of 40 s, although for 17 of them every one
    of those realisations was bad at 3. Longer rational searches from the qubit certificate (90 to 900 s per row) found
    the first rational realisation of 34 more rows, and the same search then covered them. For row 200 (qubits only at
    multiplicity 2) rational realisations were found only at κ = 2, from its qubit certificate at multiplicity 2; at
    κ = 1 the best attempts miss one cut.
  - **Outputs of R38 that R38 did not use.** R38's strategy probe had reported unimodular realisations of rows 539 and
    631 (κ = 1; exact, as checked now), but its output was not merged into the data. Row 574 is covered by two rational
    realisations from R38's runs (bad primes {163, 181, 2833, 7699} at κ = 2, and {2, 3, 5}: empty intersection). Row
    607 is covered by the matrix of R38's rejected report (Smith invariant 127, so bad primes {127}) and a rational
    realisation with bad primes {2, 3, 5}.
  - **Open: 5 rows** (587, 589, 612, 630, 635). No rational realisation is known for any of them; their best attempts,
    with the cuts they miss, are in `results/2026-10-14_r39-sandbox/open_rows.json`. The best attempts miss one to three
    cuts (row 612: ten or more), mostly with the rational rank above the target. This is what a configuration that
    exists only in characteristic 2 would produce, and also what an unfinished search produces [empirical].
  - The exploratory screen `s33_ci_oddchar.py --pairs all` (every pair, about 1,330 per row) on the 22 rows open at
    the time excludes none.
- **Errata (R38).**
  - `results/2026-10-13_r38-sandbox/probes/strategy_sample20.txt` counted kernel reports as certified. The report for
    row 574 under "ring Z[omega] from the qubit certificate" is not unimodular: bad primes 163, 181, 2833 and 7699
    (the kernel's extra checks then stopped at 97). That setting certified 2 of 20 (465, 541), not 3. The other reports
    are exact. A dated note is added to the file. Nothing in R38's data or claims rests on it.
  - The header and §5 of the R32 note said that 564 of the 641 rays lie in Stab^(Z_d)_5 for every d and that 64 are
    open, leaving out the 13 rows with covering certificates; 577 was meant, as in §4.9, the README and the R38 entry.
    Corrected.
- **`scripts/s32_unimodular.py`.** `--cover` follows the search order above. New options: `--rounds`; `--rings`
  (ring strategies for rows without a rational realisation); `--kappa2` (multiplicity 2 from two copies of the qubit
  certificate; a row whose qubit certificate is at multiplicity 2 gets this search anyway). For a newly covered row a
  smallest covering subset of its rational realisations is stored (`minimal_cover`; the rows covered in R38 keep all
  their pieces); a rational realisation without bad primes goes to the unimodular file. Certification is unchanged.
- **Data and gates.** `data/s32_stab5_unimodular.npz`: 566 rows (rows 539, 631 added).
  `data/s32_stab5_cover.npz`: 70 rows, 185 pieces (116 rational realisations, 3 of them at
  κ = 2; 69 GF(2) / GF(3) realisations). Manifest updated; gates G36, G36b and G37 unchanged (they take the
  counts from the manifest); `selftest.py --full` still has 47 gate lines.
- **Referee.** A separate agent, with its own code (its own Smith normal form, Bareiss ranks with minors, and sympy's
  Smith normal form; two implementations of the ranks over GF(p)), confirmed that every R38 entry of the two data files
  is unchanged and re-verified every new entry: the unimodular certificates of rows 539 and 631, the 94 new rational
  realisations with their exact bad primes (all 116 stored bad-prime sets agree), the 55 new GF pieces, and the
  coverage of all 70 rows. No discrepancy; every number about the data files in this entry and in §4.9 agrees with its
  computation. Its wording remarks (row 200, the count of 21 rows, row 607, the open rows, §5) are taken into the text.
- **Results.** `results/2026-10-14_r39-sandbox/` (2-core sandbox, no server time): campaign log and raw kernel output,
  the open rows with their best attempts, checks, self-tests, the screen, the calibration for p = 3.
- **Note and README.** §4.9 of the R32 note (R39 paragraph, counts, open rows), its header and §5 (corrected) and §7;
  README: status, §0, §2b, §3, §6, §7 (⓪′).
- `jobs.json`: the server job `s32-unimodular` now runs the R39 covering search on the 5 open rows (not run).
  `VERSION`: R39; `CITATION.cff` unchanged. No GitHub Release for R39 (it would mint a new Zenodo DOI).

## R38 (2026-10-13) -- five parties: the known extreme rays over every local dimension (roadmap item S5(a))
- **Result (certified; checked by an independent referee, an agent with its own code; not yet by hand; not in the paper).** 577 of the 641 known extreme-ray orbits of Stab₅
  lie in Stab^(Z_d)_5 for every d ≥ 2, hence are extreme rays of it (`docs/R32_reduction_and_local_dimension.md`
  §4.9). This is the five-party version of the four-party lift behind Theorem 2 of the paper (Theorem 7 of the note).
  - 564 rows have a *unimodular realisation*: an integer symmetric matrix W whose 31 cut blocks have Smith
    invariants 0 and 1 only and rational rank κ·r(X); W mod d is a graph state of Z_d qudits with entropy vector
    κ·r for every d at once. κ = 1 for 513 rows, κ = 2 for 50, κ = 3 for 1.
  - 13 more rows have *covering certificates*: rational realisations (right rational cut ranks) whose bad
    primes intersect inside the primes of stored qubit / qutrit realisations of the row. Every prime is covered, and
    every d by Corollary 3 of the paper.
  - The 12 BCHS ledger rays without an explicit realisation (rows 6, 7, 15, 34, 39, 40, 42, 43, 44, 46, 48, 50) now
    have one, at κ = 2: all 641 known extreme rays have explicit qubit realisations (629 before), and
    577 have explicit qutrit realisations (96 before).
  - Open: 64 rows (200, 488, 492, 493, 495, 496, 497, 501, 502, 503, 504, 506, 508, 509, 525, 539, 549, 556, 557, 559, 560, 564, 569, 572, 573, 574, 575, 576, 578, 579, 580, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 593, 594, 595, 596, 597, 598, 605, 606, 607, 609, 610, 611, 612, 613, 630, 631, 633, 634, 635, 636, 637, 638, 640). For 22 of them a rational realisation is known, and the primes
    it leaves uncovered are listed in `results/2026-10-13_r38-sandbox/cover_open.json`.
- **`cc/s32_unimodular.c` (new).** Search kernel: annealing over integer entries in {-1, 0, 1} with ranks over chosen
  primes (default GF(2), GF(3) and modulo 2^31 - 1 for Q), extra checks for every prime below 1000 and 2^61 - 1;
  neutral moves (party-local unimodular basis changes, integer local complementations); restriction of scalars from
  five quadratic orders with the unimodular trace form (κ = 2). It proves nothing by itself.
- **`scripts/s32_unimodular.py` (new).** Search strategies (plain, neutral, ring-*, qutrit-start, gf3-search-start,
  and `--cover` for covering certificates); exact certification (`--check`: diagonal forms over Z; `--full` adds the
  rational rank with t x t minors of gcd 1, ranks over GF(p) for p <= 19, and Z_4 / Z_6 state vectors of the small
  certificates; `--check-cover`: bad primes recomputed and coverage); `--self-test` (neutral moves keep every cut
  Smith form; ring expansion; the Z_d cut formula against state vectors). One report of the kernel (row 607, Smith
  invariant 127) failed the exact check and was rejected; the kernel's extra checks were then extended from primes up
  to 97 to primes below 1000.
- **Data and gates.** `data/s32_stab5_unimodular.npz` (gates **G36**, **G36b**) and `data/s32_stab5_cover.npz`
  (gate **G37**); manifest entries. `selftest.py --full` now has 47 gate lines.
- **`scripts/s33_ci_oddchar.py` (new, exploratory).** Common information of every disjoint pair plus the seven-variable
  forms that hold in odd characteristic (DFZ (65), (LB)_2, BKL (6.17) corrected; 8,400 relabelled instances) on the
  seven-element extension. It excludes none of the rows tested (floating-point verdicts).
- **Referee.** A separate agent, with its own code (checked against sympy's Smith normal form), re-verified all 564
  unimodular certificates, all 36 covering pieces and the coverage of all 13 rows, simulated 331 graph states for
  d = 2, 3, 4, 6, and checked the arguments (lemma, covering, ring construction, neutral moves). No discrepancy; it
  asked for two hypotheses to be stated (u ≠ u′ in the basis change; the column convention of M(v)), now in §4.9.
- **Results.** `results/2026-10-13_r38-sandbox/` (2-core sandbox, no server time): campaign log, summary, checks,
  self-tests, the covering run, the screen, probes.
- **Note and README.** §4.9 of the R32 note (new); §4.7, §5 and §7 updated. README: §0, §2b, §3, §6, §6a, §7 (⓪′),
  §9. The five-party pilot of R35 (`s30_dimension_profile.py --n5-pilot`, 40 of the first 60) is superseded.
- `jobs.json`: server job `s32-unimodular` (not run). No earlier data changed. `VERSION`: R38; `CITATION.cff`
  unchanged (it describes the archived release R37). No GitHub Release for R38 (it would mint a new Zenodo DOI).

## R37 (2026-10-11) -- archived release for the paper; verification scripts; licences and citation metadata
- **Paper.** R37 is the version cited by *Stabilizer entropy cones as self-dual slices of linear rank cones: a
  five-party gap and dependence on the local dimension* (X. Song and R. Cao, 2026). Its results are those of rounds
  R27-R36 (`docs/A2_stabilizer_vs_QLR.md`, `docs/R32_reduction_and_local_dimension.md`); Table VI of the paper maps
  them to the notes, the gates and the scripts in `checks/`. The proofs used in the paper and the transcriptions they
  rely on were checked line by line by hand before submission, and its literature claims were checked against the
  sources.
- **`checks/` (new).** The paper's verification scripts: ten scripts `check_*.py`, the helper module `common.py`,
  the C helper `q2_gf3.c`, `README.md` (what each script covers) and `references_log.md` (how each reference of the paper was
  checked). They read `data/` and `results/` but not the repository code, except for the transcriptions of the
  published seven-variable inequalities used by `check_sec6.py`. Run from `checks/`, every script ends with
  `ALL OK`, except `check_sec4b_facet.py` and `check_facets119.py`, which print counts (1,414 tight vectors of rank
  30; 117 facets, 119 with `--with-literature`). `check_sec4a.py` also recomputes, from state vectors, the state
  (R6) of Linden et al. (as printed, S(ab) = 1 instead of 2 on ray 6 of their Table 1, and its entropy vector is not
  an extreme ray) and the qubit state that N. Osthus gave for that ray.
- **Notes.** The status lines of both notes now point to the paper, which is authoritative where they differ. The
  R32 note §6 lists the prior work found by the literature check: Matúš 2012 (the lift and fact (c) of Theorem 1);
  Jowett-Mo-Whittle 2016; Sarvepalli-Raussendorf 2010 (the matroid case of the converse of Theorem 1); Osthus 2024
  (for prime d, the CSS entropy formula and the doubling; the question whether the local dimension matters; the
  misprint in (R6)); Peña-Macias-Sarria 2021 (inequalities for the configuration J - I for every t; the failure of
  BKL's (6.17)); Walter-Witteveen 2021 with BCHS (Theorem 7 for prime d); the conjecture of Linden et al.;
  Feng-Luo 2024. §3.2, §4.4, §4.5 and §5 of the R32 note are updated accordingly; in particular, for squarefree d
  the Z_d cone lies in the closed convex hull of the prime cones, so the question of §5 is open only when d is
  divisible by the square of a prime.
- **Conventions.** README §6a documents the mask order of subsets, the bits of the purifier and of the common
  information, the block order of realisation certificates, and the two orders of the Shannon elemental inequalities
  used by stored certificates (`epr1kit.stabcert.elemental` and `elemental7` in `scripts/s30_dimension_profile.py`).
- **Human checks.** README §7: the eleven single-source lines of `data/dfz_ref28.csv` and the Ingleton form (39)
  were checked by hand against DFZ; the novelty and literature checks of item ⓪ were done for the paper.
- **Licences and citation.** `LICENSE` is now the full MIT text (code: `epr1kit/`, `scripts/`, `cc/`, `checks/`,
  `selftest.py`); `LICENSE-DATA` states CC BY 4.0 for `data/` and `results/` and the GPL-3 licence of
  `bin/normaliz`. `CITATION.cff`: title "quantum-linear-rank-cone: Data, certificates and code for stabilizer
  entropy cones" (as cited in the paper), version R37, ORCIDs and affiliations. `requirements.txt` adds sympy (needed
  by `checks/` only). `VERSION`: R37.
- No code in `epr1kit/`, `scripts/` or `cc/` and no data changed.

## R36 (2026-10-02) -- six parties already separate qubits from every odd local dimension
- **Theorem 9 (proved).** Stab^(2)_n is not contained in Stab^(Z_d)_n for every odd d and every n >= 6
  (`docs/R32_reduction_and_local_dimension.md`, section 4.8). Before, this was known from n = 7 (Theorem 2), and the
  smallest number of parties with a dimension dependence was 5, 6 or 7; it is now 5 or 6.
  - Witness: the seven-qubit CSS state of the [7,3] simplex code (Fano matroid), one qubit per party; this is the
    logical |0> of Steane's [[7,1,3]] code. h_S = r_F + r_F^* (no coloops); checked against the 128-amplitude state
    vector.
  - Certificate: the integer form f9 (67 terms, `data/s31_six_party.npz`) with f9 >= 0 in every odd characteristic
    and for abelian groups of odd order, and f9(h_S) = -3.
  - Method (transfer): F_i^L = U_i cap (U_j + U_k) for the Fano lines L through i, Fhat_i = sum_L F_i^L. Lemmas 0,
    A-E bound every rank of the derived family by linear forms in h (Lemma A: F_i^L in F_j^L + F_k^L). Lemma 4 (LB)_2
    applied to the derived family and bounded term by term gives f9/3. All bounds are attained at h_F, so
    f9(h_F)/3 = (LB)_2(Fano) = -1.
  - g = f9 + f9 o D is the self-dual characteristic-dependent inequality that Proposition 8 asked for
    (g(h_F) = -6). So f9 is not a nonnegative combination of Shannon and the eight forms of Proposition 8 in any
    relabelling.
  - Invariant form: on Fano-invariant self-dual vectors, f9/3 = 210 beta + 30 alpha + 21 nu - mu (alpha = pairwise
    mutual information, mu = I(i; L - i), nu = I(i; non-line pair), beta = v_P - 2 mu >= 0). With DFZ (65) as the base
    form the constants become 154, 26, 5 (computed, not refereed).
  - Exact faces: the Shannon faces at h_F and at h_NF are two-dimensional, cone(r, r^*) (ranks 125 and 126,
    certified). So h_F is outside LR^(q)_7 iff no nonzero a r_F + b r_F^* is inside, and the reverse direction at
    n = 6 for the non-Fano state is the same question for cone(r_NF, r_NF^*) and LR^(2)_7 (open).
- **Referee.** An independent referee (separate agent, own code) confirmed Lemmas 0, A-E, the composition, f9
  (recomputed exactly, no difference), the values, Theorem 1(d) and Lemma 4(G) for odd composite d, and found no
  violation of f9 >= 0 in about 8.2 M arrangements each over GF(3), GF(5), GF(7), all 105,413,504 seven-tuples of
  points of PG(2,3) or zero, and about 13,000 subgroup families of odd-order abelian groups. Over GF(2), f9 < 0
  exactly on the 168 Fano frames. Wording fixes applied: one direction only, scope of Proposition 8, dependence on
  Theorem 1(d), invariant form stated on the slice, status flag.
- **BKL.** The seven-variable inequalities of Blasiak-Kleinberg-Lubetzky (arXiv:1108.2489v1, eqs. (6.8), (6.17))
  were transcribed and checked term by term against the rendered PDF (`data/bkl_char_dependent.csv`).
  - Eq. (6.17) is false over every field as printed: A, B, C, Z = <e1>, <e2>, <e3>, <e1 + e2 + e3>, W = X = Y = 0
    gives -3. The proof bounds three step-4 terms in the text but drops them from (6.16). The corrected form equals
    (LB)_2 - 6 I(A;B) and passes the exhaustive checks; only it is used. (Pitfall 21; section 4.5.)
  - Proposition 8 now covers eight forms: BKL (6.8) and (6.17) corrected also have exact Shannon certificates for
    f + f o D (supports 59 and 98, denominators 1 and 2).
  - Also noted: the proof of BKL Thm 6.4 ends with "Lambda_odd . r(F) < 0" where Lambda_even . r(N) < 0 is meant.
- **Reverse direction (exploratory, not certified).** The same transfer with (LA)_2 gives +27 on the non-Fano state;
  no family of 2-dimensional subspaces over GF(2) realises h_NF (exhaustive search). Open.
- **Code and data.**
  - New `scripts/s31_six_party.py` (`--derive [--out]`, `--witness`, `--invariant`, `--face`, `--lemmas N`,
    `--stress N`, `--exhaustive`).
  - `scripts/s30_dimension_profile.py`: eight forms, `bkl_check`, `--bkl`.
  - New data: `data/bkl_char_dependent.csv`, `data/s31_six_party.npz`; `data/s30_selfdual_shannon_certs.npz`
    regenerated with eight certificates (the first six unchanged). Manifest updated.
- **Gates.**
  - G33: now 8/8 certificates.
  - G34 (BKL values on the Fano / non-Fano configurations and on the counterexample; transcription sha) and G34b
    (`--full`: exhaustive minima 0, 0 and -3 for the printed (6.17)).
  - G35 (Theorem 9: f9 recomputed = stored, state vector, values, invariant form, exact faces, lemma checks with the
    actual subspaces, random arrangements) and G35b (`--full`: exhaustive {0,1}-point configurations of GF(2)^3,
    GF(3)^3, GF(5)^3; 200 families per field for the lemma checks).
  - 44 gate lines.
- **Docs.** R32 note: header, section 4 intro, 4.2 (BKL), 4.5 (BKL misprint), 4.7 (summary, Proposition 8), new 4.8,
  sections 5, 6, 7. README: R36 status, pipeline rows s30/s31, pitfall 21, data inventory, gate counts. Raw output:
  `results/2026-10-02_r36-sandbox/`.

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
