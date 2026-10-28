# quantum-linear-rank-cone

Machine-certified data and tools for the **five-party quantum linear rank cone**
(QLR₅): an exact H-representation, certified extreme rays, and a self-test
suite of verification gates. Self-contained: no `psitip` dependency, all input
data bundled, a static Normaliz binary included, resumable long runs.

> **Status (R48, 2026-10-20).** R38–R48 add results obtained after the paper (§0,
> first items); they are certified by the gates (the computations of R48 by the scripts in `checks/` only) and
> checked by an independent referee (an agent), not yet by hand. R47 adds a consequence for quantum secret sharing: for the access structure of Steane's code, every
> stabilizer scheme of odd prime local dimension has information ratio at least 4/3, while qubits admit an ideal
> scheme. R48 shows that the five-party inequality (★) and 118 further single-CI inequality classes are facets of
> the stabilizer cone of every local dimension, and that stabilizer schemes with information ratio 2 exist for
> both access structures of R47 (linear schemes give stabilizer schemes); its second and third addenda lower this to
> 7/5 for both: for the access structure of Steane's code and every odd local dimension, and for the other access
> structure with qubits. Its fifth addendum shows that among all Z_d only qubits admit an ideal stabilizer scheme
> for the access structure of Steane's code, and only the odd d for the other. Its seventh addendum extends the bound
> 4/3 for the access structure of Steane's code from odd prime (and squarefree) to every odd local dimension and,
> counted in qudits, to every d that is not a power of 2. Its eighth addendum changes no result: it updates the
> numbers of lemmas, propositions, remarks, tables, equations and Supplementary Notes of the revised paper quoted in
> `checks/`, the notes and this README to the version submitted in October 2026 (theorem and corollary numbers are
> unchanged).
> Release R37
> is the archived version cited by the first version of the paper *Stabilizer entropy cones as self-dual slices of
> linear rank cones: a five-party gap and dependence on the local dimension* (X. Song and R. Cao,
> 2026, [arXiv:2610.18905v1](https://arxiv.org/abs/2610.18905v1)); release R48 is the version cited by its revised
> version (October 2026), whose Supplementary Table S5 maps each result to the notes, the gates and the scripts in
> `checks/`.
> Everything in `data/` is machine-certified by the
> gates in `selftest.py`, and the numbers quoted in the paper are recomputed by the
> separate scripts in `checks/`. The proofs used in the paper and the
> transcriptions they rely on were checked line by line by hand before the paper
> was submitted, and its literature claims were checked against the sources.
> The rest of the repository remains a research notebook with reproducible
> artefacts, not a refereed publication; where the notes in `docs/` and the
> paper differ, the paper is authoritative. `CHANGELOG.md` records every round
> and all errata.

The repository grew out of an internal project (code name *EPR-1*, "complete
characterisation of the quantum entropy cone", plan A: reduce the n = 5 case to
decisions about three cones). Internal shorthand that survives in file names
and the changelog is explained in the glossary at the end.

---

## The paper and its verification scripts (R37)

- **The paper** ([arXiv:2610.18905](https://arxiv.org/abs/2610.18905)). It shows that, after a linear lift,
  the n-party stabilizer entropy cone of local dimension d is the self-dual slice of the cone of rank
  functions of n + 1 subgroups of abelian groups of exponent dividing d (for
  prime d, of subspaces over GF(d)). From this it derives a five-party gap
  between the stabilizer cone and QLR₅ and the dependence of the stabilizer cone
  on the local dimension. Its results are those of rounds R27-R36, written up in
  `docs/A2_stabilizer_vs_QLR.md` and `docs/R32_reduction_and_local_dimension.md`;
  Table VI of the paper maps each result to these notes, to the gates of
  `selftest.py` and to the scripts below. (R47, R48) A revised version of the
  paper (October 2026) also contains the results of R38–R40, R42–R44, R46
  (§4.13), R47 and R48 and cites release R48; its Supplementary Table S5 plays the
  role of Table VI. The scripts `checks/check_rev_*.py` check these results, and
  `checks/README.md` maps every script to both versions.
- **`checks/`.** The verification scripts of the paper. They read `data/` and
  `results/` but not the code of the repository, except for the transcriptions
  of the published seven-variable inequalities used by `check_sec6.py`. Run
  them from `checks/`, e.g. `cd checks && python3 check_sec3.py`. They need
  numpy, scipy and sympy; `check_sec4a.py` also needs `lrs`, and `check_sec3.py`
  compiles a small C helper with `gcc`. Each script ends with `ALL OK`, except
  two that print counts; `checks/README.md` lists what each one covers.
- **Conventions** for subsets, the purifier and the stored certificates: §6a.

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
- [x] **R48 (2026-10-20): (★) is a facet of the five-party stabilizer cone of every local dimension; secret sharing schemes with ratio 2** (`docs/A2_stabilizer_vs_QLR.md`
  §3, R48 addendum; proved; computation checked by `checks/check_rev_sec4.py`, no repository data; not yet checked by hand).
  Six-vertex graph states with 0/1 weights have cut blocks with minors in {0, ±1, ±2}, so over every odd prime field
  their entropy vectors are the rational cut-rank vectors: 835 distinct vectors of rank 31, and the same 140 that are
  tight for qubits are tight on (★) (rank 30). With Corollary 3 of the paper, (★) is a facet of Stab^(ℤ_d)_5 for every d
  and of Stab₅.
  The same holds for 119 of the 747 single-CI classes (`checks/check_rev_stab5.py`, which also re-verifies the
  realisations of R32 §4.9 for every d without Smith forms). Secret sharing (R32 §4.17): the CSS state of a linear
  scheme is a stabilizer scheme (Lemma 22), and schemes with ratio 2 exist for Γ_F (every odd d) and Γ_NF (qubits)
  (Proposition 23; `checks/check_rev_qss.py`, part (f)). Addendum: Γ_F has no ideal stabilizer scheme when d has an
  odd prime factor, Γ_NF none when d ≡ 2 (mod 4) (Proposition 24); the Steane state lies within 1/12 bit of the cones of
  odd local dimension in every entropy (R32 §4.16; `checks/check_rev_qss.py`, part (g)). Addendum 2: Γ_F has stabilizer schemes with information
  ratio 7/5 for every odd d (Proposition 25; part (h)), so its optimal ratio for odd squarefree d lies in [4/3, 7/5];
  the logical zero states of the quantum Hamming codes of length 2^m − 1 ≥ 15 lie outside the cone of every odd local
  dimension, with the seven qubits of a Fano subplane as parties (R32 §4.18, Proposition 26; `checks/check_rev_sec6.py`).
  Addendum 3: Γ_NF has stabilizer schemes with information ratio 7/5 for qubits and Galois qudits of characteristic 2
  (R32 §4.19, Lemma 27 and Proposition 28; a GF(4) code that differs from the cube code in two entries; part (i)), so
  its optimal ratio lies in [4/3, 7/5]; for Γ_F no (λ, ω)-decomposition into ideal linear schemes over a field of odd
  characteristic does better than 7/5 (Proposition 29; part (j)).
  Addendum 5: for Z_d qudits, Γ_F has an ideal stabilizer scheme if and only if d = 2, and Γ_NF if and only if d is odd
  (R32 §4.20, Proposition 30, superseding Proposition 24; a forcing lemma certified by exact linear-programming
  certificates, `data/s39_ideal_certificates.json` and `scripts/s39_ideal_schemes.py`, and a normal form over Z_{2^a};
  `checks/check_rev_qss.py`, part (k)).
  Addendum 7: every stabilizer scheme for Γ_F of odd local dimension, not only squarefree, has line sums
  Σ_{P∈L} S(P) ≥ 4 S(R), so information ratio ≥ 4/3 and average share ≥ 4/3 of the secret (R32 §4.21, Lemma 31a and
  Proposition 31: the kernel argument of Jafari and Khazaei holds in every abelian group of odd order; approximate
  schemes lose at most (11/2)ε; `checks/check_rev_qss.py`, part (m)). This supersedes the bound 43/41 of R47 for
  Z_9, Z_27, …; its certificates stay in `data/`. Counted in qudits, the bound holds for Z_d whenever d is not a power
  of 2, that for Γ_NF whenever d ≡ 2 (mod 4), both for secrets of whole qudits, and for d = 2r, r odd, Γ_NF and (r > 1)
  Γ_F have schemes with ratio 7/5 (products over prime powers; part (n)).
- [x] **R47 (2026-10-20): stabilizer secret sharing — the information ratio depends on the local dimension; distances
  of the witnesses from the cones; checks for the revision** (`docs/R32_reduction_and_local_dimension.md` §4.15,
  §4.16; Lemma 20 and Proposition 21 proved, the 4/3 bounds from Jafari–Khazaei, the 43/41 bound certified and
  computer-assisted; independent referee (agent); not yet checked by hand).
  - Access structures of the cube states with the origin as the reference system: Γ_F (qubits; Steane's code with a
    reference qubit; minimal authorised sets = the seven Fano lines; the port of AG(3, 2)) and Γ_NF (GF(p), p odd;
    WXY unauthorised, ABCZ authorised; the port of R8). Both cube states are ideal schemes.
  - Lemma 20: in normal form the lift of a stabilizer scheme of prime local dimension p is the entropy function of a
    GF(p)-linear classical scheme for the same access structure, so bounds for linear schemes apply.
  - Proposition 21: contracting a player turns Γ_F into the Fano access structure F and Γ_NF (players A, B, C, Z) into
    the non-Fano access structure N of Jafari and Khazaei (J. Cryptol. 2021), whose bounds give: every stabilizer scheme
    for Γ_F of odd squarefree local dimension (or Galois, odd characteristic), pure or mixed, has
    Σ_{P∈L} S(P) ≥ 4 S(R) on every Fano line and information ratio ≥ 4/3; every qubit (or GF(2^k)) scheme for Γ_NF has
    ratio ≥ 4/3. For every odd local dimension, Lemma 4 (LB) with t = 2 on all seven-element minors gives ≥ 43/41
    (linear programs over self-dual polymatroids with exact rational dual certificates and exact feasible points).
  - §4.16: every vector of Stab^(Z_d)_6, d odd, differs from the Steane vector by ≥ 1/952 bit in some entropy
    (‖f9‖₁ = 2856 on the 63 cut entropies); the non-Fano vector is ≥ 1/129 from the qubit cone; the seven-party cube
    states ≥ 1/45 and 2/105.
  - New `scripts/s38_secret_sharing.py` (gates G43, G43b), `data/s38_qss_certificates.json`; independent checks
    `checks/check_rev_sec5.py`, `checks/check_rev_sec6.py` (modules of the revision: Lemma 4*, Theorems 12, 14,
    Corollary 13) and `checks/check_rev_qss.py`.
- [x] **R46 (2026-10-16): fewer parties for two primes — incomparable from n = 2p + 2 on; the T₈ and T₈′ vectors lie
  in every cone** (`docs/R32_reduction_and_local_dimension.md` §4.13, §4.14, §4.5; Lemma 4\*, Theorem 17,
  Corollary 19 and Proposition 18 proved; certified; independent referee (agent); not yet checked by hand; not in the
  paper).
  - Lemma 4\*: Lemma 4 applied with the derived subspace C = (A₁ + B₁) ∩ … ∩ (A_M + B_M), its C-terms bounded by
    those of the 2M variables (modular law, dimension counts), gives two (2M)-variable inequalities, (LB\*) for
    characteristic not dividing t and (LA\*) for characteristic dividing t, also for abelian groups. For M ≥ 4 each
    fails by 1 on T(t) = [I_M | J_M − I_M] in the other characteristic.
  - Theorem 17: with doubled codes, **Stab⁽ᵖ⁾_n and Stab⁽q⁾_n are incomparable for every odd prime p, every prime
    q ≠ p and every n ≥ 2p + 2** (the paper's Theorem 2′: 2p + 3; p = 2 is Theorems 9 and 12, n ≥ 6). For p = 3 this
    is n = 8 in both directions, without DFZ. Corollary 19: Theorem 5 from 2·max(d, d′) + 2 parties.
  - Proposition 18: both seven-party candidates of R45 fail. 2·S_T8′ is the entropy vector of a qutrit state and of a
    qubit state (two qudits per system: a "twisted" two-copy family over GF(3); one copy over GF(4)); 2·S_T8 that of a
    Z_p state for every prime p ≥ 5 (an integer two-copy family whose only bad primes are 2 and 3), and 6·S_T8 that of
    a qubit state (N1 over GF(8) and AG(3, 2)). So S_T8 and S_T8′ lie in Stab^(Z_d)_7 for every d.
  - DFZ arXiv:1401.2507 (the PDF provided by the user): the misprint of Theorem 3.1 is in the final expansion of its
    proof (46, 3, 7 instead of 50, 7, 3); Theorem 4.1 checked, and its balanced form is ≥ 42 on h_T8′, as
    Proposition 18 requires.
  - Checks in gates G42 (values of the forms, Theorem 17 witnesses from the [G | G] codes, Z_N values, Proposition 18,
    DFZ) and G42b (all 8-tuples of subspaces of GF(2)², GF(3)², GF(5)²; coordinate lines with all lines of GF(2)⁴ and
    GF(3)⁴; random arrangements; M = 8); new `scripts/s37_derived_c.py`, `data/s37_derived_forms.npz`; sandbox logs
    in `results/2026-10-16_r46-sandbox/`.
- [x] **R45 (2026-10-15): six parties, one qudit per system over 13 fields — only characteristic 2 stands out**
  (`docs/R32_reduction_and_local_dimension.md` §4.10, *Other fields*; certified, Propositions 15 and 16 proved;
  two independent referees (agents); not yet checked by hand; not in the paper).
  - The entropy vectors of GF(q)-linear stabilizer states of seven systems with one qudit each, over GF(2), …,
    GF(49) (13 fields; exhaustive up to GF(5) with the new enumerator `cc/s40_enum7q.c`, sampled beyond), fall into
    136 orbits of S₇. Over every tested field with at least seven elements: 134 common orbits plus S_NF (odd
    characteristic) or S_F (characteristic 2).
  - Consequence for two odd primes at n = 6: for any two of 3, 5, 7, 11, 13, no state with one qudit per system
    witnesses Stab⁽ᵖ⁾_6 ⊄ Stab⁽q⁾_6 (certified for p = 3, 5; empirical for p = 7, 11, 13, whose lists are samples).
    A witness needs S(P) ≥ 2 for one of the seven systems (the purifier included).
  - Proposition 15: the U₂,ₙ states need q ≥ n − 1 (why the symmetric vector is missing over GF(5) and the U₂,₆
    orbit over GF(4)). Proposition 16: eight vectors in F⁴, with v₁..v₄ a basis, no three dependent and the ten
    dependencies of T₈′ (T₈ with its circuit-hyperplane {w₁..w₄} relaxed), are [I₄ | J₄ − I₄] up to scaling; its
    matroid is T₈ in characteristic 3, T₈′ in characteristic ≥ 5, AG(3, 2) in characteristic 2. The CSS states of T₈
    and T₈′ are the candidates at seven parties (a heuristic; not settled).
  - Checks in gates G41 (the orbit table: 1,646 certificates re-checked with a separate GF(q), the structure above)
    and G41b (exhaustive GF(2)–GF(5), target searches for Proposition 15, a sampling run, Proposition 16); new
    `scripts/s36_one_qudit_fields.py`, `cc/s39_sample7q.c`, `cc/s40_enum7q.c`, `data/s36_one_qudit_fields.npz`;
    sandbox logs in `results/2026-10-15_r45-sandbox/`.
  - Also: the paper's arXiv number (2610.18905) in this README and in `CITATION.cff`; the CI time caps raised for
    G41b (§9).
- [x] **R44 (2026-10-15): six parties, powers of two — Stab^(Z_d)_n ⊆ Stab^(Z_{2^a})_n exactly when d divides 2^a,
  for every n ≥ 6** (`docs/R32_reduction_and_local_dimension.md` §4.12, Theorem 14; proved; independent referee
  (agent); not yet checked by hand; not in the paper).
  - For every a ≥ 1, f^(a) = f12 + (1 − 1/a)·U_Z ≥ 0 on Stab^(Z_{2^a})_6, with U_Z = S(Z) + S(AYZ) − S(AY) + S(BXZ)
    − S(BX) + S(CWZ) − S(CW); f^(a)(S_NF) = −1/a, and f^(a) = 1 − j/a at the seven-point state over Z_{2^j}. So odd d
    and 2^a are incomparable from n = 6 on (the question R43 left open), and Z_2 ⊊ Z_4 ⊊ Z_8 ⊊ … from n = 6 on.
  - Proof: the symplectic transfer of Theorem 12 over finite abelian groups with an alternating form, with a version
    of Lemma 4 (LA) for groups of exponent 2^a (Lemma 4(G′): coefficient 3 − 1/a instead of 2). R44 also proves the
    entropy formula of Theorem 1(d), step 1, for every d (until now cited).
  - Checks in gates G40 (exact values, the tight case step by step) and G40b (every step of the proof on random Z_4
    and Z_8 states, state vectors, Lemma 4(G′) tests, two exhaustive families of 2^21 codes); new script
    `scripts/s35_two_power_transfer.py` and annealer/enumerator `cc/s38_minform_zn.c`; sandbox logs in
    `results/2026-10-15_r44-sandbox/`.
  - Erratum: the range of the question R43 left open was too large (for d = 5 against Z_4, n ≥ 11 was already
    decided); corrected in §5 of the note.
- [x] **R43 (2026-10-15): from six parties on, the qubit cone differs from the Z_d cone of every d ≥ 3**
  (`docs/R32_reduction_and_local_dimension.md` §4.11, Corollary 13; proved; independent referee (agent); like
  Theorem 12 not yet checked by hand; not in the paper).
  - The *seven-point state* over Z_d (the CSS state of the seven nonzero vectors of {0,1}³; the Steane state for
    d = 2, the non-Fano state for odd primes d) has f12 = −log₂(d / gcd(d, 2)) bits, negative for every d ≥ 3. So
    Stab^(Z_d)_n ⊄ Stab⁽²⁾_n for every d ≥ 3 and n ≥ 6, including d = 4, 8, 16, … (Theorem 12 covered d with an odd
    prime factor). With Theorems 5 and 9: from n = 6 on, the qubit cone is strictly smaller than the Z_d cone of
    every even d ≥ 4 (for Z_4 this needed seven parties before) and incomparable with that of every odd d. Galois
    qudits GF(2^k) still give exactly the qubit cone.
  - Checks in gates G39 (exact values for d = 2..16, Smith invariants of all column sets) and G39b (state vectors
    for d = 4, 6, 8); sandbox logs in `results/2026-10-15_r43-sandbox/`.
  - Erratum: an item of §5 written in R42 called two known directions open (Z_{2^k} against qubits and against odd
    d); corrected there.
- [x] **R42 (2026-10-14): six parties, the other direction, for the closed cones** (roadmap item M3, finished;
  `docs/R32_reduction_and_local_dimension.md` §4.11; proved; independent referee (agent), not yet checked by hand;
  not in the paper).
  - *Theorem 12:* Stab^(Z_d)_6 ⊄ Stab⁽²⁾_6 for every d ≥ 3 that is not a power of 2; in particular
    Stab⁽ᵖ⁾_6 ⊄ Stab⁽²⁾_6 for every odd prime p. With Theorem 9 (R36), **at six parties the qubit cone and the cone
    of every odd local dimension are incomparable** (known before from n = 7, Theorem 2). This answers the question
    of the paper's Discussion.
  - Witness: the non-Fano state. Certificate: a new inequality f12 ≥ 0 (31 terms on the six-party cuts), valid for
    every qubit stabilizer state, with f12(S_NF) = −1; it is 0 at the Steane state.
  - Proof: a quantitative ("symplectic") version of Proposition 10. Isotropy of the stabilizer group bounds the
    dimensions of the restrictions of the line pieces (Lemmas T1–T4), and Lemma 4 (LA)₂ is applied to the subspace
    arrangement dual to their sum (Lemma T5). f12∘π is a self-dual characteristic-2 linear rank inequality on seven
    variables that is not blind on the self-dual slice.
  - Checks: gate G39 (f12 pinned in the manifest; f12 ≥ 0 on all 10,773 vectors of the one-qubit-per-system list,
    and on the qutrit list f12 < 0 only at S_NF) and G39b (every step of the proof on 40 states); an adversarial
    annealing search (`cc/s37_minform7.c`, `--transfer-anneal`) found no negative value [empirical]. The referee
    recomputed f12, checked every step and ran its own tests (26,549 states with the full chain of the proof, all
    2^21 seven-vertex qubit graph states, about 6 M random states, about 78 M annealing proposals).
  - Sandbox runs in `results/2026-10-14_r42-sandbox/`.
- [x] **R41 (2026-10-14): six parties, the other direction, without closure** (roadmap item M3;
  `docs/R32_reduction_and_local_dimension.md` §4.10; proved; independent referee (agent), not yet checked by hand;
  not in the paper).
  - *Proposition 10, strengthened:* no qubit stabilizer state has an entropy vector that is a positive multiple of
    S_NF, the entropy vector of the non-Fano state, with any number of qubits per system. New step: isotropy at the two
    other points of each line forces the line pieces at the diagonal points W, X, Y to coincide, so the non-Fano
    matroid would get a multilinear representation over GF(2), which does not exist.
  - *Corollary 11:* S_NF is not a finite nonnegative combination of qubit stabilizer entropy vectors (uses the
    certified face of §4.8). So, for every odd p, the cone generated by six-party GF(p) stabilizer entropy vectors is
    not contained in the one generated by the qubit ones. For the closed cones of the paper (the Discussion's
    question) it stays open: S_NF could be a limit of such combinations. (R42: settled, Theorem 12; see R42 above.)
  - **Erratum (R40).** R40's Proposition 10(d) said that in any qubit realisation of m·S_NF the line pieces at one of
    W, X, Y differ; there is no such realisation, so (d) was true only vacuously, and the R40 text built on it (also
    the targeted m = 2 search proposed in §7 ⓪″) is corrected.
  - New check `line_pieces_check` (Proposition 10, steps 1 and 1′, on the Fano and non-Fano states with one and two
    qudits per system) in gate G38b; sandbox log in `results/2026-10-14_r41-sandbox/`.
- [x] **R40 (2026-10-14): all 641 known five-party extreme rays over every local dimension; six parties, the other
  direction** (roadmap items S5(a), finished, and M3, steps 1–2; `docs/R32_reduction_and_local_dimension.md` §4.9,
  §4.10; independent referee (agent), not yet checked by hand).
  - **All 641 known extreme-ray orbits of Stab₅ are extreme rays of Stab^(Z_d)_5 for every d ≥ 2** [certified]: 566
    by unimodular realisations and 75 by covering certificates. R40 covered the last five rows (587, 589, 612, 630,
    635) with longer rational searches from the near misses, then searches right mod the uncovered primes. So no
    known extreme ray can witness a dimension dependence at n = 5. Every row now also has an explicit qutrit
    realisation in the data.
  - **Six parties, the other direction** (is Stab⁽ᵖ⁾_6 ⊄ Stab⁽²⁾_6 for odd p?):
    - Exhaustive survey of the states of seven systems with at most one qudit per system [certified]: 10,773 qubit and
      27,334 qutrit entropy vectors; the qubit-only ones are the 30 relabellings of the Steane state; the 16,591
      qutrit-only ones form 53 orbits, and 52 of them are qubit vectors at multiplicity 2 (explicit 14-qubit graph
      states). The remaining orbit is that of the non-Fano state h_NF, the only possible witness of this kind for
      p = 3.
    - *Proposition 10 [proved]:* no GF(2^k)-linear state has entropy vector k·S_NF (the non-Fano vector in units of
      log 2^k), while every odd q has one. ~~In any qubit realisation of m·S_NF with m qubits per system, the two line
      pieces at one of the three diagonal points restrict to different, non-Lagrangian subspaces.~~ (R41: true only
      vacuously; there is no such realisation, see R41 above.) Checked by exhaustive search over GF(2), …, GF(9)
      and, for S_NF, GF(16) (`cc/s36_target7.c`, and in Python for q ≤ 4).
    - Annealing found no qubit state with 2·S_NF or 3·S_NF [empirical]; every run stopped at the residual of m
      Steane copies.
  - New: `cc/s35_enum7.c`, `cc/s34_anneal7.c`, `cc/s36_target7.c`, `scripts/s34_six_party_survey.py`,
    `data/s34_seven_qudit_vectors.npz` (gates G38, G38b); `data/s32_stab5_cover.npz` and the manifest updated;
    sandbox runs in `results/2026-10-14_r40-sandbox/`.
- [x] **R39 (2026-10-14): 636 of the 641 known five-party extreme rays over every local dimension** (roadmap
  item S5(a), continued; `docs/R32_reduction_and_local_dimension.md` §4.9; certified; independent referee (agent), not
  yet checked by hand).
  - **636 of the 641 known extreme-ray orbits of Stab₅ are extreme rays of Stab^(Z_d)_5 for every d ≥ 2** (577 in
    R38): 566 by unimodular realisations (κ = 1: 515, 2: 50, 3: 1) and 70 by covering certificates.
  - What closed 59 of the 64 rows R38 left open: the covering search now starts its searches "rational and right
    mod q" from a rational realisation of the row and asks for all uncovered primes at once; longer rational searches
    found the first rational realisation of 35 more rows (row 200 at κ = 2). Unimodular realisations of two
    rows (539, 631) were already in the output of R38's strategy probe, which had not been merged.
  - 5 rows remained open (587, 589, 612, 630, 635); R39 knew no rational realisation for any of them (best attempts
    in `results/2026-10-14_r39-sandbox/open_rows.json`). R40 (above) covered them.
  - Errata (R38): a count in R38's strategy probe file; the header and §5 of the R32 note counted 564 where 577 was
    meant (`CHANGELOG.md`, R39).
  - `scripts/s32_unimodular.py --cover` (new search order; `--rounds`, `--rings`, `--kappa2`); data files and manifest
    updated; sandbox runs in `results/2026-10-14_r39-sandbox/`.
- [x] **R38 (2026-10-13): the known five-party extreme rays over every local dimension** (roadmap item S5(a);
  `docs/R32_reduction_and_local_dimension.md` §4.9; certified; independent referee (agent), not yet checked by hand).
  Counts as of R38; R39 and R40 (above) update them.
  - **577 of the 641 known extreme-ray orbits of Stab₅ are extreme rays of Stab^(Z_d)_5 for every d ≥ 2.**
    - 564 have a *unimodular realisation*: an integer graph state W whose 31 cut blocks have Smith invariants 0 and
      1 only and rational rank κ·r(X) (κ = 1 for 513, 2 for 50, 3 for 1); W mod d realises κ·r for every d.
    - 13 more have *covering certificates*: rational realisations whose bad primes (those dividing a Smith
      invariant of a cut block) intersect inside the primes of stored GF(2) / GF(3) realisations, so every prime p
      is covered, and every d by Corollary 3 of the paper.
  - The 12 BCHS ledger rays without an explicit realisation have one now (κ = 2): all 641 known extreme rays have
    explicit qubit realisations, and 577 have explicit qutrit realisations (96 before).
  - 64 rows remain open (200, 488, 492, 493, 495, 496, 497, 501, 502, 503, 504, 506, 508, 509, 525, 539, 549, 556, 557, 559, 560, 564, 569, 572, 573, 574, 575, 576, 578, 579, 580, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 593, 594, 595, 596, 597, 598, 605, 606, 607, 609, 610, 611, 612, 613, 630, 631, 633, 634, 635, 636, 637, 638, 640). An odd-characteristic common-information screen
    (`scripts/s33_ci_oddchar.py`, exploratory) excludes none of them.
  - Data `data/s32_stab5_unimodular.npz` (gates G36/G36b) and `data/s32_stab5_cover.npz` (gate G37);
    `cc/s32_unimodular.c`, `scripts/s32_unimodular.py`; sandbox runs in `results/2026-10-13_r38-sandbox/`.
- [x] **R36 (2026-10-02): six parties already separate qubits from every odd local dimension**
  (`docs/R32_reduction_and_local_dimension.md` §4.8).
  - *Theorem 9 (proved; independent referee):* Stab⁽²⁾_n ⊄ Stab^(Z_d)_n for every odd d and every n ≥ 6. So the
    smallest number of parties with a dimension dependence is 5 or 6 (it was 5, 6 or 7).
    - Witness: the seven-qubit simplex-code (Fano) state, i.e. the logical |0⟩ of Steane's [[7,1,3]] code, one
      qubit per party; h_S = r_F + r_F*.
    - Certificate: a new seven-variable inequality f9 (67 terms, `data/s31_six_party.npz`), valid in every odd
      characteristic and for abelian groups of odd order, with f9(h_S) = −3. It comes from a *transfer*: Lemma 4
      (LB)_2 is applied to the derived family F̂_i = Σ_L U_i ∩ (U_j + U_k) (sums over the Fano lines through i), and
      every term is bounded by a linear form in h.
    - f9 + f9∘D is the self-dual characteristic-dependent inequality that Proposition 8 asked for; f9 is not a
      combination of Shannon and the eight known forms in any relabelling.
    - Invariant form: on Fano-invariant self-dual vectors f9 ≥ 0 reads μ ≤ 210β + 30α + 21ν, which h_F (μ = 1,
      α = β = ν = 0) violates.
    - Open: the reverse direction at n = 6 (is the non-Fano state outside the qubit cone?). Exact reduction: this
      asks whether LR⁽²⁾₇ meets the two-dimensional Shannon face cone(r_NF, r_NF*) outside 0. (R42: it is outside;
      Theorem 12.)
  - *Proposition 8 extended:* the seven-variable inequalities of Blasiak–Kleinberg–Lubetzky (arXiv:1108.2489v1,
    eqs. (6.8), (6.17); transcribed and checked against the PDF, `data/bkl_char_dependent.csv`) are blind at n = 6
    as well (exact Shannon certificates; eight forms in total).
  - *Misprint:* BKL's eq. (6.17) is false over every field as printed (A, B, C, Z = e₁, e₂, e₃, e₁ + e₂ + e₃,
    W = X = Y = 0 gives −3). The proof drops three terms; the corrected form is used (pitfall 21).
  - Gates G34/G34b (BKL) and G35/G35b (Theorem 9); `scripts/s31_six_party.py`, `s30_dimension_profile.py --bkl`.
- [x] **R35 (2026-10-02): how far down the dimension dependence goes** (`docs/R32_reduction_and_local_dimension.md`
  §4.7).
  - *Theorem 7 (proved):* for n ≤ 4 the stabilizer cone does not depend on the local dimension at all:
    Stab^(Z_d)_4 = h⁻¹(Shannon + Ingleton) for every d ≥ 2.
    - Each of the 46 extreme rays has an integer graph state whose cut blocks have Smith invariants 0 and 1 only,
      so it works over every GF(p).
    - Data: `data/s30_stab4_unimodular.npz` (the qubit certificates of R32 with signs on 5 of them).
  - *Proposition 6 (proved):* the one-cut observation of R33, for every t. Over GF(ℓ), ℓ ∤ t − 1, the generalised
    cube matroid differs from the rational one exactly on the balanced sets with ℓ | t − a or ℓ | a − 1. For
    ℓ = t prime this is a single cut.
  - *Proposition 8 (proved, exact certificates):* at n = 6, the six seven-variable characteristic-dependent
    inequalities considered (DFZ (65), (91), Lemma 4 t = 2, Peña–Sarria t = 2) are blind in every relabelling.
    - For each f, f + f∘D (D = Kaced's duality) is a nonnegative combination of Shannon inequalities. So on the
      self-dual slice f follows from strong subadditivity and holds for every pure state.
    - BKL's seven-variable inequalities (added in R36) and substitution instances were not covered.
    - Data: `data/s30_selfdual_shannon_certs.npz`.
  - So the dependence is absent for n ≤ 4 and was open for n = 5, 6 (R36 settles one direction at n = 6, R42 the
    other). It is
    present pairwise: qubits against odd p from n = 7, primes p < q from 2p + 3, and Stab^(Z_d) ⊄ Stab^(Z_d′)
    (d ∤ d′) from 2·max(d, d′) + 3.
  - Pilot, not certified: 40 of the first 60 known extreme rays of Stab₅ lift to all fields (superseded by R38).
  - Checked by an independent referee (all correct; scope and text fixes applied).
  - Gates G33/G33b; `scripts/s30_dimension_profile.py`.
- [x] **R34 (2026-10-02): Z_d qudits -- the stabilizer entropy cone remembers the local dimension**
  (`docs/R32_reduction_and_local_dimension.md` §4.6).
  - *Theorem 5 (proved):* for d, d′ ≥ 2 and n ≥ 2·max(d, d′) + 3, Stab^(Z_d)_n ⊆ Stab^(Z_d′)_n exactly when d
    divides d′ (entropies in bits). For d | d′ the inclusion holds for every n.
  - Examples: qubits ⊊ Z_4 from 7 parties, Z_4 ⊊ Z_8 from 11. Galois qudits GF(4) give exactly the qubit cone,
    so the two standard theories of local dimension 4 differ.
  - Tools: Theorem 1(d), the reduction of R32 for every d (subgroups of abelian groups of exponent dividing d in
    place of subspaces), and Lemma 4(G), Lemma 4 for abelian groups ((LB) by a counting argument).
  - Witnesses: the generalised cube states of R33 over Z_N, N = p^{j+1} (`scripts/s29_zd_stabilizer.py`,
    `data/s29_zd_witnesses.npz`, gates G32/G32b).
  - The Z_4 cube state (8 qudits) lies outside every prime cone Stab⁽ᵖ⁾_7. Its entropy vector in bits is the
    qubit cube vector plus the qutrit cube vector, so it is inside their convex hull.
  - Whether Z_d states can ever leave the convex hull of the prime cones is open. Up to one party it is the open
    question Γ^MixL = Γ^Abl of Khazaei, arXiv:2608.09543.
  - The A2 exclusion certificates use only Shannon, Ingleton and one common information, so they hold for every
    local dimension.
  - Checked by an independent referee (all items correct, minor fixes applied).
- [x] **R33 (2026-10-02): the stabilizer cones of any two distinct primes are incomparable.** For primes p < q
  and every n ≥ 2p + 3, Stab⁽ᵖ⁾_n and Stab⁽q⁾_n are incomparable (Theorem 2′,
  `docs/R32_reduction_and_local_dimension.md` §4.4).
  - The proof is self-contained through Lemma 4: characteristic-dependent inequalities for the configuration
    J − I, proved in the note and checked by an independent referee.
  - Peña–Sarria arXiv:1905.00003v3, Example 6 (transcribed from the PDF) serves as a cross-check.
  - Witnesses: the "generalised cube states" on the 2p + 4 points e_i, c − e_i, c, 0 of F^{p+1}, with the code of
    affine functions over GF(p) resp. GF(q) (`scripts/s28_generalized_cubes.py`, `data/s28_gencube_t3.npz`,
    `data/s28_gencube_t5.npz`, gates G31/G31b).
  - They are self-dual by explicit weights and violate the pulled-back Lemma 4 (LB) resp. (LA) by −2.
  - In all computed cases (t = 2, 3, 5) the two vectors differ on a single cut.
  - For p = 2 they are the cube states of R32.
  - Side result: Theorem 3.1 of DFZ arXiv:1401.2507v1 is misprinted. As printed, it is violated over every field
    by C = X = a line and all other subspaces zero (pitfall 20); it is not used.
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
  45,292 catalogue orbits (reproduced by `checks/check_dfz6.py`, which carries the
  transcription of DFZ eqs. (44)–(63)). Write-up with full proofs:
  `docs/A2_stabilizer_vs_QLR.md` (the proofs used in the paper were checked by hand
  before submission, R37).
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

- Python ≥ 3.8, `pip install -r requirements.txt` (numpy, scipy; sympy is needed
  only by the scripts in `checks/`).
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
python3 selftest.py --full            # 49 gate lines + ALL PASS, ~8 min; must be all green first
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
| `s32-unimodular` | S / X | done in the sandbox, R38–R40 (2026-10-13/14, 2 cores: all 641 rows certified); no open row is left, so the X variant finishes at once; kept for rows added to the catalogue later (the R39 covering search with rounds of 1,200 s) | rows 587, 589: one covering round at 3 s (~25 s) |

**Rule:** every tier-X job is handed over only after its smoke variant has
passed in the sandbox (and, optionally, in the manual CI workflow
`.github/workflows/smoke.yml`). The analysis side then pulls `main`, reads
`results/`, and refreshes the catalogue assets in `data/` at milestones.

## 3. Pipeline stages

| Stage | Command essentials | Purpose / expectation |
| --- | --- | --- |
| selftest | `selftest.py [--full] [--workers N]` | 59 regression gate lines (G0–G43 with G2b and G10b, G29b, G31b–G36b, G38b–G43b) + ALL PASS; run first on any machine |
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
| s28 | `s28_generalized_cubes.py --t T --fields p q … [--random N] [--exhaustive01] [--exhaustive2] [--out …]` ; `--dfz1401` | R33: Theorem 2′ -- generalised cube states over GF(p) and GF(q), explicit self-duality weights, contraction to the configuration J − I, pulled-back Lemma 4 (LA)/(LB) with value −2 (Peña–Sarria Example 6 as cross-check); sanity tests; `--dfz1401` reproduces the misprint of pitfall 20 |
| s29 | `s29_zd_stabilizer.py [--cases t:N …] [--state-vector] [--lemma-tests K] [--out …]` | R34: Theorem 5 -- generalised cube states over Z_N (N = p^k; exact subgroup arithmetic via Smith forms): freeness, unit self-duality weights, h_S = 2r, contraction to J − I over Z_N, pulled-back Lemma 4(G) values; `--state-vector`: the Z_4 cube state against its 4^8 amplitudes and the identity Z_4 cube = qubit cube + qutrit cube; `--lemma-tests`: Lemma 4(G) on point configurations of Z_4^k and Z_9^3 |
| s30 | `s30_dimension_profile.py [--n4] [--one-cut T …] [--six] [--bkl] [--probe6] [--n5-pilot K]` ; `--n4-search --out …` ; `--solve --out-certs …` | R35: Theorem 7 -- the 46 four-party extreme rays over every field (unimodular integer graph states); Proposition 6 -- one-cut rank changes over GF(ℓ) vs Q; Proposition 8 -- exact Shannon certificates for the self-dual symmetrisations of the eight seven-variable characteristic-dependent inequalities considered (BKL added in R36); `--bkl` BKL transcription checks and the (6.17) misprint (R36); `--probe6` encoding check; `--n5-pilot` five-party lifting pilot (not certified) |
| s31 | `s31_six_party.py [--derive [--out …]] [--witness] [--invariant] [--face] [--lemmas N] [--stress N] [--exhaustive]` | R36: Theorem 9 -- the transfer form f9 recomputed exactly from Lemmas 0, A–E and Lemma 4 (LB)_2; the simplex-code state vector (h_S = r_F + r_F*, f9 = −3); Fano-invariant form; exact Shannon faces of h_F and h_NF; lemma-by-lemma checks with actual subspaces over GF(2), GF(3), GF(5); random and exhaustive validity tests |
| s32 | `s32_unimodular.py --check [--full]` ; `--search --seconds S [--resume] [--rows …] [--strategies …] --out …` ; `--cover --out … --cover-out … [--rows …] [--seconds S] [--rounds R] [--rings …] [--kappa2]` ; `--check-cover` ; `--self-test` | R38–R40: unimodular realisations and covering certificates of the known extreme rays of Stab₅ (kernel `cc/s32_unimodular.c`: annealing over {-1, 0, 1}, neutral moves, restriction of scalars from quadratic orders); exact certification by diagonal forms over Z, rational rank + minors with gcd 1, GF(p) ranks and Z_d state vectors; covering certificates by recomputed bad primes |
| s33 | `s33_ci_oddchar.py --rows … [--pairs disjoint/all] [--out …]` | R38: exploratory screen -- common information + the seven-variable forms valid in odd characteristic (DFZ (65), (LB)_2, BKL (6.17) corrected) on the extension; floating-point verdicts |
| s34 | `s34_six_party_survey.py --check [--full]` ; `--self-test` ; `--enumerate --out …` ; `--survey --lists … --out … --seconds S` ; `--nonfano --lams 2,3 --seconds S` ; `--fields 2,3,4,5,7,8,9 [--controls K]` ; `--fields-python 2,3,4` ; `--transfer [--full --trials N]` ; `--transfer-anneal "n0,…,n6;…" --seconds S [--pull-nf μ]` | R40: six parties, the other direction -- every entropy vector of seven-system qubit and qutrit stabilizer states with at most one qudit per system (`cc/s35_enum7.c`, weighted graphs up to scaling), the qutrit-only orbits and their qubit realisations at multiplicity 2 (`cc/s34_anneal7.c`), and graphs over GF(q) with a prescribed cut-rank vector (`cc/s36_target7.c`; Proposition 10); exact checks over GF(2), GF(3) and GF(q). R42: Theorem 12 -- the inequality f12, every step of its proof on given states (gates G39/G39b), and the adversarial search `cc/s37_minform7.c` |
| s35 | `s35_two_power_transfer.py --check [--full --trials N]` ; `--anneal "a:rows:n0,…,n6;…" --seconds S [--pull μ [--singles]] [--weakened c]` | R44: Theorem 14 -- the inequality f^(a) = f12 + (1 − 1/a)·U_Z on Z_{2^a} stabilizer states: exact linear algebra over Z_{2^a}, stabilizer states as Lagrangian subgroups (graph, CSS of non-free codes, non-free Lagrangians after random symplectic maps), every step of the proof on given states without compressing them (gates G40/G40b), Lemma 4(G′) tests; annealing over CSS codes and exhaustive families of codes B + 2^s·P (`cc/s38_minform_zn.c`) |
| s36 | `s36_one_qudit_fields.py --check [--full]` ; `--enumerate q [--all] [--slice s --nslices n] [--out F]` ; `--sample q --samples N --seed S --out P` ; `--build P --out F` ; `--prop16` | R45: six parties, one qudit per system over 13 fields -- the S₇-orbit table of the entropy vectors of GF(q)-linear states with a realising graph per field and orbit (gates G41/G41b), exhaustive lists with `cc/s40_enum7q.c` (scaling and relabelling reductions; GF(2)–GF(5)), the sampler `cc/s39_sample7q.c`, target searches for Proposition 15, and Proposition 16 (minors, brute-force representability of T₈, T₈′, AG(3,2), S_T8 against S_T8′) |
| s37 | `s37_derived_c.py --check [--full]` ; `--forms` ; `--write-data` | R46: Lemma 4\* (Lemma 4 with the derived subspace C = ∩(A_m + B_m)): the forms (LB\*), (LA\*) for M = 4, 6, their values on T(t) over eight fields and over Z_N, the Theorem 17 witnesses (doubled codes, entropies from the [G \| G] code), validity tests (exhaustive over GF(q)², coordinate lines, random); Proposition 18 (the T₈ and T₈′ vectors in every cone: the two-copy families over GF(3) and over ℚ, and the other constructions); DFZ arXiv:1401.2507 Theorems 3.1 (re-expanded) and 4.1 (gates G42/G42b) |
| s38 | `s38_secret_sharing.py --verify` ; `--solve` ; `--certify --out F` | R47: stabilizer quantum secret sharing (§4.15 of the R32 note) -- the access structures Γ_F and Γ_NF of the cube states and their contractions (the Fano and non-Fano access structures of Jafari–Khazaei), the linear programs over self-dual polymatroids with Lemma 4 (LB), t = 2, on all seven-element minors (pure, and mixed with an environment), their exact rational dual certificates (`--verify`, gate G43; `--certify` regenerates them from HiGHS duals with SymPy), the re-solved optima with exact feasible points, two controls and the numerical Jafari–Khazaei programs (`--solve`, gate G43b), and the distances of the witnesses from the cones (§4.16) |
| s39 | `s39_ideal_schemes.py --verify` ; `--certify --out F` | R48 addendum 5: ideal stabilizer schemes for Γ_F and Γ_NF (§4.20 of the R32 note, Proposition 30) -- the forcing lemma (a polymatroid on R and the seven players with g(R) = g(P) = 1 and the access equalities is the rank function of AG(3,2), resp. R8), certified by exact rational multipliers for the minimum and the maximum of each of the 255 values (`--verify`, about a second; `--certify` regenerates them from HiGHS duals with an exact repair, about 20 s), and the orders of the projections of the normal form over Z_{2^a}, a ≤ 3 |
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
20. **A misprinted theorem (R33).** Theorem 3.1 of DFZ arXiv:1401.2507v1 (an
    eight-variable inequality claimed for characteristic ≠ 3) is false as
    printed: C = X = one line, all other subspaces zero, gives left side 0 and
    right side −4, over every field. An exhaustive search over all subspace
    arrangements of GF(2)² and GF(3)² finds 125 and 254 violations. The theorem
    is not used; the generalisation to all primes (R33) rests on Lemma 4,
    which is proved in the note (`s28_generalized_cubes.py --dfz1401` reproduces
    the counterexample). Moral: before citing a linear rank inequality, test it on all small
    arrangements (`s28_generalized_cubes.py --exhaustive2` style); a published
    statement can be wrong as printed. (R46: the error sits in the last expansion
    of the proof, p. 25: the coefficients of H(X|A,C,D), H(C|B,X,Y), H(C|A,W,Y)
    should be 46, 3, 7, not 50, 7, 3; `s37_derived_c.py` recomputes them from the
    proof's own bounds. A second moral: re-expanding a proof's last display is
    cheap and found the misprint at once.)
21. **A second misprinted inequality (R36).** Eq. (6.17) of
    Blasiak–Kleinberg–Lubetzky arXiv:1108.2489v1 (Theorem 6.6, characteristic ≠
    2) is false as printed: A, B, C, Z spanned by e₁, e₂, e₃, e₁ + e₂ + e₃ and
    W = X = Y = 0 give −3 over every field. The proof bounds three step-4 terms
    in the text but drops them from (6.16); `data/bkl_char_dependent.csv`
    carries both the printed form and the corrected one, and only the corrected
    one is used (Proposition 8). The same exhaustive point-configuration test
    as in pitfall 20 caught it at once (`s30_dimension_profile.py --bkl`).

## 6. Data inventory (`data/`)

| File | Rows | Meaning |
| --- | :-: | --- |
| `M_LR.npy` | 1,790 | LR instance family extracted from psitip (27 classes; historical) |
| `CLR_H.npy` / `CLR_H_fixed.npy` | 1,875 / **1,905** | classical cone H-representation (fixed = + the 30 instances of the Ingleton(39) orbit) |
| `REF28.npy` + `dfz_ref28.csv` | 28 | reference family (24 DFZ inequalities + 4 Ingleton forms); the CSV carries source tags and names the 11 single-source lines, which were checked by hand against DFZ in R37 |
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
| `s28_gencube_t3.npz` / `s28_gencube_t5.npz` | 7 / 6 | R33: entropy vectors of the generalised cube states (t = 3 over GF(3), GF(5), GF(7); t = 5 over GF(5), GF(7); 511 resp. 8,191 coordinates, purifier = origin), the pulled-back Lemma 4 forms FA (valid for char dividing t) and FB (valid otherwise), and Peña–Sarria's Fa, Fb (cross-check), gate G31 |
| `s29_zd_witnesses.npz` | 4 | R34: entropy vectors (units of log p) of the generalised cube states over Z_4 (t = 2, 127 coordinates), Z_8 (t = 4, 2,047), Z_9 (t = 3, 511) and Z_4 (t = 6, 32,767); purifier = origin; gates G32/G32b |
| `s30_stab4_unimodular.npz` | 46 | R35: for each extreme ray of Stab₄ (rays, λ, sizes, N as in `s26_stab4_certs.npz`), an integer symmetric matrix Wt ≡ W mod 2 whose visible cut blocks have Smith invariants 0 and 1 only (Theorem 7), gate G33 |
| `s30_selfdual_shannon_certs.npz` | 8 | R35, R36: for DFZ (65), (91), (LA)_2, (LB)_2, Peña–Sarria (a), (b) (t = 2), BKL (6.8) and (6.17) corrected, exact certificates den·(f + f∘D) = Σ num_i·e_i over the Shannon elemental inequalities on seven variables (fixed order of `elemental7`; Proposition 8), gate G33 |
| `bkl_char_dependent.csv` | 3 | R36: the seven-variable characteristic-dependent inequalities (6.8) [characteristic 2] and (6.17) [characteristic ≠ 2] of Blasiak–Kleinberg–Lubetzky arXiv:1108.2489v1, transcribed and checked against the PDF pages (recorded in the header); (6.17) as printed (`bkl617_odd`, false) and corrected (`bkl617_odd_corrected`); gates G34/G34b |
| `s31_six_party.npz` | 1 | R36: the transfer form f9 (127 coefficients, mask order over A, B, C, W, X, Y, Z; f9 · h ≥ 0 in odd characteristic) and h_F = r_F + r_F*, the h-vector of the seven-qubit simplex-code state (Theorem 9), gates G35/G35b |
| `s32_stab5_unimodular.npz` | 566 | R38, R39: unimodular realisations of known extreme rays of Stab₅: row (in `stab5_extreme_reps.npy`), rep, kappa, sizes, N, W (int8, padded to 64 × 64), method, detail; every cut block of W has Smith invariants 0 and 1 and rational rank kappa · r(X), so W mod d realises kappa · r for every d (§4.9 of the R32 note), gates G36/G36b |
| `s32_stab5_cover.npz` | 125 + 74 | R38–R40: covering certificates for the other 75 rows: pieces (row, kind Q = rational realisation with its bad primes, or GF = realisation over GF(prime), kappa, sizes, N, W, bad, detail); for each row the bad-prime sets of the Q pieces intersect inside the primes of its GF pieces, so the ray lies in Stab^(Z_d)_5 for every d (§4.9), gate G37 |
| `s34_seven_qudit_vectors.npz` | 10,773 + 27,334 | R40: every entropy vector (63 cuts) of seven-system stabilizer states with at most one qudit per system over GF(2) (`gf2`) and GF(3) (`gf3`); one representative of each of the 53 S₇-orbits of qutrit-only vectors (`orbit_vector`, `orbit_size`) with a GF(3) graph realising it (`orbit_weights`), and a 14-qubit graph state realising twice the vector (`qubit2`, `has_qubit2`: all but the non-Fano orbit, `nonfano`), §4.10 of the R32 note, gates G38/G38b; the manifest entry also records the counts of `cc/s36_target7.c` for S_F and S_NF over GF(2), …, GF(9) |
| `s36_one_qudit_fields.npz` | 136 orbits | R45: the S₇ orbits (36,430 vectors) of the entropy vectors of GF(q)-linear stabilizer states of seven systems with one qudit each, over GF(2), 3, 4, 5, 7, 8, 9, 11, 13, 16, 25, 27, 49 (`reps`, `sizes`, `fields`), for every field and orbit met over it a graph realising the representative (`cert`, 255 = none; edge order (0,1), (0,2), …, (5,6)), and the Steane and non-Fano orbits (`steane`, `nonfano`); §4.10 of the R32 note, gates G41/G41b; the manifest entry also records the counts of the exhaustive runs of `cc/s40_enum7q.c` |
| `s37_derived_forms.npz` | 7 forms | R46: (LB\*) and (LA\*) of Lemma 4\* for M = 4 (`LB4`, `LA4`; variables A, B, C, D, W, X, Y, Z = bits 0..7) and M = 6 (`LB6`, `LA6`; A₁..A₆, B₁..B₆ = bits 0..11), DFZ arXiv:1401.2507 Theorem 4.1 (`DFZ41`) and Theorem 3.1 as printed and re-expanded from its proof (`DFZ31_printed`, `DFZ31_reexpanded`); int64, index = subset mask, f·h ≥ 0; §4.13 and §4.5 of the R32 note, gates G42/G42b |
| `s38_qss_certificates.json` | 2 certificates | R47: rational dual certificates of the two secret-sharing linear programs (`F_pure_LB2`, `F_mixed_LB2`; bound 43/41): for each, the bound and lists of [row descriptor, multiplier] for the inequality and equality rows used; §4.15 of the R32 note, gates G43/G43b |
| `s39_ideal_certificates.json` | 1020 certificates | R48 addendum 5: for Γ_F (AG(3,2)) and Γ_NF (R8) and each nonempty set X of the eight systems, a lower and an upper certificate of the forcing lemma, lists of [row descriptor, multiplier] for the elemental inequalities and the equalities used; §4.20 of the R32 note, `checks/check_rev_qss.py` part (k) (no gate) |
| `s27_cube_witnesses.npz` | 4 | R32: entropy vectors of the qubit and odd-p cube states (127 coordinates, parties A,B,C,W,X,Y,Z = vertices 100,…,111, purifier 000) and the pulled-back seven-party inequalities F65 (odd p) and F91 (p = 2), gate G30 |
| `sixvar_template.csv` | — | transcription template for s9 |
| `manifest_shas.json` | — | canonical sha256 of every row family (`sha_rows`; `sha_rows_wide` for wide coordinates) |

Build products (`QLR_H_*.npy`) are not committed; s1 regenerates them in
seconds and checks the sha.

## 6a. Conventions in data and certificate files

- **Subsets (mask order).** A vector on the nonempty subsets of k named elements
  is stored with the subset X at index (mask of X) − 1, where the j-th element
  in the stated order is bit j (least significant bit first). For five parties
  A, …, E are bits 0-4 (31 coordinates); for four parties (`s26_*`, `s30_*`) the
  parties are bits 0-3 (15 coordinates); seven-variable forms and vectors list
  their variables in the order given in the table above (for example A, B, C, W,
  X, Y, Z = bits 0-6 in `s31_six_party.npz`; 127 coordinates).
- **The purifier.** Where a purifier is needed it is the next bit (F = bit 5 for
  five parties, bit 4 for four parties). Entropy vectors are stored on the
  subsets of the parties only; the value on a set that contains the purifier is
  that of its complement (purity), as in `epr1kit.stabcert.entropy6`.
- **Exclusion certificates** (`s17_exclusions.json`, `s18_exclusions.npz`,
  `s22_*.json`). The seven elements are A, …, F = bits 0-5 and the common
  information Z = bit 6. X and Y are masks over A, …, F. A Shannon multiplier
  with index k (or a typed term `["E", k]`) refers to the k-th row of
  `epr1kit.stabcert.elemental(7)`: first the seven rows h(N) − h(N∖i) ≥ 0,
  i = 0, …, 6, then, for each pair i < j in lexicographic order, the rows
  I(i; j | K) ≥ 0 with K running over the subsets of the other five elements in
  decreasing numerical order, from all five down to the empty set. Typed terms
  `["I", a, b, c, d, K]` are the conditional Ingleton rows of
  `epr1kit.stabcert.ingleton_row`. The stored S-form F′ has 31 integer
  coefficients in mask order.
- **Self-dual Shannon certificates** (`s30_selfdual_shannon_certs.npz`) index the
  rows of `elemental7` in `scripts/s30_dimension_profile.py`. The order is the
  same except that K runs in increasing numerical order, from the empty set up.
- **Realisation certificates** (`s16_*`, `s21_*`, `s24_*`, `s26_*`, `s30_stab4_unimodular.npz`,
  `s32_stab5_unimodular.npz`, `s32_stab5_cover.npz`).
  The qudits of each party form one consecutive block, in the order of the
  parties, with the purifier last; `sizes` gives the block sizes, and λ is the
  multiple of the ray that the state realises. In `s16_realisations.npz`, row i
  of the adjacency matrix is a uint64 whose bit j is the entry (i, j). In
  `s24_qubit_l2_certs.npz` the GF(4) matrices encode 0, 1, a, 1 + a as 0, 1, 2, 3,
  with a² = a + 1. In `s32_stab5_unimodular.npz`, kappa plays the role of λ, `row`
  indexes `stab5_extreme_reps.npy`, W is an integer matrix (int8) with zero
  entries inside each party, and `method` names the search strategy of
  `scripts/s32_unimodular.py` that found it.

## 7. What next, and how

**⓪ A2 after R30, R32.** (a) Done for the paper (R37): the proofs it uses were
read line by line, the transcriptions it relies on were checked against the
sources, and a prior-work search was made; the prior work found is credited in
the paper and listed in `docs/R32_reduction_and_local_dimension.md` §6. Still
optional: write to the BCHS authors; tell the authors of DFZ arXiv:1401.2507,
whose IEEE version still prints the false Theorem 3.1 (pitfall 20); check the
published version of BKL's eq. (6.17) (pitfall 21; Peña-Macias and Sarria noted
the error in 2021); read the "technical assumption on the local dimension" of
Majenz's Cor. 5.12 (arXiv:1810.12845) by hand. (b) The 1,155 undecided orbits
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

**⓪′ (R38–R40) The known extreme rays of Stab₅ over every local dimension: done.** All
641 rows are certified for every d: 566 by unimodular realisations
(`data/s32_stab5_unimodular.npz`) and 75 by covering certificates
(`data/s32_stab5_cover.npz`). A five-party dimension dependence would now need an extreme ray of Stab₅ that is not
yet known (or, for d divisible by a square, a vector of Stab^(Z_d)_5 outside Stab₅). The server job
`s32-unimodular` has nothing left to do.

**⓪″ (R40–R42) Six parties, the other direction (roadmap item M3): done.** Among the states with at most one qudit
per system, only the non-Fano state can witness Stab⁽³⁾_6 ⊄ Stab⁽²⁾_6 (R40, §4.10 of the R32 note); no qubit
stabilizer state has an entropy vector that is a positive multiple of S_NF, and S_NF is not a finite nonnegative
combination of qubit entropy vectors (R41, Proposition 10 and Corollary 11); S_NF lies outside the closed qubit cone,
so Stab^(Z_d)_6 ⊄ Stab⁽²⁾_6 for every d that is not a power of 2 (R42, Theorem 12, §4.11; a quantitative
"symplectic" version of the proof of Proposition 10), and for every d ≥ 3 by the seven-point state over Z_d (R43,
Corollary 13). Not yet checked by hand; for the paper, a revision could state that from six parties on the qubit cone
differs from the Z_d cone of every d ≥ 3. The survey for p = 5 is no longer needed for this question. R44 settles the
powers of two: for n ≥ 6, Stab^(Z_d)_n ⊆ Stab^(Z_{2^a})_n exactly when d | 2^a (Theorem 14, §4.12). Still open at
n = 6: two odd primes, and Z_d against Z_d′ when d′ has an odd prime factor (for example Z_4 against Z_6). R45 (§4.10,
*Other fields*): for two odd primes among 3, 5, 7, 11, 13 there is no witness among the states with one qudit per
system (certified for p = 3, 5, empirical for 7, 11, 13), and no seven-variable inequality that separates two odd
characteristics is known to us; so this question is not within reach of the present methods. The next candidates are
at seven parties: the CSS states of T₈ (characteristic 3) and T₈′ (characteristic ≥ 5), Proposition 16; a
certificate would be an eight-variable characteristic-dependent inequality not blind on the self-dual slice (DFZ's
T₈-based ones, arXiv:1401.2507, are the starting point). R46: both fail, the T₈ and T₈′ vectors lie in every
Stab^(Z_d)_7 (Proposition 18, §4.14), so no candidate is left at seven parties; at eight parties, 3 and every q ≠ 3
are incomparable (Theorem 17, §4.13), and in general p < q from n = 2p + 2. Making the GF(7) list exhaustive would take about 15 to 30
CPU-hours with `cc/s40_enum7q.c` (815,738,564,501 graphs; not registered as a server job: it would only move p = 7
from empirical to certified).

**⓪‴ (R47) Stabilizer secret sharing.** (R48 addendum 7: 4/3 for Γ_F in every odd local dimension, R32 §4.21.) Done
for the bounds in §4.15 of the R32 note: 4/3 for Γ_F in odd squarefree
local dimension and for Γ_NF with qubits (through Lemma 20 and the results of Jafari and Khazaei), 43/41 for Γ_F in
every odd local dimension (certified linear program). Open: the optimal information ratios of stabilizer schemes for
Γ_F and Γ_NF in the other characteristic (linear schemes for the six-player ports attain 4/3; whether stabilizer
schemes for the seven-player structures do is not known); Z_{2^a} schemes for Γ_NF (Lemma 4(G′) is the natural
tool); the access structures of the generalized cube states of Theorem 2′ for other pairs of primes.

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

**③ Human checks.** (a) Done (R37): the 11 single-source lines of
`data/dfz_ref28.csv` (`dfz10, 11, 13, 14, 15, 16, 19, 20, 21, 22, 23`) and the
Ingleton form (39) were checked by hand against DFZ eqs. (10)–(23) and (39). (b) The
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
  request (all 57 gate lines must pass). Every step has a time cap (job 60 min, `selftest.py --full` 40 min;
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
- **Citing:** `CITATION.cff` holds the metadata. Tagged releases are archived at
  Zenodo; R37, the first, is the version cited by the paper, and its DOI is
  given in the release notes on GitHub.
- **Licences:** code MIT (`LICENSE`), data CC BY 4.0 (`LICENSE-DATA`),
  `bin/normaliz` is an unmodified GPL-3 binary of Normaliz (provenance in
  `bin/`).

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
- **cube state, generalised cube state** — the CSS state of the code of affine
  functions on the points e_i, c − e_i, c, 0 of F^{t+1} (t = 2: the cube
  {0,1}³), one qudit per point, purifier = origin; witnesses of Theorems 2, 2′
  and 5 (`s27`, `s28`, `s29`).
- **Z_d qudit** — a qudit with the Weyl–Heisenberg group of Z_d, d not
  necessarily prime (R34). Galois qudits over GF(p^k) are different: they give
  the cone of the prime p.
- **mask order** — coordinate index = bitmask of the subset minus one, with
  parties A..E = bits 0..4; the literature's lexicographic-by-size order is
  converted by the permutation used in `scripts/` and `shc_table3_map.json`.
