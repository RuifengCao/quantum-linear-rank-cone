# Stabilizer entropy cones as self-dual slices of linear rank cones, and their dependence on the local dimension

Technical note accompanying round R32. Status flags as in `docs/A2_stabilizer_vs_QLR.md`:
**[proved]** — proof in this note; **[cited]** — published result, statement checked against the source;
**[certified]** — exact computation re-checked by a selftest gate.

**R37 (2026-10-11).** The results of this note are presented, with full proofs, in the paper *Stabilizer entropy cones
as self-dual slices of linear rank cones: a five-party gap and dependence on the local dimension* (X. Song and R. Cao,
2026); Table VI of the paper maps its theorems to those of this note. Those proofs were checked line by line by hand
before the paper was submitted, the statements it cites were checked against the sources, and the numbers it quotes
are recomputed by the scripts in `checks/`. Where this note and the paper differ, in statements, numbering,
attributions or references, the paper is authoritative. The prior work found by the literature check is listed in §6.

Round R33 adds §4.4: the stabilizer cones of **any two** distinct primes are incomparable (Theorem 2′). The proof
is self-contained, through Lemma 4. R33 also adds §4.5, a note on a misprint in DFZ arXiv:1401.2507.

Round R34 adds §4.6: Z_d qudits for every d ≥ 2. Theorem 1 extends to subgroups of abelian groups (Theorem 1(d)),
Lemma 4 extends to abelian groups (Lemma 4(G)), and for enough parties Stab^(Z_d) ⊆ Stab^(Z_d′) holds exactly
when d divides d′ (Theorem 5).

Round R35 adds §4.7: how far down the dependence goes. There is none for n ≤ 4 and any d (Theorem 7). The one-cut
observation of §4.4 is now proved (Proposition 6). At n = 6, the characteristic-dependent seven-variable
inequalities considered here are blind in every relabelling, because they collapse to Shannon inequalities on the
self-dual slice (Proposition 8).

Round R36 adds §4.8: six parties already suffice in one direction. Stab⁽²⁾_6 ⊄ Stab^(Z_d)_6 for every odd d
(Theorem 9). The witness is the seven-qubit simplex-code (Fano) state. The certificate is a new seven-variable
inequality, obtained by applying Lemma 4 to a family derived from the subspaces ("transfer"). R36 also adds the
inequalities of Blasiak–Kleinberg–Lubetzky to Proposition 8 (now eight forms) and records a misprint in their
eq. (6.17) (§4.5).

Round R38 adds §4.9: 577 of the 641 known extreme rays of Stab₅ lie in Stab^(Z_d)_5 for every d, 564 by unimodular
realisations and 13 by covering certificates (certified; independent referee, not yet checked by hand; the paper does
not contain it). Until R39 this paragraph and §5 counted only the 564 (erratum, CHANGELOG R39).

Round R39 continues §4.9: **636 of the 641** lie in Stab^(Z_d)_5 for every d (566 by unimodular realisations,
70 by covering certificates; same status). 5 rows remain open.

Round R40 completes §4.9: **all 641** known extreme rays of Stab₅ lie in Stab^(Z_d)_5 for every d (566 by unimodular
realisations, 75 by covering certificates; same status). R40 also adds §4.10, the other direction at six parties.
Among the states with at most one qudit per system, only the non-Fano state can witness Stab⁽³⁾_6 ⊄ Stab⁽²⁾_6
(certified). No GF(2^k)-linear state has entropy vector k·S_NF (Proposition 10, proved; the
non-Fano entropy vector in units of log 2^k).

Round R41 strengthens Proposition 10 (§4.10): no qubit stabilizer state has an entropy vector that is a positive
multiple of S_NF, and S_NF is not a finite nonnegative combination of qubit stabilizer entropy vectors (Corollary 11).
So the cones generated without closure already separate in the other direction at n = 6; for the closed cones of the
paper the question stays open (limits). Proved (Corollary 11 uses the certified face ranks of §4.8); independent
referee; not yet checked by hand; not in the paper. R41 also corrects R40's Proposition 10(d), which was true only
vacuously (erratum in §4.10). (R42: the closed cones are settled as well; next paragraph.)

Round R42 adds §4.11: **Stab^(Z_d)_6 ⊄ Stab⁽²⁾_6 for every d that is not a power of 2** (Theorem 12), so at six
parties qubits and every odd local dimension are incomparable (with Theorem 9; before, from n = 7). The witness is
the non-Fano state; the certificate is a new inequality f12 ≥ 0 on qubit stabilizer states, with f12(S_NF) = −1. Its
proof is a quantitative ("symplectic") version of the proof of Proposition 10: isotropy bounds the dimensions of the
restrictions of the line pieces, and Lemma 4 (LA) is applied to a dual subspace arrangement built from them.
f12∘π is a self-dual characteristic-2 linear rank inequality on seven variables that is not blind on the self-dual
slice. Proved; independent referee; not yet checked by hand; not in the paper.

## 1. Setting

* n visible parties A₁, …, A_n and a purifier F; E is the set of all n + 1 parties. A *pure-symmetric* vector S
  assigns S(X) to every X ⊆ E with S(X) = S(E∖X) and S(∅) = S(E) = 0; equivalently it is any vector on the
  2ⁿ − 1 nonempty subsets of the visible parties, extended by purity (the coordinates used throughout the
  repository).
* h_S(X) = S(X) + Σ_{P∈X} S(P) for X ⊆ E (called h_norm in the A2 note).
* π(r)(X) = r(X) + r(E∖X) − r(E): the connectivity function of a polymatroid r on E.
* r is *self-dual* if r(E∖X) = r(E) − Σ_{P∈X} r(P) + r(X) for all X, i.e. r = r* with
  r*(X) = Σ_{P∈X} r(P) + r(E∖X) − r(E).
* LR⁽ᵖ⁾_{n+1}: the closed convex cone generated by rank functions X ↦ dim Σ_{P∈X} U_P of families of n + 1
  subspaces of vector spaces over GF(p). LR_{n+1}: the closed convex hull of their union over all primes.
* Stab⁽ᵖ⁾_n: the closed convex cone generated by the entropy vectors (units of log p) of stabilizer states of prime
  local dimension p, each party holding any number of qudits. Stab_n: closed convex hull of the union over primes
  (the Stab₅ of the A2 note).
* CSS⁽ᵖ⁾_n: the closed convex cone generated by the entropy vectors of CSS states |C⟩ = |C|^{-1/2} Σ_{c∈C} |c⟩,
  C ⊆ GF(p)^m linear, the m qudits distributed among the parties.

## 2. Theorem 1 (reduction) [proved]

For every prime p and every n ≥ 1, as sets of pure-symmetric vectors,

    Stab⁽ᵖ⁾_n  =  { S : h_S ∈ LR⁽ᵖ⁾_{n+1} }  =  π(LR⁽ᵖ⁾_{n+1})  =  CSS⁽ᵖ⁾_n ,

and the same holds with all primes together (Stab_n, LR_{n+1}). Equivalently, S ↦ h_S maps Stab⁽ᵖ⁾_n bijectively
onto the self-dual part of LR⁽ᵖ⁾_{n+1}.

*Proof.*

1. *Stab ⊆ {h_S ∈ LR}.* By Lemmas 1 and 2 of the A2 note, h_S is the rank function of n + 1 subspaces over
   GF(p) for every stabilizer state of local dimension p. The proofs there do not use the number of parties.
   h_S is linear in S and LR⁽ᵖ⁾ is closed.
2. *CSS formula.* Let M be a k × m matrix over GF(p) whose columns are distributed among the parties, r(X) the
   rank of the columns of the parties in X, and C the row space of M. For X ⊆ E write c = (c_X, c_{X^c}) and
   K_X = {c ∈ C : c_{X^c} = 0}. The reduced state of |C⟩ on X is |C|^{-1} Σ_{c,c′∈C, c_{X^c}=c′_{X^c}} |c_X⟩⟨c′_X|;
   grouping the codewords by c_{X^c}, it is proportional to the projector onto the span of the coset states
   Σ_{u ∈ v + (K_X)_X} |u⟩, v ∈ C_X, which are orthonormal up to a common factor and p^{r(X) − dim K_X} in number.
   Since dim K_X = r(E) − r(E∖X),

       S(X) = r(X) + r(E∖X) − r(E) = π(r)(X).

   Every family of subspaces is the column-group family of some matrix, so π(LR⁽ᵖ⁾) ⊆ CSS⁽ᵖ⁾ ⊆ Stab⁽ᵖ⁾.
3. *π(h_S) = 2S for pure-symmetric S:* π(h_S)(X) = S(X) + S(E∖X) − S(E) + Σ_{P∈X} S(P) + Σ_{P∉X} S(P)
   − Σ_{P∈E} S(P) = 2S(X).
4. If h_S ∈ LR⁽ᵖ⁾, then S = ½ π(h_S) ∈ π(LR⁽ᵖ⁾) ⊆ CSS⁽ᵖ⁾ ⊆ Stab⁽ᵖ⁾. With step 1 the four sets coincide.

All primes: a family over GF(p^k) becomes a family over GF(p) with k times the ranks (restriction of scalars),
and a family over a field of characteristic 0 has the same rank function as a family over a suitable finite field
(take a representation over a finitely generated ring and reduce modulo a maximal ideal that avoids the
finitely many nonzero minors). Hence LR_{n+1} is exactly the cone cut out by all (n + 1)-variable linear rank
inequalities, and the statement for all primes follows from the one for each p. ∎

**Remarks.**

* *Pure symmetry is necessary.* For n = 1 the non-pure vector S(A) = S(F) = S(AF) = 1 has h_S = (2, 2, 3), the rank
  function of two planes meeting in a line, but it is not an entropy vector of a pure state.
* *No closure needed for the doubling.* For any stabilizer state with Lagrangian L ⊂ GF(p)^{2N}, the classical code
  L, with each party owning the x- and z-coordinates of its qudits, has column-group rank function
  r_L(X) = dim π_X(L) = |X| + S(X) (proof of Lemma 2). Its CSS state therefore has entropy
  π(r_L) = 2S: **every stabilizer entropy vector, doubled, is a CSS entropy vector of the same local dimension.**
  Gate G29 checks this against state vectors for random graph states (L = rowspace [I | W]).
* *Balanced form.* Replacing each subspace V of a family by V ∩ Σ_{U≠V} U lowers every rank containing V by the
  same amount H(V | rest). So for a valid linear rank inequality f, the balanced version
  f_bal(h) := f(h) − Σ_V f(e_V) H(V | rest) is again valid; e_V is the modular unit vector of V. This is the
  linear-rank analogue of Chan's balancing theorem. For a pure-symmetric S every element of h_S has H(V | rest) =
  h_S(E) − h_S(E∖V) = S(V) − S(E∖V) = 0, so f_bal(h_S) = f(h_S). A balanced form vanishes on modular functions,
  so f_bal(h_S) = f_bal(S). Hence

      Stab_n = { S pure-symmetric : f(S) ≥ 0 for every balanced (n + 1)-variable linear rank inequality f },

  with the purifier as the (n + 1)-st variable. The inclusion ⊆ goes back to Gross–Walter (arXiv:1302.6990,
  Thm 2); for every prime it follows from Lemma 2 of the A2 note. Theorem 1 supplies the converse.
* *Symplectic structure is irrelevant at the cone level.* Any linear representation of h_S, with no isotropy
  condition, gives a stabilizer (CSS) state with entropy 2S.

## 3. Consequences

### 3.1 S7′ holds [proved]

Stab₅ is cut out by the six-variable linear rank inequalities (purifier as the sixth variable), applied to h_S, or
equivalently, by the balanced ones applied directly to S. This answers the open question S7′ of the A2 note
(§9) in the affirmative. It also reframes the 1,155 undecided catalogue orbits: each asks whether one explicit
six-element polymatroid h_r lies in the closed cone of linearly representable polymatroids. This is a classical
question with no quantum content left.

### 3.2 Why QLR fails exactly at five parties [certified for n = 4; remark for n = 5]

QLR_n (BCHS) is built from n-variable linear rank inequalities, while Theorem 1 needs (n + 1)-variable ones on the
self-dual slice.

* **n = 4.** On five elements (A, B, C, D and the purifier E):
  * the 24 five-variable inequalities of DFZ (`data/dfz_ref28.csv`, arXiv:0910.0284 eqs. (1)–(24)) pulled back
    through h_S give 820 distinct rows, and every one of them is implied by the pulled-back Shannon and Ingleton
    (eqs. (36)–(39)) rows (LP, floating point);
  * the cone h⁻¹(Shannon + Ingleton) has exactly 46 extreme rays (lrs, exact arithmetic);
  * every one of the 46 is realised by a qubit graph state, 41 at λ = 1 and 5 at λ = 2 (`data/s26_stab4_certs.npz`,
    re-verified by `stabcert.verify_realisation_n`, gate G29; exact enumeration re-run by G29b).

  Hence Stab₄ = h⁻¹(Shannon + Ingleton on five elements). This exact equality follows from the realisations alone
  and needs neither the LP nor DFZ's five-variable classification. It agrees with the four-party characterisation of
  Linden–Matúš–Ruskai–Winter (LIPIcs TQC 2013, arXiv:1302.5453, Sec. 5; checked by hand in R37). As printed,
  their state (R6) has S(ab) = 1 instead of 2 on ray 6 of their Table 1, as N. Osthus noted (M.S. thesis,
  University of Minnesota, 2024, Sec. 3.8); `checks/check_sec4a.py` recomputes both states.
* **n = 5.** The new six-variable inequalities do cut the self-dual slice: these are the exclusions of R27–R30
  (A2 note §3–§7).

## 4. Theorem 2 (dependence on the local dimension) [proved; §4.2 uses DFZ, §4.4, §4.6 and §4.8 are self-contained]

For every n ≥ 7 and every odd prime p,

    Stab⁽²⁾_n ⊄ Stab⁽ᵖ⁾_n   and   Stab⁽ᵖ⁾_n ⊄ Stab⁽²⁾_n.

The first non-inclusion already holds for n = 6 (Theorem 9, §4.8), and so does the second (Theorem 12, §4.11; R42).

### 4.1 Witnesses: the cube states

Put one qudit on each vertex x ∈ {0,1}³ of the cube. The code C_F ⊂ F⁸ consists of the evaluations of the affine
functions f(x) = a₀ + a·x; its generator matrix has rows 1, x₁, x₂, x₃.

* Over GF(2), C is the first-order Reed–Muller code RM(1,3), i.e. the extended Hamming [8,4,4] code. Its CSS state is
  local-Clifford equivalent to the graph state of the 3-cube.
* Over GF(p), p odd, the same construction gives a [8,4] code on the eight cube points only. This is not RM(1,3),
  which would use all p³ points.

The visible parties are the seven nonzero vertices, in DFZ's order A, B, C, W, X, Y, Z = 100, 010, 001, 110,
101, 011, 111; the origin 000 is the purifier. Their entropies [certified, `scripts/s27_cube_states.py`, G30]:

* S(X) = 2 r(X) − |X|, where r is the rank function of the code's matroid.
* Both states are 3-uniform: S(X) = |X| for |X| ≤ 3.
* On 4-sets, S = 2 for coplanar quadruples and 4 otherwise:
  * GF(2): 14 coplanar quadruples, namely 6 faces, 6 diagonal rectangles and 2 tetrahedra;
  * GF(p): 12 coplanar quadruples. The tetrahedra are not coplanar there, because their 4 × 4 minor is ±2.
* The two entropy vectors differ in exactly one coordinate. This is the bipartition {A, B, C, Z} | {W, X, Y, F}
  (odd- versus even-weight vertices), with 2 for qubits and 4 for odd p.
* The entropies agree with explicit state vectors (2⁸ and 3⁸ amplitudes, all 254 proper cuts).

### 4.2 Proof

1. **Identical self-duality** [certified]. Both matroids satisfy r* = r. So S = 2r − |X| and h_S = 2r.
2. **Contraction at the origin** [certified]. Quotienting by the column (1, 0, 0, 0) of the origin leaves exactly
   the seven vectors e₁, e₂, e₃, e₁+e₂, e₁+e₃, e₂+e₃, e₁+e₂+e₃ of DFZ's Theorems 8.7 and 8.9. Over GF(2) this is
   the Fano configuration; over GF(p) it is the non-Fano configuration.
   * All 4 × 4 minors of the generator lie in {0, ±1, ±2}. So for every odd p the matroid, and hence the entropy
     vector, is the same as over GF(3).
3. **DFZ, Section 8** [cited]. Dougherty–Freiling–Zeger, *Achievable Rate Regions for Network Coding*,
   arXiv:1311.4601v1. The statements were checked against the PDF: p. 35 (Thm 8.6, eq. (65)), pp. 39–40
   (Thm 8.7), p. 43 (Thm 8.8, eq. (91)) and pp. 48–49 (Thm 8.9). Coefficients are in `data/dfz_char_dependent.csv`.
   * (65) holds for subspaces over any field of odd characteristic, and fails for the Fano configuration over
     GF(2).
   * (91) holds over characteristic 2, and fails for the non-Fano configuration over GF(p), p odd.

   Blasiak–Kleinberg–Lubetzky (arXiv:1108.2489, Thms 6.4 and 6.6) give seven-variable inequalities with the same
   properties. They are transcribed in `data/bkl_char_dependent.csv` (R36, checked against the PDF; eq. (6.17) is
   misprinted, §4.5) and are not used here.
4. **Pull-back** [proved]. For a stabilizer state of local dimension p, h_S is a family of eight subspaces over
   GF(p) (Lemma 2). Passing to the quotient by the purifier's subspace gives a family of seven subspaces with
   rank function T ↦ h_S(T ∪ F) − h_S(F). Define
   * F65(S) := (65) applied to this contraction, and
   * F91(S) := (91) applied to it.

   Then F65 ≥ 0 on Stab⁽ᵖ⁾_7 for every odd p, and F91 ≥ 0 on Stab⁽²⁾_7. Both are linear in S, so they hold on the
   closed cones.
5. **Values** [certified, G30].
   * F65(S_qubit) = 2 × (value of (65) on the Fano configuration) = −2.
   * F91(S_qudit) = −2.
6. **n > 7.** Merge the extra (empty) parties into the purifier. ∎

The inequalities F65 and F91 are stored, as 127 coefficients in mask order over A..Z, in
`data/s27_cube_witnesses.npz`. Neither depends on knowing which party is the purifier beyond the contraction
used above.

**Validation, not proofs:**

* (65) has no violation on random subspace families over GF(3) and GF(5), and (91) none on random families over
  GF(2). A one-off independent re-implementation, which is not part of the repository, also checked them
  exhaustively:
  * (91) on all 7-tuples of vectors in GF(2)³ (2.1 M) and GF(2)⁴ (268 M);
  * (65) on GF(3)³ up to scalars (105 M);
  * both are violated by the corresponding Fano / non-Fano tuples.
* The pulled-back F65 is never negative on random qutrit and p = 5 graph states with 1–3 qudits per party, and F91
  never on random qubit graph states (`s27_cube_states.py --random`).

### 4.3 Doubled codes: a sandwich for the dimension dependence [proved]

**Proposition 3.** Let G = [G₁ | … | G_n] be a matrix over GF(p) whose column blocks represent a polymatroid r on
n elements. Take the CSS state of the code generated by [G | G], with party i holding the qudits of G_i and the
purifier holding all the copies. Then for every set X of visible parties

    h_S(X) = r(X) + Σ_{i∈X} r(i).

*Proof.* The columns outside X include all the copies, which span the whole column space. So
S(X) = r(X) + r(E) − r(E) = r(X), and in particular S(i) = r(i). ∎

**Corollary (sandwich).** For primes p ≠ q and every n:

    LR⁽ᵖ⁾_n ⊄ LR⁽q⁾_n   ⟹   Stab⁽ᵖ⁾_n ⊄ Stab⁽q⁾_n   ⟹   LR⁽ᵖ⁾_{n+1} ⊄ LR⁽q⁾_{n+1}.

*Proof.* The right implication is Theorem 1. For the left one:

1. Take a GF(p)-representable polymatroid r and an inequality f, valid on LR⁽q⁾_n, with f(r) < 0. Such an r
   exists, since LR⁽ᵖ⁾_n is generated by representable polymatroids.
2. The modular units e_V lie in LR⁽q⁾_n, so f(e_V) ≥ 0. Hence the balanced version f_bal (Remark in §2) is valid
   on LR⁽q⁾_n and satisfies f_bal(r) ≤ f(r) < 0.
3. For a stabilizer state of local dimension q, the restriction of h_S to the visible parties is a GF(q) family
   (Lemma 2), so f_bal(h_S) ≥ 0 there. Because f_bal is balanced, this is a linear inequality on S.
4. The doubled-code state of r has h_S = r + (modular) on the visible parties, so it gives f_bal(r) < 0. ∎

So the dependence of Stab_n on the local dimension is sandwiched between that of the n- and the (n + 1)-variable
linear rank cones.

* For n = 7, the corollary gives a second proof of Theorem 2, through the 14-qudit doubled Fano and non-Fano
  states. On them the balanced (65) and (91) take the value −1; see `s27_cube_states.py` and G30. The cube states
  are the more economical witnesses.
* For n = 5, a dimension dependence would require a characteristic-dependent six-variable linear rank
  inequality.

### 4.4 Theorem 2′ (any two primes) [proved] (R33)

For primes p < q and every n ≥ 2p + 3,

    Stab⁽ᵖ⁾_n ⊄ Stab⁽q⁾_n   and   Stab⁽q⁾_n ⊄ Stab⁽ᵖ⁾_n.

The proof is self-contained: it uses Lemma 2 of the A2 note and Lemma 4 below. Peña–Sarria's inequalities
(arXiv:1905.00003v3, Theorem 3 / Example 6) are of the same kind and serve as a cross-check, but they are not needed.
(R37) For prime fields, Theorem 2′ also follows from Theorems 7 and 9 of Peña-Macias and Sarria (J. Inf. Optim.
Sci. 42, 345 (2021), arXiv:1903.11587v2), which treat the same configuration for every t; Lemma 4 adds that both
inequalities hold for finite abelian groups, which Theorems 5 and 9 need (see the paper, Sec. V B).

**Lemma 4 (characteristic-dependent inequalities for the configuration J − I) [proved].** Let M = t + 1 ≥ 3. Let
A₁, …, A_M, B₁, …, B_M and C be subspaces of a finite-dimensional vector space over a field F. Write A = Σ_k A_k
and A_{≠i} = Σ_{j≠i} A_j, and define the error terms (all ≥ 0):

* α = Σ_k H(A_k) − H(A),
* γ = H(C | A),
* δ_i = I(C; A_{≠i}),
* ε_i = H(C | A_i, B_i),
* ζ_i = H(B_i | A_{≠i}),
* η_i = H(B_i | A_i, C).

Then

* (LB) if char F does not divide t:  M·H(C) ≤ H(B₁, …, B_M) + M·[γ + M·α + Σ_i (δ_i + ε_i + ζ_i)];
* (LA) if char F divides t:  H(B₁, …, B_M) ≤ (M − 1)·H(C) + M·α + M·γ + Σ_i (ζ_i + η_i).

On the configuration A_i = ⟨e_i⟩, B_i = ⟨c − e_i⟩, C = ⟨c⟩ in F^M every error term vanishes, H(C) = 1, and
H(B₁, …, B_M) = rank(J − I) = M − 1 if char F | t and M otherwise. The determinant of J − I is (−1)^t·t. When
char F | t the rank is exactly M − 1, because an (M − 1)-minor equals ±(t − 1) ≡ ∓1. So (LB) fails by 1 when
char F | t, and (LA) fails by 1 when char F ∤ t.

*Proof.* Facts used:

* (F1) codim_X(X ∩ W₁ ∩ … ∩ W_k) ≤ Σ codim_X(X ∩ W_j).
* (F2) For W ⊆ W′: codim_X(X ∩ W) ≤ codim_X(X ∩ W′) + codim_{W′}(W). Indeed
  codim_X(X ∩ W) = codim_X(X ∩ W′) + codim_{X∩W′}(X ∩ W), and (X ∩ W′)/(X ∩ W) embeds into W′/W.
* (F3) codim_X(X ∩ (Y + Z)) = H(X | Y, Z).

*Step 0 (making the A_i independent).*

* Let O_i = A_i ∩ A_{≠i}. Then dim O_i = I(A_i; A_{≠i}) ≤ α. Let O = Σ_i O_i, so dim O ≤ M·α and O ⊆ A_{≠i} for
  every i.
* Let q: V → V/O. The images Ā_i = q(A_i) are independent. Indeed, if a_i − Σ_{j≠i} a_j = Σ_k o_k with o_k ∈ O_k,
  then a_i − o_i ∈ A_i ∩ A_{≠i} = O_i, so q(a_i) = 0.
* Write π_i for the projection of Ā = ⊕_k Ā_k onto Ā_i. Then q(A_{≠i}) = ⊕_{j≠i} Ā_j = ker π_i.

*(LB).*

1. Let B′_i = B_i ∩ A_{≠i} and C₁ = C ∩ ⋂_i (A_i + B′_i) ⊆ A.
   * We have codim_{A_i+B_i}(A_i + B′_i) ≤ dim B_i − dim B′_i = ζ_i.
   * By (F2) and (F3), codim_C(C ∩ (A_i + B′_i)) ≤ ε_i + ζ_i. Also codim_C(C ∩ A) = γ.
   * By (F1), dim C₁ ≥ H(C) − γ − Σ_i (ε_i + ζ_i).
2. For x ∈ C₁ write x = a + b with a ∈ A_i and b ∈ B′_i ⊆ A_{≠i}. Uniqueness of the decomposition in ⊕Ā_k gives
   π_i(q x) = q a, so q x − π_i(q x) = q b ∈ q(B_i).
3. Let K_i = q(C₁) ∩ ker π_i. Its preimage is (C₁ + O) ∩ A_{≠i} = (C₁ ∩ A_{≠i}) + O (modular law, O ⊆ A_{≠i}).
   Hence dim K_i ≤ dim(C ∩ A_{≠i}) = δ_i.
4. Let D be a complement of Σ_i K_i in q(C₁). Then
   dim D ≥ dim C₁ − dim O − Σ_i δ_i, and D ∩ ker π_i = 0 for all i.
5. Put β_i(x) = x − π_i x for x ∈ D; by step 2, β_i(D) ⊆ q(B_i).
   * Each β_i is injective: β_i x = 0 forces x ∈ Ā_i ⊆ ker π_k for k ≠ i.
   * The β_i(D) are independent. Suppose Σ_i β_i(x_i) = 0. Applying π_k gives π_k(Σ_{i≠k} x_i) = 0, so
     Σ_{i≠k} x_i = 0 for every k. Summing over k gives t·Σ_i x_i = 0. When char F ∤ t this means Σ_i x_i = 0, and
     then every x_k = 0.
   * Hence H(B₁, …, B_M) ≥ dim q(Σ B_i) ≥ M·dim D. Substituting dim O ≤ M·α gives (LB).

*(LA).*

1. Let C_A = C ∩ A, with codim γ, and B″_i = B_i ∩ A_{≠i} ∩ (A_i + C_A).
   * By (F1)–(F3), codim_{B_i}(B″_i) ≤ ζ_i + η_i + γ. This uses codim_{A_i+C}(A_i + C_A) ≤ γ.
2. For b ∈ B″_i write b = a + x with a ∈ A_i and x ∈ C_A.
   * Since q b ∈ ker π_i, we get q a = −π_i(q x).
   * So q b = β_i(q x), where β_i(y) = y − π_i y on q(C_A) ⊆ Ā.
3. The map Φ(y₁, …, y_M) = Σ_i β_i(y_i) on q(C_A)^M contains Σ_i q(B″_i) in its image.
   * Its kernel contains the diagonal, because Σ_i β_i(y) = M·y − y = t·y = 0 when char F | t.
   * Hence dim q(Σ_i B″_i) ≤ (M − 1)·dim q(C_A) ≤ (M − 1)·H(C).
4. Finally:
   * H(B₁, …, B_M) ≤ dim(Σ_i B″_i) + Σ_i codim_{B_i}(B″_i);
   * dim(Σ_i B″_i) ≤ dim q(Σ_i B″_i) + dim O ≤ (M − 1)·H(C) + M·α.

   Together these give (LA). ∎

*Remarks.*

* The γ terms are not needed: C₁ ⊆ A holds automatically, and B_i ∩ A_{≠i} ∩ (A_i + C) ⊆ A_i + C_A.
* Also dim O ≤ (M − 1)·α.
* The proof works for M ≥ 2; (LA) is non-empty only for M ≥ 3.
* We keep the stated (slightly weaker) forms, which are the ones encoded and tested.
* An independent referee re-derived every step and stress-tested both forms. Billions of configurations over GF(2),
  GF(3), GF(5) and GF(7) for t = 2, 3 were checked, exhaustive up to GL in several families, together with
  adversarial searches. The valid form was never violated, while the same searches do find the wrong-characteristic
  violations and the violations of deliberately weakened forms.

`scripts/s28_generalized_cubes.py:lemma_forms` encodes (LA) and (LB). They are sanity-checked, which is not a
proof:

* on all {0,1}-point configurations of GF(2)³ and GF(3)³ (t = 2; 2.1 M each). The valid form is never violated,
  and the wrong-characteristic form is violated there (gate G31b);
* on random {0,1}-point configurations of GF(q)⁴ (t = 3) and on random subspace arrangements.

Peña–Sarria's Example 6 (a) and (b) (typed from the rendered page 4 of arXiv:1905.00003v3, `ps_forms`) behave in the
same way and give the same values below.

**Generalised cube states.** Let t ≥ 2, M = t + 1 and c = (1, …, 1) ∈ F^M. Take the 2t + 4 points

* A_i = e_i (i = 1, …, M),
* B_i = c − e_i (i = 1, …, M),
* C = c,
* the origin 0,

one qudit on each, and the code of affine functions on them (generator rows 1, x₁, …, x_M). The parties
A_1, …, A_M, B_1, …, B_M, C are visible and the origin is the purifier. For t = 2 these are the eight vertices of
the cube {0,1}³, i.e. the cube states of §4.1. For Theorem 2′ take t = p, over F = GF(p) and over F = GF(q).

*Proof of Theorem 2′.*

1. **Self-duality** [proved]. Give the points the weights d = −1 on each e_i, +1 on each c − e_i, −(t − 1) on c and
   t − 1 on 0. The identities
   * Σ d = 0,
   * Σ d·x_j = −1 + (M − 1) − (M − 2) = 0,
   * Σ d·x_j x_k = (M − 2) − (M − 2) = 0 (j ≠ k)

   hold over ℤ, so G diag(d) Gᵀ = 0. If char F ∤ (t − 1), all weights are nonzero, so C^⊥ = C·diag(d). Since the
   code has half the length as dimension, its matroid is identically self-dual: S = 2r − |X| and h_S = 2r.
   * For t = p this holds over GF(p), since p ∤ p − 1, and over GF(q) for q > p, since 0 < p − 1 < q.
   * Over GF(2) with t odd it fails; that pair is covered by Theorem 2.
2. **Contraction** [proved]. Contracting the origin leaves exactly the configuration of Lemma 4. This is also the
   counterexample configuration of Peña–Sarria, Theorem 3, with M(n, t) = t + 1.
3. **Pull-back** as in §4.2, with Lemma 4 in place of DFZ. The values below are certified for t = 3 (GF(3) against
   GF(5), GF(7)) and t = 5 (GF(5) against GF(7)) by gate G31; t = 2 coincides with §4.1. Peña–Sarria's (b)/(a) give
   the same values.
   * The GF(p) state takes the value −2 on the pulled-back (LB), which holds for every q-qudit stabilizer state
     since q ∤ p.
   * The GF(q) state takes the value −2 on the pulled-back (LA), which holds for every p-qudit stabilizer state.
4. n > 2p + 3: merge empty parties into the purifier. ∎

Notes:

* For t = 2, Lemma 4 also re-proves Theorem 2 without DFZ.
* The choice t = min(p, q) always works, since self-duality over GF(q) needs q ∤ t − 1. With t = max(p, q) one
  would need p ∤ q − 1, which fails for instance whenever p = 2.

**One cut again.** In every case computed (t = p = 2, 3, 5 against primes q > p up to 11; for p = 2 every odd q gives
the same vector, §4.1) the two entropy vectors differ in exactly one
coordinate. It is the bipartition {A_1, …, A_M, C} | {B_1, …, B_M, 0}, with entropy t over GF(p) and t + 2 over
GF(q). Proposition 6 (§4.7) proves this for all t and q.

### 4.5 Misprints in DFZ arXiv:1401.2507 and BKL arXiv:1108.2489 (not used)

**DFZ (R33).** Theorem 3.1 of arXiv:1401.2507v1 (p. 11) states an eight-variable inequality claimed to hold over every finite
field of characteristic ≠ 3; Theorem 3.2 shows that it fails for T₈ over GF(3). As printed, the inequality is false
over every field.

* *Smallest counterexample.* Take C = X = one line and all other subspaces zero. The left side H(A) is 0. On the
  right side, 3H(X) − 17H(C) + 7H(C|D,Y,Z) + 3H(C|A,W,Y) = 3 − 17 + 7 + 3 = −4, and every other term vanishes.
* *Search.* An exhaustive search over all subspace arrangements of GF(2)² and of GF(3)² finds 125 and 254 violating
  arrangements respectively.

The final expansion of the proof (p. 25) probably carries a coefficient error; the published IEEE version was not
checked. Theorem 4.1 of the same paper (valid for characteristic 3) passes the same exhaustive check, but it is not
used either. Section 4.4 needs neither. (R37) Sec. V A of the paper uses DFZ's Theorems 4.1 and 4.2, as stated
there, together with Proposition 3 of §4.3, for one direction at eight parties; that inequality was not checked
independently beyond the exhaustive test above.

**BKL (R36).** Blasiak–Kleinberg–Lubetzky, arXiv:1108.2489v1, Section 6, give two seven-variable inequalities on the
points of GF(2)³ (transcribed term by term from the rendered PDF, p. 19 eq. (6.8) and p. 22 eq. (6.17), into
`data/bkl_char_dependent.csv`, variables named as in DFZ):

* (6.8), Theorem 6.4: valid in characteristic 2, −1 on the non-Fano configuration. It passes our checks: minimum 0
  over all {0,1}-point configurations of GF(2)³ (G34b).
* (6.17), Theorem 6.6: claimed valid in characteristic ≠ 2, −1 on the Fano configuration. **As printed it is false
  over every field.** Take A = ⟨e₁⟩, B = ⟨e₂⟩, C = ⟨e₃⟩, Z = ⟨e₁ + e₂ + e₃⟩ and W = X = Y = 0: the value is −3 over
  GF(2), GF(3) and GF(5) (G34). Over GF(3)³ the printed form reaches −3 on {0,1}-point configurations.
  * (R37) Peña-Macias and Sarria (arXiv:1903.11587v2, before their Theorem 9) had already noted that (6.17)
    fails, with the same counterexample, and gave a different inequality in its place (their Theorem 9). This
    note did not know their paper when it was written.
  * The proof of Theorem 6.6 bounds the step-4 terms Σ_{i+j=111} δ[111 | i, j] in the text, but (6.16) drops them.
    Restoring them adds 3·[H(Z,A,Y) − H(A,Y) + H(Z,B,X) − H(B,X) + H(Z,C,W) − H(C,W)] to (6.17).
  * The corrected form (`bkl617_odd_corrected`) has minimum 0 over all {0,1}-point configurations of GF(3)³ (G34b)
    and is −1 on the Fano configuration. It equals Lemma 4 (LB)₂ (§4.8 labelling) minus 6·I(A; B), so it is slightly
    stronger than (LB)₂. Its validity rests on BKL's proof with the dropped terms restored; it is not re-proved here.
* The proof of Theorem 6.4 ends with "Consequently Λ_odd · r(F) < 0"; from the statement it should read
  Λ_even · r(N) < 0 (a typo without consequences).

Both forms enter Proposition 8 (§4.7), (6.17) in its corrected form. Neither is needed elsewhere.

### 4.6 Composite local dimensions: the cone remembers d [proved] (R34)

**Setting.**

* A Z_d qudit carries the Weyl–Heisenberg group of Z_d, for any d ≥ 2 (even d included). The stabilizer group of a
  pure stabilizer state on m qudits maps injectively onto a Lagrangian subgroup L ⊆ Z_d^{2m}, with |L| = d^m.
* Entropies are measured in one fixed unit (bits), so that the cones of all local dimensions live in one space.
  Stab^(Z_d)_n is the closed convex cone generated by the entropy vectors of pure Z_d stabilizer states (any number
  of qudits per party; n visible parties plus a purifier). For a prime p it is Stab⁽ᵖ⁾_n.
* For a finite abelian group V and subgroups U_P, the *subgroup rank function* is H(X) = log |Σ_{P∈X} U_P|.
  Γ^(d)_m is the closed convex cone generated by subgroup rank functions on m elements with V of exponent dividing
  d. For a prime p, Γ^(p)_m = LR⁽ᵖ⁾_m.

**Theorem 1(d) [proved].** For every d ≥ 2 and every n,

    Stab^(Z_d)_n  =  { S : h_S ∈ Γ^(d)_{n+1} }  =  π(Γ^(d)_{n+1})  =  closed cone of CSS states over Z_d.

*Proof.*

1. *Entropy* [cited]. S(X) = log |proj_X L| − n_X·log d, where n_X is the number of qudits held by X
   (Gross–Walter arXiv:1302.6990, Thm 2, valid for every d; LMRW arXiv:1302.5453, Prop. 9, is equivalent).
2. *h_S is a subgroup rank function.*
   * Let L_P ⊆ L consist of the elements supported on party P. Since L = L^⊥, proj_P L = L_P^⊥ inside
     Z_d^{2n_P}, so log |L_P| = 2n_P·log d − log |proj_P L|. Hence h_S(X) = log |proj_X L| − Σ_{P∈X} log |L_P|.
   * Let Γ be the image of L in ⊕_P proj_P(L)/L_P. The kernel of proj_X L → ⊕_{P∈X} proj_P(L)/L_P is exactly
     ⊕_{P∈X} L_P, so h_S(X) = log |proj_X Γ|.
   * Dualising (characters of Γ), log |proj_X Γ| = log |Σ_{P∈X} U_P|, where U_P is the image of the dual of the P-th
     projection. The exponent of the character group of Γ divides d.
   * Contraction by the purifier (passing to V/U_F) stays in the class.
3. *CSS states.* For a subgroup C ⊆ ⊕_P Z_d^{m_P}, the state |C⟩ = |C|^{-1/2} Σ_{c∈C} |c⟩ is a stabilizer state
   (X-shifts by C, Z-phases by C^⊥, |C|·|C^⊥| = d^m). As in step 2 of Theorem 1 (cosets of the subcode supported
   on X), S(X) = π(r_C)(X) with r_C(X) = log |proj_X C|. Every subgroup rank function arises in this way: embed V
   into some Z_d^k, let the columns of a block G_P generate U_P, and take for C the row space of [G_1 | … | G_{n+1}];
   row space and column space of an integer matrix mod d have the same size (Smith form).
4. π(h_S) = 2S for pure-symmetric S (§2, step 3). Steps 1–4 give the equalities as in Theorem 1. The doubling
   remark of §2 also holds: the code L, with each party owning the x- and z-coordinates of its qudits, has
   r_L(X) = S(X) + n_X·log d, so 2S is exactly a CSS entropy vector over Z_d. ∎

**Lemma 4(G) [proved].** Lemma 4 (§4.4) holds verbatim for subgroups A_i, B_i, C of a finite abelian group V,
with H = log |·|:

* (LB) if multiplication by t is injective on V, i.e. gcd(t, |V|) = 1;
* (LA) if t·V = 0.

*Proof.* Facts (F1)–(F3) and Step 0 hold for subgroups, with indices in place of codimensions. The modular law holds
in the subgroup lattice of an abelian group, and t stays injective on V/O.

* (LA) uses only these facts, the decomposition Ā = ⊕ Ā_k with its projections π_i, and the map Φ, whose kernel
  contains the diagonal when t·V = 0. All of this is valid for groups.
* (LB): the complement D of step 4 need not exist in a group, so steps 4–5 are replaced by counting.
  * Put Y = q(C₁), K_k = Y ∩ ker π_k (so log |K_k| ≤ δ_k by step 3) and Φ(y₁, …, y_M) = Σ_i β_i(y_i) on Y^M.
    Its image lies in q(Σ_i B_i).
  * If Φ(y) = 0 and s = Σ_i y_i, then s − y_k ∈ K_k for every k, and t·s = Σ_k (s − y_k).
  * The map y ↦ (s, (s − y_k)_k) is injective, and for given κ ∈ Π_k K_k at most one s satisfies t·s = Σ_k κ_k.
    Hence |ker Φ| ≤ Π_k |K_k|.
  * So H(B₁, …, B_M) ≥ log |q(Σ_i B_i)| ≥ M·log |Y| − Σ_k δ_k ≥ M·(H(C) − Σ_i (ε_i + ζ_i) − M·α) − Σ_k δ_k, which
    implies (LB). ∎

For vector spaces this is also a shorter proof of (LB) in Lemma 4; the bound obtained needs neither γ nor the
factor M in front of Σ δ_i.

An independent referee re-derived Theorem 1(d), Lemma 4(G), the witnesses below and Theorem 5. It also
stress-tested Lemma 4(G) with about 1.7 M evaluations, using its own encoding of the forms (identical to
`lemma_forms`):

* groups such as Z_4^5, Z_8^4, Z_9^4, Z_2 × Z_4² × Z_8 and Z_12^6;
* near-J − I families with non-unit scalings, random subgroups, and hill-climbing.

No claimed form was violated. Every form claimed only for other groups did go negative, and so did weakened forms
(dropping ε, δ or ζ from (LB), or α, η or ζ from (LA)).

**Witnesses over Z_N.** Let N = p^{j+1}, and let t ≥ 2 with t − 1 a unit mod p. Take the generalised cube code of
§4.4 over Z_N (rows 1, x₁, …, x_M on the points e_i, c − e_i, c, 0).

* The code is free: the points 0, e₁, …, e_M give a unimodular maximal minor. So |C| = N^{M+1} = |C^⊥|.
* The weights of §4.4 are units and G diag(d) Gᵀ = 0, so C^⊥ = C·diag(d). Hence r is identically self-dual and
  h_S = 2r.
* Contracting the origin gives the configuration ⟨e_i⟩, ⟨c − e_i⟩, ⟨c⟩ over Z_N. All error terms of Lemma 4
  vanish, H(C) = log N, and H(B₁, …, B_M) = (M − 1)·log N + log(N / gcd(t, N)), because the Smith form of J − I is
  diag(1, …, 1, t).
* If gcd(t, N) = p^j, then (LA)_t takes the value −log p on the configuration and −2·log p on the state.

Certified by gate G32 (units of log p; `scripts/s29_zd_stabilizer.py`, `data/s29_zd_witnesses.npz`):

| N, t | role | pulled-back value on the state |
| --- | --- | --- |
| Z_4, t = 2 (the "Z_4 cube state", 8 qudits) | Z_4 against Z_2 and against every odd p | (LA)_2 = (LB)_2 = −2 bits |
| Z_8, t = 4 | Z_8 against Z_4 | (LA)_4 = −2 bits |
| Z_9, t = 3 | Z_9 against Z_3 | (LA)_3 = −2 trits |
| Z_4, t = 6 (G32b) | Z_4 against Z_6 | (LA)_6 = −2 bits |

The entropies of the Z_4 cube state also agree with its state vector (4^8 amplitudes, all 254 cuts, deviation
9·10⁻¹⁶). The forms are valid on the smaller group in each row: (LA)_4 on the Z_4 state with t = 4 gives 0.

**Theorem 5 (divisibility) [proved].** Let d, d′ ≥ 2.

* If d | d′, then Stab^(Z_d)_n ⊆ Stab^(Z_d′)_n for every n.
* If d ∤ d′, then Stab^(Z_d)_n ⊄ Stab^(Z_d′)_n for every n ≥ 2·max(d, d′) + 3.

So for n ≥ 2·max(d, d′) + 3, Stab^(Z_d)_n ⊆ Stab^(Z_d′)_n holds exactly when d divides d′. In particular, the cones
of the local dimensions 2 ≤ d ≤ (n − 3)/2 are pairwise distinct. For two primes, Theorem 2′ gives the better bound
n ≥ 2·min(p, q) + 3.

*Proof.* If d | d′, then Γ^(d) ⊆ Γ^(d′), and Theorem 1(d) gives the inclusion. If d ∤ d′, choose a prime p with
v_p(d) > v_p(d′) =: j.

* j = 0. The GF(p) generalised cube state with t = p (§4.4; 2p + 3 ≤ 2d + 3 parties) lies in
  Stab⁽ᵖ⁾ ⊆ Stab^(Z_d), since p | d. It violates (LB)_p by −2·log p. Since p is a unit mod d′, (LB)_p is valid on
  Stab^(Z_d′) by Lemma 4(G) and Theorem 1(d).
* j ≥ 1. Take N = p^{j+1}, which divides d, and t = d′. Then t − 1 ≡ −1 is a unit mod p and gcd(t, N) = p^j. So the
  generalised cube state over Z_N (2d′ + 3 parties) lies in Stab^(Z_N) ⊆ Stab^(Z_d) and violates (LA)_{d′} by
  −2·log p. Since d′·V = 0 for every group of exponent dividing d′, (LA)_{d′} is valid on Stab^(Z_d′).

For more parties, merge empty parties into the purifier. ∎

**Consequences and remarks.**

* *Prime powers.* For j ≥ 1 and n ≥ 2p^j + 3, Stab^(Z_{p^j})_n ⊊ Stab^(Z_{p^{j+1}})_n (t = p^j, witness over
  Z_{p^{j+1}}). For example, qubits ⊊ Z_4 from n = 7 (the Z_4 cube state), and Z_4 ⊊ Z_8 from n = 11.
* *Two theories of local dimension 4.* Galois qudits give exactly the qubit cone. Restriction of scalars gives
  Stab^(GF(4)) ⊆ Stab⁽²⁾ (R31, field reduction), and extension of scalars on CSS codes gives ⊇. Z_4 qudits give a
  strictly larger cone from n = 7 on.
* *The Z_4 cube state* violates both (LA)_2 and (LB)_2, so it lies outside Stab⁽ᵖ⁾_7 for every prime p.
  * All minors of its generator lie in {0, ±1, ±2}, so log₂ |·| over Z_4 equals rank over GF(2) plus rank over
    GF(3).
  * Its entropy vector in bits is therefore the qubit cube vector (bits) plus the qutrit cube vector (trits), and
    it lies in the convex hull of the prime cones (checked exactly by G32).
* *Four parties.* Subgroups have common information: Z = U_X ∩ U_Y, and the extended family is again a subgroup
  configuration. So every inequality derived from Shannon inequalities and common information holds on every
  Γ^(d), in particular Ingleton. Hence Stab^(Z_d)_4 ⊆ h⁻¹(Shannon + Ingleton) = Stab⁽²⁾_4 (§3.2), with equality for
  even d. Theorem 7 (§4.7) gives equality for odd d as well.
* *The A2 results hold for every local dimension.* Every exclusion certificate of the A2 catalogue (A2 note
  §3–§7, including (★)) combines Shannon inequalities (and, for the CI+Ingleton certificates, Ingleton
  inequalities) on the family extended by one common information. Both hold for subgroup configurations, since the
  extended family is again one. So these certificates are valid on Stab^(Z_d)_5 for every d. In particular
  Stab₅ ⊊ QLR₅ also holds for the cone over all local dimensions, which according to the AI pre-check of §6 is the
  cone LMRW consider.

### 4.7 How far down does the dependence go? [proved] (R35)

In short:

* n ≤ 4: no dependence on the local dimension at all, for any d (Theorem 7).
* Dependence, pairwise:
  * qubits against every odd p: Stab⁽²⁾_n ⊄ Stab⁽ᵖ⁾_n from n = 6 (Theorem 9, §4.8; also against Z_d for every odd
    d), and Stab⁽ᵖ⁾_n ⊄ Stab⁽²⁾_n from n = 7 (Theorem 2; R42: from n = 6, Theorem 12, §4.11, also for Z_d with d
    not a power of 2);
  * primes p < q from n = 2p + 3 (Theorem 2′);
  * Stab^(Z_d) ⊄ Stab^(Z_d′) for d ∤ d′ from n = 2·max(d, d′) + 3 (Theorem 5).
  * GF(4) and GF(2) never differ (§4.6).
* Open: n = 5; at n = 6 all pairs of odd primes. (The direction Stab⁽ᵖ⁾_6 ⊄ Stab⁽²⁾_6, p odd, listed here until R41,
  is Theorem 12, §4.11.) The eight characteristic-dependent seven-variable inequalities considered here are blind at
  n = 6 (Proposition 8); Theorem 9 goes through a derived family of subspaces, Theorem 12 through the isotropy of the
  stabilizer group.

**Theorem 7 (four parties) [proved; computer-assisted: lrs enumeration (G29b) and stored lifts (G33)].** For every
d ≥ 2,

    Stab^(Z_d)_4  =  h⁻¹(Shannon + Ingleton on five elements)  =  Stab⁽²⁾_4 .

In particular Stab⁽ᵖ⁾_4 is the same cone for every prime p. For n ≤ 4 the cones do not depend on d.

*Proof.*

* ⊆: subgroup rank functions satisfy Ingleton (common information, §4.6), so h_S ∈ Γ^(d)_5 satisfies Shannon and
  Ingleton by Theorem 1(d).
* ⊇: the right-hand cone has exactly 46 extreme rays (§3.2; lrs, gate G29b).
  * For each of them, `data/s30_stab4_unimodular.npz` stores an integer symmetric matrix W̃ with zero diagonal. It is
    the qubit certificate of §3.2 with signs on its edges, so W̃ ≡ W mod 2.
  * On every visible cut, the block W̃[X, X^c] has Smith invariants in {0, 1} and rational rank λ·r(X).
  * So its rank over GF(p) equals its rational rank for every prime p, and W̃ mod p is a GF(p) graph state realising
    λ·r. This was also re-verified with `stabcert` for p ≤ 13.
  * Hence every extreme ray lies in Stab⁽ᵖ⁾_4 for every p, and in Stab^(Z_d)_4 by Theorem 5 (take p | d).
* For n < 4, Stab_n is the face S(D) = 0 of Stab_4 (an empty party). If S(D) = 0, then Araki–Lieb gives
  S(XD) = S(X), so the face is a copy of the n = 3 cone. It is generated by the extreme rays of Stab_4 that lie on
  it, and their realisations have D isolated (zero-rank block), so they are product states. ∎

Most rays lift with all signs +1. Rays 36 and 39 needed an exhaustive search over signs, and rays 33, 41 and 45 a
random one (`s30_dimension_profile.py --n4-search`).

**Proposition 6 (one cut) [proved].** Let t ≥ 2, M = t + 1, and let F be a field whose characteristic ℓ does not
divide t − 1 (ℓ = 0 allowed). Let r_F be the rank function of the vectors (1, x) for the 2M + 2 points
x ∈ {e_j, c − e_j, c, 0} of F^M.

* Call X *balanced* if it contains exactly one of e_j, c − e_j for every j, and put a = #{j : c − e_j ∈ X}.
* Then r_F(X) = r_ℚ(X) − 1 if X is balanced and either
  * c ∈ X, 0 ∉ X, ℓ | t − a and a ≠ t, or
  * 0 ∈ X, c ∉ X, ℓ | a − 1 and a ≠ 1.
* For every other X, r_F(X) = r_ℚ(X).

For the generalised cube states (identically self-dual when ℓ ∤ t − 1, so S = 2r − |X|), S_F = S_ℚ − 2 exactly on
the cuts of these sets. Here S_ℚ is the formal vector 2r_ℚ − |X|; it equals S_F for every prime ℓ > t. In particular:

* S_F = S_ℚ if ℓ > t;
* if ℓ = t (t prime), S_F and S_ℚ differ exactly on the single cut {A_1, …, A_M, C} | {B_1, …, B_M, 0}, by 2.

So for primes p < q, the two witnesses of Theorem 2′ (t = p) differ in exactly one coordinate, the "one cut" of §4.4.
The same computation explains why the GF(3) state with t = 5 differs from the GF(7) state in 15 coordinates
(a = 2: C(6, 2) = 15).

*Proof.* A dependency λ ∈ F^X of the vectors (1, x), x ∈ X (λ_x = 0 for x ∉ X), satisfies Σ_x λ_x = 0, and for
every coordinate j

    λ_{e_j} + Σ_{i≠j} λ_{c−e_i} + λ_c = 0.

Write u_j = λ_{c−e_j}, β = Σ_j u_j and γ = −β − λ_c. The coordinate equations say λ_{e_j} = u_j + γ, and then the
first equation says λ_0 = (t − 1)·β + t·λ_c. So a dependency is determined by the u_j (for c − e_j ∈ X) and λ_c (if
c ∈ X), subject to:

* (i) u_j + γ = 0 for every j with e_j ∉ X (with u_j = 0 if c − e_j ∉ X);
* (ii) (t − 1)·β + t·λ_c = 0 if 0 ∉ X;
* (iii) λ_c = 0 if c ∉ X.

The nullity |X| − r_F(X) depends on F only through the rank of this system. Let K = {j : e_j, c − e_j ∈ X},
L = {j : only c − e_j ∈ X}, Z = {j : neither}, a = |L| and s = Σ_{j∈K} u_j.

* Z ≠ ∅. Then (i) forces γ = 0 and u_j = 0 on L, hence β = s = −λ_c, and (ii) becomes λ_c = 0. All coefficients
  are 0 or ±1 and the rank is the same over every field.
* Z = ∅. Then (i) gives u_j = −γ on L, so β = s − a·γ. The definition of γ becomes E1: (a − 1)·γ − s − λ_c = 0, and
  (ii) becomes E2: −a(t − 1)·γ + (t − 1)·s + t·λ_c = 0.
  * c, 0 ∈ X: only E1, rank 1.
  * c ∈ X, 0 ∉ X: if K ≠ ∅, E1 and E2 are independent over every field (proportionality would need
    −μ = t − 1 and −μ = t). If K = ∅ (X balanced), the unknowns are γ and λ_c, with determinant
    (a − 1)·t − a·(t − 1) = a − t. So the rank drops iff ℓ | t − a; over ℚ iff a = t.
  * 0 ∈ X, c ∉ X: λ_c = 0 and E1 reads (a − 1)·γ = s. This has rank 1 if K ≠ ∅. If K = ∅ (balanced), the rank drops
    iff ℓ | a − 1; over ℚ iff a = 1.
  * c, 0 ∉ X: λ_c = 0.
    * If K ≠ ∅, E2 + (t − 1)·E1 = −(t − 1)·γ, so the rank is 2 because ℓ ∤ t − 1.
    * If K = ∅, the only unknown is γ, with coefficients a − 1 and a(t − 1), which never vanish together (ℓ ∤ t − 1,
      and ℓ cannot divide both a and a − 1).

So the rank over F differs from the rank over ℚ exactly in the two balanced cases listed, where the nullity grows by
one. ∎

`s30_dimension_profile.py --one-cut` compares the predicted sets with the exact rank functions over GF(ℓ) and ℚ for
t = 2, …, 7 and ℓ ≤ 13 (gates G33/G33b); they agree in all 30 cases.

**Proposition 8 (six parties: eight characteristic-dependent inequalities are blind) [proved; computer-assisted:
exact certificates (G33)].**

Let D be the duality map h ↦ h*, with h*(X) = Σ_{P∈X} h(P) + h(E∖X) − h(E). This is the dual polymatroid of Kaced
(arXiv:1611.04109), who shows that the polymatroid and Ingleton regions are closed under D. D maps LR⁽ᵖ⁾ into itself:
the dual of a linear polymatroid is linear over the same field.

The eight seven-variable characteristic-dependent inequalities considered here are DFZ (65) and (91), Lemma 4 (LA)_2
and (LB)_2, Peña–Sarria (a) and (b) with t = 2, and BKL (6.8) and (6.17), the latter in its corrected form (§4.5;
the BKL pair was added in R36). For each of them, f + f∘D is a nonnegative rational combination of the Shannon
elemental inequalities on seven variables. The exact certificates are in `data/s30_selfdual_shannon_certs.npz`, with
supports of 54–114 inequalities and denominators 1–46 (BKL: 59 and 98, denominators 1 and 2). Since the Shannon
cone is invariant under relabelling and D commutes with relabellings, the same holds for every bijective assignment of
the seven variables to the seven parties.

* Hence f(h) = ½·(f + f∘D)(h) ≥ 0 for every self-dual polymatroid h on seven elements.
* For every pure quantum state on six parties plus a purifier, h_S is a self-dual polymatroid. It is submodular by
  strong subadditivity, monotone by Araki–Lieb (strong subadditivity alone suffices, since h_S(E) = h_S(E∖i)), and
  self-dual by purity.
  * So, pulled back to six parties in any relabelling, these eight inequalities follow from strong subadditivity.
    They hold for all states and cannot distinguish local dimensions.
  * Substitution instances, in which variables are unions of parties, are not covered by this argument. A random
    probe of 120,000 substitutions per form on 29 self-dual states found no violation (independent referee, R35; the
    six forms of R35).
  * Inequalities obtained by applying these forms to a *derived* family of subspaces are not covered either. Theorem 9
    (§4.8) is of this kind.
* A dimension dependence at n = 6 needs a characteristic-dependent seven-variable inequality g with g∘D = g. Indeed,
  if f is valid over GF(q) and f(h_S) < 0 for a GF(p) state, then g = f + f∘D is valid over GF(q) (D preserves
  LR⁽q⁾) and g(h_S) = 2f(h_S) < 0. Section 4.8 constructs one (g = f9 + f9∘D).

Encoding check (implied by the certificates, not needed for the proof): f(r + r*) never goes negative on all
{0,1}-point configurations of GF(p)³, or on 2.2 M random ones of GF(p)⁴ (p = 2 for (65) and (LB)_2; p = 3, 5 for
(91) and (LA)_2), although f itself reaches −1 there (`--probe6`). The eight-variable forms used at n = 7 are different:
the cube states violate them, so their symmetrisations are not Shannon.

**Five parties.**

* By Theorem 1, a dependence at n = 5 is equivalent to a characteristic dependence of the six-variable linear rank
  cone on its self-dual slice. No characteristic-dependent six-variable inequality is known.
* *Pilot (R35), superseded by §4.9.* The lifting method of Theorem 7, run on the qubit certificates of the first 60
  known extreme rays of Stab₅, lifted 40 of them (`s30_dimension_profile.py --n5-pilot 60`, no lifts stored).
  Section 4.9 (R38–R40) certifies all 641 for every d.

### 4.8 Six parties: qubits are not inside any odd local dimension [proved; refereed] (R36)

**Theorem 9 [proved; the expansion of f9 below is computer-assisted (G35)].** For every odd d ≥ 3 and every n ≥ 6,

    Stab⁽²⁾_n ⊄ Stab^(Z_d)_n .

In particular Stab⁽²⁾_6 ⊄ Stab⁽ᵖ⁾_6 for every odd prime p. Since Stab⁽²⁾ ⊆ Stab^(Z_d′) for even d′ (Theorem 5), also
Stab^(Z_d′)_n ⊄ Stab^(Z_d)_n for d′ even, d odd and n ≥ 6. Only this direction is shown: whether some state of odd
local dimension leaves Stab⁽²⁾_6 is open (see the last remark). (R42: it does; the non-Fano state, Theorem 12,
§4.11.)

*Labelling.* A, B, C, W, X, Y, Z = 100, 010, 001, 110, 101, 011, 111 (as in §4.1), with the seven Fano lines ABW,
ACX, BCY, AYZ, BXZ, CWZ, WXY. Every point lies on three lines.

*Witness.* The CSS state of the [7, 3] simplex code over GF(2), whose generator has the seven nonzero vectors of GF(2)³
as columns, with one qubit per party and any one party as the purifier. This is the logical |0⟩ of Steane's
[[7, 1, 3]] code (the even-weight half of the Hamming code; Nielsen–Chuang's |0_L⟩ up to the order of the qubits).
By Theorem 1 (step 2) its entropy is S = π(r_F), where r_F is the rank function of the Fano matroid. The Fano matroid has no coloops, so

    h_S = r_F + r_F* =: h_F ,

with values 2 on points, 4 on pairs, 5 on lines, 6 on the other triples and on the line complements, 7 on the other
4-sets and on every larger set. The entropies agree with the 128-amplitude state vector (G35).

*Base inequality.* Lemma 4 (LB) with t = 2 and (A₁, A₂, A₃, B₁, B₂, B₃, C) = (A, B, C, Y, X, W, Z):

    (LB)₂ = 9h(A) + 9h(B) + 9h(C) − 12h(ABC) + 6h(Z) + 3h(ABCZ) − 3h(ABZ) − 3h(ACZ) − 3h(BCZ) + 3h(AYZ) + 3h(BXZ)
          + 3h(CWZ) − 3h(AY) − 3h(BX) − 3h(CW) + 3h(BCY) + 3h(ACX) + 3h(ABW) + h(WXY)  ≥ 0 .

It holds for subspaces over every field of characteristic ≠ 2 and for subgroups of every finite abelian group of odd
order (Lemma 4(G)). On the Fano matroid it takes the value −1. Its positive-coefficient sets are the points A, B, C,
Z, the seven lines and ABCZ; the negative ones are AY, BX, CW, ABC, ABZ, ACZ and BCZ. Applied to h_F directly it
is useless: (LB)₂ ≥ 17 on h_F in every relabelling (Proposition 8). The proof transfers it to a derived family.

*Derived family.* Let U_P (P ∈ E) be subspaces with rank function h. For a line L = {i, j, k} through i put

    F_i^L = U_i ∩ (U_j + U_k),        F̂_i = Σ_{L ∋ i} F_i^L ⊆ U_i ,

and let ĥ be the rank function of the family (F̂_P). Write m_i^L = I(i; L∖i) = h(i) + h(L∖i) − h(L).

* **Lemma 0.** dim F_i^L = m_i^L.
* **Lemma A.** F_i^L ⊆ F_j^L + F_k^L. If v = a + b ∈ U_i with a ∈ U_j and b ∈ U_k, then a = v − b ∈ U_j ∩ (U_i + U_k)
  = F_j^L, and likewise b ∈ F_k^L.
* **Lemma B (points).** Let L₁, L₂, L₃ be the lines through i and Q_ab = (L_a ∪ L_b)∖{i}. Then
  F_i^{L_a} + F_i^{L_b} ⊆ U_i ∩ U_{Q_ab}, which has dimension I(i; Q_ab). By submodularity,
  ĥ(i) ≤ I(i; Q_ab) + I(i; Q_ac) − m_i^{L_a} for each a. Averaging over a gives
  ĥ(i) ≤ u_i := ⅓·[2·Σ_{a<b} I(i; Q_ab) − Σ_a m_i^{L_a}].
* **Lemma C (lower bounds).** ĥ(i) ≥ l_i := ⅓·Σ_L m_i^L, and ĥ(X) ≥ h(X) − Σ_{i∈X} (h(i) − l_i) for every X, since
  Σ_{i∈X} U_i is spanned by Σ_{i∈X} F̂_i together with complements of the F̂_i in the U_i.
* **Lemma D (lines).** For a line L = {i, j, k}: ĥ(L) ≤ u_i + u_j + u_k − m_i^L, because F_i^L ⊆ F̂_i ∩ (F̂_j + F̂_k)
  by Lemma A. f9 uses the first point of L in the order A, B, C, W, X, Y, Z.
* **Lemma E.** ĥ(ABCZ) ≤ ĥ(ABCYZ) ≤ u_A + u_B + u_C + (u_Y − m_Y^{BCY}) + (u_Z − m_Z^{AYZ}). Add F̂_Y to
  F̂_A + F̂_B + F̂_C, then F̂_Z, using Lemma A on the lines BCY and AYZ.

These lemmas hold for subspaces over any field. With log |·| in place of dim they also hold for subgroups of any
finite abelian group, because the subgroup lattice is modular.

*The inequality f9.* The derived family lies over the same field (or in the same group). So (LB)₂(ĥ) ≥ 0 in
characteristic ≠ 2, or for groups of odd order. Bound every positive term from above (Lemmas B, D, E) and every
negative term from below (Lemma C). This gives 0 ≤ (LB)₂(ĥ) ≤ f9(h)/3 with the integer form (67 terms,
`data/s31_six_party.npz`)

    f9 = + 36·h(A) + 45·h(B) + 54·h(C) + 18·h(W) + 21·h(X) + 21·h(Y) + 45·h(Z) − 10·h(AB) − 10·h(AC) − 22·h(BC)
         − 42·h(AW) − 51·h(BW) − 36·h(CW) − 42·h(AX) − 36·h(BX) − 51·h(CX) − 13·h(WX) − 45·h(AY) − 42·h(BY)
         − 51·h(CY) − 10·h(WY) − 13·h(XY) − 13·h(AZ) − 10·h(BZ) − 10·h(CZ) − 51·h(WZ) − 51·h(XZ) − 51·h(YZ)
         − 36·h(ABC) + 103·h(ABW) + 103·h(ACX) + 115·h(BCY) + 36·h(WXY) − 9·h(ABZ) − 9·h(ACZ) − 9·h(BCZ) + 88·h(CWZ)
         + 88·h(BXZ) + 100·h(AYZ) + 98·h(BCWX) + 92·h(ACWY) + 92·h(ABXY) + 48·h(ABCZ) + 104·h(AWXZ) + 98·h(BWYZ)
         + 98·h(CXYZ) − 42·h(ABCWX) − 42·h(ABCWY) − 42·h(ABCXY) − 14·h(ABWXY) − 14·h(ACWXY) − 20·h(BCWXY)
         − 14·h(ABCWZ) − 14·h(ABCXZ) − 42·h(ABWXZ) − 42·h(ACWXZ) − 36·h(BCWXZ) − 20·h(ABCYZ) − 42·h(ABWYZ)
         − 36·h(ACWYZ) − 42·h(BCWYZ) − 36·h(ABXYZ) − 42·h(ACXYZ) − 42·h(BCXYZ) − 20·h(AWXYZ) − 14·h(BWXYZ)
         − 14·h(CWXYZ)  ≥ 0 .

f9 is linear, so it holds on the closed cones LR⁽q⁾₇ (q odd) and Γ^(d)₇ (d odd).

*Value.* f9(h_F) = −3. This can be checked by hand: at h_F every bound is attained (m_i^L = 1, u_i = l_i = 1, the
line bounds are 2, the ABCZ bound is 3, and the lower bounds equal the Fano ranks), so f9(h_F)/3 = (LB)₂(r_F) = −1.
Also f9(r_F) = −3, while on the non-Fano configuration f9 = 81 (h) and 33 (r).

*Proof of Theorem 9.* For a Z_d stabilizer state on six parties, h_S is a subgroup rank function on seven elements
in an abelian group of exponent dividing d, hence of odd order (Theorem 1(d), steps 1–2; Theorem 1 for d prime). So
f9(h_S) ≥ 0, and by linearity and closedness f9 ≥ 0 on h(Stab^(Z_d)_6). The simplex-code state has f9(h_S) = −3.
For n > 6, merge the extra (empty) parties into the purifier. ∎

The result inherits the status of Theorem 1(d), steps 1–2 (the cited entropy formula of Gross–Walter).

**Remarks.**

* *Beyond Proposition 8.* g = f9 + f9∘D is a self-dual characteristic-dependent inequality, as Proposition 8 requires:
  g(h_F) = −6 while h_F is a polymatroid, so g is not a Shannon consequence. Hence f9 is not a nonnegative
  combination of Shannon inequalities and the eight forms of Proposition 8 in any relabelling (that would make g a
  Shannon consequence). On h_F the derived family sees only the Fano part: I_{r_F*}(i; L∖i) = 0 on every line, so the
  dual half of r_F + r_F* drops out of the F_i^L.
* *Invariant form.* On Fano-invariant self-dual vectors (a four-dimensional slice; self-duality forces v_E = 7v_P/2)
  write v_P for the value on points, α = 2v_P − h(pair) (pairwise mutual information), μ = v_P + h(pair) − h(line)
  (= I(i; L∖i)), ν = v_P + h(pair) − h(non-line triple), and β = v_P − 2μ ≥ 0. There
  f9/3 = 210·v_P + 30·α + 21·ν − 421·μ = 210·β + 30·α + 21·ν − μ. So every Fano-invariant self-dual vector of
  Γ^(d)₇, d odd, satisfies μ ≤ 210·β + 30·α + 21·ν, while h_F has α = β = ν = 0 and μ = 1. With DFZ (65) as the base
  form, and every positive-coefficient set bounded by the same credit rule (Lemma A along a suitable order of a
  superset), the constants become 154, 26 and 5 (`s31_six_party.py --invariant`; computed, not refereed).
* *Sanity, not proofs.* With the actual subspaces, every bound of Lemmas 0, A–E and the whole chain are attained with
  equality on many random families over GF(2), GF(3) and GF(5); over GF(2), f9 reaches −3 on the Fano configuration
  (G35, G35b).
* *Referee.* An independent referee (separate agent, own code) re-derived Lemmas 0, A–E and the composition and
  recomputed f9 exactly (no difference). It confirmed f9(h_F) = −3 from r_F + r_F* and from the state vector, and
  checked Theorem 1(d) and Lemma 4(G) for odd composite d. Stress tests found no violation of f9 ≥ 0:
  * about 8.2 M arrangements each over GF(3), GF(5), GF(7);
  * all 105,413,504 seven-tuples of points of PG(2, 3) or zero (minimum 0);
  * about 13,000 subgroup families of Z_9², Z_3² × Z_9, Z_9³, Z_27², Z_3 × Z_9², Z_9 × Z_27, Z_5 × Z_25 and Z_15².

  Over GF(2), f9 < 0 exactly on the 168 Fano frames of GF(2)³ (minimum −3). Forms claimed only for other
  characteristics, and deliberately weakened variants of f9, do go negative in the same searches. Its wording fixes
  were applied: one direction only, the scope of Proposition 8, the dependence on Theorem 1(d), and the invariant
  form stated on the slice.
* *Faces [certified, exact ranks; `s31_six_party.py --face`, G35].* The Shannon elemental inequalities tight at h_F
  have rational rank 125, and those tight at r_F and at r_F* rank 126 (upper bounds from the vectors they annihilate,
  lower bounds from the rank over GF(1000003)). So the smallest face of Γ₇ containing h_F is two-dimensional,
  and r_F, r_F* are extreme rays of Γ₇ lying in it: the face is cone(r_F, r_F*). Every decomposition of h_F into
  polymatroids stays in this face. Since D preserves LR⁽q⁾₇ and swaps r_F and r_F*, the key step h_F ∉ LR⁽q⁾₇ of
  Theorem 9 is equivalent to: LR⁽q⁾₇ contains no nonzero vector a·r_F + b·r_F* (q odd). The same holds verbatim for
  the non-Fano pair r_NF, r_NF* (ranks 125 and 126).
* *The reverse direction at n = 6 (R41: settled without closure; R42: settled, Theorem 12, §4.11).* The natural
  candidate is the non-Fano state, the same code over GF(3), with h_NF = r_NF + r_NF*. By the previous remark,
  h_NF ∈ LR⁽²⁾₇ if and only if LR⁽²⁾₇ contains a nonzero vector a·r_NF + b·r_NF*.
  * The same construction with (LA)₂ gives +27 on h_NF. The points W, X, Y lie on only two non-Fano lines each, and
    I(W; ABCZ) = 2, so for a general subspace family their two line pieces need not coincide. (R41: the analogous
    objects of a stabilizer state with entropy vector m·S_NF, the restrictions to W of its elements supported on the
    two lines, must coincide, by isotropy at the other two points of each line (Proposition 10, step 1′, §4.10);
    over GF(2) there is no such state.)
  * No family of 2-dimensional subspaces over GF(2) realises h_NF exactly (exhaustive search, not in the
    repository). (R41: a special case of Corollary 11.)
  * An exploratory hand argument (not refereed) suggested that an exact GF(2) family with rank function
    a·r_NF + b·r_NF* must have a = b. (R41: superseded; no exact family with (a, b) ≠ (0, 0) exists.)
  * (R40) A qubit state with entropy vector m·S_NF gives the exact point a = b = m: the code r_L of the doubling remark
    (§2), with S(P) qubits per system, has rank function h_{m·S_NF} = m·h_NF. Conversely, an exact point a = b = m
    gives a qubit state with 2m·S_NF (step 2 of §2). For m = 2, 3 annealing found no such state (§4.10) [empirical].
  * (R41) There are no exact points: no qubit stabilizer state has an entropy vector that is a positive multiple of
    S_NF, and no GF(2) family has rank function a·r_NF + b·r_NF* with (a, b) ≠ (0, 0) (Proposition 10, Corollary 11,
    §4.10) [proved]. So h_NF can lie in LR⁽²⁾₇ only as a limit, and the question of this remark is now about limits
    only.
    R40 had written instead that in any qubit realisation the line pieces at one of W, X, Y differ; that was true
    only vacuously (erratum in §4.10).
  * (R42) h_NF ∉ LR⁽²⁾₇, and S_NF ∉ Stab⁽²⁾_6: the characteristic-2 inequality g12 = f12∘π is valid on LR⁽²⁾₇ and
    g12(h_NF) = −2 (Theorem 12, §4.11) [proved]. Its transfer uses the isotropy of a stabilizer group, which a
    general subspace family lacks; that is why the rank-only transfer above gives +27.
  * Among the states over ℤ_3 with at most one qudit per system, h_NF is the only remaining candidate for this
    direction (§4.10) [certified].

### 4.9 Five parties: the known extreme rays over every local dimension [certified; independent referee] (R38–R40)

Theorem 7 lifts every four-party extreme ray to all local dimensions with one integer matrix. This section does the
same for the known five-party extreme rays, as far as the search reaches.

*Definition.* Let r be a five-party vector (31 coordinates) and κ ≥ 1. A **unimodular realisation** of κ·r is an
integer symmetric matrix W on N = κ·(r(A) + … + r(F)) vertices, system q owning κ·r(q) of them, with zero entries
inside each system (in particular a zero diagonal), such that for every nonempty X ⊆ {A, …, E} the cut block
W[Q_X, Q_{N∖X}] has Smith invariants in {0, 1} and rational rank κ·r(X).

*Lemma (every d at once) [proved; independent referee].* If W is a unimodular realisation of κ·r, then for every d ≥ 2 the graph state of
ℤ_d qudits with weight matrix W mod d has entropy vector κ·r in units of log d. Hence r ∈ Stab^(ℤ_d)_5 for every d,
as an entropy vector (no closure).

*Proof.* For a graph state of ℤ_d qudits, S(X) = log |W[Q_X, Q_{N∖X}] ℤ_d^{N∖X}| (paper, Appendix A 1). Write the
block as U·D·V with U, V unimodular over ℤ and D = diag(1, …, 1, 0, …, 0) with κ·r(X) ones. U and V stay invertible
mod d, so the image is isomorphic to ℤ_d^{κ r(X)}. ∎ (For prime d this is the argument of Theorem 7; for composite d
it avoids Corollary 3 of the paper.)

*Result [certified: exact Smith forms by two independent methods, gates G36, G36b, G37; independent referee].*
`data/s32_stab5_unimodular.npz` stores a unimodular realisation for **566 of the 641** known extreme-ray orbits
of Stab₅ (`data/stab5_extreme_reps.npy`): κ = 1 for 515, κ = 2 for 50, κ = 3 for 1; at most 40 vertices.
`data/s32_stab5_cover.npz` covers the other **75** rows by covering certificates (below). So **all 641** lie in
Stab^(ℤ_d)_5 for every d ≥ 2 (R38: 577; R39 added 2 unimodular and 57 covered rows, R40 the last 5 covered rows).

* Each of these 641 rays lies in Stab^(ℤ_d)_5 for every d ≥ 2. Since Stab^(ℤ_d)_5 ⊆ QLR₅ for every d (paper,
  Sec. IV A) and the rays are extreme rays of QLR₅, they are extreme rays of Stab^(ℤ_d)_5 for every d.
* The 12 rays of the BCHS ledger that had no explicit realisation (rows 6, 7, 15, 34, 39, 40, 42, 43, 44, 46, 48, 50;
  their membership in the stabilizer cone rested on hypergraph realisations and the theorem of Walter and Witteveen,
  in the closure and up to scaling) all have one now, at κ = 2. So **every one of the 641 known extreme rays has an
  explicit qubit realisation** (W mod 2 for these twelve), and κ = 2 is the least multiple for qubits at rows 43 and
  44 (no realisation at multiplicity 1 over GF(2) or GF(3), exhaustive search of R31).
* Every row also has an explicit qutrit realisation in the data: a unimodular certificate mod 3, or a covering piece
  that is right at 3. So all 641 rows have one (96 before R38, from the stored qutrit certificates).

*Covering certificates [proved; certified by G37; independent referee].* One matrix for all primes is not needed.

* Call an integer symmetric matrix W, zero inside each system, a *rational realisation* of κ·r if every cut block has
  rational rank κ·r(X). Let bad(W) be the set of primes dividing a Smith invariant of some cut block; it is finite and
  computed exactly. For every prime p ∉ bad(W), W mod p realises κ·r over GF(p): a block of rational rank t keeps
  rank t mod p exactly when p does not divide its t-th determinantal divisor.
* Suppose rational realisations W_1, …, W_m (of multiples of r) and realisations over GF(p_1), …, GF(p_l) satisfy
  bad(W_1) ∩ … ∩ bad(W_m) ⊆ {p_1, …, p_l}. Then r ∈ Stab⁽ᵖ⁾_5 for every prime p. By Corollary 3 of the paper
  (Stab⁽ᵖ⁾ ⊆ Stab^(ℤ_d) when p divides d), r ∈ Stab^(ℤ_d)_5 for every d. A unimodular realisation is the case
  bad(W) = ∅.
* `data/s32_stab5_cover.npz` stores such certificates for the 75 rows that have no unimodular realisation here:
  125 rational realisations (3 of them at κ = 2), with their bad primes, and 74 realisations over
  GF(2) or GF(3) (the stored qubit and qutrit certificates of those rows). Gate G37 recomputes every bad-prime set and
  checks coverage.
* Examples. Row 456 (R38): two rational realisations with bad primes {2, 3} and {2, 13}, and the qubit certificate.
  Row 607 (R39): the matrix of the rejected report below (Smith invariant 127, so bad primes {127}) and a rational
  realisation with bad primes {2, 3, 5}; the intersection is empty, so no GF piece is needed.
* The search (`--cover`, R39 order) runs rounds over the rows not yet covered. A row without a rational realisation
  gets a rational search (residual modulo 2^31 − 1 only) from its qubit certificate; a row with one gets a search that
  is rational and right modulo all its uncovered primes at once, started from its rational realisation with the fewest
  uncovered primes. R38 started the latter from the qubit certificate, one prime at a time.

*Method.* `cc/s32_unimodular.c` anneals the inter-party entries over {-1, 0, 1}. The residual is the sum over the
31 cuts of |rank − κ·r(X)| over GF(2), GF(3) and modulo 2^31 − 1 (a stand-in for ℚ). When it reaches 0, the ranks are
checked for every prime below 1000 and for 2^61 − 1, and a prime that fails joins the residual. A matrix reported
this way is still not certified: one report (row 607) had the Smith invariant 127 and was rejected by the exact check.
The search proves nothing; `scripts/s32_unimodular.py` certifies.

* *Starts.* The stored qubit certificate of the row (all signs +1, then random signs), a stored qutrit certificate,
  a fresh GF(3) realisation found by a GF(3)-only search, or random matrices.
* *Neutral moves [proved; independent referee; tested by `--self-test`].* Two moves keep every cut block in its Smith class, so they never change
  the residual; they only move the search to another representative.
  * A party-local unimodular basis change, W ↦ EᵀWE with E = I ± e_u e_{u′}ᵀ, u ≠ u′ in one system.
  * An integer local complementation at a vertex c, W ↦ W − a·w_c w_cᵀ (a = ±1, w_c the column of c), followed by
    clearing the diagonal and the entries inside each system. It is induced by the symplectic shear
    (x_c, z_c) ↦ (x_c + a z_c, z_c) on c and party-local gates; W[c][c] = 0 is used.

  On a random sample of 20 hard rows, at 10 s per row, the number certified rose from 1 without these moves to 8 with
  them (one step in five) and to 11 (one in two, longer anneals).
* *Restriction of scalars [proved; independent referee].* Let R = ℤ[θ] with θ² = sθ + t. The trace form B(x, y) = Tr(xy/(2θ − s)) is
  integral and unimodular on R: by Euler's formula its Gram matrix in the basis {1, θ} is G = [[0, 1], [1, s]], of
  determinant −1. Multiplication by v is self-adjoint for B, so G·M(v) is symmetric, where M(v) is the matrix of
  multiplication by v (its columns are the coordinates of v·1 and v·θ; with the transposed convention the product is
  not symmetric in general). For a symmetric matrix V over R with zero entries inside each system, the integer matrix W
  with 2 × 2 blocks G·M(V_uv) is symmetric, and its cut blocks are (I ⊗ G)·M(V_X). If V_X has unit Smith invariants
  over R (a principal ideal domain for the five rings used), so does M(V_X) over ℤ, and the rational rank doubles.
  The search does not rely on this: it certifies W directly. Five rings are available: ℤ[ω] (Eisenstein), ℤ[i],
  ℤ[(1 + √5)/2], ℤ[√2] and ℤ[√−2].
  * ℤ[ω] fits the rays that qubits realise only at multiplicity 2: ℤ[ω]/(2) = GF(4) and ℤ[ω]/(√−3) = GF(3). Started
    from a qutrit certificate (y = 0, so the residual mod 3 is 0 at once), it certified at κ = 2, most within
    seconds, 17 of the 21 rows it was run on: rows that were open at the time and have a qutrit
    certificate at multiplicity 1 but no qubit certificate at multiplicity 1. Row 40, the entropy vector of an absolutely maximally entangled state of four
    parties (S = 1 on singles and 2 on pairs; no qubit state has it), was found at once.
  * Reduced at the prime √−3, a ℤ[ω] certificate at κ = 2 gives a qutrit realisation at multiplicity 1. For row 305
    this was the first one: GF(3)-only searches had not found it.
* *Strategies used* (`--strategies`; counts are the certificates each strategy contributed in this round):
  * `plain` (from the stored qubit certificate, one entry at a time): 397 (κ = 1: 387, 2: 9, 3: 1)
  * `neutral` (the same with neutral moves): 129 (κ = 1: 126, 2: 3)
  * `qutrit-start` (from a stored qutrit certificate (Eisenstein ring, or neutral at κ = 2)): 20 (κ = 2: 20)
  * `ring-eisenstein` (ℤ[ω] from a qubit certificate at multiplicity 1 or from random matrices): 14 (κ = 2: 14)
  * `gf3-search-start` (ℤ[ω] from a fresh GF(3) realisation): 3 (κ = 2: 3)
  * `ring-gauss` (ℤ[i] from a qubit certificate): 1 (κ = 2: 1)
  * R39 added rows 539 and 631 (`neutral`, κ = 1), found by R38's strategy probe (below).

*Round R39 [certified; independent referee].* R39 closed 59 of the 64 rows that R38 left open.

* *Search order.* The covering search now starts each search "rational and right mod q" from a rational realisation
  of the row, the one with the fewest uncovered primes, and asks for all uncovered primes at once. For 17 of the 21
  rows whose rational realisations from R38's runs left primes uncovered, every one of those was bad at 3. On four of
  them, searches over GF(3) at multiplicity 1 and over GF(9) (ℤ[i], κ = 2) from random starts found nothing in 30 s,
  while this search found a rational realisation right mod 3 for all four
  (`results/2026-10-14_r39-sandbox/probes/calibration_p3.txt`). It covered all 21 within a few rounds of 40 s.
* *First rational realisations.* Rational searches from the qubit certificate, 90 to 900 s per row (R38's covering
  run: 15 s per round), found one for 34 more rows; the search above then covered them. On the 21 rows without one at
  the time, ring searches (ℤ[ω], ℤ[i]; 60 s) and searches with entries in {−2, …, 2} (90 s) found none. Row 564
  needed 300 s: its first four rational realisations were all bad at 3. Row 200, which qubits realise only at
  multiplicity 2, has rational realisations at κ = 2 from its qubit certificate at multiplicity 2; at κ = 1 the best
  attempts miss one cut.
* *Outputs of R38 that R38 did not use.* Rows 539 and 631 have unimodular realisations in the output of R38's
  strategy probe, which was never merged into the data. Row 574 is covered by two rational realisations from R38's
  runs, at κ = 2 with bad primes {163, 181, 2833, 7699} (a report of the probe) and at κ = 1 with bad primes
  {2, 3, 5}; R38's covering run had not combined them. Row 607 is covered by the matrix of R38's rejected report
  (bad primes {127}) together with a rational realisation found in R39 (bad primes {2, 3, 5}).
* *Erratum (R38 probe file).* `results/2026-10-13_r38-sandbox/probes/strategy_sample20.txt` counted kernel reports as
  certified. The report for row 574 under "ring Z[omega] from the qubit certificate" is not unimodular (bad primes
  above; the kernel then checked primes up to 97 only), so that setting certified 2 of 20 (rows 465, 541). The other
  reports in the file are exact. A dated note is added to the file; nothing in R38's data or claims rests on it.

*Round R40 [certified; independent referee].* The five rows that R39 left open, 587, 589, 612, 630 and 635
(N = 24, 24, 30, 26, 28), are covered (`results/2026-10-14_r40-sandbox/`).

* First rational realisations: rows 587 and 589 from rational searches of 1200 s, started from the qubit
  certificate and from R39's best near miss; rows 630 and 635 from 600-s rounds started from the best near miss at
  start temperature 0.3 instead of 1.0; row 612 in the second such round of a later push.
* Searches that are rational and right mod the uncovered primes, started from these realisations, then covered all
  five (60 s each; row 612: 300 s, third round). Row 630 has a rational realisation with bad primes {2}; each of the
  others has two whose bad primes meet in {2}. The stored qubit certificate supplies GF(2).
* A separate agent, with its own code, re-verified the 14 new pieces (structure; ranks on all 31 cuts; bad primes by
  its own Smith form, by minors and by sympy) and the coverage. It found the 185 earlier pieces and the unimodular
  file unchanged, and every row 0–640 certified one way or the other. Not yet checked by hand.

*The 5 rows open after R39 (closed in R40) [empirical].*

* R39 knew no rational realisation for any of them. Its best attempts, with the cuts they miss, are in
  `results/2026-10-14_r39-sandbox/open_rows.json`. They missed one to three cuts (row 612: ten or more), mostly with
  the rational rank above the target, which is the rank over GF(2) of the qubit certificate. A configuration that
  exists only in characteristic 2 would produce such misses, and so does a search that has not finished; R39 could
  not tell the two apart. R40 found that the searches had not finished: longer runs from the near misses succeeded.
  Their qutrit realisations come from the covering pieces that are right at 3.
* (R38) At this size (N = 18–30), searches modulo a single small odd prime are hard, while searches for rational
  realisations (residual modulo 2^31 − 1 only) often succeed within seconds; the covering certificates exploit this.
  For example, a GF(3)-only search found no qutrit realisation of row 636 at multiplicity 1 (60 s), but its first
  rational realisation has bad primes {2, 5, 7}: its reduction mod 3 is one
  (`results/2026-10-13_r38-sandbox/probes/badprimes_example.txt`).
* *Odd-characteristic screen (exploratory; `scripts/s33_ci_oddchar.py`).* On the seven-element extension by the
  common information of a pair, the screen imposes the forms of Proposition 8 that hold in odd characteristic:
  DFZ (65), (LB)₂ and the corrected BKL (6.17), in all 5,040 relabellings. It excluded none of the 77 rows
  that were open before the covering run.
  * Pairs: every disjoint pair (about 280 per row).
  * Verdicts: floating-point LP; every pair was feasible.
  * A row excluded this way would lie outside every odd-characteristic cone. That would be a five-party dimension
    dependence and, by duality, a six-variable linear rank inequality that holds in odd characteristic and fails in
    characteristic 2; none is known.
  * (R39) With every pair (`--pairs all`, about 1,330 per row) on the 22 rows open at the time
    (`results/2026-10-14_r39-sandbox/ci_all_pairs.log`): no exclusion either.

*Independent referee (R38).* A separate agent, with its own code and without reading the repository's verification
code, re-verified all 564 unimodular certificates (structure, rational ranks and Smith invariants on all 31 cuts,
compared with sympy's Smith normal form; ranks over GF(p) for p ≤ 13), all 36 covering pieces and the coverage of all
13 rows; simulated the graph states of 331 (certificate, d) pairs for d = 2, 3, 4, 6 (entropies equal to κ·r·log d
within 3·10⁻¹⁴); and checked the lemma, the covering argument, the ring construction and the neutral moves (their
hypotheses u ≠ u′ and W[c][c] = 0, and the column convention of M(v), are now stated). No discrepancy. Not yet
checked by hand.

*Independent referee (R39).* A separate agent, with its own code and without reading the repository's verification
code, checked that the R38 entries of both data files are unchanged and re-verified every new entry: the unimodular
certificates of rows 539 and 631, the 94 new rational realisations with their exact bad primes (its own Smith form,
Bareiss ranks with minors, and sympy's Smith normal form agree), the 55 new GF(2) or GF(3) pieces, and the coverage
of all 70 rows. No discrepancy; the numbers in this section agree with the data. Its wording remarks are taken into
the text. Not yet checked by hand.

*What this shows.* Every certified ray is independent of the local dimension. This supports, but does not prove, that
Stab^(ℤ_d)_5 is independent of d (the main open question of the paper):

* Stab₅ itself is not known: the catalogue has 1,155 undecided orbits and the enumeration of QLR₅ is incomplete.
* Since R40 no known ray remains open.

### 4.10 Six parties, the other direction: states with at most one qudit per system [certified; Proposition 10 and Corollary 11 proved; empirical for h_NF] (R40, R41; the closed cones: §4.11, R42)

Question (§5): is Stab⁽ᵖ⁾_6 ⊄ Stab⁽²⁾_6 for some odd p? The only candidate so far is the non-Fano state
h_NF = r_NF + r_NF* (§4.8). This section surveys all stabilizer states of seven systems (six parties and the
purifier) with at most one qudit per system, for p = 2 and p = 3. It also shows that no qubit stabilizer state has an
entropy vector that is a positive multiple of that of the non-Fano state, and that S_NF is not a finite nonnegative
combination of qubit stabilizer entropy vectors (Proposition 10, Corollary 11; R41).

*Enumeration [certified: exhaustive; the GF(2) list recomputed independently, gates G38, G38b].*

* Every stabilizer state over ℤ_p (p prime) is local-Clifford equivalent to a weighted graph state. For qubits this is
  Theorem 1 of Van den Nest, Dehaene and De Moor (quant-ph/0308151), who attribute the d-level case to Schlingemann
  (quant-ph/0111080). Their constructive proof works verbatim over GF(p): local Fourier transforms make the X-part of
  a generator matrix invertible, the generator matrix then becomes (I | A) with A symmetric, and local phase gates
  clear the diagonal.
* A party-local Clifford leaves S(P) qudits per system. The stabilizer elements supported on P form an isotropic
  subspace of dimension n_P − S(P); a Clifford on P maps it to Z-type elements on n_P − S(P) of its qudits, which then
  factor out.
* Scaling a vertex by a unit (row and column times u) changes no cut rank.
* So the entropy vectors of the states with every S(P) ≤ 1 (units of log p) are exactly the cut-rank vectors of the
  weighted graphs over GF(p) on seven vertices, one per system.
* `cc/s35_enum7.c` enumerates these graphs up to scaling. Normal form: the edge from each vertex to its smallest
  earlier neighbour has weight 1; every scaling orbit contains such a graph.
  * p = 2: all 2^21 graphs, 10,773 distinct vectors.
  * p = 3: 255,602,200 graphs in normal form, 27,334 distinct vectors.
  * 10,743 vectors are in both lists.
* Checks: the GF(2) list is recomputed in numpy from all 2^21 graphs, independently of the C code (G38b); 20,000
  random GF(3) graphs, with arbitrary weights, all land in the GF(3) list (G38b).

*Result.*

* **Qubit-only** (in the GF(2) list, not in the GF(3) one): 30 vectors, the relabellings of the Steane state
  [certified]. By Theorem 9 they lie outside Stab^(ℤ_d)_6 for every odd d.
* **Qutrit-only**: 16,591 vectors in 53 orbits of S₇ (the seven systems relabelled, the purifier among them)
  [certified]. In 41 orbits every S(P) = 1; in the other 12 some system is unentangled, so they are states of at most
  five parties.
* For **52 of the 53 orbits**, `cc/s34_anneal7.c` found a qubit graph state with two qubits per qutrit realising 2·S,
  most within seconds (60-s budget per orbit). These vectors lie in Stab⁽²⁾_6 [certified: exact ranks over GF(2),
  G38; state vectors in the self-test]. None of these states is a tensor product of two states with one qubit per
  system: for no orbit is 2·S a sum of two vectors of the GF(2) list [certified: all pairs, G38b].
* The remaining orbit (210 vectors) is that of **h_NF**. No qubit state with entropy vector m·S_NF was found for
  m = 2 or 3 [empirical]. (R41: none exists, for any m; Proposition 10.)
  * Search class: m qubits per system, graph form, no edges inside a system. This loses no generality: the
    reductions above, and party-local controlled-Z gates remove the edges inside a system. The class has 2^84 graphs
    for m = 2 and 2^189 for m = 3.
  * The *residual* of a run is the smallest L1 distance, Σ |cut rank − target| over the 63 cuts, that it reached.
  * m = 2: best residual 2 in each of seven runs: the 60-s survey run, two runs of 120 s, and four of 600 s (two from
    random starts, two from two copies of the best attempt at m = 1).
  * m = 3: best residual 3 in each of three runs: two of 120 s and one of 600 s from three copies of the best attempt.
  * Baseline [certified, G38b]: S_F and S_NF differ on one cut (2 against 3), and S_F is the only vector of the GF(2)
    list at distance 1 from S_NF; no sum of two GF(2) vectors is closer than 2 to 2·S_NF. So m copies of a Steane
    state already have residual m, and no run got below this baseline.
  * Controls: 2·S_Steane was found at once from a random start, and 3·S_Steane in 600 s from a random start and from
    copies (one 120-s run ended at residual 2). Copies of the Steane state solve these targets, so they are weak
    controls. The stronger control is the survey above: its 52 realisations are not products of two states with one
    qubit per system (other tensor splits are not excluded).

*Field-linear states, and combinations of qubit states.*

**Proposition 10 [proved] (R40; strengthened in R41).** Use the labels of §4.8: the non-Fano lines are ABW, CWZ, ACX,
BXZ, AYZ, BCY. The points A, B, C, Z lie on three lines each, the points W, X, Y on two (WXY is the seventh line of
the Fano plane, not a line here). Every line contains two of A, B, C, Z. Let F = GF(q), and let L be an F-linear
stabilizer state of the seven systems with m F-qudits each (a Lagrangian subspace of ⊕_P F^{2m}) whose entropy
vector, in units of log q, is m·S_NF. For a line ℓ let L_ℓ be the elements of L supported in ℓ, and π_P the
restriction to system P.

* (a) At every point P, the restrictions π_P(L_ℓ), ℓ ∋ P, are all equal to one Lagrangian subspace D_P of F^{2m}.
  (R40 proved this at A, B, C, Z; R41 adds W, X, Y.)
* (b) Hence the non-Fano matroid has an m-fold (multilinear) representation over F, so q is odd.
* (c) A GF(q)-linear stabilizer state whose entropy vector (units of log q, any number of qudits per system) is a
  positive multiple of S_NF exists if and only if q is odd. In particular, no qubit stabilizer state has an entropy
  vector that is a positive multiple of S_NF, and no GF(2^k)-linear state has entropy vector k·S_NF in units of
  log 2.

*Proof.* Let L_X be the elements of L supported in X ⊆ E. The doubling remark of §2 gives
dim π_X(L) = m|X| + S(X); over GF(q) this follows from GF(p) by restriction of scalars. With purity,
dim L_X = dim L − dim π_{E∖X}(L) = m|X| − S(X). For S = m·S_NF this means:

* L_P = L_{PP′} = 0 for points P, P′;
* dim L_ℓ = m for each line ℓ, and π_P is injective on L_ℓ for P ∈ ℓ (its kernel is L_{ℓ∖P} = 0).

1. If lines ℓ ≠ ℓ′ meet in P, then for s ∈ L_ℓ and s′ ∈ L_ℓ′ isotropy gives ω_P(π_P s, π_P s′) = ω(s, s′) = 0,
   since the supports meet only in P. So π_P L_ℓ′ lies in the symplectic complement (π_P L_ℓ)^⊥, and equals it by
   dimension (m each). At a point on three lines the three restrictions are pairwise complements of each other, so
   they are all equal to one subspace D_P with D_P = D_P^⊥. This is (a) at A, B, C, Z.
1′. (R41) Let ℓ = {P, P′, D} with P, P′ among A, B, C, Z. For s, s′ ∈ L_ℓ, isotropy gives
   0 = ω(s, s′) = ω_P(π_P s, π_P s′) + ω_{P′}(π_{P′} s, π_{P′} s′) + ω_D(π_D s, π_D s′). The first two terms vanish,
   because π_P s, π_P s′ lie in the Lagrangian D_P, and likewise at P′. So π_D L_ℓ is isotropic of dimension m, hence
   Lagrangian. By step 1 the restriction of the other line through D is its complement, that is, π_D L_ℓ itself. This
   is (a) at W, X, Y.
2. By (a), π_P L_ℓ = D_P for every point P and every line ℓ ∋ P. Let Z₀ = Σ_ℓ L_ℓ, and identify ⊕_P D_P with F^{7m}
   through bases of the D_P; then Z₀ ⊆ L ∩ ⊕_P D_P.
   * Let V_P be the span of the columns of P in a generator matrix of Z₀^⊥, the complement for the standard dot
     product. The linear dependencies among the columns are exactly Z₀. So for every set S of points,
     r(S) := dim Σ_{P∈S} V_P = m|S| − dim(Z₀ ∩ F^{mS}) ≥ m|S| − dim L_S = m·S_NF(S), with equality on the lines
     (L_ℓ ⊆ Z₀).
   * dim Z₀ ≥ 4m: L_ABW, L_ACX, L_AYZ, L_BCY are independent. In a relation t₁ + t₂ + t₃ + t₄ = 0, restriction to W
     gives t₁ = 0, then X gives t₂ = 0, then Z gives t₃ = 0, so t₄ = 0. Hence r(E) = 7m − dim Z₀ ≤ 3m.
   * A non-line triple T has r(T) ≥ 3m, so r(E) = 3m and r(S) = 3m for every S that contains a non-line triple, which
     is every set of four or more points. Since S_NF = r_NF on sets of at most three points, the rank function of
     (V_P) is m·r_NF.
3. In characteristic 2 this is impossible. The rank values give V_A ⊕ V_B ⊕ V_C = F^{3m}, and V_Z meets each of
   V_A + V_B, V_A + V_C, V_B + V_C only in 0. So V_Z is the graph of isomorphisms V_C → V_A and V_C → V_B; using
   them to identify V_A and V_B with V_C gives V_Z = {(z, z, z)}.
   * V_W ⊆ (V_A + V_B) ∩ (V_C + V_Z), which has dimension 2m + 2m − 3m = m. So V_W = {(z, z, 0)}; likewise
     V_X = {(z, 0, z)} and V_Y = {(0, z, z)}.
   * V_W + V_X + V_Y is the image of (z₁, z₂, z₃) ↦ (z₁ + z₂, z₁ + z₃, z₂ + z₃). Its kernel {(z, z, z) : 2z = 0} has
     dimension m in characteristic 2, so dim(V_W + V_X + V_Y) = 2m ≠ 3m = m·r_NF(WXY). This proves (b).
4. (c), "only if": entropies of a GF(q)-linear state are integers in units of log q and S_NF(P) = 1, so the multiple is
   an integer m; an F-linear party-local Clifford (Witt's theorem over F) reduces each system to S(P) = m qudits, as
   above, and (b) applies. "If": the non-Fano matroid is representable over GF(p) for odd p, hence over GF(q), and the
   CSS state of the code has entropy vector S_NF (§2, step 2). A GF(2^k)-qudit is k qubits, and a GF(2^k)-linear
   state is a qubit state whose entropies in units of log 2 are k times those in units of log 2^k (restriction of
   scalars). ∎

**Corollary 11 [proved; uses the certified face ranks of §4.8, G35] (R41).** S_NF is not a nonnegative combination
of finitely many entropy vectors of qubit stabilizer states, and h_NF is not a nonnegative combination of finitely
many rank functions of subspace families over GF(2). Since S_NF is a GF(p) stabilizer entropy vector for every odd
prime p, the convex cone generated by the entropy vectors of six-party GF(p) stabilizer states is not contained in the
one generated by the qubit ones, for every odd p. For the closed cones Stab⁽ᵖ⁾_6 and Stab⁽²⁾_6 of the paper this
remains open: S_NF could be a limit of such combinations. (R42: it is not; Theorem 12, §4.11.)

*Proof.* Let h_NF = Σ c_i h_i with c_i > 0, where each h_i is the rank function of a nonzero family of subspaces over
GF(2) (for qubit states, h_i = h_{S_i}; §2, step 1). Each h_i is a polymatroid, and by the face remark of §4.8
[certified, G35] the smallest face of Γ₇ containing h_NF is cone(r_NF, r_NF*). So h_i = a_i·r_NF + b_i·r_NF* with
a_i, b_i ≥ 0, not both 0 since h_i ≠ 0. The CSS state of the family (§2, step 2) is a qubit stabilizer state with
entropy vector π(h_i) = (a_i + b_i)·S_NF, because π(r*) = π(r); and a_i + b_i = h_i(A) is a positive integer. This
contradicts Proposition 10(c). For S_NF = Σ c_i S_i, apply this to h_NF = Σ c_i h_{S_i}. ∎

*Erratum (R40).* R40 stated (c) for m = 1 only, and added (d): for q = 2, in every qubit state with entropy vector
m·S_NF and m qubits per system, the two line pieces at one of W, X, Y restrict to different, non-Lagrangian
subspaces. Step 1′ shows that they always coincide; with (b) there is no such state at all, and (d) was true only
vacuously. The R40 text built on (d) (§4.8, "What this shows" below, §5, README ⓪″ with its targeted search at
m = 2, CHANGELOG R40) described a structure that cannot occur. Neither R40 nor its referee noticed step 1′. Corrected
in R41 (the CHANGELOG entry of R40 is left as it was).

*Numerical check of steps 1 and 1′ (G38b, `line_pieces_check` in the s34 self-test).* On states that exist (the CSS
states of the non-Fano code over GF(3) and of the Fano code over GF(2), with one and with two qudits per system, after
random symplectic maps on every system): every L_ℓ has dimension m, and at every point the restrictions of the line
pieces through it are one Lagrangian subspace. This checks the conventions, not the proposition.

*Check [certified: exhaustive searches; the counts for q ≤ 4 also in gate G38b, all of them in
`results/2026-10-14_r40-sandbox/m3_fields.log`].* `cc/s36_target7.c` lists every graph over GF(q) on seven
vertices, in the normal form above, with a given cut-rank vector. Each vertex in turn imposes conditions on its edge
weights: a weight vector must lie in, or must not lie in, the row space of a block of the graph so far. The first
kind is linear and is solved first; the second kind is tested on the solutions.

* S_NF: no graph over GF(2), GF(4), GF(8) or GF(16) (GF(16): 910,604,736 six-vertex partial graphs examined, 75 min
  on one core). Over GF(3), GF(5), GF(7) and GF(9) there are 5,400, 94,200, 621,320 and 2,564,136 normal forms.
* S_F: 532, 27,504 and 1,342,144 normal forms over GF(2), GF(4) and GF(8), and none over the odd fields, as
  Theorem 9 requires.
* Random targets are found over every field, and S_F over GF(16).
* A separate Python implementation (G38b) reproduces the counts over GF(2), GF(3) and GF(4). The counts over GF(2) and
  GF(3) agree with those of `cc/s35_enum7.c` (results/2026-10-14_r40-sandbox/m3_fields.log).

*Independent referee (R40).* A separate agent, with its own code and without reading the repository's code for this
section, re-enumerated both lists (GF(2): all 2^21 graphs; GF(3): another normal form, also 255,602,200 graphs) and
found them equal to the stored ones. It confirmed the orbits, their sizes and representatives, the 52 qubit
certificates, the Steane and non-Fano orbits, the tensor-product and distance statements, and, with two exhaustive
programs of its own, every count of the Check paragraph from GF(2) to GF(9). It checked Proposition 10 line by line;
step 2 follows its shorter argument, and its wording points are taken into the text. Not yet checked by hand.

*Independent referee (R41).* A separate agent, without reading the repository's code, checked Proposition 10 (with
step 1′), Corollary 11 and the erratum line by line, recomputed the face ranks (125 at h_NF, 126 at r_NF and r_NF*),
and confirmed steps 1, 1′ and 2 numerically on the non-Fano CSS state over GF(3) with up to three qutrits per system
and on the Fano state over GF(2) with up to two qubits per system, after random local symplectic maps (dim Z₀ = 4m,
rank function m·r). No error; its wording points (the status of Corollary 11, "positive multiple", the closure
caveats, the reason for the +27) are taken into the text. Not yet checked by hand.

*What this shows.*

* Among the states over ℤ_3 with at most one qudit per system, a witness of Stab⁽³⁾_6 ⊄ Stab⁽²⁾_6 can only be the
  non-Fano state. Every other qutrit-only vector is a qubit vector at multiplicity 2.
* (R41) For the cones generated without closure, the non-Fano state is a witness: no qubit stabilizer state has an
  entropy vector that is a positive multiple of S_NF, and no finite nonnegative combination of qubit entropy vectors
  equals S_NF (Proposition 10, Corollary 11). The annealing runs at m = 2, 3 are consistent with this.
* For the closed cones of the paper the question remained open in R41: S_NF could still be a limit of such
  combinations (whose terms are never multiples of S_NF). A linear inequality valid on Stab⁽²⁾_6 and violated by
  S_NF would settle it. By Theorem 1 it would be a characteristic-2 linear rank inequality on seven variables that is
  not blind on the self-dual slice; none was known to us (Proposition 8 shows that the eight forms considered there
  are blind). (R42: §4.11 gives one, f12, and settles the question: Theorem 12.)
* Where to look: steps 1 and 1′ use the isotropy of the stabilizer group, while the transfer of §4.8 (which gave
  +27 with (LA)₂) works with ranks of subspaces derived from the family. Every family has a stabilizer state, its
  doubled CSS state (§2). A quantitative version of steps 1, 1′ and 2 for that state, with every defect bounded by a
  linear form in the entropies as in the transfer, would give such an inequality, if the bounds close. (R42: they
  close; Lemmas T1–T5 of §4.11.)
* Not covered: p ≥ 5 with one qudit per system (for p = 5 the normal forms number about 2·10¹¹), and systems with
  more than one qudit.

### 4.11 Six parties, the other direction for the closed cones [proved; independent referee] (R42)

**Theorem 12 [proved].** For every d ≥ 3 that is not a power of 2,

    Stab^(Z_d)_6 ⊄ Stab⁽²⁾_6 ,

and hence Stab^(Z_d)_n ⊄ Stab⁽²⁾_n for every n ≥ 6. For n ≥ 7 this was known: Theorem 2 gives
Stab⁽ᵖ⁾_n ⊄ Stab⁽²⁾_n for an odd prime p dividing d, and Stab⁽ᵖ⁾ ⊆ Stab^(Z_d) by Theorem 5. The new case is n = 6.
In particular:

* Stab⁽ᵖ⁾_6 ⊄ Stab⁽²⁾_6 for every odd prime p;
* with Theorem 9, Stab⁽²⁾_n and Stab^(Z_d)_n are incomparable for every odd d ≥ 3 from n = 6 on (before: from
  n = 7);
* for even d that is not a power of 2, Stab⁽²⁾_n ⊊ Stab^(Z_d)_n from n = 6 on (before: from n = 7).

The witness is the non-Fano state, and the certificate is the inequality f12 ≥ 0 below. This answers the question
left open in §4.7, §4.8, §4.10 and §5 (and in the paper's Discussion). The display of f12 as 31 terms is
computer-checked (G39); the proof uses only the structured form.

*Labels.* As in §4.8 and §4.10: A, B, C, W, X, Y, Z = 100, 010, 001, 110, 101, 011, 111, Z is the purifier, and
the non-Fano lines are ABW, ACX, BCY, AYZ, BXZ, CWZ. The *corners* A, B, C, Z lie on three lines each, the
*diagonal points* W, X, Y on two, and every line consists of two corners and one diagonal point. The six *corner
pairs* AB, AC, BC, AZ, BZ, CZ each lie on exactly one line. Two distinct lines meet in exactly one point.

*The form.* For a pure-symmetric vector S on E = {A, B, C, W, X, Y, Z} put

* k(T) = Σ_{P∈T} S(P) − S(T) for T ⊆ E;
* d(P, ℓ) = k(ℓ) − k(ℓ∖P) for a line ℓ through P;
* e_P(ℓ′, ℓ″) = 2S(P) − d(P, ℓ′) − d(P, ℓ″) at a corner P, for two of the lines through P;
* U_Z = 4S(Z) − d(Z, AYZ) − d(Z, BXZ) − d(Z, CWZ);
* U_W = 2S(W) − d(W, ABW) + e_A(ACX, AYZ) + e_B(BCY, BXZ);
* U_X = 2S(X) − d(X, ACX) + e_A(ABW, AYZ) + e_C(BCY, CWZ);
* U_Y = 2S(Y) − d(Y, BCY) + e_B(ABW, BXZ) + e_C(ACX, CWZ);

and

    f12 = 5·U_Z + 4·U_W + U_X + U_Y + k(WXY) + 6·k(ABC) + Σ_{corner pairs} k − Σ_{lines} k − 3·k(ABW) − 3·k(CWZ).

Expanded on the 63 cuts of the visible parties (purity folds the sets containing Z onto their complements;
`transfer_terms`, manifest entry `s34_transfer_form`):

    f12 = + 3·S(A) + 3·S(B) + 3·S(C) − 3·S(W) − 5·S(AB) − 2·S(AC) − S(AW) − S(AX) − 5·S(AY) − 2·S(BC) − S(BW)
          − 5·S(BX) − S(BY) − 5·S(CW) − 4·S(CX) − 4·S(CY) − 6·S(ABC) + 10·S(ABW) + 7·S(ACX) + 7·S(BCY) − S(WXY)
          + 11·S(ABXY) + 11·S(ACWY) + 11·S(BCWX) − 5·S(ABCWX) − 5·S(ABCWY) − 2·S(ABCXY) − S(ABWXY) − S(ACWXY)
          − S(BCWXY) + 2·S(ABCWXY)  ≥ 0 .

(By purity S(ABXY) = S(CWZ), S(ACWY) = S(BXZ), S(BCWX) = S(AYZ) and S(ABCWXY) = S(Z). For instance, the
coefficient 2 of S(Z) is 20 − 15 = 5 from 5·U_Z, + 3 from the pairs AZ, BZ, CZ, − 3 from the lines through Z and
− 3 from −3·k(CWZ); the coefficient 11 of S(CWZ) is 5 from 5·U_Z, 1 each from U_X and U_Y, 1 from the lines and 3
from −3·k(CWZ).)

*Values [by hand].* For S_NF (one qutrit per system; S = π(r_NF)): S = 1 on points; k = 0 on pairs and on the
non-line triples (ABC and WXY among them); k = 1 on the six lines. So every d(P, ℓ) = 1, every e_P = 0,
U_Z = U_W = U_X = U_Y = 1, and

    f12(S_NF) = 5 + 4 + 1 + 1 − 6 − 3 − 3 = −1 .

For the Steane state S_F everything is the same except k(WXY) = 1 (WXY is a Fano line), so f12(S_F) = 0.

*Setting of the proof.* By the reductions of §4.10 (a party-local Clifford factors out n_P − S(P) qubits of each
system) every qubit stabilizer entropy vector is that of a state in which every system P holds exactly S(P) qubits.
Let L ⊆ ⊕_P GF(2)^{2S(P)} be its stabilizer Lagrangian, ω = ⊕ ω_P, π_P the restriction to system P, and for X ⊆ E
let L_X be the elements of L supported in X. As in the proof of Proposition 10 (with S(P) qubits on P),

    dim L_X = Σ_{P∈X} S(P) − S(X) = k(X),

in particular L_P = 0 for every system P. For a line ℓ and P ∈ ℓ the kernel of π_P on L_ℓ is L_{ℓ∖P}, so
dim π_P(L_ℓ) = d(P, ℓ). If two lines ℓ, ℓ′ meet in P, then for s ∈ L_ℓ and s′ ∈ L_ℓ′ isotropy gives
ω_P(π_P s, π_P s′) = ω(s, s′) = 0, so π_P(L_ℓ) ⊥ π_P(L_ℓ′) (step 1 of Proposition 10). For a subspace V of a
symplectic space write rk(V) for the rank of the form restricted to V; then dim(V ∩ V^⊥) = dim V − rk(V).

**Lemma T1 (corners) [proved].** Let V₁, V₂, V₃ be pairwise orthogonal subspaces of a symplectic space of dimension
2s, of dimensions d₁, d₂, d₃. Then rk(V_i) ≤ 2s − d_j − d_k for {i, j, k} = {1, 2, 3}, and
dim(V₁ + V₂ + V₃) ≤ 4s − d₁ − d₂ − d₃.

*Proof.* Let U = V_j + V_k. Since U ⊆ V_i^⊥, V_i ∩ U ⊆ V_i ∩ V_i^⊥, so dim(V_i + U) ≥ d_i + dim U − (d_i − rk V_i)
= dim U + rk V_i. The intersection V_j ∩ V_k is orthogonal to V_i, V_j and V_k, so it lies in (V_i + U)^⊥. Hence
d_j + d_k = dim U + dim(V_j ∩ V_k) ≤ dim U + 2s − dim(V_i + U) ≤ 2s − rk V_i. For V_Σ = V₁ + V₂ + V₃: because the
V_i are pairwise orthogonal, the pull-back of the form on V_Σ along the surjection V₁ ⊕ V₂ ⊕ V₃ → V_Σ is the direct
sum of the forms on the V_i, so rk(V_Σ) = Σ_i rk(V_i) (pull-back along a surjection keeps the rank). Since
V_Σ ∩ V_Σ^⊥ ⊆ V_Σ^⊥, dim V_Σ − rk V_Σ ≤ 2s − dim V_Σ, so
dim V_Σ ≤ s + ½·Σ_i rk(V_i) ≤ s + ½·Σ_i (2s − d_j − d_k) = 4s − Σ_i d_i. ∎

**Lemma T2 (a line transports the rank) [proved].** For a line ℓ = {P, P′, Q}:
rk(π_Q L_ℓ) ≤ rk(π_P L_ℓ) + rk(π_P′ L_ℓ).

*Proof.* For s, s′ ∈ L_ℓ, isotropy gives ω_P(π_P s, π_P s′) + ω_P′(π_P′ s, π_P′ s′) + ω_Q(π_Q s, π_Q s′) = 0. So the
pull-back of the form on π_Q L_ℓ along the surjection π_Q: L_ℓ → π_Q L_ℓ is minus the sum of the pull-backs of the
forms on π_P L_ℓ and π_P′ L_ℓ. Pull-back along a surjection keeps the rank, and the rank of a sum of forms is at most
the sum of the ranks. ∎

**Lemma T3 (the points) [proved].** Let Z₀ = Σ_ℓ L_ℓ ⊆ L (the six lines) and D_P = π_P(Z₀) = Σ_{ℓ∋P} π_P(L_ℓ).
Then dim D_Z ≤ U_Z, dim D_W ≤ U_W, dim D_X ≤ U_X and dim D_Y ≤ U_Y.

*Proof.* At Z the restrictions of L_AYZ, L_BXZ, L_CWZ are pairwise orthogonal, of dimensions d(Z, ·), and
Lemma T1 gives dim D_Z ≤ U_Z. At W, D_W = V + V′ with V = π_W L_ABW and V′ = π_W L_CWZ ⊆ V^⊥ (the two lines meet
only in W). So dim D_W ≤ dim(V + V^⊥) = 2S(W) − dim(V ∩ V^⊥) = 2S(W) − d(W, ABW) + rk(V). By Lemma T2,
rk(V) ≤ rk(π_A L_ABW) + rk(π_B L_ABW), and Lemma T1 at the corner A (lines ABW, ACX, AYZ) and at the corner B
(lines ABW, BCY, BXZ) bounds these by e_A(ACX, AYZ) and e_B(BCY, BXZ). For X use the line ACX (corners A and C,
other lines ABW, AYZ and BCY, CWZ), for Y the line BCY (corners B and C, other lines ABW, BXZ and ACX, CWZ). ∎

**Lemma T4 (supports in Z₀) [proved].** Let z(T) be the dimension of the space of elements of Z₀ supported in T.
Then (i) z(T) ≤ k(T), so z(P) = 0 for every system P; (ii) z(ℓ) = k(ℓ) for every line ℓ; (iii)
z(ABCZ) ≥ k(ABW) + k(CWZ) − dim D_W.

*Proof.* (i) Z₀ ⊆ L. (ii) L_ℓ ⊆ Z₀. (iii) The map L_ABW ⊕ L_CWZ → GF(2)^{2S(W)}, (s, s′) ↦ π_W s + π_W s′, has image
D_W (only these two lines contain W), so its kernel K has dimension k(ABW) + k(CWZ) − dim D_W. On K, the map
(s, s′) ↦ s + s′ takes values in the elements of Z₀ supported in ABCZ, and it is injective: s + s′ = 0 means
s = −s′ ∈ L_ABW ∩ L_CWZ = L_W = 0. ∎

**Lemma T5 (the dual arrangement) [proved].** There are subspaces R_P (P ∈ E) of one vector space over GF(2) with
dim Σ_{P∈T} R_P = Σ_{P∈T} dim D_P − z(T) for every T ⊆ E.

*Proof.* As in step 2 of the proof of Proposition 10. Z₀ ⊆ ⊕_P D_P; identify each D_P with GF(2)^{dim D_P} through a
basis, let G be a generator matrix of the complement of Z₀ for the standard dot product, and let R_P be the span of
the columns of P. The linear dependencies among the columns of G are exactly the vectors of Z₀, so the columns of T
have rank Σ_{P∈T} dim D_P − dim(Z₀ ∩ ⊕_{P∈T} D_P), and Z₀ ∩ ⊕_{P∈T} D_P consists of the elements of Z₀ supported
in T. ∎

Lemmas T1–T5 and their proofs hold verbatim over any field (with that field in place of GF(2)). For the GF(q)-linear
states with entropy vector m·S_NF, which exist only for odd q (Proposition 10), every bound of Lemmas T1 and T3 is an
equality, every e_P vanishes and each D_P is Lagrangian: Lemmas T1–T3 are the quantitative form of steps 1 and 1′
of Proposition 10. Only the last step below needs characteristic 2.

*Proof of Theorem 12.* Write δ_P = dim D_P (unrelated to the δ_i of Lemma 4 (LB)), and r(T) = Σ_{P∈T} δ_P − z(T)
for the rank function of Lemma T5. Apply Lemma 4 (LA) with t = 2 (characteristic 2) and
(A₁, A₂, A₃, B₁, B₂, B₃, C) = (A, B, C, Y, X, W, Z), the labelling of §4.8, to the family (R_P). By Lemma T4,
r(P) = δ_P and r(ℓ) = Σ_{P∈ℓ} δ_P − k(ℓ), so

* H(B₁, B₂, B₃) = r(WXY) = δ_W + δ_X + δ_Y − z(WXY), H(C) = δ_Z;
* α = z(ABC), γ = H(Z | A, B, C) = δ_Z − z(ABCZ) + z(ABC);
* ζ₁ = H(Y | B, C) = δ_Y − k(BCY) + z(BC), ζ₂ = H(X | A, C) = δ_X − k(ACX) + z(AC),
  ζ₃ = H(W | A, B) = δ_W − k(ABW) + z(AB);
* η₁ = H(Y | A, Z) = δ_Y − k(AYZ) + z(AZ), η₂ = H(X | B, Z) = δ_X − k(BXZ) + z(BZ),
  η₃ = H(W | C, Z) = δ_W − k(CWZ) + z(CZ).

The terms δ_A, δ_B, δ_C cancel, and (LA) becomes

    0 ≤ 5·δ_Z + δ_W + δ_X + δ_Y + z(WXY) + 6·z(ABC) + Σ_{corner pairs} z − 3·z(ABCZ) − Σ_{lines} k .

Bound z(WXY), z(ABC) and the z of the corner pairs from above by k (Lemma T4(i)), and −3·z(ABCZ) from above by
−3·k(ABW) − 3·k(CWZ) + 3·δ_W (Lemma T4(iii)). Then δ_Z, δ_W, δ_X, δ_Y have the coefficients 5, 4, 1, 1, all
positive, and Lemma T3 bounds them from above. The result is 0 ≤ f12(S). So f12 ≥ 0 on the entropy vector of every
qubit stabilizer state of seven systems, and since f12 is linear, on the closed cone Stab⁽²⁾_6 (Z the purifier).

The non-Fano matroid is representable over GF(p) for every odd prime p, and the CSS state of its code C_NF has
entropy vector π(r_NF) = S_NF (Theorem 1, step 2), so S_NF ∈ Stab⁽ᵖ⁾_6 while f12(S_NF) = −1. If d has an odd prime
factor p, the Z_d CSS state of the image of C_NF under Z_p ≅ (d/p)·Z_d ⊆ Z_d has entropy vector (log₂ p)·S_NF in
bits (Theorem 1(d), step 3), so S_NF ∈ Stab^(Z_d)_6. For n > 6, merge the extra (empty) parties into the
purifier. ∎

Status: Lemmas T1–T5 and the assembly are elementary linear algebra; the characteristic enters only through
Lemma 4 (LA) [proved; independent referee, R33]. Theorem 12 also uses Theorem 1 (step 2), for composite d step 3 of
Theorem 1(d) (the CSS formula, proved here), and the reductions of §4.10. It does not use the cited entropy formula
of Gross–Walter. The incomparability for composite odd d goes through Theorem 9 and inherits its status
(Theorem 1(d), steps 1–2).

**Remarks.**

* *A self-dual characteristic-2 inequality.* With π the connectivity function of §1, g12 = f12∘π (63 terms on the
  nonempty subsets of seven elements, h(E) included) is a linear rank inequality valid over every field of
  characteristic 2: a family over GF(2) with rank function h has a CSS state with entropy vector π(h) (Theorem 1,
  step 2); GF(2^k) reduces to GF(2) by restriction of scalars, and any other field of characteristic 2 to a finite
  one by the reduction at the end of the proof of Theorem 1 (the residue fields are then finite of characteristic 2).
  Since π(h*) = π(h), g12∘D = g12. g12(r_NF) = g12(r_NF*) = −1, so g12(h_NF) = −2, and g12(r_F) = 0. So g12 is
  the inequality that Proposition 8 asks for in this direction, a
  characteristic-2 seven-variable inequality that is not blind on the self-dual slice. It is not a nonnegative
  combination of Shannon inequalities and the eight forms of Proposition 8 in any relabelling, since each of those is
  ≥ 0 at the self-dual polymatroid h_NF. We know of no earlier inequality of this kind; no literature check was made.
* *Why this transfer works where that of §4.8 gave +27.* The derived family of §4.8 sees only ranks, and for a
  general subspace family the two line pieces at W, X, Y need not coincide. A stabilizer state carries isotropy:
  Lemma T2 transports the (small) rank of the form from the two corners of a line to its diagonal point, which bounds
  dim D_W, dim D_X, dim D_Y, and Lemma T1 bounds dim D_Z. The terms that cannot be bounded this way (δ_A, δ_B, δ_C)
  cancel in (LA).
* *Sharpness [certified, G39].* f12 = 0 at S_F; on the GF(2) list of §4.10 (10,773 vectors, at most one qubit per
  system) the minimum of f12 is 0, attained 38 times; on the GF(3) list (27,334 vectors) f12 < 0 only at S_NF
  itself. The other 209 vectors of the non-Fano orbit are separated by the relabelled forms f12∘σ.
* *Probably not optimal [empirical].* f12 − Σ_{corner pairs} k and f12 − k(ABC), which are stronger, were not violated
  on the GF(2) list or on random states (the referee's tests and ours), while f12 − Σ_{corner pairs} k − k(ABC) is
  violated (−1 on the GF(2) list). The proof gives neither.
* *Numerical check of the proof (G39b; checks, not proofs).* `transfer_chain` evaluates every step on a given
  state: the dimensions d(P, ℓ) of the restrictions, the bounds of Lemmas T1 and T3, Lemma T4(iii), (LA) on the
  rank function of Lemma T5, and f12 ≥ 0. It holds on random compressed qubit graph states with one to three qubits
  per system, on compressed CSS states of random GF(2) families pushed towards the non-Fano dependencies, and on one
  and two copies of the Steane state after random local symplectic maps (where f12 = 0); Lemma T5 is also checked by
  building the arrangement explicitly, and the compression (symplectic reduction) against uncompressed states. The
  gate runs 40 states; a sandbox run of 400 found no failure (`results/2026-10-14_r42-sandbox/`).
* *Adversarial search [empirical].* `cc/s37_minform7.c` (`--transfer-anneal`) minimises f12 by simulated annealing
  over qubit graph states with a given number of qubits per system, optionally pulled towards m·S_NF; the wrapper
  recomputes the entropies of the best graph exactly and runs `transfer_chain` on it. Nine sandbox runs of 60–240 s
  (1, 2 or 3 qubits per system, also unequal; two of them pulled towards 2·S_NF and 3·S_NF), about 6·10⁸ moves:
  no run went below 0, six reached 0, and every reported minimiser passed every step of the proof
  (`results/2026-10-14_r42-sandbox/transfer_anneal.log`).
* *Not covered.* The proof needs a field of characteristic 2. For Z_{2^k} with k ≥ 2 (the Z_4 cube state violates
  (LA)₂), whether Stab^(Z_d)_6 ⊄ Stab^(Z_{2^k})_6 for odd d is open, and so is Stab⁽²⁾_6 ⊊ Stab^(Z_{2^k})_6
  (strict from n = 7, by the Z_4 cube state). Two odd primes at n = 6 remain open as well.
* *Referee.* A separate agent, which read this note and the data but not the code of this section, checked
  Theorem 12 line by line: the 31-term expansion recomputed from the structured form (no difference); the values at
  S_NF and S_F (both vectors also from state vectors; f12(S_NF) = −1 from the non-Fano code over GF(3), GF(5), GF(7)
  and GF(11)); the setting; Lemmas T1–T5; every term of (LA) and the direction of every bound; the consequences; and
  the remarks on g12 (63 terms, g12∘D = g12, the four values). Verdict: proved as stated. With its own code it
  checked every step on 26,549 qubit stabilizer states (about 19,900 of them needing the reduction to S(P) qubits;
  among them Steane copies, CSS states of families pushed towards the non-Fano dependencies with m = 1, 2, 3, tensor
  products and the 38 vectors of the GF(2) list with f12 = 0), enumerated all 2^21 seven-vertex qubit graph states
  (minimum 0, on 686 graphs and 38 vectors), and evaluated 5.85 M random graph states and about 78 M annealing
  proposals: f12 never went below 0. Weakened forms (without the e-terms, with 3·S(Z) in U_Z, or with U_W instead of
  4·U_W) do go negative on the same families, so the tests can detect a violation. Its corrections are in the text:
  for n ≥ 7 the statement was already known (Theorems 2 and 5), so the new case is n = 6; Theorem 12 does not depend
  on the Gross–Walter formula; g12 over infinite fields of characteristic 2; notation (V_Σ, Q); and the statements
  elsewhere in this note that still called the question open. Not yet checked by hand.

## 5. Open questions

* **The smallest number of parties with a dimension dependence** is 5 or 6: for n ≤ 4 there is none, for any d
  (Theorem 7), and at n = 6 there is one (Theorem 9); at n = 6 qubits and every odd local dimension are incomparable
  (Theorems 9 and 12, R42).
  * n = 5: by Theorem 1 this is equivalent to a characteristic dependence of the six-variable linear rank cone on
    its self-dual slice. No characteristic-dependent six-variable inequality is known; the smallest known ones
    have seven variables. All 641 known extreme rays of Stab₅ lie in Stab^(Z_d)_5 for every d (§4.9, R38–R40),
    so none of them can witness a dependence at n = 5. A witness would need an extreme ray of Stab₅ that is not yet
    known (the 1,155 undecided orbits of the catalogue are candidates, not the only ones), or, for d divisible by a
    square, a vector of Stab^(Z_d)_5 outside Stab₅ (the ℤ_d question below).
  * ~~n = 6, the other direction: is Stab⁽ᵖ⁾_6 ⊄ Stab⁽²⁾_6 for some odd p?~~ Settled by Theorem 12 (§4.11, R42):
    yes, for every odd p, and Stab^(Z_d)_6 ⊄ Stab⁽²⁾_6 for every d that is not a power of 2. The witness is the
    non-Fano state; the certificate f12 is a quantitative ("symplectic") version of the proof of Proposition 10.
    History: among the states over ℤ_3 with at most one qudit per system only the non-Fano state could be a witness
    (§4.10, R40); for the cones generated without closure it was one (Proposition 10, Corollary 11, R41); the direct
    probes had been negative (they are explained by Proposition 8): symmetrised Fano and non-Fano polymatroids
    r + r* never violate (65) or (91) in any of the 5,040 relabellings (minimum 19 and 16), and a 240-second
    annealing probe over seven-party qutrit graph states found nothing either (best value 4 ≥ 0;
    `s27_cube_states.py --probe6`).
  * n = 6 for odd primes p < q: open; Theorem 2′ needs 2p + 3 parties.
  * n = 6 for Z_{2^k}, k ≥ 2: open in both directions against qubits and against odd d (§4.11, "Not covered").
* **The transfer for other primes.** Can the derived-family argument of §4.8 lower the bound 2p + 3 of Theorem 2′
  for odd p, for example by recovering the configuration J − I of Lemma 4 from a self-dual sum r + r*?
* ~~**Different odd primes, and all pairs of primes.**~~ Settled by Theorem 2′ (R33): for primes p < q the cones
  are incomparable from n = 2p + 3 on. Whether fewer parties suffice for odd p is open. For example, at n = 8
  one direction for p = 3, Stab⁽q⁾_8 ⊄ Stab⁽³⁾_8, follows from the eight-variable characteristic-3 inequality of
  DFZ arXiv:1401.2507 (Theorems 4.1 and 4.2, not checked independently here) and doubled codes (§4.3); the other
  direction would need a correct inequality for characteristic ≠ 3, and the printed Theorem 3.1 is not usable
  (§4.5).
* **The 1,155 undecided orbits of the A2 catalogue**: see §3.1.
* **Are Z_d stabilizer states ever beyond all prime ones?** Is Stab^(Z_d)_n always contained in the closed convex
  hull of the prime cones ∪_p Stab⁽ᵖ⁾_n? (The Z_4 cube state is inside.)
  * (R37) For squarefree d the answer is yes: a group of exponent dividing d is the direct sum of its Sylow
    subgroups, which are vector spaces over prime fields, and every subgroup is the direct sum of its
    intersections with them, so every subgroup rank function is a sum of linear rank functions over prime fields
    (the paper, Sec. VII). The question is open only when d is divisible by the square of a prime.
  * By Theorem 1(d), this asks whether the self-dual slices of the mixed-linear cone Γ^MixL and of the abelian
    cone Γ^Abl agree.
  * The doubled-code sandwich of §4.3 holds verbatim for subgroups. So the n-party question lies between the
    questions "Γ^MixL_m = Γ^Abl_m?" for m = n and m = n + 1 variables.
  * Khazaei (arXiv:2608.09543) proves Γ^MixL ⊆ Γ^Abl and writes "We do not know whether the first inclusion can be
    strict."
* **Smaller party counts for Theorem 5**, e.g. Z_4 against Z_6 below n = 15. (R42: for d not a power of 2 against
  d′ = 2, from n = 6, Theorem 12.)
* ~~**Five parties over every field**: lift or exclude the known extreme rays.~~ All 641 are certified for every d
  (§4.9, R38–R40). The n = 5 question now concerns extreme rays not yet known (first item).
* ~~**Six parties**: is there a characteristic-dependent seven-variable inequality g with g∘D = g?~~ Yes: g = f9 + f9∘D
  (Theorem 9, R36), and for characteristic 2, g12 = f12∘π (Theorem 12, R42). BKL's Theorems 6.4 and 6.6 were checked
  and are blind (Proposition 8).

## 6. Literature status

**R37 (2026-10-11).** For the paper, every cited statement was checked against its source by hand, and a
prior-work search was made. The paper gives the final attributions. The main prior work that the earlier record
below did not know:

* F. Matúš, *Polymatroids and polyquantoids* (WUPES 2012, arXiv:1210.7931): the lift h_S, and the bijection
  between pure-symmetric submodular vectors and (tight) self-dual polymatroids, which is fact (c) of Theorem 1.
  Jowett, Mo and Whittle (Adv. Appl. Math. 81, 1 (2016), arXiv:1605.01455) found the same correspondence for
  connectivity functions.
* P. Sarvepalli and R. Raussendorf (Phys. Rev. A 81, 052333 (2010)): an identically self-dual matroid that is
  representable over a finite field gives a pure-state quantum secret sharing scheme, built on a CSS state whose
  lift is twice the rank function. This is the matroid case of the converse direction of Theorem 1.
* N. Osthus (M.S. thesis, University of Minnesota, 2024, https://hdl.handle.net/11299/269174): for prime d, the
  entropies of CSS states, and twice those of general stabilizer states, as ranks of submatrices, so that the
  doubling 2S = λ(r_L) follows from his formulas; the question whether the local dimension matters; and the
  misprint in the LMRW state (R6) (§3.2).
* V. Peña-Macias and H. Sarria (J. Inf. Optim. Sci. 42, 345 (2021), arXiv:1903.11587v2): characteristic-dependent
  inequalities for the configuration J − I for every t over fields (their Theorems 7 and 9; §4.4), and the
  failure of BKL's (6.17) (§4.5).
* M. Walter and F. Witteveen (J. Math. Phys. 62, 092203 (2021)) together with Bao et al. (arXiv:2002.05317,
  §3.3.1): for prime d, Theorem 7 also follows from the literature.
* Linden et al. (LIPIcs TQC 2013, Sec. 5) conjectured that the inequalities obtained from common information,
  added to the classical ones, describe the stabilizer cone; the paper (Sec. VII) discusses what its results say
  about this conjecture.
* Feng and Luo (J. Phys. A 57, 395304 (2024)): the two usual definitions of qudit stabilizer states agree exactly
  when d is squarefree. This repository uses the definition through stabilizer groups, which allows non-free
  Lagrangian subgroups.

**R42 (2026-10-14).** No literature check was made for §4.11 (Theorem 12, the inequalities f12 and g12, and the
symplectic transfer). Before any of them is presented as new, a human check in the communities listed below is
required, as for the earlier results.

Earlier record (R32–R36):

So far there is only an AI pre-check. It found no statement of the converse direction of Theorem 1, or of
Theorem 2, in LMRW (arXiv:1302.5453), Gross–Walter (arXiv:1302.6990), BCHS (arXiv:2006.16292) or Majenz's thesis
(arXiv:1810.12845, Ch. 5); the search used web summaries, not full texts. Two pieces are standard:

* the entanglement of CSS and graph states via matroid connectivity and cut-rank functions;
* the use of characteristic-dependent linear rank inequalities in network coding.

A human check in at least three communities (quantum entropy cones, network coding and linear rank inequalities,
matroid theory) is required before any of Theorems 1, 1(d), 2, 2′, 5, 7 and 9 is presented as new. Cited inputs and
how they were checked:

* DFZ arXiv:1311.4601v1, Thms 8.6–8.9: checked term by term against the rendered PDF pages (R32).
* Peña–Sarria arXiv:1905.00003v3, Thm 3 / Example 6: transcribed from the rendered page 4 and tested as in §4.4
  (R33). It is used only as a cross-check, since Lemma 4 is proved here. Its publication status was not checked.
* BKL arXiv:1108.2489v1, eqs. (6.8) and (6.17): checked term by term against the rendered PDF pages 19 and 22 (R36);
  (6.17) is misprinted (§4.5). Used only in Proposition 8.

For §4.6 (R34), from an AI pre-check of web texts (not yet checked by hand):

* Gross–Walter (arXiv:1302.6990), Thm 2: the entropy formula used in step 1 of Theorem 1(d), for every d.
* LMRW (arXiv:1302.5453): arbitrary local dimension through abstract error groups (Example 7 is the Weyl–Heisenberg
  group of Z_d), entropy through abelian quotients (Prop. 9). Their cone is the closure over all stabilizer states,
  without fixing d; no dependence on d is discussed.
* Majenz's thesis (arXiv:1810.12845): Cor. 5.12 holds "under a technical assumption on the local Hilbert space
  dimension". The assumption could not be read from the web text; it must be checked by hand.
* BCHS (arXiv:2006.16292) treats qubits only.
* Khazaei, *Separating Abelian and Homomorphic Entropy Cones* (arXiv:2608.09543): the open inclusion
  Γ^MixL ⊆ Γ^Abl quoted in §5.
* No statement was found that stabilizer entropy cones depend on d, or that compares Z_4 with qubits or GF(4).

For §4.7 (R35), again from an AI pre-check:

* LMRW, Thm 13 characterises the four-party cone over all local dimensions jointly. Its realisations are not
  qubit-only: two ray families are realised by qutrit states. No statement was found that the cone is the same for
  each single d (Theorem 7).
* Kaced (arXiv:1611.04109): the duality operator D, closure of the polymatroid and Ingleton regions under D, and the
  self-duality of Shannon and Ingleton (checked against the PDF text). Nothing was found on characteristic-dependent
  inequalities restricted to self-dual polymatroids.

For §4.8 (R36), from an AI pre-check (web search, abstracts only):

* No statement was found that separates stabilizer cones of different local dimensions at six parties, or that
  derives a characteristic-dependent inequality by applying a known one to a family built from intersections
  U_i ∩ (U_j + U_k). Closest topics seen: Peña-Macias, *Access structures for finding characteristic-dependent linear
  rank inequalities* (Kybernetika 2023; inequalities from n × n binary matrices, no seven-variable or self-dual
  results); Bérczi et al. (arXiv:2507.10709; skew-representability, tensor products, rank inequalities).
* Peña-Macias–Sarria-Zapata, arXiv:1903.11587 (characteristic-dependent inequalities via complementary vector
  spaces) was not transcribed or checked here.

## 7. Reproduce

    python3 scripts/s26_reduction.py --self-test                     # CSS formula and doubling vs state vectors
    python3 scripts/s26_reduction.py --n4 --lrs                      # n = 4: LP + exact rays (lrs) + certificates
    python3 scripts/s27_cube_states.py --state-vector --random 400   # Theorem 2 witnesses and checks, ~10 s
    python3 scripts/s27_cube_states.py --probe6 240                  # exploratory six-party probe
    python3 scripts/s28_generalized_cubes.py --t 2 --fields 2 3 5 --exhaustive01             # Lemma 4, t = 2, ~25 s
    python3 scripts/s28_generalized_cubes.py --t 3 --fields 3 5 7 11 --random 300 --exhaustive01   # Theorem 2', p = 3
    python3 scripts/s28_generalized_cubes.py --t 5 --fields 3 5 7 11                         # Theorem 2', p = 5
    python3 scripts/s28_generalized_cubes.py --dfz1401                                       # the misprint of §4.5
    python3 scripts/s29_zd_stabilizer.py --cases 2:4 4:8 3:9 6:4 --state-vector --lemma-tests 2200000   # §4.6
    python3 scripts/s30_dimension_profile.py --n4 --one-cut 2 3 4 5 6 7 --six                 # §4.7
    python3 scripts/s30_dimension_profile.py --bkl                                           # BKL checks, §4.5
    python3 scripts/s31_six_party.py --derive --witness --invariant --face --lemmas 200 --stress 400 --exhaustive  # §4.8
    python3 scripts/s32_unimodular.py --check --full                 # §4.9: exact checks of the stored certificates
    python3 scripts/s32_unimodular.py --self-test                    # §4.9: neutral moves, ring expansion, Z_d formula
    python3 scripts/s32_unimodular.py --check-cover                  # §4.9: covering certificates
    cp data/s32_stab5_unimodular.npz u.npz       # §4.9: the covering search, on a copy and into a new cover file
    python3 scripts/s32_unimodular.py --cover --rows 587,589 --seconds 300 --out u.npz --cover-out c_new.npz
    python3 scripts/s33_ci_oddchar.py --rows 456,488                 # §4.9: odd-characteristic CI screen (exploratory)
    python3 scripts/s34_six_party_survey.py --check --full           # §4.10: survey data, GF(2) list recomputed
    python3 scripts/s34_six_party_survey.py --self-test              # §4.10: conventions, cut ranks, line pieces (R41)
    python3 scripts/s34_six_party_survey.py --enumerate --out lists.npz   # §4.10: cc/s35_enum7.c, p = 2, 3 (~15 min)
    python3 scripts/s34_six_party_survey.py --fields 2,3,4,5,7,8,9   # §4.10: cc/s36_target7.c, S_F and S_NF (~4 min)
    python3 scripts/s34_six_party_survey.py --fields-python 2,3,4    # §4.10: the same counts in Python (~2 min)
    python3 scripts/s34_six_party_survey.py --transfer --full --trials 400   # §4.11: f12 and every step of its proof
    python3 scripts/s34_six_party_survey.py --transfer-anneal "2,2,2,2,2,2,2;3,3,3,3,3,3,3" --seconds 120   # §4.11
    python3 selftest.py                                              # gates G29-G39 (G29b, G31b-G39b with --full)

Raw output of the sandbox runs: `results/2026-10-02_r32-sandbox/`, `results/2026-10-02_r33-sandbox/`,
`results/2026-10-02_r34-sandbox/`, `results/2026-10-02_r35-sandbox/`, `results/2026-10-02_r36-sandbox/`,
`results/2026-10-13_r38-sandbox/`, `results/2026-10-14_r39-sandbox/`, `results/2026-10-14_r40-sandbox/`,
`results/2026-10-14_r41-sandbox/`, `results/2026-10-14_r42-sandbox/`.
