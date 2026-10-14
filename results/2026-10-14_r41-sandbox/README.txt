R41 sandbox runs (2026-10-14 afternoon; 2-core container; no server time used).  Roadmap item M3: six parties, the
other direction (docs/R32_reduction_and_local_dimension.md section 4.10, Proposition 10 strengthened and
Corollary 11; both are proofs, so this folder holds only checks).

line_pieces_check.json  scripts/s34_six_party_survey.line_pieces_check(): the CSS states of the non-Fano code over
                        GF(3) and of the Fano code over GF(2), with m = 1 and m = 2 qudits per system, after random
                        symplectic maps on every system (30 random transvections each).  Per state: the stabilizer is
                        isotropic; every line piece L_l has dimension m; at every point the restrictions of the line
                        pieces through it are one Lagrangian subspace (Proposition 10, steps 1 and 1'); the elements
                        supported on WXY have dimension m for Fano and 0 for non-Fano; moved_points counts the points
                        of the line ABW whose restriction is no longer the coordinate z-plane (the maps act).
s34_selftest.log        python3 scripts/s34_six_party_survey.py --self-test (now with the line-piece check)
selftest_full.log       python3 selftest.py --full (49 gate lines)
