# AUDIT 7: KQ winding fallback, second-order enclosures, D2E2 and the D2 rules, continuation runs, checkers, documentation; final-pass scripts

Independent audit, 2026-10-23, 07:56–09:55. The audit covers the round-7 request and its three addenda, and verifies
the fixes made during the round.

**Scripts and outputs:** `audit/round7/`. Hashes of the audited files are in `round7/hashes_round7.txt` (08:22) and
`round7/hashes_round7_addenda.txt` (09:46). `nb_cert2.py` was patched again at 09:41; its re-tests are listed below.

**Rules kept.** At most 2 cores were used, and nothing was written into `code/` or the project documents:
- numba caches go to `audit/numba_cache`, and Python runs with `-B`;
- the project's `check_all.py` was run only through `g12_check_all_guarded.py`, which redirects `RESULTS.md` into `round7/test_g12*/`.

Two slips:
- at 08:35 I copied `bulk.json.records.jsonl` (27 MB) to the system temp directory and deleted it a minute later;
- for about a minute three of my processes overlapped.

## Verdict

**No soundness issue was found.** This covers:
- the round-7 changes (KQ winding fallback, bb_replay/check_all/runio, the c3_run resume);
- the addenda: the second-order enclosures for the f-chart, a-chart and C3/KQ polynomial charts (`nb_cert2`, `nbcoef_f2`, `nbcoef_a2`, `nbeng_c3_{KBQe,Ca}_2`), the chart D2E2, the D2Tp coverage rules (i) and (ii), and the continuation-level run structure.

Every finished or partial run I checked passes my own replay, re-certification with the final code, and 100-digit
spot checks.

**All round-7 issues are closed:**
- the two checker problems found early in the round (R7-1, R7-2) were fixed during the round, and I verified the fixes;
- the unaudited enclosure change (R7-3) has now been audited: code review plus independent 50- and 100-digit tests, 0 violations;
- the documentation items were corrected.

The proof is complete only when the open continuation levels (bulk, KQ*, Ca) have finished. The final pass is ready
for that (`round7/final_pass.sh`).

| id | severity | finding | status |
|---|---|---|---|
| R7-1 | major (checking tool) | **`bb_replay --recert` gave false problems on every pre-parity run.** The final engines stop at total ≤ 11; `bb_replay` required ≤ the recorded bound (≤ 10). Measured (60 leaves each): L2P 19, B2 22 and D2K 12 false "problems". My own re-certification of **every** leaf of the 19 finished bb_run families with the final code: 182,592 leaves, 0 failures; 18,627 now certify at total 11 rather than their recorded bound. | **Fixed** 09:40 (≤ max(bound, 11) for parity engines). Verified: L2P, B2, D2K, F0 and D2E2, 60 leaves each, 0 problems (`g14b`). |
| R7-2 | major (checking tool) | `check_all` skipped the initial-box grid test for the C3 families (it required `'exclude' in H`). | **Fixed** 09:40 (test runs whenever `domain` and `n0` exist). Verified by reading and by a guarded run. My `g2` finds every finished grid complete. |
| R7-3 | major (audit coverage) | `nb_cert2` (second-order enclosure) entered KQ at 08:40, unaudited, and was later extended to the a-chart, C3 Ca and bulk (`--order2`). | **Audited (§6):** code review and independent tests, 0 violations. FLOAT_RIGOUR §5 now covers it. |
| R7-4 | minor (checking tool) | `check_all`'s completeness did not enforce per-family bound limits (KA ≤ 2; C3 targets; zone ≤ 5). | **Fixed** (`code_limits`). Verified by reading. `g2` enforces the same limits independently. |
| R7-5 | minor (float write-up) | FLOAT_RIGOUR's winding paragraph kept the R4-3 errors, and its libm list omitted cosh/sinh (and atan, atan2, tan). | **Fixed** 09:44, with my measurements (`g8`: cosh/sinh ≤ 1.06 ulp, atan2 ≤ 0.68 ulp). |
| R7-6 | minor (docs) | The e ≤ 0.2 exhaustiveness step said "l² ≥ 1.96e ⇒ ε'' ≤ 0.51". | **Fixed** (l² > 2e ⇒ ε'' < 0.5). |
| R7-7 | minor (docs) | "Exactly quasi-homogeneous" was wrong (the correct statement is weighted order 10, degrees 10–28), and λ was undefined in the paper. | **Fixed.** |
| R7-8 | minor (docs) | Stale or inexact statements (§5 table). | **Fixed** (spot-checked: driver order, D2b extremes, μ sliver, RESULTS/failure wording, C3_CORE winding-first, CORNER_E1, c3_run docstring, "seven rounds"). |
| R7-9 | minor (hygiene) | `runio.resolve` silently preferred `.json` over `.json.gz`. | **Fixed** (raises if both exist). |
| R7-10 | info | The C3 chart engines used by the runs (01:21) predate min-over-rounds. | README now notes the engine versions. |
| R7-11 | info | The record files mix engine generations. The KQ files hold count-only first-order, winding-first first-order, transitional 08:40/08:58, and final (09:09) records; bulk holds 3,825 first-order records plus `--order2`; Ca holds first- and second-order records. Headers describe only the last process. | Acceptable, because every generation is sound (first-order audited in rounds 1–6, second-order in §6) and the final pass re-certifies every winding leaf and samples the rest, with the 07:50 code first and the final code second (`g1`, `g9`). Recommend per-record provenance in future runs. |
| R7-12 | info | My checker's pruning predicate in float flagged one Ca pruned leaf. | False alarm (float rounding of l² − 1). With exact arithmetic, all 21,177 pruned Ca leaves pass (`g22`); `g2` now uses exact predicates. |

## (1) KQ winding fallback (kq_run.py 07:50): OK

**(a) Variables.**
- KQ boxes (e, l, ω, t), with t = κ²/(1+κ²), are exactly nb_cert chart 1.
- `wind()` uses chart-1 data, and `wind_box` is affine in t.

**(b) Positive multiple.**
- The f-chart polynomial is t·P(f), with P = μ²k²sin²fσ² − D²B² and μ² = (1−t)/t.
- P(f) = c(f)·Pol(u(f)) with c = D⁶(du/df)²/(μ²(1+e)²r²) > 0, up to normalisation constants.
- σ² cancels against ρ² = r²σ², so axis crossings need no exception.
- The identity holds for every physical point (e < 1, 0 < t < 1). The faces e = 1 and t = 1 contain no physical point.

**(c) Mixed records.**
- `--resume` keeps completed records only.
- Old records lack `wL` and `nwind`, and every reader uses `.get` or digit codes.
- The header is written by the last process.
- Evidence: KQ1 snapshot (427 records; 68 winding leaves). 0 replay problems. All 68 winding leaves re-certified, plus 60 count leaves. 100-digit checks on 180 points: W reproduced everywhere, all N_mult even (`g1_kq_wind_snap.out`).

**Stress tests at 100 digits** (`g5`, 07:50 code). 1,500 random boxes, biased to the faces, of which 557 were winding-certified.
- Each was checked at 8 points.
- 0 violations; minimum margin 0.123.

## (2) bb_replay.py: OK

- PARITY_ENGINES whitelist; KQ via `kq_run.init()` and `certify`.
- Engines a, f and KQ re-certify with `nb_cert2` (`certify_box_deep2`, `wind_box2`).
- **Functional test:** 234 leaves of new KQ1 records, `--recert -1`: 0 problems (`g10`).
- **Cross-engine check** (`g23`, bulk). The first-order leaves made before the 09:26 resume re-certify with the second-order engine that `bb_replay` uses: 300 regular and 300 winding leaves, 0 failures. The converse fails, as expected (265/300 and 300/300): the new bulk leaves depend on `nb_cert2`.

## (3) check_all.py (09:40): OK

- **R6-1, R7-2 and R7-4 fixes** are in place.
- **Continuations** are replayed with each file's own header scale, and their roots must be owed by earlier files.
- **D2 rules (i) and (ii):** see §7.
- **Guarded run on the current data** (`test_g12b/RESULTS.md`):
  - every finished family is complete, including C1p, D2Tp and D2E2;
  - the base files of KQ2, KQB, KQK1, KQK2 and Ca are INCOMPLETE until their continuation levels finish;
  - bulk and KQ1 have no final file yet.

## (4) c3/c3_run.py: OK

- **Certification** is unchanged except for the engine module: the second-order `nbeng_c3_Ca_2` when it exists. The header records `engine_module`.
- **Resume** keeps completed records only.
- **The `--boxes`/`--scale` continuation mode** passes the scale to the workers and records it in the header. `certify` does not depend on the scale.
- **Ca:**
  - Base file: 34 done, 478 open; 48,981 leaves, **39,980 at near count 9**, which rely on the parity composite 9 + 0 + 2 ≤ 11. Phi: 570 leaves, all 0. KA: 64 leaves, all exactly 2.
  - Ca.r1 partial (42 records, 1,180 leaves): 0 replay problems with its own scale, 30 spot points 0 violations, 60 leaves re-certified 0 problems (`g4_car1_snap.out`).

## (5) Documentation against code

**Verified correct** (own computations):
- **Sharpness example** (`g6`): Pol has 10 simple roots; the 2-D count is 10 (3 min, 5 saddles, 2 max; Euler characteristic 0); the complex pair has |log|z|| = 1.81.
- **D1b leading term** (`g7`): r₁ = (3 sin 3u − sin u)/4, rank 4.
- **C3 forms** (`g11`, from my own P(f)): P̂₀, B̂, the E = m = 0 form, the level-2 leading part; only even powers ν¹⁰…ν²⁸.
- **FLOAT_RIGOUR §1–§4 numbers.**
- **FLOAT_RIGOUR §5** (new): term and degree counts recomputed. `nbcoef_c3_Ca_2` has 49 terms, total degree 26 and power chains up to 18, so N ≤ 69 ≤ 78.
- **Atlas table and seams against `launch_runs.sh`** (`g2` header comparison).
- **The KQ "aphelion pair at τ = O(ε²)":** holds for μ bounded below.

**Mismatches found at 08:22, all corrected at 09:44** (R7-5–R7-8), among them:
- FLOAT_RIGOUR's winding weights and libm list;
- "every issue raised was fixed";
- the 1.96e step, "exactly quasi-homogeneous" and the missing λ;
- the driver order, the D2b extremes and the μ sliver;
- "failures (must be 0)" and "continuation runs";
- C3_CORE §5, CORNER_E1 "(open)", the c3_run docstring and zone comment, "six rounds".

**Not checked:** the literature statements (paper §1, README "Prior work"), and STATEMENT.md, LOG.md and LIT_REVIEW.md (unchanged).

## (6) Second-order enclosures (addendum A and addendum 2; R7-3)

**Code review.**
- **f-chart (`box_coeffs2`).** The Hessian in (e, l, ω) via the chain rule through p = l sin ω, q = l cos ω: all six entries were derived by hand.
- **a-chart (`box_coeffs2a`).** With e = sin η and b = cos η:
  - H_ηη = C_ee b² − 2C_eb eb + C_bb e² − C_e e − C_b b;
  - H_ηl and H_ηω as stated;
  - the path bounds |e'|, |b'| ≤ r_η, |e''| ≤ r_η², |e'''| ≤ r_η³, and |sin x − sin y| ≤ |x − y| for the abs-poly arguments.
- **Remainder (Faà di Bruno).** g''' = ΣC_ijk g1g1g1 + 3ΣC_ij g2g1 + ΣC_i g3:
  - ordered-triple multiplicities 1/3/6;
  - sorted-pair and sorted-triple slot orders match the generators (`gen_nb2`, `gen_nb2a`, `gen_c3_2`);
  - for C3 there is no chain rule, and the remainder is (1/6)Σ mult·sup|∂³c|·r r r.
- **Assembly.** Σ_v r_v|Σ_k ∂_v c_k T_k| + Σ_m q_m|Σ_k H_{k,m} T_k| + Σ_k R_k|T_k| is valid because d is common to all k.
- **Shared parts.** `certify_box2`/`certify_box_deep2` are line-for-line `nb_cert`'s loops with `arc_signs2`; `wind_box2` is `wind_box` with `labels2`; `away_count_t2` is `c3_away.away_count_t` with `arc_signs2`.
- **Wiring.** `kq_run` (`away_count_t2`, `wind_box2`, `nbeng_c3_KBQe_2.count`; wind 8/16384), `c3_run` (`engine_module`) and `bb_run --order2` (a, f) were all read.

**09:41 patch.** The radius is now min(t_old, t_new), two complete radii:
- t_old repeats `nb_cert.arc_signs` / `wind_cert.comp_bounds` exactly, with the same aab and pads;
- t_new carries its own aab_new: the magnitudes entering s, rd, g, H and R.

Both are complete upper bounds, so the minimum is one too, and the final radius is ≤ the first-order radius arc by arc.
Checked in `g24`:
- 120,142 signs certified by `nb_cert.arc_signs`: 0 lost or changed by `arc_signs2`, 6,703 additional;
- 21,506 labels of `wind_cert.labels`: 0 lost.

**Float argument** (agrees with the new FLOAT_RIGOUR §5):
- Centre data: N ≤ 78, so ≤ γ₇₈·Pabs < 10⁻¹³·Pabs.
- Chain-rule and evaluation errors of gradient and Hessian are inside the e1/e2 pads, which enter R.
- Sums and t-combinations are inside 10⁻¹²·pabs, taken over both t-parts.
- Trig-argument errors (≤ 1.6e-14·k^j) are inside 2e-14·aab_new.
- Winding weights: (|kL| + 3)u ≤ 1.03e-14 fits the spare of the 2e-14 term, as in R4-3.

**Independent tests** (exact values at 50 or 100 digits; 0 violations everywhere):

| test | scope | result |
|---|---|---|
| `g15` (f-chart data) | coef_val2 = nbcoef_f; second derivatives; abs-polys and third derivatives; per-coefficient Taylor inequality | 0/33,600; worst 0.926 |
| `g19` (pre-patch, 500 boxes) | enclose_point2 modes 1/2/0, 96,000 each; comp_bounds2, 36,000; Hessian 2e-15 | 0 violations; tightest slack 4.5e-8 relative, the expected corner attainment |
| `g19` (patched) | 23,040 per mode; comp_bounds2 8,640 | 0 |
| `g21` (a-chart, pre-patch) | (A)–(C) data; Hessian 6e-16; Taylor 0/3,360 (worst 0.49); enclose_point2 11,520 per mode (worst 0.9999988); end-to-end 150 bulk boxes: certify_box2 139 certified (first-order 120), N_mult ≤ total at 100 digits | 0; all N_mult even |
| `g21` (patched) | same, 80 boxes; end-to-end 120 bulk boxes: 111 certified (first-order 98) | 0; Taylor 0/2,240 (worst 0.52); enclose_point2 7,680 per mode (worst 0.9999976); all N_mult even |
| `g17` (C3 KBQe_2, Ca_2) | data identical to the first-order modules; second derivatives 4e-30; abs and third derivatives; Taylor 6,240 per chart; point_enclosure2 7,680 per mode and chart | 0; smallest slack 1.1e-7 (KBQe), 1.3e-7 (Ca) |
| `g5b` (KQ winding, final engine) | 1,000 random boxes, 425 certified, 8 points each | 0; margin ≥ 0.122 |
| `g16` (KQ count, `away_count_t2` + KBQe_2) | 400 random boxes, 253 certified, 8 points each | 0; slack 0 attained (valid) |
| `g1` on KQ1, KQK2, KQB snapshots (09:00) | 506 winding leaves (**438 certify only with nb_cert2**); 150 count leaves per file (6, 26 and 3 certify only with the current code) | 0 problems; 100-digit 0 violations |
| `g13` (KQ1 at 09:30) | 10,894 winding leaves (10,813 labelled only by the second-order enclosure) | W reproduced; max \|step\| ≤ π − 2.30; wrap 6e-15 |
| `g1` on KQK1.r1 partial | 890 leaves, continuation scale | 0 problems; 80 count leaves certify only with the current code; 60 points 0 violations |

So a large part of the KQ, bulk and Ca results now rests on the second-order engines. They were checked as above.

## (7) D2E2 and the D2 coverage rules (addendum B)

**Generation.**
- `nbcoef_D2E2` against my physics, at 50 digits: max relative error 1.0e-18. Abs-polys: 2,498 terms correct. D2E checked the same way (`g20`).

**Rule (i).**
- It is sound because the D2 charts are only needed on C1p's owed boxes. Every other point of the D2 cube lies in C1p's domain (l ≥ 0.995) and is certified there.
- Exact check: 20,333 of D2Tp's 20,334 owed boxes have SD2·ν_lo·E_lo > 6.107e-6. So does every one of the 20,245 failed leaves (all at e ≥ 1.25e-3).

**Rule (ii).**
- The map (ν, E, P, K) ↦ (Eν, P/E, K/E, T = 1/E) sends a D2Tp point to the same physical point (e, p, κ, τ'). I checked this exactly to 1e-45.
- The charts differ only by positive normalisations.
- The remaining box maps into D2E for T ≤ 1.05 and into D2E2 for T ≥ 1.04; check_all's float version uses outward factors 1 ± 1e-12.
- 49 sampled points of that box (corners included) all lie in certified D2E or D2E2 leaves, with 100-digit counts ≤ the leaf bounds.
- My exact implementation in `g2` agrees: D2Tp and C1p are complete.

**D2E2 run:**
- 256 leaves, own replay OK;
- **all** leaves re-certified (232 regular, 24 winding), 0 mismatches;
- 135 spot points, 0 violations.

## (8) Continuation levels (addenda 3 and 5)

- **`kq_run.py` and `c3_run.py --boxes … --scale …`:** the scale goes to the workers (`init`) and is recorded in the header. Under spawn, workers re-import the module and get the scale through `init`.
- **Closing a level** (TMAX = −1, `--resume`) stores each unprocessed box as an incomplete record, with tree '' and `left` = [box]. This replays exactly.
- **`next_level.py`** lists the open boxes (left stacks and failed leaves) of a closed level.
- **`check_all` and my `g2`** require every continuation root to be owed by an earlier file, replay each file with its own scale, and require no owed box to remain.
- **Partial chains checked now:**
  - KQ2, KQB, KQK1 and KQK2 base files: 0 problems, owed 321, 115, 440 and 4;
  - Ca base: 0 problems, owed 478;
  - KQK1.r1 and Ca.r1: as above.
- `code/c3/rect_cert.py` is unused (not part of any certificate).

## Runs checked in this round (current data)

| run(s) | result |
|---|---|
| all finished families (C1, C1p, C1pp, B1, B2, B3±, L2·, D2E/E2/P/K/Tp/Tm, F0, KA, 13 C3 charts) | `g2`: own replay, grids, codes and bounds, owed boxes with rules (i)/(ii), exact pruning predicates: **all OK** |
| the 19 earlier bb_run families | **every leaf** re-certified with the final code (`g9`), 0 failures; D2E2: every leaf; D2Tp: 400 leaves |
| 23 families, then D2Tp, D2E2 and partial bulk | `g3`, 100 digits, own maps: 0 violations; N_mult always even; max 10 (KQK1, KQK2, D2Tp) |
| bulk partial (3,495 records) | 5,734,059 leaves, 709,555 winding: 0 replay problems, 135 spot points 0 violations, 3,000 winding leaves W reproduced |

## Final-pass scripts (`round7/final_pass.sh`)

Run `sh final_pass.sh` from `audit/round7/` once every continuation level has finished. It uses 2 cores and writes
to `round7/final/`; alerts go to `final/ALERTS.txt`. The scripts read `.json` or `.json.gz` and pick up continuation
files (`NAME.rK.json`) automatically. They were tested on the current data.

| step | script | checks |
|---|---|---|
| structure | `g2_replay_all.py` | Own replay of every family and continuation (each with its own scale): grid and duplicates, headers against `launch_runs.sh`, allowed codes and per-family bounds, −15 ≤ L < 0, owed boxes (C1p→D2, rules (i)/(ii) exact), pruned boxes (exact). |
| project checker | `g12_check_all_guarded.py` | `check_all.py` unchanged, with RESULTS.md redirected. |
| re-certification | `g9_recert_all.py` | **Every** leaf of every family except bulk and KQ*. Bulk, KQ* and continuations: 300 regular/deep leaves per family plus **every** winding leaf. Engines a and f: first-order first, then `nb_cert2` (counted). |
| KQ | `g1_kq_wind.py` | Every winding leaf and 300 count leaves of KQ* and continuations, first with the audited 07:50 code, then the final code (counted separately); 100-digit W and N_mult. |
| 100-digit spots | `g3_spot_all.py` | Own chart maps (including D2E2), all whole-circle families, KA, continuations. |
| D2 | `g20_d2e2_checks.py` | D2E2 against physics; rules (i)/(ii) with 100-digit checks inside D2E/D2E2 leaves. |
| C3 | `g4_c3_spot.py` | Own replay including zone codes; 100-digit spots; 200 leaves re-certified (current engine module). |
| winding margins | `g13_wind_margins.py` | All winding leaves of F0 and KQ* (first- or second-order labels), and 20,000 sampled from bulk. |

## Pending

- **Running:** bulk (`--order2`), KQ1.r1, KQ2/KQK1/KQK2/KQB.r2, Ca.r2, and any further levels. Then run the final pass and check `final/ALERTS.txt`.
- **Not reproduced:** the authors' A/B timings, and the "3,509 → 159 leaves" KQ1 measurement.
