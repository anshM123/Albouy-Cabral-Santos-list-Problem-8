# AUDIT 6: engine KQ (final e ≥ 0.95 architecture), parity everywhere, min over rounds

Independent audit, 2026-10-23, 07:31–07:50.

**Scripts and outputs:** `audit/round6/`. Hashes of the audited code are in `round6/hashes_round6.txt`.

**Rules kept.** At most 2 cores were used, and nothing was written outside `audit/`.

## Verdict

**No soundness issue was found.** KQ is a clean certificate:
- both parts count zeros of the same P with multiplicity;
- the near and away parts overlap;
- μ, ε, P_z and Q_z are rounded outward;
- T_MIN is enforced.

Parity is correctly applied to whole-circle totals and to the C3 composite (≤ 9 + 0 + 2). The final atlas covers e ≥ 0.95 (500,000 samples, 0 uncovered). Every finished or partial run I checked passes replay, re-certification and 100-digit spot checks.

| id | severity | finding |
|---|---|---|
| R6-1 | minor (hygiene) | `check_all.PARITY_ENGINES` includes **KA**. KA is a partial (away-arc) certificate with target 2, so 'l' can never occur for it; but accepting "total 11" there would be logically wrong. K and KB are listed too, although they are no longer in the atlas. Remove KA (and K, KB). |
| R6-2 | minor (hygiene) | `bb_replay.py` accepts 'l'/'L' for **any** engine. It should use the same whitelist as `check_all`. `bb_replay` also has no `--recert` dispatch for engine 'KQ'; I re-certified KQ leaves with `kq_run.certify` directly (f5). |
| R6-3 | info | The round-5 tooling issues are fixed: 'l'/'L' are accepted (R5-1), and `bb_replay`'s winding dispatch includes D2E/D2P/D2K (R5-2). K and KB were never completed and are not in the final atlas, so their encoding issue is moot. |

## (a) KQ certificate: OK

- **Same function in both parts.**
  - **near** |φ| ≤ 0.9: Rolle runs in τ ∈ [−tan 0.45 − 1e-9, tan 0.45 + 1e-9] on Q(τ) = P(f)(1+τ²)⁶ (chart KBQe in ε = 1 − e).
  - **away** |φ| ≥ 0.9 − 1e-9 (minus pad): the f-chart polynomial on arcs, via `c3_away.away_count_t`.
  - P = μ²k²sin²fσ² − D²B² = F₋F₊, and the f-chart polynomial (1−t)(PhiF − κ²PsiF) = (1−t)κ²·P is a positive multiple of it for t < 1. Both parts therefore count zeros of P with multiplicity, and Q(τ) preserves multiplicities.
  - The near τ-range reaches |φ| = 0.9 + 1.6e-9, so the two parts overlap.
- **Parity.** total = aw + q ≥ N_mult(P), which is even, so ≤ 11 gives ≤ 10. KQ uses no branch functions, so it needs no σ > 0 condition: axis crossings are simply double roots of P.
- **Rounding.**
  - **μ:** μ(t) = √((1−t)/t) is decreasing in t. mu1 = up(√((1−t₀)/t₀)) and mu0 = dn(√max((1−t₁)/t₁, 0)), with t₁ clamped to 1 and mu0 ≥ 0. 1 − t is exact for t ≥ 0.5 (Sterbenz), and division plus sqrt cost ≤ 1.5 ulp, inside the 2-ulp padding.
  - **ε:** [dn(1 − e₁) ∨ 0, up(1 − e₀)].
  - **P_z, Q_z:** 2-ulp-padded libm sin/cos (monotone on [0, π/2]) with outward products. The box in (P_z, Q_z) is a superset of the image of (l, ω), which is safe.
- **T_MIN.** `certify` refuses t₀ < 0.5. The KQ1 domain starts at t = 0.5 and F0 ends at t = 0.5 (shared edge), so μ ≤ 1 throughout.
- **KBQe polynomial.** Its coefficients were checked against my own physics in AUDIT_3 (≤ 5.5e-101). The `count` engine is AUDIT_3's `root_signs` and `dp_linear`, now with min over rounds.
- **Re-certification** (`f5_kq_recert.py`, `kq_run.certify`) on completed-initial-box snapshots of KQ1, KQ2, KQK1 and KQK2: **569 leaves, 0 problems.**
- **100-digit spot checks** (`f1_kq_spot.py`, own formula for P, N_mult(P) on the whole circle):

  | run (partial) | leaves | spot points | violations | max N_mult |
  |---|---|---|---|---|
  | KQ1 | 7,325 | 120 | 0 | 8 |
  | KQ2 | 428 | 120 | 0 | 6 |
  | KQK1 | 2,138 | 120 | 0 | 6 |
  | KQK2 | 119 | 120 | 0 | **10** (the bound 11 → 10 is tight) |

  Every N_mult was even.

## (b) Coverage of e ≥ 0.95: OK (`f4_coverage_final.py`)

**Domains** (from `launch_runs.sh`):

| region | covering run(s) |
|---|---|
| e ≤ sin 1.25355… = 0.9501 | bulk |
| t ≤ 0.5 | F0 |
| 0.5 ≤ t ≤ 0.9975 | KQ1 |
| t ≥ 0.9975, e ≤ 0.99 | KQ2 |
| t ≥ 0.997456, e ≥ 0.99, l ≤ 0.95 | KQK1 |
| t ≥ 0.997456, e ≥ 0.99, l ≥ 0.95, ω ≤ 3π/20 | KQK2 |
| t ≥ 0.997456, e ≥ 0.99, l ≥ 0.95, ω ∈ [3π/20, 9π/20] | KQB |
| e ≥ 0.99, l ≥ 0.95, ω ≥ 9π/20, μ ≤ 0.0505 | C3 + KA |

- **Seams use identical float literals:** 0.5; 0.9975; 0.997456; 0.95; 0.47123889803846897 (3π/20); 1.413716694115407 (9π/20, also KA's ω₀).
- **C3 threshold.** C3's μ ≤ 0.0505 means t ≥ 0.99745624, which is ≥ 0.997456. The strip t ∈ [0.997456, 0.99745624) lies inside KQ1 (t ≤ 0.9975).
- **Sampling.** 500,000 points with e ≥ 0.94, stressing every seam to ±1e-7 and t to within 1e-14 of 1: **0 uncovered.**

## (c) C3 parity: OK

- **Targets.** Ca now has near target 9: non-zone leaves count y ∈ [−3, 3] ≤ 9, and zone leaves count the outer part ≤ ZONE_OUTER = 5, plus the inner part ≤ 4 from the L2 charts. CK, CE and Cm used target 8.
- **Total over the whole circle.** For every p in the cell: near ≤ 9 (|τ| ≤ 3δ, Rolle on G(y) with multiplicity), intermediate 0 (Phi charts), away ≤ 2 (KA, Rolle on the f-chart polynomial with multiplicity).
- **Parity.** These three parts cover the circle, f = 0 included in the away arc, and count zeros of the same real trigonometric polynomial P with multiplicity. So the total ≥ N_mult(P), which is even: ≤ 11 gives ≤ 10.
- **Ca** had not finished. Rerun `round3/c2_c3_replay_spot.py` on it; that script checks zone leaves (outer ≤ c, inner ≤ 4).

## (d) Parity acceptance in the checkers: OK, except R6-1 and R6-2

| checker | behaviour |
|---|---|
| `check_all` | Accepts 'l'/'L' only if the run's `engine` is in PARITY_ENGINES: all whole-circle Rolle and winding engines, D2*, KQ; also KA, K, KB (R6-1). C3 files carry `chart`, not `engine`, so 'l' is rejected for them, and they never produce it anyway. |
| `verify_mp` | Treats 'l'/'L' as an effective bound of 10. |
| `bb_replay` | Accepts 'l'/'L' as 11 everywhere (R6-2). |

## (3) Min over rounds: OK

`c3_away.away_count_t` and `c3/make_c3.count` now keep best = min(best, dp) over rounds. Each round recomputes every sign on a complete overlapping cover, so each round's dp bound is valid on its own, and so is the minimum.

## (e) Replays and spot checks of finished runs: OK

| run | result |
|---|---|
| F0 (f-chart, t ≤ 0.5) | `bb_replay --recert 200`: 0 replay problems, 0 re-certification problems. 400 records, 12,354 leaves (3,925 winding, 1,761 at total 11), 0 fails. 100-digit spot check: 120 points, 0 violations, max N_mult 6, all even. |
| D2Tm (re-run with parity) | `bb_replay --recert 200`: 0 problems. Own full re-certification of all 256 leaves (10 at total 11): 0 mismatches, 0 spot violations (`f3_d2tm.out`). |
| KQ1, KQ2, KQK1, KQK2 (partial) | see (a) |

**Earlier finished runs** are unchanged since AUDIT_5: C1, C1p, C1pp, B1, B2, B3±, L2*, D2E/P/K, KA, and the C3 charts except Ca.

## Pending

- **Unfinished runs:**
  - **KQ1, KQ2, KQK1, KQK2, KQB:** rerun `f1_kq_spot.py` and `f5_kq_recert.py` on the final JSONs.
  - **bulk (fixed wind_cert):** `bb_replay --recert`, plus `round4/d5_step_margins.py` or `round5/e3_r4_fixes.py` for the winding leaves.
  - **Ca:** `round3/c2_c3_replay_spot.py`.
  - **D2Tp:** `round6/f3_d2_recert_parity.py`.
- **Final `check_all` report.** Once everything is finished, run it once over all families.
