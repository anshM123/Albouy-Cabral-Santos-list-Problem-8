# AUDIT 5: R4 fixes, the D2b corner, the parity lemma in K/KB, splits, replays

Independent audit, 2026-10-23, 03:13–03:40.

**Scripts and outputs:** `audit/round5/`. Hashes of the audited code are in `round5/hashes_round5.txt`.

**Rules kept.** At most 2 cores were used, and nothing was written outside `audit/`.

## Verdict

**No soundness issue was found.** The R4 fixes are correct and verified on every existing winding leaf. The D2b derivation, divisibility, cone coverage, pruning and the C1p-fails ⊂ cube claim are all correct. The parity argument for K and KB is valid.

**Two tooling gaps** would make the final checks fail or crash. Both need fixing before the K, KB and D2 runs can be closed.

| id | severity | finding |
|---|---|---|
| R5-1 | **major (tooling, blocks verification)** | Under the parity lemma, K and KB now certify leaves with **total 11**. `bb_run` (engine K) and `kb_run` encode such a leaf as `chr(ord('a') + 11) = 'l'`. But `check_all.py` accepts only 'a'..'k' and reports **"bad code l"**, and `bb_replay.py` raises `SystemExit('bad code l')`. The K1, K2 and KB runs had not finished, so no 'l' code was observed yet, but such leaves are the purpose of the change. **Fix:** accept 'l' as "total 11, hence ≤ 10 by parity" **only** for the K and KB families, and keep it invalid for every other engine. |
| R5-2 | minor (tooling) | `bb_replay.py --recert` crashes on D2E, D2P and D2K (`TypeError: 'NoneType' object is not callable`). `bb_run` uses the winding certificate for these engines, but the `bcf` dispatch in `bb_replay` lists only a, f, C1*, B1, B2 and L2E/P/A±. **Fix:** add 'D2E', 'D2P' and 'D2K' to that list. The structural replay works (0 problems). I re-certified every D2 leaf with my own script (below). |
| R5-3 | minor (hygiene) | `check_all` drops C1p's owed boxes once the five D2 files merely **exist** (`d2_ok`). This is sound only because the D2 families are checked in the same report. Better to require that D2 is complete: no fails, nothing incomplete, nothing owed. |

## (a) D2b corner: OK

- **Symbolic check** (`e4_d2_checks.py`, sympy, exact):
  - **Base value:** C1p's function F = (1−τ)Φ/e² − τ(Φ + Ψ) at e = p = κ = 0 equals 4s²c²(1 − 2τ), so it vanishes identically at τ = ½.
  - **Linear part:** the first-order part in (e, p, κ, τ′) is exactly −8τ′s²c² + 2κ(2s⁴ + 4s²c²) − 4p·s·c³ − 4e·s²c³. I also derived each partial derivative by hand.
  - **Independence:** the four functions have rank 4. The linear part therefore vanishes only at the origin, and one homogeneous blow-up resolves the corner.
- **Divisibility.** `gen_bu.py` asserts that F_C1p(ν·x) vanishes at ν = 0, with b reduced mod b² = 1 − e². At the base point G0 = Φ/e² − Φ − Ψ = 0 and G1 = 8s²c², which confirms it.
- **Generated code** (50 digits, own chart maps): nbcoef_D2E/P/K/Tp/Tm match the physics f·(1−τ)/(e²ν) to ≤ 1.5e-18. That gap is the non-dyadic constants (5e-3, 1.05), the R2-3 effect, which the exact mean-value argument covers. The 10,742 abs-polynomial terms are exact.
- **Cone coverage.** The D2 runs use ν ∈ [0, 1.01] and ratios ∈ [0, 1.05], with T ∈ [−1.05, 1.05] for the affine charts. 200,000 points of the cube {e, p, κ ≤ SD2, |τ′| ≤ SD2} all lie in a chart domain.
- **Pruning (`bb_run.prune`).** For each D2 chart it uses p_lo = SD2 · (min ratio or 1) · ν_lo and κ_hi = SD2 · (max ratio or 1) · ν_hi. The ratios and ν are nonnegative and κ ≤ 1, so min over the box of p² + (1 − κ)² ≥ p_lo² + (1 − κ_hi)². "> 1 + 1e-12" therefore implies l > 1 everywhere, as required. D2E, D2P and D2K produced no 'u' leaves.
- **C1p failures ⊂ cube.** C1p has 659 failed leaves, 0 incomplete boxes, and none of them owed beyond D2. Mapped through e = ε, p ≤ sin ω_hi, κ ≤ 1 − l_lo·cos ω_hi and τ − ½, they satisfy e ≤ 6.11e-6, p ≤ 6.71e-4, κ ≤ 2.01e-3 and |τ′| ≤ 2.14e-3. That is inside the cube with a factor ≥ 2.3 to spare. `check_all`'s drop test uses exactly these bounds (with p ≤ sin ω_hi since l ≤ 1). See R5-3.
- **Same function.** The D2 charts certify the same function as C1p up to the positive factor ν, so the root counts transfer.

## (b) Parity lemma in K and KB: OK (argument valid; see R5-1 for the encoding)

**Every part bounds zeros of the same real trigonometric polynomial P = F₋F₊ (degree ≤ 6), with multiplicity:**
- **near F₋:** Rolle runs on [π − φ₀ − pad, π + φ₀ + pad]. Since σ² > 0 is checked on every near arc, F∓ are smooth there and mult_P = mult_{F₋} + mult_{F₊}.
- **near F₊ (Lemma N2):** at most one zero, and it is **simple**, because F₊′ < 0 on the interval between 0 and φ_B, which contains every zero.
- **outer F₊:** Rolle runs, with multiplicity.
- **middle (KB):** Q(τ) = P(f)(1+τ²)⁶ has the same roots and multiplicities.
- **away:** the f-chart polynomial (1−t)(PhiF − κ²PsiF) is a positive multiple of P, with t = 1/(1+μ²).

**Conclusion.** Overlapping arcs only overcount, so the certified total is ≥ N_mult(P) on the circle. That number is even for a real trigonometric polynomial ≢ 0. So total ≤ 11 gives N_mult(P) ≤ 10, and N ≤ N_mult(P) by the reduction lemma, since P is a positive multiple of Pol composed with the anomaly map.

**Code.** `nb_corner.certify_corner_box` uses budget = 11 − nm − 1 − po and returns tot ≤ 11. `kb_run` uses LIMIT_P = 11, and its two-stage evaluation only changes refinement targets. The C3 partial counts (near ≤ 8) correctly do **not** use parity, since parity holds only for full-circle totals.

## (c) R4 fixes: OK

- **`arg_in_branch`** measures each arc's two endpoints in the label's continuous branch: Re>0 (atan(im/re)), Im>0 (atan2 with im > 0), Re<0 (π + atan(im/re)) and Im<0 (atan2 with im < 0). There is no normalisation step.
  - **Rejection:** values outside the half-plane are rejected (sentinel 99).
  - **Exactness:** at each shared node the computed and true values lie in the convex intersection of the two half-planes, so the per-arc branch errors cancel exactly.
  - **Wrap-around:** the last arc reuses the value at left[0] (R4-2).
- **Crafted antipodal case** (AUDIT_4 d4) through the new `arg_in_branch`: total −2π, winding −1. Correct.
- **Guards:** L = 0.1 and L = −20 are rejected, and W ∈ [1, 6] is enforced. F ≡ 1 gives W = 6, correct since Q = z⁶. `choose_L` clips at −15.
- **Every winding leaf re-certified with the fixed `wind_box`:** W reproduced in all of them (0 mismatches).

  | run | winding leaves |
  |---|---|
  | C1p (certified with v1) | 24,550 |
  | C1 | 34,471 |
  | C1pp | 6,138 |
  | B1 | 5,171 |
  | L2E | 290 |
  | D2E, D2P, D2K | 587 |

- **Synthetic stress test** on the fixed code: 4,000 polynomials, 1,579 certified, 0 wrong. Q ≡ 0 is rejected.

## (d) Replays and re-certification: OK

| run | result |
|---|---|
| D2E, D2P, D2K | Own replay (`e2_d2_recert.py`). **Every** leaf re-certified (256 + 1,130 + 464): regular leaves reproduce their bound exactly, 0 deep-leaf mismatches, and the 587 winding leaves reproduce W. 180 spot points at 100 digits (own physics, real roots ≤ bound and roots inside \|z\| < e^L = W): 0 violations. No pruned leaves. |
| C1p | `bb_replay --recert 200`: 0 replay problems, 0 re-certification problems (winding leaves re-certified with the fixed code). 44,140 leaves (24,550 winding, 2,121 deep), 659 fails, all owed to D2. |
| C1, C1pp, B1, L2E | All winding leaves reproduced with the fixed code. Their regular and deep leaves were already replayed and re-certified in AUDIT_4. |
| C3 CE, Cm | Own replay plus 100-digit spot checks (`e5_c3_CE_Cm.out`). CE: 171,636 leaves; Cm: 2,861 leaves; 0 replay problems; 150 points with 0 violations. Pruning: about 1.06 M samples, 0 violations. |

**Splits (item 4): OK.** The pieces tile the old domains on identical float edges, and the K exclusion is unchanged in K1 and K2:
- F1 = F1a (l ≤ 0.8) ∪ F1b1..4, with e-slabs 0.95 / 0.9625 / 0.975 / 0.9875 / 1.0;
- F2 = F2a ∪ F2b (split at l = 0.8);
- K = K1 (e ≤ 0.995) ∪ K2 (e ≥ 0.995).

**Chart 'fe' (item 5)** is unused in every run, so it was not audited.

## Not checked / pending

- **Unfinished runs:** D2Tp, D2Tm, F1a, F1b1..4, F2a, F2b, K1, K2, KB, bulk and Ca (C3). When they finish:
  - D2Tp/D2Tm: rerun `round5/e2_d2_recert.py`;
  - the f-chart and bulk runs: `bb_replay --recert` and `round4/d2_wind_spot.py` (the bulk winding leaves need the fixed code);
  - K, KB: their checkers first need R5-1;
  - Ca: rerun `round3/c2_c3_replay_spot.py`.
- **'l' codes:** K and KB will need R5-1 fixed before `check_all` and `bb_replay` can pass.
