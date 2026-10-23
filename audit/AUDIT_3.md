# AUDIT 3: the C3 core, the l ≈ 1 band (KB) and the float-rigour patches

Independent audit, 2026-10-23, 01:27–01:55.

**Specification:** `C3_CORE.md`.

**Scripts and outputs:** `audit/round3/`. Hashes of the audited code are in `round3/hashes_round3.txt`. The audited `nb_cert.py` is the 01:21 version.

**Rules kept.** At most 2 cores were used, and nothing was written outside `audit/`: numba caches went to `audit/numba_cache` and Python ran with `-B`.

## Verdict

**No critical or major issue was found.** Items (a)–(g) all check out. Every finished C3, KA and runs/ certificate I tested is confirmed by my own 100-digit computations.

| id | severity | item | finding |
|---|---|---|---|
| R3-1 | minor (float write-up) | trig-argument padding | The premise "\|k·u + jπ/2\| < 50" is **false for the away arcs**. In K, KB and nb_corner the arcs reach u = 3π − 0.9 (KA: 3π − PHI_A), so \|arg\| ≤ 6·8.52 + 7π/2 ≈ 62.2. With the stated 3u·\|arg\| bound that would give 2.1e-14 > 2e-14. **The 2e-14·k^j padding is nevertheless sufficient.** Count the roundings separately (k·u; j·fl(π)/2; the sum; \|fl(π) − π\|·j/2; libm 1 ulp): the error is ≤ 2u·\|arg\| + 4.3e-16 + 1.1e-16 ≈ 1.43e-14 < 2e-14. Correct the written bound. |
| R3-2 | minor (float write-up) | non-dyadic constants and sum rounding | The γ-argument holds; see "Float rigour" below. The constants argument (R2-3) is valid only if phrased through the **exact** mean-value theorem for the exact polynomial. The write-up should also cover the rounding of the trig and power sums s, g_v, w_v. It is covered by 1e-12(mag + tot) and by the spare in `dwid`'s 1e-13·a1 pad. |
| R3-3 | minor (code hygiene) | `c3_run.phys_bounds` | 'L2aYm' and 'L2KYm' end in 'm', so they take the m-chart branch. Their μ "lower bound" ν₀³·M2·η₀² ignores that mp ∈ [0, R2²] is the root variable, so it is not a lower bound. **It cannot fire:** it is ≤ 0.37³·0.01·1.0002 = 5.1e-4 < MUMAX = 0.0505, and 'o' never triggers from it. Fix: test `ch[3:] == 'm'`. |
| R3-4 | minor (docs) | C3_CORE.md, STATEMENT.md, working notes | **Missing reference:** `PROOF.md` is cited but does not exist yet.<br>**Wrong factor:** C3_CORE.md §1 has "F∓ = (2σD³/(κ(1+e)))·g∓′". The verified identity (AUDIT_1 a9) is F∓ = g∓′·D³σ/(μ(1+e)) = κD³σ g∓′/(1+e), so the stated factor is off by 2μ². Only positivity matters.<br>**Stale location:** LOG.md line 53 still puts the √2 parabola gap "at k ≈ 2"; it is the limit k ↓ 1. |
| R3-5 | info (completeness) | runs | The C3 core is **not complete**: Ca, CE and Cm have not finished. Ca is the only chart that needs the level-2 zone. KB and the bulk, B1, C1, C1p, C1pp, F1, F2, K and L2E runs are still running. **Near-curve leaves (R2-2):** the fallback now calls deep, but deep only certifies at h ≈ 1e-5, i.e. normalised width ≈ 2e-5 = `minw`. So a few boxes near ω = π/2, e = l², t → 1 may still be recorded as failed. Check bulk's fails, or lower `minw` there. |

---

## (a) C3 coverage logic: OK

- **Cell radius:** on the cell, δ = max(Q_z, √κ, √ε, μ^{1/3}) ≤ max(0.15643, 0.24839, 0.1, 0.0505^{1/3}) = **0.369626 < NU1 = 0.37**.
- **Sampling** (200,000 physical points, max δ 0.369626; `c4_coverage.py`):
  - **near** |τ| ≤ 3δ: the dominant level-1 chart (Ca/CK/CE/Cm, ν = δ, ratios ≤ 1) contains the point;
  - **intermediate** 3δ ≤ |τ| ≤ 1.11: PhiP/PhiM with ν = |τ|/3 ∈ [δ, 0.37] contains the point, and the root m = μ/ν³ ≤ 1;
  - **away** |τ| ≥ 1.11: the point lies in KA's domain.
  - **0 membership failures.**
- **No gaps in τ.** The three ranges partition the τ-line. Since ν is a continuous box variable on [0, 0.37], PhiP/PhiM cover |τ| ∈ [3δ, 1.11] without gaps.
- **Away arc.** PHI_A = 2·atan(3·0.37) − 1e-9 = 1.6749674242 < 2·atan(1.11) = 1.6749674252, so the near/intermediate part and KA overlap by 1e-9.
- **Total.** near ≤ 8, intermediate 0, away ≤ 2. All three count zeros of the same function: Q(τ) = P(f)(1+τ²)⁶ and the f-chart polynomial with t = 1/(1 + μ²) are positive multiples of F₋F₊ (AUDIT_1 a9).
- **Superset in (Q_z, κ).** The C3 charts cover the cell's (Q_z, κ) superset; KA covers the (l, ω) cell exactly. A point in the superset but outside the (l, ω) cell (e.g. ω = 1.41 < 9π/20, l = 0.951) is covered entirely by KB, so no gap arises.

## (b) Level-2 zone logic: OK

- **Exact identities** (sympy, by hand):
  - at E = m = 0, P̂₀ = −(x⁴/4)·x²·(x² − 3ax + 4K)²;
  - the leading level-2 part is 32K·y²[4m² − 2K(E + 2y²)²], so at most 4 roots near y = 0 (K > 0 on physical Ca/CK points, because the leading-order l ≤ 1 condition is a² ≤ 2K).
- **Inner roots for a zone point p** (E ≤ E2, m ≤ M2, so η_p ≤ 1):
  - roots with |y| ≤ Y2·η_p lie in L2·E or L2·m (η = η_p, ratio ≤ 1, Y ∈ [−1, 1]), target 4;
  - roots with Y2·η_p ≤ |y| ≤ Y2 lie in L2·Y± (η = |y|/Y2 ∈ [η_p, 1] ⊂ [0, ETA1], Ep, mp ≤ 1), which certify no zero.
  - With the outer count on |y| ∈ [Y2, Y0] (≤ 4), the near count is ≤ 8.
- **Physical points have η_p > 0** (ε > 0), so the 6-fold limit point E = m = 0 is never needed.
- **Sampling:** 61,735 zone samples, 0 failures (c4).
- **Cone ranges:** η ≤ 1.0001 and ratios ≤ 1.001 (squared or cubed where weighted) contain all exact ratios ≤ 1.
- **CK needs no zone.** It produced no 'Z' leaves. That is consistent: in CK, K = 1, and the double roots are real only for a² ≥ 16/9, while a ≤ 1.001.

## (c) Pruning in `c3_run.phys_bounds`: OK (except the harmless R3-3)

- **'u' test, Ca and L2a.** The test is "1 ≤ 2K − ν²K² fails everywhere". On the domain (ν²K ≤ 0.137 < 1), 2K − ν²K² is increasing in K and decreasing in ν. So its maximum is at (ν₀, K₁), as coded.
- **'u' test, CE, Cm, Phi.** These skip only if a₀² > max(2K − ν²K²).
- **'u' test, CK and L2K.** Never unphysical, since a² ≤ 1.002 < 2 − ν².
- **'o' test.** It uses valid lower bounds for each chart, with a relative margin of 1e-9.
- **Physical points are never at ν = 0.** Every physical point has ν = δ(p) > 0 (ε > 0, μ > 0).
- **Empirical check** (`c3_prune_check.py`): sampling the corners, random points and root values of every pruned leaf in the 11 finished C3 runs gave about 514,000 checks: **0 'u' leaves contain an l ≤ 1 point, and 0 'o' leaves contain a cell point.**

## (d) `gen_c3` against C3_CORE.md: OK

- **Substitutions.** Read line by line, they are exactly the ones in §2–§4:
  - Ca, CK, CE and Cm (dominant variable = ν);
  - PhiP and PhiM (τ = ±3ν, root m);
  - the level-2 E, m and Y± charts (E = E2η², m = M2η², y = ±Y2η).
  - The exact ν¹⁰ and η⁶ divisions are asserted.
- **Self-test maps.** The physical maps in `phys()` agree with my own maps.
- **Generated code against my own physics** (`c5_c3_coeffs.py`, 100 digits; values of Σ c_k r^k against Q(τ)/norm built from the trig formula):
  - Ca, CK, CE, Cm, PhiP, PhiM, KBQe: ≤ 3.4e-99;
  - the 8 level-2 charts: ≤ 3.4e-17, which is the non-dyadic constants 0.01 (R2-3), covered by R3-2;
  - gradients are consistent;
  - 21,831 abs-polynomial terms are exact.
- **c3_sym claim confirmed.** G(y; ν) − P̂₀(2y) shrinks by exactly 100× per 10× in ν. So it is O(ν²), with only even powers.

## (e) KB engine: OK

- **`corner2.check_plus_sub` is the same Lemma N2.**
  - β₀ = −max over 32 sub-arcs of B′ (upper bound) and |B′|max = max over the sub-arcs; these are valid, tighter bounds.
  - The sub-arcs are padded and overlap, covering [π − φ₀ − pad, π + φ₀ + pad], and σ² > 0 is checked on each.
  - Conditions (2)–(4) are identical to `nb_corner.check_plus`, which AUDIT_1 and AUDIT_2 accepted.
- **Arcs overlap:**
  - near |φ| ≤ 0.1 + pad (F₋ Rolle; F₊ by N2 or Rolle);
  - middle τ ∈ ±[tan 0.05 − 1e-9, tan 0.45 + 1e-9], i.e. φ ∈ [0.1 − 2e-9, 0.9 + 2e-9] (Q(τ) in ε, monomial basis);
  - away |φ| ≥ 0.9 − 1e-9 − pad (f-chart).
  - σ > 0 on the near arc: an l = 1 orbit crosses the axis at |φ| = π/2 − ω ≥ π/20 > 0.1, and the code also checks σ² > 0.
- **Outward rounding.** ε = 1 − e is outward-rounded (exact by Sterbenz, then padded). P_z = l sin ω and Q_z = l cos ω use 2-ulp-padded libm values (sin increasing, cos decreasing on [3π/20, 9π/20]) and products rounded outward.
- **Empirical check** (`c7_kb_check.py`, 30 random band boxes):
  - 14 certified (totals 6–10), and no sampled point exceeds its total (100-digit count of the real roots of P on the circle);
  - Lemma N2 accepted 26 boxes, and F₊ changed sign at most once on |φ| ≤ 0.1 at every sampled point.
- **Band width.** KB covers ω ∈ [3π/20, 9π/20]. That contains the whole range ω ≥ π/2 − 0.9 = 0.671 where the l = 1 axis crossing enters K's near arc.

## (f) Coverage of e ≥ 0.95: OK, by design (runs pending)

- **For e ∈ [0.99, 1] and t ≥ 0.9975**, i.e. μ ≤ 1/√399 = 0.0500626 ≤ 0.0505:
  - K covers all (l, ω) except its 28 skipped initial boxes (l ∈ [0.95, 1], ω ∈ [3π/20, π/2]);
  - KB covers l ∈ [0.95, 1], ω ∈ [3π/20, 9π/20];
  - the C3 charts plus KA cover l ∈ [0.95, 1], ω ∈ [9π/20, π/2].
- **For e ∈ [0.95, 0.99]:** F1 (t ≤ 0.9975) and F2 (t ≥ 0.9975).
- **For e ∈ [0.9499, 0.95]:** the bulk (e ≤ 0.9501) and F2.
- **Float seams:**
  - K's skip starts at l = 0.9500000000000001 ≥ KB's 0.95 (overlap);
  - ω = 3π/20 is the same float in K's grid and KB's domain;
  - 9π/20 is the same float in KB, KA and the C3 literal 1.413716694115407;
  - F2 and K meet at e = 0.99, the same float.
- **Sampling** (400,000 points with e ≥ 0.94, including points within 1e-6 of each seam): **0 uncovered** (c4).

## (g) Replays and spot checks: OK

- **`bb_replay.py --recert 40` (fixed R2-1 dispatch)** on runs/B2, B3p, B3m, L2P, L2Ap, L2Am, L2Tp, L2Tm and KA: **0 replay problems and 0 re-certification problems** (`c1_replay_recert.out`).
- **C3 runs** (my own replay of the `c3_run` tree format, `c2_c3_replay_spot.py`): PhiP, PhiM, CK and all 8 level-2 charts have 0 replay problems and 0 fails.
- **100-digit spot checks** (25 leaves × 3 points per run = 825 points).
  - Method: own chart maps, then Q(τ) from the trig formula, then interpolation in the root variable (residual ≤ 2.4e-93), then real roots with multiplicity on the certified interval.
  - **0 violations.** The largest count was 4 on CK and the L2·E/m charts, and 0 roots on Phi and L2·Y±.
- **KA leaves** (`c8_ka_spot.py`, 64 leaves × 3 points): **0 violations**. The maximum is 2 real roots on the away arc, so the bound of 2 is attained.

## Float rigour (authors' item 1)

| piece | verdict |
|---|---|
| Evaluation error ≤ γ_n·Pabs, crad = 1e-13·Pabs(\|x_c\|) | **OK.** Operation counts (`c6_opcounts.py`, all 31 generated modules): ≤ 265 terms, total degree ≤ 40, power chains ≤ 17 multiplications, so n = deg + T − 1 ≤ 300 and γ₃₀₀ + u = 3.34e-14 < 1e-13. |
| Non-dyadic constants (R2-3) | **OK with the right argument.** Apply the mean-value form to the exact polynomial: P(x) = P(x_c) + ∇P(ξ)·(x − x_c). Then \|P(x_c) − cmid\| ≤ (γ₃₀₀ + u)·Pabs(\|x_c\|) ≤ crad. For the gradient, \|∇P(ξ) − dcen\| ≤ [Pabs_∇(\|x_c\| + rr) − Pabs_∇(\|x_c\|)] + (γ₃₀₀ + u)(Pabs_∇ terms) ≤ `dwid`, using its 1e-13·a1 pad. The constants' rounding never needs to be covered at box points, only at the centre. |
| Trig argument, 2e-14·Σ k^j(…) | **Sufficient**, but the written bound needs correcting (R3-1). |
| Final sums, 1e-15(\|vc\| + tail) | **OK.** About 4 operations, error ≤ 4.4e-16 relative. |
| Mean-value radii padding | **OK.** It is present in nb_cert, make_engines and make_engines_bu (01:21), and in the c3 engines. |
| `root_signs` (c3) | **OK.** No trig. The power sums have relative error γ_≈14, covered by 1e-12·mag and the `dwid` spare. The Taylor bound uses (\|r_m\| + h)^(k−j−2). |

## AUDIT_2 items

| item | status |
|---|---|
| R2-1 `--recert` dispatch | **fixed**, verified on 9 runs |
| R2-2 deep fallback at tot = 999 | **fixed** (c9: deep is now called and certifies at h ≤ 1e-5); see the `minw` caveat in R3-5 |
| R2-3 constants | **resolved** by the argument in R3-2 |
| R2-4 wording | **fixed** in STATEMENT, LOG, CORNER_E1 and DEGENERACY_TREE; R3-4 leftovers |

## Not checked / pending

- **Unfinished runs:** Ca, CE and Cm (C3), KB, bulk, B1, C1, C1p, C1pp, F1, F2, K and L2E.
- **PROOF.md** (not yet written).
- **C3 zone leaves.** The full level-2 zone accounting has only been exercised by CK, which needs no zone. Ca is the decisive run. When it finishes, rerun `round3/c2_c3_replay_spot.py` on `c3/runs/Ca.json`. Its spot check handles zone ('A'+c) leaves: outer ≤ c and inner ≤ 4.

## Scripts (`audit/round3/`)

| script | checks |
|---|---|
| `c1_replay_recert.out` | output of `bb_replay --recert` (9 runs) |
| `c2_c3_replay_spot.py` | C3 replay and 100-digit leaf checks |
| `c3_prune_check.py` | 'u'/'o' pruning |
| `c4_coverage.py` | (a)(b)(f) membership and coverage |
| `c5_c3_coeffs.py` | C3 coefficients against physics; c3_sym leading order |
| `c6_opcounts.py` | γ_n counts |
| `c7_kb_check.py` | KB and N2 sub-arc check |
| `c8_ka_spot.py` | KA leaves |
| `c9_deep_fallback_v2.py` | R2-2 fix |
