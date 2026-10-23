# AUDIT 1: Problem 8, Q2 (circle–ellipse, at most 10 stationary points)

Independent adversarial audit, 2026-10-22, 23:05–23:55.

**Scope:** the items 1–7 of the audit request, plus two follow-up questions from the authors received during the audit:
- the run change to bulk_a2 and the axis-aligned claim;
- the `certify_box` stall rule, `fails.jsonl`, the C1/C1p/C1pp charts, the corner engine, and the H1 claim.

**Method.** Every check below uses my own code in `audit/`. I re-derived and re-computed independently and did not reuse the project's numerics. The only exception is where the object under test *is* project code, such as `certify_box`.

**Resources.** At most 2 cores were used. Nothing was written outside `audit/`: numba caches went to `audit/numba_cache` and Python ran with `-B`.

**The code changed under me.** `nb_cert.py` changed twice during the audit:

| version | hash prefix | time | change |
|---|---|---|---|
| v1 | `a0c726` | 22:20 | the version I read and hashed (`hashes_audited.txt`) |
| v2 | `d29ba8` | 23:36 | stall rule and 4096-arc cap in `certify_box` |
| v3 | `fbf634` | 23:45 | adds `certify_box_deep`; `certify_box` identical to v2 |

`box_coeffs`, `arc_signs` and `dp_count` are byte-identical in all three versions.
- **Tests on v1:** a4, a5, a6, the first a6b, and the first a6c (`a6c_curve6.out`).
- **Repeated on v2:** a6 and a6b (`*_newcode.out`), plus a6c (second run), a6d, a6e, a13 and a14.
- **No mix-up:** a13 shows the engine is deterministic. The changed results came from the edit, not from noise.

---

## Verdict

- **No critical issue was found.**
  - No certificate returned a bound that a sampled configuration exceeded.
  - Every generated formula matches an independent derivation exactly.
  - The reduction lemma, its multiplicity accounting, the symmetry reductions and the CORNER_E1 formulas are all correct.
- **Major issues (4):**
  - (M1) a coverage hole at the seam e = 0.2;
  - (M2) the certified leaves are not recorded, although the docs say they are;
  - (M3) floating-point rigour is not yet established;
  - (M4) the δ ≥ 0.836 margin claim is false.
- **Minor issues:** listed in the issue table below.
- **Not a proof yet.** Nothing is complete:
  - **Bulk:** the bulk run is in progress and already has structural failure clusters.
  - **Missing regions:** e > 0.95 with κ > 3, the C3 sub-corner, and the seams are not covered.
  - **Rigour:** the float rigour is not established.

### Issue list

| id | severity | item | one-line description |
|---|---|---|---|
| M1 | **major** | coverage | Seam hole: e ∈ (0.2, 0.2000020], l ∈ [0.408, 1] lies in no chart. Float-level seam slivers exist too (θ_max = fl(arctan 2) against epp ≤ 0.5; the bulk ends at e = 0.9499982). |
| M2 | **major** | nb_run*.py | "Every certified leaf is recorded (box + certified count)" is false: only counts, volumes and histograms are stored. The remaining stack of an incomplete box is truncated to 50 entries (20 in `fails.jsonl`), so `recheck_fails.py` cannot complete it. |
| M3 | **major** | nb_cert, nbeng, nb_corner | Floating-point rigour is not established (LOG acknowledges this). Arcs leave float gaps (up to 8.9e-16 measured), and the last arc ends at fl(2π) < 2π. Trig arguments, MVT radii and a1−a0 differences are padded heuristically. All effects are ~1e-15, far below the observed margins. |
| M4 | **major** | STATEMENT / LOG claim | "Off-circle gap δ ≥ 0.836 … a robust margin" is false. δ = 0.5557 at a genuine N = 10 configuration with e = 0.94998 (inside the bulk range; 60-digit check, 2-D count 10). δ → 0 as e → 1 at small s (0.28 at e ≈ 0.99). |
| m1 | minor | LOG claim | Parabola–circle smallest complex-pair gap is √2 ≈ 1.414 (l = 1, ω = 0, k ↓ 1), not √3 ≈ 1.732. Exact factorisation below. "At most 8" is unaffected. |
| m2 | minor | CORNER_E1 (C3 paragraph) | "B = 2 sin f · D at P_z = 1, Q_z = 0" holds only at e = 1. In general B = 2 sin f (cos f + e). The stated limit P = −4 sin²f D⁴ (e = 1, μ = 0) is correct. |
| m3 | minor | new certify_box (stall rule) | Sound, but it loses 30 of 267 certifiable boxes (11%, face-biased sample). It makes the 6-fold curve uncertifiable for e ≥ 0.9 at every box size down to h = 1e-8; the old loop certified at h ≤ 1e-5. Expect permanent failed leaves there unless `recheck_fails.py` is run. |
| m4 | minor | authors' H1·H2 formula | The exact identity is P_f = (1 − t)·4 sin²f·H1·H2, with a positive factor. Harmless. |
| m5 | minor | nb_run / nb_run_chart | `OUT.fails.jsonl` is opened in append mode, so reusing an output name mixes stale entries. The `nb_run_chart.py` docstring is copy-pasted ("CHART: a … or f"). |
| m6 | minor | structural failures | Three failure clusters are near-degeneracies, not counterexamples; they need deep re-checks or local charts. True maximum counts: 8 on curve6, 10 in the near-triple cluster, ≤ 8 in the bulk_a3 cluster. |

---

## 1. Statement fidelity: OK

- **Canonical text.** Checked against acs2012.txt lines 206–217. Q2: one conic is a circle, the other an ellipse with e > 0, and they share a focus. For a circle the focus is its centre.
  - The excluded cases (coinciding conics; concentric coplanar circles) cannot occur when e > 0.
  - So Theorem 1 as stated in STATEMENT.md is exactly Q2.
- **Parametrisation.**
  - **Orientation.** Rotating about the circle's axis removes Ω. (P, Q) is the z-x-z Euler frame with I ∈ [0, π] and ω ∈ [0, 2π), so every relative orientation is covered.
  - **Scale.** Scaling to a = 1 leaves the stationary points unchanged, so s = R/a ∈ (0, ∞) is the right parameter.
  - The project's domain e ∈ [0, 1) is a superset of Q2's e ∈ (0, 1). That is harmless.
- **Dependence on (e, s, P_z, Q_z) only: OK.** X_z = (cos u − e)P_z + β sin u Q_z and ρ² = r² − X_z², so Pol depends only on these four parameters. At a point where ρ = 0, the stationary points over that u are the solutions v of A/2 = s(X'_x cos v + X'_y sin v). Their number depends only on X'_x² + X'_y² = |X'|² − X'_z², which again depends only on (e, P_z, Q_z).
- **Symmetries: OK.**
  - **Proof.** ω → ω + π maps (P_z, Q_z) to −(P_z, Q_z), and Pol is invariant. ω → −ω gives Pol_new(u) = Pol_old(−u) (checked by hand). I → π − I leaves (P_z, Q_z) unchanged.
  - **Numerics.** `a10_spotchecks.py` S1: in 40 random configurations, the 2-D Newton count is the same under Ω ≠ 0, I → π − I, ω → −ω and ω → ω + π (40/40).
  - **Domain.** These give l = sin I ∈ [0, 1] and ω ∈ [0, π/2], with (P_z, Q_z) filling the quarter disc.

## 2. Reduction lemma: OK

**Proof.**
- **Points with ρ(u) > 0.** ∂_v d² = 2sρ sin(v − φ(u)), so for fixed u there are exactly 2 critical v. By the envelope argument, ∂_u d² there equals g∓'(u). Also Pol = ρ²·g₋'·g₊'.
  - One branch vanishing at u gives multiplicity ≥ 1.
  - Both vanishing at u requires A = (ρ²)' = 0, which gives multiplicity ≥ 2.
- **Points with ρ(u₀) = 0.** d² is constant in v. The condition A/2 = s(X'_x cos v + X'_y sin v) has at most 2 solutions, because X'_xy ≠ 0. Indeed X × X' = β(1 − e cos u)(P × Q) ≠ 0, so X' is never parallel to X, which lies on the axis. Pol(u₀) = Pol'(u₀) = 0, so the multiplicity is ≥ 2.
- **Points with Pol(u) ≠ 0.** No stationary point lies over such a u.
- **Conclusion.** N ≤ #real roots of Pol with multiplicity, whenever Pol ≢ 0.

**Pol ≢ 0 for e ∈ (0, 1): OK** (`a9_corner_symbolic.py`, item 5).
- **Short proof.** The extreme coefficient is c₀ = e⁴(p² − e² − b²q²)/16 + i·b e⁴ p q/8, independent of S.
  - If c₀ = 0 then q = 0 and p = e.
  - At q = 0, p = e the next coefficient is c₁ = e⁵(1 − e²)/4 ≠ 0.
  - So Pol ≢ 0 for e ∈ (0, 1). For e = 0, l > 0 we have Pol = −S((ρ²)')² ≢ 0.

**Numerical attempt to break it.** `a3_reduction.py` and `a3b_reduction_extra.py` use a damped 2-D Newton from 56×56 starts, built from the raw 3-D vectors.

| family | configs | result |
|---|---|---|
| R (random) | 150 | N_2D = N_Pol in 150/150 |
| L1 (l = 1 exactly: ρ = 0 at two points) | 120 | N_2D ≤ N_Pol always; strict in 37 cases, e.g. (4, 8), (6, 8), (8, 10). These are the expected over-counts at ρ = 0. |
| L1⁻ (l = 1 − δ) | 100 | equal in 100/100 |
| W0 (ω ∈ {0, π/2}: A = (ρ²)' = 0 at u = 0, π) | 120 | 0 violations; 3 "ambiguous" (6, 4) cases (below) |
| AX (P_z = 1) | 120 | N_2D = 8 in 120/120 |
| T (search-maximised, N_Pol = 10) | 40 | N_2D = 10 in 40/40 |

- **The ambiguous W0 cases are exact double roots** at u = 0, π that split under rounding. At 50 digits |Pol|, |Pol'| ≤ 1.1e-13, and this residual comes only from cos(π/2) = 6e-17 in the inputs.
- **No violation (N_2D > N_Pol) in any family.** In every case each 2-D critical point's u is a root of Pol.

## 3. Generated coefficients: OK (exact)

**`a1_symbolic.py`** (system python + sympy). This is an independent pipeline: real trig-polynomial algebra with product-to-sum formulas, with no z = e^{iu} substitution and no re/im. Every generated function was run on sympy *symbols* and compared monomial by monomial. The maximum coefficient difference is **0.000** for all 16 objects:
- `coeffs_ce` (Re, Im of Q(z), including S);
- `coeffs2_ce` (phi_ab, psi_ab);
- `coeffs3_ce` (V and all partials D in e, b, p, q);
- `coeffsF` (V and partials in e, p, q);
- `nbcoef_a` and `nbcoef_f`: values, partials, |value| and |partial| polynomials, with layout out[d+1, f, part, k].

**Negative control (`a1_negcontrol.py`).** A planted 1e-6 perturbation is reported as 1e-6, and swapping Φ ↔ Ψ is reported as 4.0.

**Also verified:**
- **Planar identity:** Pol = A²(r² − S) at l = 0.
- **Concentric coplanar circles:** Pol ≡ 0 at e = l = 0.

**Raw geometry (`a2_geometry_dft.py`).** Φ and Ψ were built from the 3-D vectors (I, ω), with derivatives from mpmath.diff at 40 digits and a 32-point DFT. The worst relative error against `coeffs3_ce` is 3.8e-16 (a-chart) and 4.6e-16 (f-chart, D⁸/(1+e)⁴ scaling).

**f-chart derivation: OK.**
- **Hand derivation.** r = (1+e)/D, (r²)' = 2(1+e)²e sin f/D³, (ρ²)' = (1+e)²[(σ²)'D + 2e sin f σ²]/D³. Multiplying by D⁸/(1+e)⁴ gives PhiF − κ²PsiF, and B = (σ²)'D + 2e sin f σ² simplifies via ζ' cos f + ζ sin f = Q_z (also asserted symbolically).
- **Root counts (a2).** The a-chart and f-chart agree (κ = s/(1 − e), anomaly map) in 399/400 cases.
  - The single mismatch (e = 0.9436, l = 0.006) is a tolerance artefact of *my* harness: a near-real complex pair sits within 1e-7 in one chart's metric.
  - A 60-digit sign scan (`a2b_mismatch.py`) gives the same structure in both charts.

## 4. Soundness of the certificate (`nb_cert.py`): OK, apart from float rigour (M3)

### Line-by-line review

| component | checks | status |
|---|---|---|
| `box_coeffs` (chart a) | chain rule ∂_η = c_e b − c_b e, ∂_l = c_p sin ω + c_q cos ω, ∂_ω = l(c_p cos ω − c_q sin ω); variation bounds such as vDe·b_max + \|De\|·db (derived by hand; the ω bound is slightly conservative); image box \|p − p_c\| ≤ rl + l_max·rw (≤ dp, likewise q, e, b); abs-polynomial inequality \|P(x) − P(x₀)\| ≤ Pabs(\|x₀\| + r) − Pabs(\|x₀\|) | OK |
| `box_coeffs` (chart f) | same pieces with e, P_z, Q_z as variables | OK |
| `arc_signs` | k^j cos(ku + jπ/2) derivative formula; k = 0 skipped for j > 0; Taylor remainder h\|f^(j+1)(u_m)\| + h²/2·sup\|f^(j+2)\| with sup ≤ Σ k^(j+2)(\|a_k\| + \|s_k\|) | OK |
| t-endpoint trick | f_t = (1 − t)Φ − tΨ is affine in t for fixed (u, x); a strict sign at both ends holds on the whole interval | OK |
| `dp_count` (proof) | runs are closed arcs, and arc st (no zeros) is excluded, so no run wraps; a zero at a shared endpoint is counted in both runs (over-count only); the DP is optimal because dp is monotone and maximal runs are used; with no zero-free arc it returns 999 | OK |
| `dp_count` (`a4_dpcount.py`) | the numba and Python versions equal an exhaustive brute force on 52,500 sign patterns (M ≤ 7, J ≤ 2); on 3,000 random degree-6 trig polynomials with exact sign tables the true root count never exceeds dp_count (bound attained 2,455 times) | OK |
| control flow | the "hopeless" exit, the v2/v3 stall rule, the 4096/arc_cap break and `certify_box_deep` are all `break` paths that return False; a True requires total ≤ 10 from `dp_count` | OK |

**Answer to the authors' item 1.** The stall change can never produce a false True. Code reading shows it only adds early exits. Empirically (`a14_regression.py`), on 400 boxes v2 never certified a box the old loop rejected (new_only = 0), and the end-to-end test on v2 has 0 violations (below).

### Empirical soundness

Ground truth is an independent evaluator (closed forms + FFT; `nbport.py`).

| test | scale | violations |
|---|---|---|
| **a5** port vs numba signs | 400 random boxes, 4 size scales, both charts | 0 mismatches (port is exact) |
| **a5** pointwise parameter enclosure (f^(j)(u_m; x) ∈ vc ± rad) | 22.1 M checks | **0**; max \|f − vc\|/rad = **0.996** (nearly tight, so the test is sensitive) |
| **a5** arc × box enclosure [lo, hi] | 29.0 M checks | **0** |
| **a5** numba sign claims at sampled (u, x, t) | 7.9 M checks | **0** |
| **a6** end-to-end, v1 | 460 boxes (376 certified), 60 points × 3 t per box | **0**; bound attained (slack 0) in 215 boxes; the 26 certified boxes centred on N = 10 configurations all have B = 10 |
| **a6** end-to-end, v2 | 340 boxes (276 certified) | **0** |
| **a6b** faces ω = π/2, ω = 0, l = 0, l = 1, t = 0, t = 1, v1 and v2 | 2 × 480 boxes, sampling on the faces | **0** |
| **a6d** python v5 engine (`certify5_ce`) | 160 boxes (137 certified) | **0** |

The refined final arcs were also checked: 0 sign violations (a6).

### Float rigour (M3)

The design is sound; these items are not a proof yet.
- **Arc gaps.** Arc endpoints are computed in floats. Consecutive arcs leave gaps of up to 8.9e-16 (measured in a6), and the last arc ends at fl(2π) < 2π.
- **Trig arguments.** cos and sin of (ku + jπ/2) use the float math.pi. The 1e-12·mag pad does not provably cover the argument error for terms whose own |c| is missing from mag.
- **MVT radii.** The radii `rr` are the unpadded half-widths about a rounded centre.
- **Abs-polynomial differences.** The rounding of a1 − a0 in the l/ω directions is padded relative to |Dp| + |Dq| + vD, not relative to a1.

All of these are ~1e-15 effects, many orders below the margins. A rigorous version needs:
- outward-rounded (or `mpmath.iv` / arb) evaluation;
- dyadic arc endpoints with an enclosure of π;
- radii padded for the centre rounding.

### Completeness (not soundness)

**Cost-map failures.** `a7_costmap_fails.py`: 67% of the 331 failed leaves in `costmap.json` lie on the curve **ω = π/2, e = l², t → 1**.
- Along it, (ρ²)' = 2l²(1 − l²) sin u (1 + cos u), so f₁ = −Ψ has a **6-fold root at u = π**.
- The other 33% form a cluster at e ≈ 0.436, l ≈ 0.99, ω ≈ 0.26, t ≈ 0.197 (`a7b_cluster2.py`): two roots 0.01 apart with a third nearby (a near-triple root). The true count is ≤ 10 there, and the off-circle pair is far from the circle (|log|z|| = 2.26).
- With v1, both certify only for tiny boxes (`a6c_curve6.out`): h ≤ 1e-4 for e ≤ 0.6, h ≤ 1e-5 for e = 0.9.
- With v2, curve6 at e ≥ 0.9 fails for every h down to 1e-8 (`a13_determinism.out`, `a6e_curve6_high_e.out`); this is m3.

**bulk_a3 so far.** `a12_analyze_bulk.py` on the running `fails.jsonl`: 69 failed leaves, all at e = 0.20000–0.20003, l ≈ 0.2, ω ≈ 0.19, t ≈ 0.509 (certified total 11). `a7c_cluster3.py` shows the true count in that hull is 4, 6 or 8, never 10 (two double-root transitions coincide).

## 5. Branch-and-bound coverage (`nb_run.py`): OK, except M1 and M2

- **Incomplete and failed boxes.** A box that hits the time cap is marked `done = False`. Its certified leaves still count towards the printed "volume fraction", but completeness must be judged by `done` and `fails` only.
- **Recording.** Failures are recorded with their totals. v2 also streams them to `OUT.fails.jsonl` (append mode; m5).
- **Tiling.** Initial boxes come from `linspace` edges and splits at shared float midpoints, so the domain is tiled exactly: leaves + fails + unprocessed stack = initial box.
- **M2: leaves are not recorded.**
  - Neither `nb_run.py`, `nb_run_chart.py` nor `nb_run_corner.py` records certified leaves. They store `leaves = [count, volume]` and a histogram only, contrary to their docstrings and to LOG ("records every certified leaf for later independent re-checking").
  - `left` is truncated to 50 entries (JSON) and 20 (`fails.jsonl`), and `recheck_fails.py` reads only that list. Incomplete initial boxes must therefore be re-run in full.
- **M1: the domain ends do not meet the neighbouring charts.**
  - The bulk runs (bulk_a2, bulk_a3) use η ∈ [0.20136, 1.25323], i.e. **e ∈ [0.2000020, 0.9499982]**. The e < 0.2 charts end at e = 0.2 exactly, and the e → 1 work starts at 0.95 or above.
  - Coverage sampling in `a15_charts.py`: 5,714 of 400,000 points are uncovered, all with e ∈ (0.2, 0.2000020] and l ∈ [0.408, 1].
  - **Fix:** use η_lo ≤ 0.20135 and η_hi ≥ 1.25324, or overlap the neighbouring charts. Also widen the float seams C1 θ_max = fl(arctan 2) (about 1e-16 short) against C1pp epp ≤ 0.5, and write an explicit coverage lemma.
- **Runs I could not analyse:**
  - `bulk_a.json` was killed.
  - `bulk_a2.log` ends with `PermissionError`.
  - `nb_f_wide.json` was never written and the run is not active.
  - `bulk_a3`, `run_C1`, `run_C1p` and `run_C1pp` are still running. At 23:49 bulk_a3 stood at 475/8000 boxes with 69 failed leaves, and run_C1pp at 200/256 with 0 failed. `a12_analyze_bulk.py FILE` is ready for the final JSONs.

## 6. Numerical claims in LOG / STATEMENT

| claim | verdict | evidence |
|---|---|---|
| max N = 10, never 12 | **OK** | a3 T-family (40 configurations, N_2D = N_Pol = 10). The a8 searches never found 12 in the interior; the gap → 0 only towards e → 1. |
| off-circle gap δ ≥ 0.836 | **FALSE (M4)** | a8: 0.554 (e ≤ 0.95), 0.784 (e ≤ 0.90), 0.275 (e ≤ 0.999). Minima are always at the top e-edge with small t (the C2 regime). Confirmed in a11 (60 digits): e = 0.94998, l = 0.999996, ω = 0.000364, t = 0.009032 has 10 roots on the circle, gap 0.5557, 2-D count 10. |
| parabola–circle ≤ 8 | **OK** | a10 S3: 30,000 random configurations, maximum 8. |
| parabola–circle gap ≈ 1.732 | **inaccurate (m1)** | At l = 1, ω = 0: PolP = 16w²(1 − w²)²((1 + w²)² − k²), gap √(1 + k) → √2 as k ↓ 1, with 8 real roots (a11). At nearby interior points the count drops to 2. |
| max #crit(g₋) = 6, #crit(g₊) = 4; joint counts never (6,6) or (8,·) | **OK** | a10 S4: 12,000 random + 1,500 × 60 targeted configurations. Maxima 6 and 4; joint histogram {(2,2), (2,4), (4,2), (4,4), (6,2), (6,4)}. |
| corner probes (C1, C1′ ≤ 8; C2 tight at 10; C3 ≤ 10 rarely) | **OK** | a10 S5 with an independent evaluator: C1 max 8, C1′ max 8, C2 max 10, C2b max 10, C3 max 10 (1 of 8,000). Odd counts near e → 1 are float artefacts, as LOG says. |
| planar Pol = A²(r² − S), N ≤ 6 | **OK** | a1 (symbolic) |
| e = 1 a-chart: Pol = 4c²r² sin²u (r² − Sc²), exactly 10 real roots (6 at u = 0) | **OK** | checked by hand: b = 0 gives Z = −p r and ρ² = c²r²; the 10 real roots need 0 < sc < 2 |

## 7. CORNER_E1.md: OK, except m2

`a9_corner_symbolic.py` (sympy, exact):
- **F∓ formula.** F∓ = 2μe(1+e) sin f σ ∓ D·B equals g∓'·D³σ/(μ(1+e)), a positive multiple where σ > 0 and D > 0. Also F₋F₊ = μ²(PhiF − κ²PsiF) with κ = 1/μ.
- **Aphelion values.** B(π) = −2P_zQ_z(1 − e), B'(π) = 2P_z² − 2e + 2(e − 1)Q_z² (= −2c² at e = 1), and D(π) = 1 − e.
- **AM–GM threshold.** min over φ of (2c²ε + 3c²φ² − 2|P_zQ_z|εφ) = ε(6c⁴ − εk²)/(3c²) with k = |P_zQ_z|. This is ≥ 0 iff ε ≤ 6c⁴/k², which is ≥ 24c⁴ because k ≤ 1/2.
- **C3 limit.** At e = 1, μ = 0, P_z = 1: F₋F₊ = −4 sin²f D⁴, a 10-fold zero at π plus a double zero at 0. The intermediate claim "B = 2 sin f · D" holds only at e = 1 (m2).
- **Lemma N2 logic (hand check): OK.**
  - **Step 3 (sign argument).** For φ > max(0, φ_B), μk sin f σ ≤ 0 and DB < 0 (D ≥ 1 − e > 0). For φ < min(0, φ_B), both terms are ≥ 0 and DB > 0.
  - **Step 4 (monotonicity).** On the interval between 0 and φ_B: |φ| ≤ |φ_B| ≤ 2|P_zQ_z|ε/β₀ = Cε, |D'| ≤ |φ|, |B| ≤ |B'|_max·|φ − φ_B| ≤ |B'|_max·Cε, and DB' ≤ −β₀ε. Hence F₊' ≤ C²|B'|_max ε² − β₀ε < 0 whenever C²|B'|_max ε < β₀.
- **Mixed certificate.** N ≤ zeros(P on |φ| ≥ φ₀) + zeros(F₋ on |φ| < φ₀) + zeros(F₊ on |φ| < φ₀) is valid. It needs σ > 0 on the inner arc (ensured by c ≥ c₀) and the reduction lemma's multiplicity accounting on the outer arcs. The "6" for the outer arcs is a branch-and-bound target, not yet proved.

## 8. The authors' later items

| item | verdict | evidence |
|---|---|---|
| **Axis-aligned case (P_z = 1, Q_z = 0), N = 8** | **OK** | (a) P_f = (1 − t)·4 sin²f·H1·H2 and H2(f) = −H1(−f), exact (a9). (b) Geometry: the orbit lies in the xz-plane with the apse line along z; the evolute cusp sits at distance a e² from the centre, while (±R, 0) sits at height a e above the centre. So (±R, 0) is outside the evolute, giving exactly 2 non-degenerate normals each, none at the apses (R > 0). The apses give v = ±π/2. Total 4 + 2 + 2 = 8. (c) a10 S2: H1 has exactly 2 zeros in 20,000/20,000 (e, κ), including e → 1. (d) a3 AX: N_2D = 8 in 120/120. |
| **1. stall rule / 4096 cap: no false True** | **OK** | Code reading (only `break`s returning False) and a14 (new_only = 0). Completeness regression: m3. |
| **2. `OUT.fails.jsonl`** | **OK, with m5 and M2** | Failures and incomplete boxes are logged and flushed. Append mode can mix runs, and `left` is truncated to 20. |
| **3. C1 / C1p / C1pp algebra** | **OK** | Hand check: C1 gives f/ε² = (1 − t)Φ/ε² − tΨ/ε², with A/ε = 2 cos θ sin u r and Ψ/ε² = (A/ε − sin θ (zt²)')². For C1p and C1pp, with t = ε²τ/(1 − τ) (resp. epp²τ/(1 − τ)): (1 − τ)f_t/e² = (1 − τ)F_A − τF_B, F_B = ε²Φ/e² + Ψ (resp. epp²Φ/e² + Ψ/l⁴). Numerically (a15) the identity holds to 1.1e-15. Every DER entry was checked by hand (b, ib through e = ε cos θ or epp·l²). τ ∈ [0, 1] covers t ∈ [0, ∞), which includes non-physical t > 1; that is harmless. |
| **3. generated nbcoef_C1\*** | **OK** | a15: values match the a-chart rescaled to 2.0e-15, 1.9e-15 and 1.5e-15. Box-variable gradients match complex-step derivatives to 3.5e-14 relative or better. Abs-polynomials are exact (difference 0). |
| **3. nbeng engines** | **OK** | The atom enclosures are monotone and correct: cos and sin are monotone on [0, π/2], and b, 1/b are monotone in e with the extreme e taken at the right corners. The MVT through atoms is valid; it needs no sign assumption, only \|a\| ≤ \|a_c\| + Da. a16: C1, C1p and C1pp each certified 150/150 boxes, including boxes near the blow-up faces, with **0 violations**. |
| **3. coverage of {e < 0.2}** | **OK, except M1** | The argument holds: if l ≥ 0.4, C1p applies; otherwise √(e² + l⁴) < √0.0656 = 0.2561 < 0.26, and θ ≤ arccot 0.5 gives C1 while e/l² ≤ 0.5 gives C1pp, with l < 0.4 ≤ 0.52. Two seam problems: the hole up to 0.2000020 against the bulk, and float-level θ/epp seams. |
| **4. corner engine `nb_corner.py`** | **OK in logic; float rigour as M3** | `near_derivs`: ζ derivatives, (σ²)⁽ʲ⁾ for j ≤ 4, the σ⁽ʲ⁾ recursion (σ' = (σ²)'/2σ, σ'' = ((σ²)'' − 2σ'²)/2σ, σ''' = ((σ²)''' − 6σ'σ'')/2σ, σ'''' = ((σ²)'''' − 8σ'σ''' − 6σ''²)/2σ), Leibniz sums for (sin f σ)⁽ʲ⁾ and (DB)⁽ʲ⁾, and B⁽ʲ⁾ = (σ²)⁽ʲ⁺¹⁾ + 2e(sin⁽ʲ⁾ − Q_z ζ⁽ʲ⁾) are all correct. `trig_arc` is valid on ⊂ [π/2, 3π/2], and φ₀ = 0.9 < π/2. μ enters affinely, so the endpoint check is valid. `check_plus` implements N2 exactly as in §7 and is μ-independent (needs only μ ≥ 0); \|φ_B\| ≤ 2\|P_zQ_z\|ε/β₀ uses the exact B(π); the condition "2C²\|B'\|max ε₁ < β₀" is sufficient, with a factor 2 to spare. Coverage: F₋ on [π − φ₀, π + φ₀]; F₊ on [π − φ_in, π + φ_in] (N2) together with the two outer pieces; P on [π + φ₀, 3π − φ₀]. Together these cover the circle when φ_in ≤ φ₀ (0.12 ≤ 0.9); only endpoints overlap, and float gaps are as in M3. The P_z/Q_z bounds use unpadded math.sin/cos (M3). a16 part 2: 60 corner boxes (57 certified); every component bound held at sampled points (dense sign counts): F₋ near, F₊ inner (≤ 1), F₊ outer, P away and the total, all **0 violations**. |
| **5. H1 has 2 zeros** | **OK** | Proved by the evolute argument, and confirmed independently (20,000 samples). |

## 9. Not checked / unclear

- **Final results of the running jobs.** bulk_a3, run_C1, run_C1p, run_C1pp and the full corner run were still running. nb_f_wide.json does not exist. Run `a12_analyze_bulk.py` on each final JSON.
- **Coverage of e ∈ [0.95, 1) as a whole.**
  - The f-chart covers κ ≤ 3, i.e. μ ≥ 1/3.
  - The corner smoke run covers e ≥ 0.98, μ ≤ 0.05, l ≤ 0.9.
  - Not covered by any run or chart yet: μ ∈ (0.05, 1/3), e ∈ [0.95, 0.98) at small μ, and the C3 sub-corner (l > 0.9, e → 1, μ → 0).
- **`c3_blowup.py`:** new, and not in this audit's scope.
- **Not attempted:** a fully rigorous (directed-rounding) re-implementation, and a check of the bulk's certified leaves (they are not recorded; M2).

## Appendix: scripts (in `audit/`) and their outputs

- **Coefficients and reduction:** `a1_symbolic.py` and `a1_negcontrol.py` (exact coefficients), `a2_geometry_dft.py` and `a2b_mismatch.py` (geometry, chart correspondence), `a3_reduction.py` and `a3b_reduction_extra.py` (reduction lemma).
- **Certificate soundness:** `a4_dpcount.py` (dp_count), `a5_enclosures.py` with `nbport.py` (enclosures), `a6_certify_bruteforce.py`, `a6b_faces.py` and `a6d_certify5.py` (end-to-end), `a6c_curve6.py`, `a6c_helpers.py` and `a6e_curve6_high_e.py` (curve6), `a13_determinism.py`, `a14_regression.py` with `oldcert_port.py` (stall rule).
- **Failures and gap claims:** `a7_costmap_fails.py`, `a7b_cluster2.py` and `a7c_cluster3.py` (failure clusters), `a8_gap_search.py` (δ claim), `a11_confirm.py` (60-digit confirmations), `a12_analyze_bulk.py` (run analysis, ready for final JSONs).
- **Corner and LOG checks:** `a9_corner_symbolic.py` (CORNER_E1, H1·H2, Pol ≢ 0, evolute), `a10_spotchecks.py` (symmetries, H1, parabola, branches, corner probes).
- **New engines:** `a15_charts.py` (C1/C1p/C1pp algebra, generated code, coverage), `a16_engines.py` (new engines end-to-end).
- **Outputs and hashes:** every script's output is in the matching `*.out` file. `hashes_audited.txt` holds the v1 hashes and `a6_new_version_hash.txt` the v2 hash.
