# AUDIT 2: fixes for AUDIT_1 and the near-coincident blow-up charts (D1a/D1b)

Independent audit, 2026-10-23, 00:52–01:06.

**Scripts and outputs** are in `audit/round2/`. Hashes of the audited code are in `round2/hashes_round2.txt`.

**Rules kept.** At most 2 cores were used, and nothing was written outside `audit/`: numba caches went to `audit/numba_cache` and Python ran with `-B`.

**The code changed again during the audit.**
- **00:58:45:** nb_cert, make_engines and make_engines_bu, and every `nbeng_*`, were regenerated. The change pads the mean-value radii (`rr`) for the float box centre. This is AUDIT_1 M3(iii); it is padding only, so it can only turn True into False.
- **00:59:** all runs were restarted into `code/runs/`, and the earlier ones moved to `code/runs_old/`.
- **Which files I tested:**
  - the old L2P, L2A±, L2T± JSONs (runs_old), checked with the *current* engines;
  - the new B2, B3p and B3m JSONs (runs/).

---

## Verdict

**No critical issue was found.** Every chart identity, generated polynomial, cone-coverage claim and certificate I could test is correct.

| id | severity | item | finding |
|---|---|---|---|
| R2-1 | **major** (checking tool) | `bb_replay.py --recert` | The engine dispatch handles `a`, `f`, `C1`, `C1p` and `C1pp`; **everything else falls through to `nb_corner.certify_corner_box`**, including B1, B2, B3±, L2E, L2P, L2A± and L2T±. On `runs_old/L2Tp.json --recert 5` it reports **5/5 false "re-certification problems"** (`b7_replay_recert_bug.out`). A wrong engine could also pass a leaf spuriously. **Fix:** send every engine other than a/f/K to `nbeng_<engine>.certify_box_deep`. I did the requested re-certification myself with the right engines (b2, below). |
| R2-2 | minor (completeness) | m3 deep fallback in `bb_run.py` | The fallback runs only when the regular certifier fails with **tot < 999**. On the 6-fold curve (ω = π/2, e = l², t → 1) the stall rule stops at tot = 999, so the fallback never runs. Example: e = 0.9 and 0.93 at h = 1e-5 and 2e-6 give regular (False, 999), the fallback is not called, and `certify_box_deep` alone gives (True, 10) (`b6_deep_fallback.out`). Expect failed leaves along that curve for e ≳ 0.9 in `runs/bulk` (e ≤ 0.9501). **Fix:** also call deep when tot == 999 for small boxes, or stop the stall counter while unresolved arcs keep decreasing. |
| R2-3 | minor (float rigour; part of M3) | generated constants | B2 and all L2 charts emit **double-rounded non-dyadic constants**: products with 0.6825, 0.01, 1.05 and their powers. Against a 50-digit evaluation of the physics, the polynomials differ by 1.2e-18 (B2), 3.4e-16 (L2E) and ≤ 3.6e-19 (other L2), even with the binary-float range constants. B1 and B3±, whose constants are all dyadic, match to ~1e-49. So the certified polynomial is the exact chart function up to ~1e-16 relative coefficient rounding. `crad` covers this at the box centre, but not provably over the box when Pabs(\|a_c\| + Da) ≫ Pabs(\|a_c\|). The written float-error analysis must include it. Alternatively, use dyadic scales (e.g. 1/128, 11/16, 17/16) or exact rational coefficients with outward rounding. |
| R2-4 | minor (documentation) | STATEMENT.md, DEGENERACY_TREE.md | Still uncorrected: δ ≥ 0.836 (M4) and an outdated proof plan in STATEMENT.md (Rouché; "s = t/(1−t)" should read S = s² = t/(1−t)). In DEGENERACY_TREE.md, "gap √3" should be √2 (m1). "Hessian determinant 4s² sin²I > 0" at e = 0 should read ±4s² sin²I: 4 of the 8 points are saddles. All 8 are still nondegenerate. |

---

## Item 1. G1/G2 identities, divisibility and the b-reduction: OK

`b4_symbolic.py` uses exact sympy and independent trig expressions; `b3_bu_coeffs.py` compares numerically at 50 digits.

- **B1** (p = e·ph, q = e·qh): Φ − Ψ ≡ e³G1 and Φ + Ψ ≡ e²G2 exactly, so F = f/e³ = G1/2 − τ̂G2.
  - By hand, using r² − 1 = e(−2cos u + e cos²u): G1 = At²(−2cos u + e cos²u − e zh²) + 2At(zh²)' − e((zh²)')², and G2 = At²(r² − e²zh²) + (At − e(zh²)')².
- **B2 and B3±:** Φ − Ψ ≡ e²θ·G1 (B2) and ≡ e²η·G1 (B3), Φ + Ψ ≡ e²G2, with l² = eθ. All exact.
- **Leading term:** F₀ = (G1/2 − τ̂G2) at e = 0, b = 1 equals 2 sin u [(qh² − ph² − 1) sin 2u + 2ph·qh cos 2u − 4τ̂ sin u], exactly.
  - It vanishes identically **only at (ph, qh, τ̂) = (0, 1, 0)** when ph, qh ≥ 0. That is the D1b direction l = e, ω = 0, s = 1.
- **CT expansion:** f/e² = 2 sin u [θ sin 2(u+ω) − 4τ sin u − e sin 2u] + O(2) holds exactly at first order.
- **b-reduction (level 2):** reduce G1 and G2 mod b² = 1 − e² to R₀ + b·R₁.
  - **R₁ is divisible by ph for both G1 and G2**, as claimed. b enters only through the cross term 2ph·qh·b(cos u − e) sin u of zh²; its square reduces to b² = 1 − e².
  - **R₀(G1) = 0 at the D1b point.** Since τ̂ ∝ ν, the whole of F_B1(ν·x) is then divisible by ν as a polynomial, with b kept as an independent symbol. That is exactly what `gen_bu.py` asserts.
  - The reduction changes no values on b = √(1 − e²).
- **Generated code (b3):** for all 10 charts I compared coef_val with my own chart maps and the a-chart physics at 50 digits.
  - **Values:** B1, B3± agree to ≤ 1e-49; B2 and the L2 charts agree up to R2-3.
  - **Gradients** (out[1+v], chain rule through DER) agree with mpmath central differences of the value polynomials to ≤ 1.4e-18.
  - **Absolute-value polynomials:** 16,876 terms checked one by one, **0 mismatches**.
- **The authors' L2 self-test error (~1e-10)** is cancellation in the double-precision reference f/(e³ν), as suspected. At 50 digits the L2 chart polynomials agree with the physics to ≤ 3.4e-16 (L2E) and ≤ 3.6e-19 for the others.

## Item 2. Cone coverage of the C1- and B1-excluded blocks, including float edges: OK

`b5_coverage.py` uses the actual run domains, grids and `--exclude` arguments from the running `bb_run` command lines.

**C1 skip:** 64 initial boxes.
- They are exactly ε ≤ 0.065, θ_C1 ≤ 0.1845247863 (a grid edge), t ∈ [0.375, 0.625].
- So e ≤ 0.065, θ = tan θ_C1 ≤ 0.186648 and |τ| ≤ 0.125. **This lies inside the cube {e ≤ 0.07, θ ≤ 0.2, |τ| ≤ 0.13}.**

**Cones (R = 1.05):** each chart bound is ≥ R times the needed ratio, with a 1.05 margin:
- e-dominant: θ ≤ 3.0e, so ph, qh ≤ √3 = 1.732 < 1.75, and |τ̂| ≤ 1.95 < 2;
- θ-dominant: eh ≤ 0.3675 and |τ̂| ≤ 0.6825;
- τ-dominant: eh ≤ 0.5654 and thh ≤ 1.6154.

**B1 skip:** 4 initial boxes, exactly {e ≤ 0.01, ph ≤ 0.25, qh ∈ [0.75, 1.25], x3 ∈ [0.375, 0.625]}, i.e. |τ̂| ≤ 0.5. In scaled variables that is max\|X_i\| ≤ 1 < 1.01, so it lies inside the union of the L2 cones.

**Sampling** (800,000 points; e > 0, l ∈ [0, 1], ω ∈ [0, π/2], t ∈ [0, 1]): **0 uncovered**, apart from the declared C3 core of K (l ≥ 0.95, ω ≥ 9π/20, still pending).

| family | samples | uncovered |
|---|---|---|
| uniform | 100k | 0 |
| small e (down to 1e-7) | 200k | 0 |
| C1-skip block | 100k | 0: all in B1, B2 or B3 |
| B1-skip block | 100k | 0: all in L2 |
| e → 1 | 200k | 0, plus 918 in the declared C3 core |
| seams 0.199/0.2/0.2001, 0.9499/0.95/0.9501, 0.99 | 100k | 0 |

**Seams (AUDIT_1 M1): fixed.**
- **Low end:** the bulk now runs from asin 0.199, overlapping C1p (e ≤ 0.2001).
- **C1/C1pp seam:** C1pp extends to e/l² ≤ 0.51, overlapping C1 (θ_C1 ≤ arccot 0.5).
- **High end:** F1 covers e ∈ [0.95, 1] with t_f ≤ 0.9975 (κ ≤ 19.97). F2 covers e ∈ [0.9499, 0.99] with t_f ≥ 0.9975, overlapping the bulk (e ≤ 0.9501). K covers e ∈ [0.99, 1] with μ ≤ 0.0505 (κ ≥ 19.80), overlapping F1. F2 and K share the closed float edge e = 0.99.

**Float edges.**
- **Excluded blocks:** these consist of whole grid boxes, and `bb_replay` checks grid = processed ∪ skipped. So the processed region is exactly the complement, with no float slivers.
- **Chart-to-chart seams:** all have ≥ 0.05% (seams) to 5% (cones) slack.

## Item 3. The MV4 `signs` function against `nb_cert.arc_signs`: OK

- **Code comparison.** The same formulas are used for:
  - the derivative trig sums (k = 0 skipped for j > 0);
  - the pads: `rd` from `crad`, prad = Σ_v rr_v(\|g_v\| + w_v) over 4 directions, and tot += 1e-12(mag + tot);
  - the Taylor remainder h·m1 + h²/2·B_{j+2}, with B_{j+2} = Σ k^{j+2}(cm₀ + cm₁)(1 + 1e-12);
  - the hopeless flag and the sign tests.
- **The only differences:** one function instead of two t-ends; four mean-value directions instead of three; and `box_coeffs` reads `val[1 + v, 0, …]` for v ≤ 3. Arcs are padded by ARC_PAD = 1e-12 in `_run`.
- **Empirical test (`b8_mv4_signs.py`).** 40 random boxes each in B3p, L2Tp (MV4) and B2, L2P (affine). Every claimed sign of f^(j) (63,022 signed arc/derivative pairs) was checked at sampled points of box × padded arc: **5.67 M checks, 0 violations**.

## Item 4. `bb_run.py` and `bb_replay.py`: OK, except R2-1 and R2-2

**`bb_run.py`.**
- **Tree encoding.** It writes one code per popped box in pre-order: '+' for a split, '!' for a failure, and a leaf letter that carries the bound and the method. The split rule is deterministic.
- **Leaf codes are valid:** certified means tot ≤ 10, so the letter falls in 'a'..'k' or 'A'..'K'.
- **Time-out:** it stores the complete stack, so `left` is no longer truncated.
- **Failures file:** `fails.jsonl` is opened in 'w' mode, which fixes m5.
- **Exclusions:** `--exclude` skips only initial boxes lying *entirely* inside the exclusion box, and records them.

**`bb_replay.py`.** The replay is correct.
- **Stack mechanics:** it pops in the same order and pushes children (lower, upper).
- **Accounting checks:** fails are compared in order, `left` is checked for incomplete records, the leaf count must match, and volume must be conserved to 1e-9.
- **Grid check:** grid = processed ∪ skipped, with the two sets disjoint.

**Replay results** (`b1_replay.out`): 0 replay problems on all 11 JSONs.
- runs_old: B3p, B3m, B2, L2P, L2Ap, L2Am, L2Tp, L2Tm (all 256/256 complete, 0 fails).
- runs: B2, B3p, B3m (new engines; identical bound histograms).

**Re-certification with the correct engines** (`b2_recert.py`; replaces the buggy `--recert`):

| run | leaves checked | regular certifier reproduces (True, recorded bound) exactly | deep certifier worse | brute-force violations |
|---|---|---|---|---|
| L2Tp | 256 (all) | 256/256 | 0 | 0 |
| L2Tm | 256 (all) | 256/256 | 0 | 0 |
| B3p, B3m, B2 (runs/) | 300 each, random | all | 0 | 0 |
| L2P, L2Ap, L2Am (runs_old) | 300 each, random | all | 0 | 0 |

- **Brute force:** 16 points per leaf, with the real roots of the chart function counted (biased upwards).
- **Engine change:** the runs_old leaves were certified with the 00:27 engines, before the radius padding. They are reproduced exactly by the current engines.
- **Deep-certifier claim confirmed:** in bb_run's docstring, the claim that every regular leaf is reproduced by `certify_box_deep` (same refinement sequence, larger cap) held in every sampled case.

## Item 5. Fidelity to the canonical statement: OK

- **Theorem 1 is still exactly Q2.** The charts only reparametrise (e, l, ω, t) and divide f by positive normalisations: e³ (B1), e²θ (B2), e²η (B3±), e³ν (L2*). On the domain interior (e > 0, θ > 0, η > 0, ν > 0) this preserves the real roots and their multiplicities.
- **Affine variables.** In the affine charts, τ̂ ↔ x3 and T ↔ x3 are affine bijections, and F is affine in that variable, so (1 − x3)F_A − x3F_B = F exactly.
- **Faces where the normalisation vanishes** (e = 0, ν = 0) are limit functions certified as part of the boxes. They are not physical Q2 points, and certifying them is harmless.
- **No mismatch** between the new material and STATEMENT.md's target. The stale parts of STATEMENT.md are listed under R2-4.

## Item 6. e = 0 and Lemma E0

- **Not needed for Q2.** "If only one of the eccentricities is 0" means the ellipse has e > 0, so e = 0 (two concentric inclined circles) lies outside Q2.
- **Every e > 0 is reached for all t.** C1 takes t directly. C1p and C1pp use t = ε²τ/(1 − τ) (resp. ε''²τ/(1 − τ)), so for every e > 0 the values τ ∈ [0, 1/(1 + ε²)] give all t ∈ [0, 1]. The sampling found 0 uncovered points down to e = 1e-7.
- **Lemma E0 is correct** (b4, exact): at e = 0, Pol = −S l⁴ sin² 2(u + ω), so it has 8 real roots with multiplicity. In the circle–circle normal form x̂·ŷ = cos u cos v + cos I sin u sin v, the 8 critical points have Hessian determinants ±(cos² I − 1) ≠ 0 for 0 < l = sin I ≤ 1, **including l = 1**. So N = 8, all nondegenerate.
- **Recommendation:** include E0 only if the theorem is also stated for two circles; it is optional for Q2.

## Status of the AUDIT_1 items

| AUDIT_1 | status |
|---|---|
| M1 seams | **fixed** (item 2) |
| M2 leaf recording | **fixed** (bb_run trees; replay verified). Old nb_run* superseded. |
| M3 float rigour | **partly fixed:** ARC_PAD on all arcs (overlap covers [−pad, 2π + pad]), padded MV radii (00:58), nb_corner PAD on every sub-range, 2-ulp-padded P_z/Q_z. **Still open:** trig-argument rounding, a1 − a0 rounding, coefficient-constant rounding (R2-3). These are for the promised float-error analysis. |
| M4, m1, m2, m4 wording | pending (R2-4) |
| m3 stall rule | deep fallback added; gap R2-2 |
| m5 fails file | **fixed** ('w' mode) |

## Not checked / pending

- **Running jobs (01:05, 0 failures logged so far):** bulk, B1, C1, C1p, C1pp, F1, F2, K, L2E. Re-runs of L2P, L2A± and L2T± on the 00:58 engines are not yet in runs/.
- **C3 core chart** (`code/c3/`, being built) and the K exclusion it must cover.
- **K boxes near l → 1:** the corner engine needs σ² > 0 on the near arc |f − π| ≤ 0.9. For l close to 1 and ω ≳ 0.67, σ = 0 can occur on that arc, and K returns failure there. Check K's fails once it finishes.
- **The written float-error analysis.**

## Scripts (`audit/round2/`)

| script | checks | output |
|---|---|---|
| `b3_bu_coeffs.py`, `b3b_bu_coeffs_floatconst.py` | values against 50-digit physics, gradients, abs-polynomials | `*.out` |
| `b4_symbolic.py` | G1/G2 identities, b-reduction, F₀ and D1b, CT expansion, Lemma E0 | `b4_symbolic.out` |
| `b5_coverage.py` | coverage with the actual run parameters | `b5_coverage.out` |
| `b1_replay.out` | output of `bb_replay.py` (no recert) on 11 JSONs | — |
| `b2_recert.py` | correct-engine re-certification plus brute force | `b2_recert.out` |
| `b6_deep_fallback.py` | R2-2 | `b6_deep_fallback.out` |
| `b7_replay_recert_bug.out` | R2-1 demonstration | — |
| `b8_mv4_signs.py` | sign-enclosure soundness for MV4 and affine charts | `b8_mv4_signs.out` |
