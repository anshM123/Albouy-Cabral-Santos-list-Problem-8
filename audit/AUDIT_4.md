# AUDIT 4: the winding-number certificate (`code/wind_cert.py`)

Independent adversarial audit, 2026-10-23, 02:15–02:40.

**What was audited:** `wind_cert.py`, its integration in `bb_run.py` / `bb_replay.py`, PROOF.md §2 (winding lemma) and FLOAT_RIGOUR.md §3 (winding paragraph).

**Scripts and outputs:** `audit/round4/`. Hashes of the audited code are in `round4/hashes_round4.txt`.

**Rules kept.** At most 2 cores were used, and nothing was written outside `audit/`.

## Verdict

**The lemma, the G formula, the enclosures and the driver integration are correct.** All 46,070 existing winding leaves were re-checked:
- B1: 5,171; C1pp: 6,138; L2E: 290; C1: 34,471.

**One genuine logic gap.** The floating-point winding evaluation can, in principle, pick the wrong branch for arcs labelled Re<0 and overcount W by 1 (item 4). It is **not triggered** by any existing leaf: every step is ≤ π − 2.34, with node values ≥ 6.6% inside their half-planes. The fix is a few lines.

| id | severity | finding |
|---|---|---|
| R4-1 | **major** (soundness gap in principle; not triggered) | The winding step uses principal-value normalisation with strict `while d > π` / `while d < −π`. For an arc labelled **Re<0**, the atan2 cut lies inside the half-plane. If both node values sit within rounding of the imaginary axis on opposite sides, the raw difference is exactly ±fl(π) and is **not** renormalised. The step is then off by 2π, and the total still passes the "within 0.5 of 2π·integer" test. In one direction this **overcounts W by 1**, giving a leaf bound 2 too small. The sanity checks (node values in the arc's half-plane) do not exclude it. Demonstration on crafted node values (`d4_branch_demo.py`): the code returns winding 0 where the true winding is −1. **Fix:** measure both endpoints in the label's own continuous branch (Re>0: (−π/2, π/2); Im>0: (0, π); Re<0: (π/2, 3π/2); Im<0: (−π, 0)) and drop the normalisation; or reject the box if any \|step\| > π − 1e-6. |
| R4-2 | minor (rigour hygiene) | **Wrap-around node:** the first node value is computed at θ = left[0] and the last at left[0] + 2π, separately. Telescoping is exact only if the same value is reused. Observed discrepancy ≤ 4.3e-15 relative; reuse the first value, checked against both adjacent labels.<br>**L guard:** there is no explicit check that L < 0 in `wind_box`, `bb_replay` or `check_all`. Validity requires L ≤ 0; only `choose_L` enforces L < −1e-3, so a stored or replayed L is trusted. Add `assert L < 0`, and 1 ≤ W ≤ 6, in the certifier and the checkers. |
| R4-3 | minor (float write-up) | FLOAT_RIGOUR §3 says the weights are accurate to "2 ulps", and charges their error to "the 10⁻¹²·mag term". Both statements need correcting:<br>• cosh(fl(k·L)) amplifies the rounding of k·L, so the relative weight error is ≤ (\|kL\| + 3)u.<br>• 10⁻¹²·(mag + tot) does not bound Σ\|w_k\|·\|A_k(x)\|·\|trig\| over the box, because tot contains \|g_v\| with cancellation.<br>The error **is** covered by the spare in the 2e-14·aab term (aab bounds that absolute sum). The trigonometric error with phase and \|arg\| ≤ 6·2π + π ≈ 40.9 is ≤ 9.5e-15, leaving a spare ≥ 1.05e-14 ≥ (\|kL\| + 3)u while \|kL\| ≤ 90. Observed max \|kL\| = 55.6; a guard \|L\| ≤ 15 makes it unconditional. The common (1 + 10⁻¹²) factor on all weights is harmless (a positive rescaling of G). |

## (1) The lemma: OK

- **Self-inversive.** z¹²·conj(Q(1/z̄)) = Σ conj(c_k) z^{6−k} = Σ c_j z^{6+j} = Q(z), using c_{−k} = conj(c_k).
- **Roots at 0 and degree drop.** If Q = z^{m₀}R with R(0) ≠ 0, the coefficient symmetry forces deg R = 12 − 2m₀, and R is self-inversive. Its roots pair z ↔ 1/z̄. So:
  - the roots on the unit circle number (12 − 2m₀) − 2·#R_in;
  - the roots of Q in |z| < 1 number m₀ + #R_in;
  - hence N_circle = 12 − 2·#Q_in, exactly. Roots at 0 count as inside, matched by the degree drop.
- **Real zeros ↔ circle roots.** Q(e^{iu}) = e^{6iu}F(u), and d/du = iz·d/dz. So multiplicities agree.
- **Counting inside a smaller circle.** For ρ < 1 with no root on |z| = ρ, #Q_in ≥ #{|z| < ρ} = W, so N ≤ 12 − 2W.
- **Q ≢ 0** is implied by certified labels (G ≠ 0 on the circle).
- **L ≤ 0 is essential:** for L > 0 the bound is invalid. See R4-2.
- **Synthetic check** (`d3_wind_stress.py`): 4,000 self-inversive polynomials with known roots, covering:
  - double roots on |z| = 1;
  - m₀ = 1..3 zeros at 0 with matching degree drop;
  - a root 1e-12 to 1e-3 from the test circle;
  - W up to 6.
- **Results:** 2,089 were certified, **all with the exact W**. Near-circle roots were certified only once in 800 (safe). Q ≡ 0 gives (False, −1).

## (2) The G formula, the phase −π/2 and W = 6 + wind(G): OK

- **Derivation.** Q(ρe^{iθ}) = ρ⁶e^{6iθ}·Σ c_k ρ^k e^{ikθ}. Expanding c_{±k} = (a_k ∓ i b_k)/2 gives exactly:
  - Re G = Σ cosh(kL)(a_k cos kθ + b_k sin kθ);
  - Im G = Σ sinh(kL)(a_k sin kθ − b_k cos kθ).
- **Phase.** a_k sin x − b_k cos x = a_k cos(x − π/2) + b_k sin(x − π/2), so the phase is −π/2. For k = 0, sinh(0) = 0 exactly, so the cos(−π/2) ≈ 6e-17 artefact vanishes.
- **Winding.** wind(Q) = 6 + wind(G), by the argument principle. This also holds with roots at 0.
- **Confirmation at 100 digits** (`d2_wind_spot.py`): 360 random points in 120 winding leaves (B1, C1pp, L2E, C1).
  - The count of roots of Q in |z| < e^L **equals the recorded W in every case**.
  - No root lies within 1e-30 of the test circle; the minimum relative separation is 0.25.
  - The unit-circle count is ≤ 12 − 2W.

## (3) `comp_bounds` / `labels`: OK (write-up correction R4-3)

- **Same scheme as `arc_signs`.** Value and first derivative at the arc centre (derivative phase + jπ/2, factor w_k·k^j), Taylor remainder with B₂ = Σ\|w_k\|k²(cmax₀ + cmax₁), the mean-value radius, and the pads 1e-12(mag + tot) + 2e-14·aab and 1e-15 on the final sum.
- **t-ends.** Labels are certified at both t-ends with the same t-combination as `arc_signs`. G is affine in t (the weights do not depend on t), so they hold on all of [t0, t1].
- **Empirical check** (`d6_comp_bounds.py`; charts a, B1, C1pp; L ∈ [−6, −0.05]; independent physics):
  - 967,680 enclosure checks of Re G and Im G at arc ends and inside arcs, at box corners and random points: **0 violations**;
  - 37,618 checks of certified labels at random (x, t, θ): **0 violations**.

## (4) Winding evaluation: **not provably exact as written** (R4-1)

**What holds.** At interior nodes the same computed value is reused for both adjacent arcs. Whenever every computed node value lies in both adjacent half-planes (the sanity check), the computed and true argument errors telescope exactly. So the result is exact up to atan2 rounding (observed \|tot − 2πw\| ≤ 3.2e-14), with two exceptions:

1. **the branch flip of R4-1.** It can occur only for label Re<0, the only label whose atan2 range wraps. For the other labels the raw difference lies in [−fl(π), fl(π)], so it is never renormalised.
2. **the wrap-around node of R4-2.**

**Evaluation at the box centre is legitimate.** The winding is evaluated at the float centre coefficients `cmid`. Those lie within `crad` of c(x_c), and the label enclosures cover that ball, so by homotopy the winding at `cmid` equals the winding for every box point.

**Existing leaves are safe** (`d5_step_margins.py`). Re-running the exact labelling and refinement for **every** winding leaf:

| run | leaves | W reproduced | max \|step\| | min node \|label comp\|/\|G\| | max \|kL\| | L < 0 |
|---|---|---|---|---|---|---|
| L2E | 290 | 290 | π − 2.61 | 0.123 | 34.7 | all |
| B1 | 5,171 | 5,171 | π − 2.52 | 0.087 | 25.1 | all |
| C1pp | 6,138 | 6,138 | π − 2.52 | 0.285 | 55.6 | all |
| C1 | 34,471 | 34,471 | π − 2.34 | 0.066 | 34.4 | all |

The branch flip needs a step within about 1e-15 of π and node values within rounding of the label boundary. The nearest observed margins are 2.34 rad and 6.6% of |G|.

## (5) Driver and replay: OK

- **Driver order** in `bb_run.work`: regular Rolle, then winding (only when regular fails), then deep Rolle.
- **Winding leaves:** coded '1'..'6' (`make_wind` requires ok and W ≥ 1), with L appended to `wL` in tree order. The histogram uses 12 − 2W.
- **Engines.** It is used only for 3-variable affine engines, which matches `comp_bounds`' 3 mean-value directions. MV4, K, KB and C3 are excluded.
- **`bb_replay.py --recert 150`** (`d1_replay_recert.out`), 0 replay problems and 0 re-certification problems (winding leaves included) on:
  - B1 (8,379 leaves, 5,171 winding);
  - C1pp (11,356 leaves, 6,138 winding);
  - L2E (436 leaves, 290 winding);
  - C1 (77,436 leaves, 34,471 winding, 8 deep, 0 fails; finished).
- **Wind parameters.** The header stores the wind parameters (64, 3, 1024), and the replay uses them.

## (6) Ways W ≥ 1 could be wrong: covered, except R4-1 and R4-2

| case | status |
|---|---|
| Q ≡ 0 | labels impossible, so no certificate (tested) |
| Leading coefficient vanishing, roots at 0 | counted correctly by the argument principle (lemma); 598 synthetic cases correct |
| Roots on the test circle | excluded by the labels |
| Roots near the test circle | either certified correctly or rejected (synthetic family 2) |
| L ≥ 0 | would be invalid; currently prevented only by `choose_L` (R4-2) |
| Wrong branch | R4-1 |

## Notes

- **Parity lemma.** It is correct: odd-multiplicity zeros are exactly the sign changes, which are even in number on the circle. It is not used in code (LIMIT = 10), as the authors said.
- **Runs not audited here:** bulk, C1p, F1 and F2 are still running. Rerun `d1`/`d2`/`d5` on them when they finish.

## Scripts (`audit/round4/`)

| script | checks |
|---|---|
| `d1_replay_recert.out` | `bb_replay --recert` on 4 runs |
| `d2_wind_spot.py` | 100-digit root moduli against W |
| `d3_wind_stress.py` | synthetic self-inversive stress test |
| `d4_branch_demo.py` | R4-1 demonstration and fixed branch logic |
| `d5_step_margins.py` | all winding leaves: W, step margins, L, \|kL\| |
| `d6_comp_bounds.py` | enclosure and label soundness |
