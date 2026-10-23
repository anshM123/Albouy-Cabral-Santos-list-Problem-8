# AUDIT 13: the tau-winding certificate (TW, leaf code '7')

Independent audit, 2026-10-23, 16:58–17:20. At most 2 Python processes. Scripts and outputs are in `audit/round13/`.

I audited `code/tau_wind.py` (16:58 version, which adds TW on the C3 charts), its integration in `kq_run.py`,
`c3/c3_run.py`, `check_all.py`, `bb_replay.py` and `verify_mp.py`, and every TW leaf written so far by the eight
levels running with `--tw`: KQ0T.r4, KQ2.r5, KQB.r7, KQF.r2, KQF2.r3, KQK1.r7, KQS.r2 and c3 Ca.r8.

## Verdict

**TW is sound.** The argument checks out and the code implements it. My independent 60-digit checks of 244 TW leaves
from all eight levels, and my adversarial tests, found no failure.

One checking-tool gap matters for the final verification (R13-1): no project tool re-certifies the C3 TW leaves
(Ca.r8).

| # | item | status |
|---|---|---|
| (a) | counting identity N_P = 12 − #nonreal(Q), including degree drops at either end | OK |
| (b) | homotopy invariance, including a vanishing leading coefficient | OK |
| (c) | disc bounds, T at complex z, float pads | OK |
| (d) | the argument-step logic | OK |
| (e) | check_tw re-certifies | OK (244 of 244) |
| (f) | independent spot checks on new run files | OK |
| C3 | TW on the C3 charts (Ca.r8; added 16:58) | OK as a certificate; see R13-1 and R13-2 |
| R13-1 | C3 TW leaves are not re-certified by any project tool (bb_replay handles KQ/KQF only) | **ISSUE (major, checker)** |
| R13-2 | check_all accepts which = 1 on C3 contours (only 0 is meaningful) | ISSUE (minor) |
| R13-3 | TW_CHARTS lists CK, CE and Cm, but they have no second-order engines | note |

## (a) Counting identity: OK

**The algebra.** With z = e^{if}:
- th = tan((f − π)/2) = −i(z + 1)/(z − 1);
- tp = tan(f/2) = −i(z − 1)/(z + 1).

Q(th) = P(f)(1 + th²)⁶ = Σ_k c_k (−1)^k (1 + i th)^{6+k} (1 − i th)^{6−k}. This is a real polynomial of formal degree 12.

- The real zeros of P on the circle other than f = 0 correspond one to one, with multiplicity, to the real zeros of Q.
- The multiplicity at f = 0 equals the order of the zero of G(tp) = tp¹² Q(−1/tp) at tp = 0. That is the number of
  vanishing top coefficients of Q, i.e. 12 − deg Q.
- So N_P = #real(Q) + 12 − deg Q = 12 − #nonreal(Q).
- The same holds for G, with the degree drop at f = π.
- A zero strictly inside a closed disc in the open upper half-plane is non-real, and its conjugate is another one.
  Hence N_P ≤ 10.
- Q ≢ 0 is implied, because Q has no zero on the contour.

**Numerical check (j2 a).** I built Q and G independently, by exact interpolation of P(f(s))(1 + s²)⁶ at 13 nodes.
The identity held at 41 points with 0 failures:
- 25 random points;
- 8 points with Q_z = 0, where f = 0 is a double zero: deg Q = deg G = 10;
- 8 points with e = 1, where f = π is a double zero: deg G = 10.

## (b) Homotopy: OK

The coefficients are polynomials in the box variables, and the box is connected. If F(·; x) ≠ 0 on the contour for
every x in the box, then (1/2πi)∮F′/F is a continuous integer, hence constant.

Zeros that go to infinity when the leading coefficient vanishes stay outside the bounded contour, so they don't affect
the count.

For the C3 charts, ν = 0 (the blow-up face) is part of the box. It is unphysical (ε = ν²E = 0, so e = 1), and the
homotopy simply includes it.

## (c) Disc bounds and pads: OK

The data come from `box_coeffs2`, the round-7-audited second-order enclosure:
- cmid and crad (= 1e-13 · |a|_abs);
- dcen, and dwid (the mean-value width plus 1e-13 · a1);
- hcen;
- rem (third-order Faà di Bruno remainder plus 1e-13 pads);
- qw and rr.

Every one of these is **coefficient-wise**. So |Σ_k δa_k w^k| ≤ Σ_k |δa_k| ρ^k with ρ = |z_j| + D ≥ |w| holds verbatim
for complex w.

Each term of `arc_radius`:
- **old** = Σ_v rr_v (sup_disc |∂_v F(·; x_c)| + Σ_k dwid ρ^k): the mean-value form. The disc maximum Σ_m |b_m| D^m of
  each derivative polynomial is exact.
- **new** = Σ_v rr_v sup|∂_v F| + Σ_m qw_m sup|H_m| + Σ_k rem_k ρ^k: the second-order form.
- **T** = rd + min(old, new), with relative pads.
- **dev** = Σ_{m≥1} |b_m| D^m: the exact Taylor expansion of F~(·; x_c) at z_j.

**Float errors.**
- The Taylor coefficients b_m are computed in complex floats. Their total absolute-term mass is Σ_k |cmid_k| ρ^k = mag,
  so the errors are about 30u·mag, covered by 1e-12·mag.
- The disc maxima inside 'new' are covered by 1e-12·pabs.
- Those inside 'old' are covered by dwid's own 1e-13·a1 ≥ 1e-13 |dcen| term (≫ 30u). So 'old' needs no extra pad.
- v0 > R in floats implies |b0| > R exactly, since the abs() error of about 1e-16·|b0| ≤ 1e-16·mag is far below the
  1e-12 pad.

**Arc coverage.** D = r h (1 + 1e-12) + 1e-12 (|c| + r) covers the arc (chord ≤ r h). It also covers float gaps of
about 1e-16 rad between adjacent arcs and the error in z_j.

**Empirical test of the inequality |F(w; x) − b0| ≤ R**, with my own F at 60 digits; its normalisation agrees with the
engines' polynomials to float-evaluation level, 6e-16 to 7e-11:
- **240 TW leaves** (j1): 19,200 checks at box corners and random points, w on and inside the disc boundary, arcs of
  half-width π/64 and π/16. Largest ratio 0.996.
- **The 4 G-chart leaves:** 480 checks, largest 0.848.
- **Adversarial (j2 c):** 12 random parameter boxes of width up to 0.05 (much larger T) with discs up to D = 0.4, 360
  checks. Largest 0.999.

**No violations.** Ratios near 1 are expected: on the disc boundary |b_1|·D ≈ dev is attained.

## (d) Argument steps: OK

The arcs are contiguous: the initial grid uses centres 2π(j + ½)/M0 with half-width π/M0, and a refinement splits an
arc into its two halves, in order. The shared end point p_j of arcs j and j+1 therefore lies in both discs.

Within disc j, F stays in the sector of half-angle φ_j = arcsin(R_j/|b0_j|) < π/2 about arg b0_j. So the lifted
difference arg b0_{j+1} − arg b0_j, measured through F(p_j), lies within (−(φ_j + φ_{j+1}), φ_j + φ_{j+1}) ⊂ (−π, π).
It therefore equals the principal Arg(b0_{j+1}/b0_j), and the sum telescopes to 2πW.

Float margin: v0 > R in floats gives R/v0 ≤ 1 − 2⁻⁵², hence φ ≤ π/2 − 2.1e-8. Each true step is thus at least 4.2e-8
away from ±π, while the division and atan2 errors are about 1e-15. The rounding of tot/2π is safe up to 65,536 arcs.

The wrap from the last arc to the first is included (vals[(j+1) % n]).

## (e) check_tw re-certifies: OK

- The function is deterministic, with more rounds and a larger cap than certify_tw, so any certified contour
  re-certifies along the same refinement path.
- 244 of 244 sampled leaves re-certified (j1), including 30 Ca chart leaves via check_tw_chart.
- bb_replay calls check_tw for every '7' leaf of KQ/KQF files. The R11-2 header init stays in place.

## (f) Independent spot checks on the new run files: OK

All TW leaves so far replay exactly, with '7' count = len(tw) and every contour satisfying c_im > r > 0:

| level | TW leaves at the time |
|---|---|
| KQK1.r7 | 3,044 |
| KQB.r7 | 1,831 |
| KQ2.r5 | 2,801 |
| KQ0T.r4 | 1,547 |
| KQF.r2 | 1,627 |
| KQF2.r3 | 618 |
| KQS.r2 | 709 |
| Ca.r8 | 2,365 |

Of about 20,000 contours, 4 use G (all in KQK1.r7); the rest use Q.

**At 60 digits with my own P** (FFT, degree-12 roots, Möbius map, and th = νy for Ca), on 30 leaves per level plus all
4 G-chart leaves:
- the box centre and 4–6 random points per leaf;
- exactly 1 zero strictly inside the stored contour at every point (constant over the box, as the homotopy says);
- N_P ≤ 10 everywhere;
- the smallest relative distance of any zero to the contour is 0.286 (Ca) and at least 0.48 for the KQ files.

**Negative tests (j2 b), 0 wrongly certified:**
- certify_tw (10 rounds, 6 candidate contours) fails on all 8 boxes that contain a degenerate corner e = 1, μ = 0 with
  N_P = 12 (l ≥ l*(ω) at ω = 0.22, 0.62, 1.05, 1.31).
- check_tw rejects a contour moved off the zero (W = 0) and a mirrored one in the lower half-plane.

## TW on the C3 charts (c3_run --tw, Ca.r8)

**Validity.** Chart Ca has G(y) = Q(νy)/ν¹⁰ with (Q_z, κ, ε, μ) = (ν, ν²K, ν²E, ν³M). This is my g4 map, and the
engine's polynomial agrees with it to float-evaluation level (6.7e-11 relative).

A non-real zero y of G at ν > 0 gives the non-real zero th = νy of Q. Every physical point of the chart has ν > 0,
since ε > 0 forces ν > 0. So a Ca TW leaf certifies N_P ≤ 10 on the whole circle for every physical point of its
image, independently of the near / intermediate / away (Ca + Phi + KA) decomposition.

The DEG of Ca's engine is 12, consistent with BINOM (13 × 13).

### R13-1 (major, checker)

No project tool re-certifies the c3 '7' leaves:
- bb_replay handles '7' only for KQ/KQF engines;
- final_verify.sh has no c3 re-certification step;
- check_all checks only structure and contour validity;
- verify_mp only spot-checks N ≤ 10 at sample points.

So a wrong stored contour in a c3 file would pass the project's verification. My g4 now re-certifies a sample with
`tau_wind.check_tw_chart` and spot-checks the zeros inside the contour (j3).

Fix: re-certify every c3 '7' leaf with `check_tw_chart` in the final verification (for example a c3 branch in
bb_replay or a small c3 replay).

### R13-2 (minor)

check_all accepts which ∈ {0, 1} for contours of the C3 engines. Only 0 (the chart variable y) is meaningful, and
check_tw_chart ignores 'which', so this is harmless. Tighten it to which == 0 for Ca, CK, CE and Cm.

### R13-3 (note)

TW_CHARTS = (Ca, CK, CE, Cm), but nbeng_c3_{CK,CE,Cm}_2 do not exist. `c3_run --tw` on those charts would fail
loudly at import. Only Ca uses TW now.

## Tool updates (for the final pass)

**g2:**
- '7' gives bound 10 for engines KQ, KQF, Ca, CK, CE and Cm;
- the contour is consumed in tree order and validated (which ∈ {0, 1} for KQ/KQF, 0 for C3; c_im > r > 0);
- the count is checked against len(tw);
- the record-hist cross-check uses kq_run's key 't';
- tw_leaves is reported per family.

**g1:**
- TW leaves are re-certified with check_tw;
- 100-digit spot checks: ≥ 1 zero inside the contour, the same count at every point, N ≤ 10.

**g3:** '7' leaves are checked as whole-circle leaves (≤ 10), also in KQF files.

**g9:** a new kind 'tw' re-certifies with check_tw. The final-pass scripts now run `--kinds reg,deep,tw`.

**g4:**
- '7' leaves of the C3 charts are re-certified with check_tw_chart;
- 100-digit spot checks of the zeros inside the contour (y = th/ν) and N ≤ 10;
- AUDIT_SCALE for header-less continuation records.

**g3:** the engine lookup now also works for continuation names (NAME.rK).

All of these were tested on live records (j3: KQK1.r7 with g1, KQF2.r3 with g3, Ca.r8 with g4). The round-12
versions are kept in `round13/*_r12.py`.

## Not checked

- The levels are still running; I sampled the TW leaves present by about 17:10.
- C3_CORE §8 is not written yet.
- The engine polynomial data (`box_coeffs2`) rest on AUDIT_7's audit of make_c3_2 and gen_c3_2; TW reuses them
  unchanged.
