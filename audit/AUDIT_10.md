# AUDIT 10: the perihelion chart KBPe, `kq_run --peri` (four-part count), the split-angle change

Independent audit, 2026-10-23, 12:14–12:45. Scripts and outputs are in `audit/round10/`. At most 2 Python processes
ran at once.

## Verdict

**No soundness issue was found.**
- The perihelion chart KBPe is the exact reversal of the aphelion polynomial, and it preserves zeros with multiplicity.
- The four parts of `count_peri` cover the circle with overlaps of about 2e-9 rad at every seam. Each part's bound holds against 100-digit zero counts, and so do the totals, in face mode too.
- `certify` always trying `phi_split` and then PHI_0, and re-certification trying every option, are both valid.
- `bb_replay`'s KQ and KQF paths accept leaves made with either option.

The round-9 follow-ups are done:
- R9-1: `launch_runs.sh` now has F1, KQ0, KQ0T, KQF, FE/FM/FIP/FIM and the chain levels, and C3_CORE §6 and §7 are written;
- R9-2 is fixed;
- R9-3 is left as is (harmless).

## (1) Chart KBPe: OK

**The identity.**
- τ = tan((f − π)/2) = −cot(f/2) = −1/tp.
- So tp¹²·Q(−1/tp) = P(f)·tp¹²(1 + 1/tp²)⁶ = P(f)(1 + tp²)⁶.
- The coefficient of tp^j is (−1)^j q_{12−j}. This needs deg Q ≤ 12, which is asserted in `gen_c3.py`.
- f ↦ tp = tan(f/2) is a diffeomorphism of (−π, π) onto ℝ, and (1 + tp²)⁶ > 0, so the real zeros of G are exactly the zeros of P with f ≠ π, with the same multiplicities.
- f = 0 corresponds to tp = 0; aphelion, tp = ∞, is never used (|tp| ≤ tan 0.6 + 1e-9).

**Data** (`h15`):
- The generated coefficients of nbcoef_c3_KBPe and nbcoef_c3_KBPe_2 are exactly the KBQe_2 coefficients reversed with signs. All slots match (value, first and second derivatives), 200 random points, difference 0.
- G(tp) against my P(f)(1 + tp²)⁶ at 50 digits, for e ∈ [0.83, 1], μ ∈ [0, 1] and |tp| ≤ 0.8: relative error 3e-50.
- The derivative slots match numerical derivatives (7e-30).
- The engine `nbeng_c3_KBPe_2` is the make_c3_2 engine audited in round 7 (`g17`), applied to these data.

## (2) `count_peri`: OK

**Code review.**
- Two stages, with arc target 2 and then 0. Each part is a valid upper bound for any rounds, cap and target, because min over rounds is kept.
  - aphelion: KBQe_2 on [−t_hi, t_hi]; in face mode, on [T_X, t_hi] ∪ [−t_hi, −T_X];
  - perihelion: KBPe_2 on [−tp_hi, tp_hi];
  - two arcs: `away_count_arc2` (start, length) = (π + φ₀ − 1e-9, arcl) and (φ_P − 1e-9, arcl), with arcl = π − φ₀ − φ_P + 2e-9.
- The returned total is their sum, compared with LIMIT_P = 11 (parity), or with LIMIT_OUTER = 6 in face mode.
- `away_count_arc2` is line-for-line `away_count_t2` with the arc [a₀ − PAD, a₀ + L₀ + PAD] and `dp_linear`, a linear count, which is correct for an arc.
- The angle conventions agree: the f-chart variable is the true anomaly, and aphelion is at f = π, as in KBQe.

**Coverage** (`h16`, the code's float end points, ARC_PAD on every root interval and arc, PAD inside `away_count_arc2`, exact seam arithmetic):
- The smallest overlap at any of the four seams is 1.97e-9 rad over 2,002 pairs (φ₀, φ_P), including (0.9, 1.2) and (1.25, 1.2).
- At the actual parameters the overlaps are 2.3–2.6e-9.
- For φ₀ + φ_P ≥ π the arcs are dropped, and the two polynomial parts overlap, since 2·atan(tan(x/2) + 1e-9) > x.
- In face mode the aphelion part starts at |τ| = T_X − ARC_PAD, and the face charts reach τ = 3·fl(0.02) > fl(0.06). So that seam also overlaps.

**Per-part soundness at 100 digits** (`h17`). 40 random boxes (e ∈ [0.95, 1), t ∈ [0.5, 1)), parts computed exactly as `count_peri` calls the engines, 4 points per box:
- every part ≥ the zeros of my P in its unpadded region;
- every certified total ≥ N_mult;
- **0 violations** at 160 points. Most part bounds are attained: aphelion 160, perihelion 157, arcs 158 and 160 times.

**Face mode** (`h18`): 25 random KQF boxes outside the C3 cell, 21 certified with total ≤ 6. At 84 points the 100-digit count of zeros with |τ| ≥ T_X never exceeds the total: **0 violations**, maximum 6.

**Documentation.** The ridge example of C3_CORE §7 (e = 0.9925, l = 0.714, ω = 0.001, μ = 0.0103) has N_mult = 10 at 100 digits: 4 zeros within 0.9 of aphelion and 6 within 1.2 of perihelion, as stated. The ω = 0 factorisation B = 2 sin f·[(1 − l²) − l² cos f], with its triple zero at l = 1/√2, was checked by hand.

## (3) Split angle and re-certification: OK

- `certify` tries the peri options in order ([PERI, None]) and, for each, `phi_split(bx)` and then PHI_0. Every combination is a valid count, and the first that certifies is returned; otherwise the minimum.
- **Re-certification (DEEP_ALWAYS).**
  - It adds the other peri option and allows the deep pass.
  - `bb_replay` KQ: `kq_run.init()` plus DEEP_ALWAYS tries peri ∈ {None, 1.2} × φ₀ ∈ {phi_split, 0.9}.
  - `bb_replay` KQF: `init(scale, face=True)` plus DEEP_ALWAYS makes `certify_outer` try {None, 1.2}. KQF is correctly outside PARITY_ENGINES (partial count, so total ≤ recorded bound).
  - Leaves made with or without `--peri` therefore re-certify.
- My final-pass tools (`g1`, `g9`) call the same `certify`/`certify_outer` under DEEP_ALWAYS, so they accept these leaves too. No change was needed.

## Not checked

- The authors' efficiency claims (ridge boxes about 100 times wider; the band counts with φ₀ = 1.25).
- The KQK1 r4+ files written with `--peri`: none was complete at audit time. They are covered by the final pass.

## Round 10b (12:38–12:50): weighted split rule, face 'o' pruning: OK

**Split rule 'face/K'.**
- `is_face` and `face_k` parse the rule. In rule (2) the ε and μ widths are divided by K (s2 is multiplied by K). Case (1), the face cuts, is unchanged.
- Any tiling split is a valid certificate.
- Check (`h20`, 20,000 random boxes near the face):
  - 0 invalid tilings for 'face/30', 'face' and None;
  - 'face' and None are bit-identical to `splitrule_v1_backup.py`;
  - 'face/30' cuts l or ω where 'face' cut e or t in 4,691 boxes, as intended.
- My `g2` replays each file with its own header rule and checks the tiling of every split, so it needs no change.
- KQB r5 and KQF r1 were restarted with fresh records files, and the old ones moved to `runs/superseded/`. This keeps one rule per file, which is what the replay assumes; `g2` only picks up `NAME.rK.json` in `runs/`.

**Face charts: 'o' when s_lo > SF.**
- This is sound by the same argument as the excuse rule (AUDIT_9 §1): a box whose float s_lo exceeds fl(0.02) contains no needed point.
- `check_all.prune_ok` now verifies the face-chart 'u' and 'o' leaves (outward factors 1 ∓ 1e-12 on l and ω).
- One small slack (info): `prune_ok` tests 'u' as P_z,lo² + Q_z,lo² > 1.0 without the 1e-12 margin that `face_bounds` uses. The verifier is therefore laxer than the engine by about 1 ulp. No leaf is affected, because the engine's margin is stricter, and my `g2` checks 'u' in exact arithmetic.
- My `g2` face 'o' predicate now accepts C3-cell boxes (exact) or s_lo > 0.02.
- Replay of the current base files: FE, FM, FIP, FIM, KQF and KQB all have 0 problems and 0 pruning violations.
  - FIM is complete.
  - The others still have open boxes for their running levels: FE 73, FM 67, FIP 90 (after excusing 124,446 boxes with s_lo > 0.02), KQF 309, KQB 218.

**kq_run** (`--splitrule=` parsing, docstring): no certificate code changed.
