# AUDIT 9: the face blow-up at e = 1, μ = 0 (R8-1), the face split rule, the restricted bulk, and engine KQ away from e = 1

Independent audit, 2026-10-23, 11:40–12:30. Scripts and outputs are in `audit/round9/`. The final-pass tools in
`audit/round7/` were extended (see the end).

**Rules kept.** At most 2 Python processes ran at once, and scripts were kept light. bulk.json (59 MB) was not loaded;
it is part of the final pass.

## Verdict

**No soundness issue was found.**
- **The face construction** resolves R8-1 correctly. It consists of the charts FE/FM (near), FIP/FIM (gap) and the outer run KQF, together with the KQ pruning code 'p'.
- **The excuse rule** for face-chart boxes with s_lo > SF is sound.
- **The 'face' split rule** always produces a valid tiling.
- **Engine KQ** has no hidden assumption of e near 1 or ε small, so KQ0 and KQ0T are valid replacements for the excluded bulk slab.
- **Completeness** still depends on the running levels: FE, FM, FIP, KQF, KQ0, KQ0T and the KQK1/KQK2/KQB continuations.
- **Docs:** C3_CORE §6, the D3 row of DEGENERACY_TREE, the PROOF atlas table and `launch_runs.sh` do not yet describe the round-9 runs.

| id | severity | finding |
|---|---|---|
| R9-1 | minor (docs, reproducibility) | `launch_runs.sh` (09:48) has no lines for F1, KQ0, KQ0T, KQF, FE/FM/FIP/FIM, the `--facesplit` levels, or the bulk `--exclude`. The atlas documents still list D3 as resolved by KQ. C3_CORE §6 is announced but not written. |
| R9-2 | info | The `splitrule` docstring says the face-side child "then lies in the face neighbourhood". When both e* and t* fall inside the box it takes two cuts (e* first, then t*). Tested: every corner child is prunable within two cuts. |
| R9-3 | info | FIP has many failed leaves: 96,142 at 12:00, all with bound 1. All have σ_lo ≥ 0.020007 > SF, so the excuse rule covers them. The thin margin slab σ ∈ (SF, SF·R1] costs time but is never needed. Ending the chart domain at SF would avoid this; the overlap with KQF is supplied by the arc padding anyway. |

## (1) The face blow-up (R8-1): OK

**Exactness** (`h9`, my own Q(τ) = P(f)(1+τ²)⁶ from P = μ²k²sin²fσ² − D²B², sympy):
- The lowest power of s is exactly 6 for FE (e = 1 − s², μ = s²m, τ = sy), FM (e = 1 − s²E, μ = s², τ = sy) and FIP/FIM (τ = ±3s, e = 1 − s²E, μ = s²m). So G = Q/s⁶ is a polynomial, and a positive multiple of Q for s > 0.
- With c² = 1 − P_z², the leading parts are:
  - FE: 16c²y²[4m² − c²(1+2y²)²];
  - FM: 16c²y²[4 − c²(E+2y²)²];
  - FIP/FIM: 144c²[4m² − c²(E+18)²].
- So at leading order the near part has at most 4 real zeros: the double zero y = 0, plus at most one real pair. The gap part has no zero whenever 324c² > 4·R1², i.e. P_z < 0.9938. Outside the C3 cell P_z ≤ max(0.95, sin 9π/20) = 0.988, so this holds with margin.

**Generated modules.** Both versions of nbcoef_c3_{FE,FM,FIP,FIM} (first- and second-order) agree with my G = Q/s⁶ at 50 digits, to relative error ≤ 5e-40.

**Coverage argument**, for a physical point p with ε, μ ≤ SF² (s = max(√ε, √μ) ≤ SF = 0.02), on the τ-line:
- **near, |τ| ≤ 3s:** FE (ε ≥ μ: s = √ε, m = μ/ε ≤ 1) or FM (μ ≥ ε: E = ε/μ ≤ 1), y = τ/s ∈ [−3, 3]. Target 5.
- **gap, 3s ≤ |τ| ≤ T_X = 0.06:** σ = |τ|/3 ∈ [s, SF], E = ε/σ² ≤ 1, m = μ/σ² ≤ 1. These lie in FIP (τ > 0) or FIM (τ < 0), which certify "no zero" for every box point and every m ∈ [0, R1].
- **outer, |τ| ≥ 0.06:** KQF counts the near polynomial on [T_X, T_HI] and [−T_HI, −T_X], plus the f-chart away arc from φ₀ − 1e-9. Target 6.
- **Total.** All three parts count zeros of the same P with multiplicity and overlap at their seams, through ARC_PAD at both interval ends. That gives ≤ 5 + 0 + 6 = 11, hence ≤ 10 by parity.
- **Float seam.** The needed σ ≤ fl(0.06)/3 < fl(0.02), and KQF starts at fl(0.06) − 1e-12.

**Pruning predicates** (code review; my replay re-checks every pruned leaf in exact arithmetic):
- **KQ 'p'** (`in_face`): 1 − e_lo ≤ SF²(1 − 1e-9) and μ(t_lo) ≤ SF²(1 − 1e-9), so the whole box is in the face neighbourhood. KQF's domain covers it: e ≥ 0.9996 and t ≥ 0.99999983 hold for these points, and KQF has l ∈ [0, 1], ω ∈ [0, π/2] (read from its records).
- **KQF 'o':** l_lo ≥ 0.95 and ω_lo ≥ 9π/20. The face-chart 'o' is l_lo ≥ 0.95 and atan2(P_z,lo, Q_z,hi) ≥ 9π/20, and 'u' is P_z,lo² + Q_z,lo² > 1. These are C3-cell or unphysical boxes; the face neighbourhood lies inside the C3 cell's e and μ ranges, so the C3 core plus KA cover them.
- **check_all** requires FE, FM, FIP, FIM and KQF complete before accepting any family that uses 'p'. My `g2` enforces the same.

**Excuse rule** (open or failed face-chart boxes with s_lo > SF are dropped).
- Needed points have s ≤ SF·√(1 − 1e-9) < SF in FE/FM, and σ ≤ fl(0.06)/3 < fl(0.02) in FIP/FIM.
- A box whose float lower bound exceeds fl(0.02) therefore contains no needed point.
- The rule is sound, and `next_level.py` applies the same test.

**Data so far:**
- FIM is complete: 512 records, 57,097 certified leaves (all 0 zeros) and 57,976 pruned leaves; my replay with exact pruning predicates OK. `g4`: 40 points at 100 digits, no zero; 40 leaves re-certified with `c3_run.certify`, 0 problems.
- FIP (12:00): 96,142 failed leaves, all with s_lo ≥ 0.020007. Of its open boxes, 67 have s_lo ≤ 0.02; they are still running.
- KQF (partial): leaves with totals 2–6, mostly 6. This matches the three double zeros of B at the corner.

**Empirical check of the three regions** (`h14`): 150 random face-neighbourhood points, ε and μ log-uniform in [1e-12, 4e-4], outside the C3 cell, counted with my P at 100 digits:
- max near 4 (target 5);
- max gap 0;
- max outer 6 (target 6, attained);
- max total 10, all even.

## (2) Split rule 'face': OK

- `splitrule.split_box`:
  - case (1): cut at e* = 1 − S2 or t* = 1/(1 + S2²), strictly inside the box, when the box meets {ε < S2, μ < S2};
  - case (2): bisect the largest blow-up-scaled width.
- **Unit test** (`h13`, 20,000 random boxes near the face):
  - 0 invalid tilings (exactly one coordinate cut, strictly inside, children sharing the cut);
  - 4,498 face cuts, each leaving a corner child that is 'p'-prunable within two cuts (exact check);
  - rule None reproduces the old bisection exactly (5,000 of 5,000).
- **Validity.** Any tiling split gives a valid certificate. `check_all.replay` and `bb_replay` call `split_box` when the header has `split_rule`. My `g2` now replays each file with its own rule and checks the tiling of every split.
- No finished file uses the rule yet (KQK1/KQK2 r4, KQB r5 are running).

## (3) Restricted bulk and the new runs: OK

**bulk.json** (header only, read without loading the file):
- exclude = [[0.9902515840381316, 1.2535563099909826] × [0,1] × [0,π/2] × [0,1]]. These are the last two η-slices of the 8-slice grid, i.e. e ≥ sin 0.9902516 = 0.836164.
- 2,000 skipped boxes; 6,000 completed records; 0 incomplete; 0 failed.
- `bb_run --resume` keeps one completed record per non-excluded initial box. The grid check in my final pass (roots ∪ skipped = grid, roots ∩ skipped = ∅) will confirm that the slab records were dropped.

**Replacement over e ∈ [0.835, 0.9501]:**
- F1 (f-chart, `--order2`, t ≤ 0.5), KQ0 (t ∈ [0.5, 0.9975]) and KQ0T (t ≥ 0.9975).
- Each family covers its whole t-range in its own parametrisation. The seams are in e only: 0.835 < 0.836164, and 0.9501 > 0.9499 (F0/KQ1/KQ2).
- **F1 is complete:** 400 records, 1,093 leaves (329 winding), 0 fails, own replay OK. 160 leaves re-certified, 0 mismatches; 152 of them need the second-order engine, as expected for an `--order2` run. 90 spot points at 100 digits: 0 violations, max 6.

**Engine KQ away from e = 1** (the main question): **no hidden assumption.**
- The near polynomial KBQe is `corner_Q()` with e = 1 − ε substituted. It is an exact polynomial identity, with no division and no truncation.
- My check (`h9`): nbcoef_c3_KBQe_2 against my Q(τ) at e ∈ [0.83, 0.96], μ ∈ [0, 1], |τ| ≤ 0.75 gives relative error 3e-50.
- The away arc and the winding certificate use the exact f-chart polynomial.
- The box conversion (outward-rounded ε, P_z, Q_z, μ) is general. T_MIN = 0.5 keeps μ ≤ 1. 'p' never fires (ε ≥ 0.0499).
- At t = 1 with e ≤ 0.9501, D ≥ 0.0499 > 0, so P = −D²B² has at most 8 real zeros, and KQ0T has no corner problem.
- **Stress test** (`h11`, 100 digits): 60 random boxes in the KQ0/KQ0T domains, 51 certified. Each was checked at 8 points: **0 violations**.

## Final-pass tools (updated)

| script | change |
|---|---|
| `g2` | F1, KQ0, KQ0T (header domain), KQF (codes 'a'..'g', 'o'), FE/FM/FIP/FIM (target 5/5/0/0, exact 'u'/'o'); bulk exclusion; per-file `split_rule` with a tiling check; exact 'p'; the s_lo > SF excuse; the face-families-complete requirement for 'p' |
| `g1` | 'p'/'o' codes (exact); KQF re-certification via `certify_outer`, plus a 100-digit outer-count check |
| `g3` | KQF outer count; split rule |
| `g4` | FE/FM/FIP/FIM physical maps (near count on y ∈ [−3, 3]; FIP/FIM no zero for m ∈ [0, R1]) |
| `g9` | KQF; split rule |
| `g13` | split rule; KQF skipped (no winding) |
| `final_pass.sh` | F1, KQ0, KQ0T, KQF and the face charts added; continuation files of every family are picked up automatically |

## Pending

- **Running:** FE, FM, FIP (levels), KQF, KQ0, KQ0T, and the KQK1/KQK2/KQB continuations with the face split. Then the final pass, including the bulk replay, and the final verdict.
- **Docs:** R9-1.
